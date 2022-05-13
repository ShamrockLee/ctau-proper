{ runCommand
, lib
, stdenv
, storeDir ? builtins.storeDir
, writeScript
, singularity
, writeReferencesToFile
, bash
, vmTools
, gawk
, util-linux
, runtimeShell
, e2fsprogs
, writeText
, writeShellScriptBin
, coreutils
}:

let
  # May be "apptainer" instead of "singularity"
  projectName = singularity.projectName or "singularity";
in
rec {
  shellScript = name: text:
    writeScript name ''
      #!${runtimeShell}
      set -e
      ${text}
    '';

  writeMultipleReferencesToFile = paths: runCommand "runtime-deps-multiple" {
    referencesFiles = map writeReferencesToFile paths;
  } ''
    touch "$out"
    declare -a paths=();
    for refFile in $referencesFiles; do
      while read path; do
        isPathIncluded=0
        for pathIncluded in "''${paths[@]}"; do
          if [[ "$path" == "$pathIncluded" ]]; then
            isPathIncluded=1
            break
          fi
        done
        if (( ! isPathIncluded )); then
          echo "$path" >> "$out"
          paths+=( "$path" )
        fi
      done < "$refFile"
    done
  '';

  buildImage =
    { name
    , contents ? [ ]
    , diskSize ? 1024
    , memSize ? 512
    , runScript ? "#!${stdenv.shell}\nexec /bin/sh"
    , runAsRoot ? null
    , executableFlags ? [ ]
    , buildImageFlags ? [ ]
    , singularityPackage ? singularity
    }:
    let
      # May be "apptainer" instead of "singularity"
      projectName = singularityPackage.projectName or "singularity";
      runAsRootFile = shellScript "run-as-root.sh" runAsRoot;
      runScriptFile = shellScript "run-script.sh" runScript;
      result = vmTools.runInLinuxVM (
        runCommand "${projectName}-image-${name}.img"
          {
            buildInputs = [ singularity e2fsprogs util-linux gawk ];
            layerClosure = writeMultipleReferencesToFile (contents ++ [ bash runScriptFile ]);
            preVM = vmTools.createEmptyImage {
              size = diskSize;
              fullName = "${projectName}-run-disk";
            };
            inherit memSize executableFlags buildImageFlags;
          }
          ''
            rm -rf $out
            mkdir disk
            mkfs -t ext3 -b 4096 /dev/${vmTools.hd}
            mount /dev/${vmTools.hd} disk
            mkdir -p disk/img
            cd disk/img
            mkdir proc sys dev

            # Run root script
            ${lib.optionalString (runAsRoot != null) ''
              mkdir -p ./${storeDir}
              mount --rbind ${storeDir} ./${storeDir}
              unshare -imnpuf --mount-proc chroot ./ ${runAsRootFile}
              umount -R ./${storeDir}
            ''}

            # Build /bin and copy across closure
            mkdir -p bin ./${storeDir}
            for f in $(cat $layerClosure) ; do
              cp -ar $f ./$f
            done

            for c in ${toString contents} ; do
              for f in $c/bin/* ; do
                if [ ! -e bin/$(basename $f) ] ; then
                  ln -s $f bin/
                fi
              done
            done

            # Create runScript and link shell
            if [ ! -e bin/sh ]; then
              ln -s ${runtimeShell} bin/sh
            fi
            mkdir -p .${projectName}.d
            ln -s ${runScriptFile} .${projectName}.d/runscript

            # Fill out .${projectName}.d
            mkdir -p .${projectName}.d/env
            touch .${projectName}.d/env/94-appsbase.sh

            cd ..
            mkdir -p /var/${projectName}/mnt/{container,final,overlay,session,source}
            echo "root:x:0:0:System administrator:/root:/bin/sh" > /etc/passwd
            echo > /etc/resolv.conf
            TMPDIR=$(pwd -P) ${singularityPackage}/bin/singularity $executableFlags build $buildImageFlags $out ./img
          '');

    in
    result;

  inherit (import ./definition-generator.nix { inherit lib; })
    allPrimarySectionNames
    allAppSectionNames
    toDefinitionString
    ;

  buildImageFromDef =
    args@{ name
    , contents ? [ ]
    , definitionOverrider ? null
    , executableFlags ? [ ]
    , buildImageFlags ? [ ]
    , singularityPackage ? singularity
    , ...
    }:
    let
      # May be "apptainer" instead of "singularity"
      projectName = singularityPackage.projectName or "singularity";
      layerClosure = writeMultipleReferencesToFile (contents ++ [ bash coreutils ]);
      definition = if (args.definition or null != null) then args.definition else
      (
        if lib.isFunction definitionOverrider then
          definitionOverrider
        else if builtins.isAttrs definitionOverrider then
          (d: lib.recursiveUpdate d definitionOverrider)
        else
          lib.id
      ) {
        header.Bootstrap = "scratch";
        setup = ''
          mkdir -p ''${SINGULARITY_ROOTFS}/${storeDir}
          for f in $(cat ${layerClosure}) ; do
            cp -ar "$f" "''${SINGULARITY_ROOTFS}/${storeDir}"
          done
          mkdir -p "''${SINGULARITY_ROOTFS}/bin"
          "${coreutils}/bin/ln" -s "${runtimeShell}" "''${SINGULARITY_ROOTFS}/bin/sh"
          mkdir -p "''${SINGULARITY_ROOTFS}/usr/bin"
          "${coreutils}/bin/ln" -s "${coreutils}/bin/env" "''${SINGULARITY_ROOTFS}/usr/bin/env"
        '';
        environment = {
          PATH = "${lib.makeBinPath contents}:\${PATH:-}";
        };
      };
      definitionFile = writeText "${name}.def" (toDefinitionString definition);
      # Pass for users who want to build from the command line instead of inside a VM.
      buildscriptPackage = (writeShellScriptBin "build-image" ''
        if [ "$#" -lt 1 ]; then
          echo "Expect IMAGE_PATH" >&2
          exit 1
        fi
        pathSIF="$1"
        shift
        ${lib.toUpper projectName}ENV_PATH="$PATH" "${singularityPackage}/bin/singularity" ${toString executableFlags} build ${toString buildImageFlags} "$@" "$pathSIF" "${definitionFile}"
      '') // {
        meta.mainProgram = "build-image";
      };
    in
    (runCommand "${projectName}-image-${name}.img" (removeAttrs args [ "name" "definition" "definitionOverrider" ] // {
      inherit executableFlags buildImageFlags;
      passthru = args.passthru or { } // {
        inherit singularity layerClosure definition definitionFile buildscriptPackage;
      };
    }) ''
      ${lib.toUpper projectName}ENV_PATH="$PATH" "${singularityPackage}/bin/singularity" $executableFlags build $buildImageFlags "''${buildImageFlagsArray[@]}" "$out" "${definitionFile}"
    '');

  buildImage' =
    args@{ name
    , diskSize ? 1024
    , memSize ? 512
    , passwdString ? "root:x:0:0:System administrator:/root:/bin/sh"
    , resolvString ? ""
    , ...
    }:
    vmTools.runInLinuxVM ((buildImageFromDef (removeAttrs args [ "diskSize" "passwdString" "resolvString" ] // {
      buildInputs = args.buildInputs or [ ] ++ [ singularity e2fsprogs util-linux gawk ];
      preVM = vmTools.createEmptyImage {
        size = diskSize;
        fullName = "${projectName}-run-disk";
        destination = "./${projectName}-disk-image";
      };
      inherit memSize;
    })).overrideAttrs (oldAttrs: {
      buildCommand = ''
        mkdir -p /var/${projectName}/mnt/{container,final,overlay,session,source}
        echo ${lib.escapeShellArg passwdString} > /etc/passwd
        echo ${lib.escapeShellArg resolvString} > /etc/resolv.conf
        rm -rf "$out"
        export SINGULARITY_TMPDIR="$(mktemp -d -p "$TMPDIR" "build-temp-XXXXXX")"
        mkfs "/dev/${vmTools.hd}"
        mount "/dev/${vmTools.hd}" "$SINGULARITY_TMPDIR"
      '' + oldAttrs.buildCommand;
    }));
}
