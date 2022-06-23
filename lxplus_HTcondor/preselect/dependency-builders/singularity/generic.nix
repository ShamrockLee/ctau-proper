{ lib
, buildGoModule
  # Native build inputs
, makeWrapper
, util-linux
, which
  # Build inputs
, bash
, coreutils
, cryptsetup
, fakeroot
, go
, gpgme
, libuuid
  # This is for nvidia-container-cli
, nvidia-docker
, openssl
, squashfsTools
  # Configurations
, pname
, version
, src
, projectName # "apptainer" or "singularity"
  # These cannot be correctly override with `overrideAttrs`
, vendorSha256 ? null
, deleteVendor ? false
, proxyVendor ? false
, doCheck ? true
, enableNvidiaContainerCli ? true
  # Whether the configure script treat SUID support as default
, defaultToSuid ? true
  # Whether to compile with SUID support
, enableSuid ? false
, starterSuidPath ? null
  # Remove the symlinks to `singularity*` when projectName != "singularity"
, removeCompat ? true
, extraDescription ? ""
, extraMeta ? { }
}:

buildGoModule {
  inherit pname version src;

  # All three projects provide vendored source tarball
  # so we set vendorSha256 = null;
  # One who use unvendored source can override this
  # with `override`
  inherit vendorSha256 deleteVendor proxyVendor;

  # go is used to compile extensions when building container images
  allowGoReference = true;

  strictDeps = true;

  passthru = {
    inherit
      enableSuid
      projectName
      removeCompat
      ;
  };

  nativeBuildInputs = [
    makeWrapper
    util-linux
    which
  ];

  buildInputs = [
    bash # To patch /bin/sh shebangs.
    cryptsetup
    gpgme
    libuuid
    openssl
  ];

  configureScript = "./mconfig";

  configureFlags = [
    "--localstatedir=/var/lib"
    "--runstatedir=/var/run"
  ]
  ++ lib.optional (defaultToSuid && !enableSuid) "--without-suid"
  ++ lib.optional (!defaultToSuid && enableSuid) "--with-suid"
  ;

  postPatch = ''
    if ! [ -e .git ] || ! [ -e VERSION ]; then
      echo "${version}" > VERSION
    fi
    # Patch shebangs for script run during build
    patchShebangs --build "$configureScript"
    patchShebangs --build makeit
    patchShebangs --build e2e
    patchShebangs --build scripts
    patchShebangs --build mlocal/scripts
    # Patching the hard-coded [Dd]efaultPath
    for subPath in cmd/internal/cli/actions.go internal/pkg/util/env/env.go e2e/env/env.go; do
      if [ -f "$subPath" ]; then
        echo "Patching [Dd]efaultPath in \"$subPath\""
        sed -i 's|^\([ \t]*[Dd]efaultPath[ \t]*[:]\?=[ \t]*\)"[^"]*"|\1"${lib.makeBinPath [ bash coreutils squashfsTools ]}"|' "$subPath"
      fi
    done
  '';

  postConfigure = ''
    # Code borrowed from pkgs/stdenv/generic/setup.sh configurePhase()

    # set to empty if unset
    : ''${configureFlags=}

    # shellcheck disable=SC2086
    $configureScript -V ${version} "''${prefixKey:---prefix=}$prefix" $configureFlags "''${configureFlagsArray[@]}"

    # End of the code from pkgs/stdenv/generic/setup.sh configurPhase()
  '';

  buildPhase = ''
    runHook preBuild
    make -C builddir -j"$NIX_BUILD_CORES"
    runHook postBuild
  '';

  installPhase = ''
    runHook preInstall
    make -C builddir install LOCALSTATEDIR="$out/var/lib"
    runHook postInstall
  '';

  postFixup = ''
    substituteInPlace "$out/bin/run-singularity" \
      --replace "/usr/bin/env ${projectName}" "$out/bin/${projectName}"
    wrapProgram "$out/bin/${projectName}" \
      --prefix PATH : "${lib.makeBinPath [ fakeroot ]}"
    # Explicitly configure paths in the config file
    sed -i 's|^# mksquashfs path =.*$|mksquashfs path = ${lib.makeBinPath [squashfsTools]}/mksquashfs|' $out/etc/${projectName}/${projectName}.conf
    sed -i 's|^# unsquashfs path =.*$|unsquashfs path = ${lib.makeBinPath [squashfsTools]}/unsquashfs|' $out/etc/${projectName}/${projectName}.conf
    sed -i 's|^# cryptsetup path =.*$|cryptsetup path = ${lib.makeBinPath [cryptsetup]}/cryptsetup|' $out/etc/${projectName}/${projectName}.conf
    sed -i 's|^# go path =.*$|go path = ${lib.makeBinPath [go]}/go|' $out/etc/${projectName}/${projectName}.conf
    ${lib.optionalString enableNvidiaContainerCli ''
      sed -i 's|^# nvidia-container-cli path =.*$|nvidia-container-cli path = ${lib.makeBinPath [nvidia-docker]}/nvidia-container-cli|' $out/etc/${projectName}/${projectName}.conf
    ''}
    ${lib.optionalString (removeCompat && (projectName != "singularity")) ''
      while IFS= read -r -d $'\0' filePath; do
        if [ -f "$(realpath "$filePath")" ]; then
          fileDir="$(dirname "$filePath")"
          fileName="$(basename "$filePath")"
          fileNameReplaced="''${fileName/singularity/${projectName}}"
          if [ -e "$fileDir/$fileNameReplaced" ] && [ "$(realpath "$filePath")" == "$(realpath "$fileDir/$fileNameReplaced")" ]; then
            unlink "$filePath"
          fi
        fi
      done < <(find "$out" -mindepth 1 -type l -name "*singularity*" -print0)
      for filePath in "$out"/share/*-completion/completions/singularity; do
        if [ -e "$filePath" ]; then
          rm "$filePath"
        fi
      done
    ''}
    ${lib.optionalString enableSuid ''
      chmod +x $out/libexec/${projectName}/bin/starter-suid
    ''}
    ${lib.optionalString (enableSuid && starterSuidPath != null) ''
      mv "$out"/libexec/${projectName}/bin/starter-suid{,.orig}
      ln -s ${lib.escapeShellArg starterSuidPath} "$out/libexec/${projectName}/bin/starter-suid"
    ''}
  '';

  meta = with lib; {
    description = "Application containers for linux" + extraDescription;
    longDescription = ''
      Singularity (the upstream) renamed themselves to Apptainer
      to distinguish themselves from a fork made by Sylabs Inc.. See

      https://sylabs.io/2021/05/singularity-community-edition
      https://apptainer.org/news/community-announcement-20211130
    '';
    license = licenses.bsd3;
    platforms = platforms.linux;
    maintainers = with maintainers; [ jbedo ShamrockLee ];
    mainProgram = projectName;
  } // extraMeta;
}
