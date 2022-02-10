{ lib
, writeShellScriptBin
, writeTextFile
, bash
, coreutils
, runtimeShell
, singularity
, package
, pname ? (lib.getName package) + "-singularity-buildscript"
, version ? (lib.getVersion package)
, linkBinSh ? true
, linkUsrBinEnv ? true
}:

let
  definitionString = ''
    Bootstrap: scratch

    %setup
        nix-store --export $( \
          nix-store -qR "${package}"; \
          ${lib.optionalString linkBinSh "nix-store -qR \"${runtimeShell}\";"} \
          ${lib.optionalString linkUsrBinEnv "nix-store -qR \"${coreutils}\";"} \
        ) \
        | nix-store --store "''${SINGULARITY_ROOTFS}" --import
  '' + lib.optionalString linkBinSh ''
        mkdir -p "''${SINGULARITY_ROOTFS}/bin/sh"
        "${coreutils}/bin/ln" -s "${runtimeShell}" "''${SINGULARITY_ROOTFS}/bin/sh"
  '' + lib.optionalString linkUsrBinEnv ''
        mkdir -p "''${SINGULARITY_ROOTFS}/usr/bin/env"
        "${coreutils}/bin/ln" -s "${coreutils}/bin/env" "''${SINGULARITY_ROOTFS}/usr/bin/env;"
  '' + ''

    %environment
        export PATH="${package}/bin:/bin:$PATH"
  '';

  definitionFile = writeTextFile {
    name = "${pname}-singularity.def";
    text = definitionString;
  };

  scriptString = ''
    if [ "$#" -lt 1 ]; then
      echo "Expect IMAGE_PATH" >&2
      exit 1
    fi
    ${lib.toUpper singularity.projectName}ENV_PATH="$PATH" "${singularity}/bin/${singularity.mainProgram or "singularity"}" --verbose build "$1" "${definitionFile}"
  '';

in (writeShellScriptBin "build-singularity-image" scriptString).overrideAttrs (oldAttrs: {
  inherit pname version;
  passthru = {
    inherit definitionString definitionFile;
  };
  meta = (oldAttrs.meta or { }) // {
    mainProgram = "build-singularity-image";
  };
})
