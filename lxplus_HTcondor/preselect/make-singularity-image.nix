{ lib
, stdenvNoCC
, writeTextFile
, singularity
, package
, pname ? (lib.getName package) + "-singularity-image"
, version ? (lib.getVersion package)
}:

let
  definitionString = ''
    Boostrap: scratch

    %setup
        nix-store --export $(nix-env -qR "${package}")
        | nix-store --store-path "''${SINGULARITY_ROOTFS}" --import

    %environment
        export PATH="${package}/bin:$PATH"
  '';
  definitionFile = writeTextFile {
    name = "${pname}-singularity.def";
    text = definitionString;
  };
in stdenvNoCC.mkDerivation {
  passthru = {
    inherit definitionString definitionFile;
  };
  inherit pname version;
  nativeBuildInputs = [
    singularity
    package
  ];
  buildInputs = [
    package
  ];
  dontUnpack = true;
  buildPhase = ''
    runHook preBuild
    SINGULARITYENV_PATH="$PATH" singularity build --no-check-root "$out" "${definitionFile}"
    runHook postBuild
  '';
  dontInstall = true;
  dontFixup = true;
}
