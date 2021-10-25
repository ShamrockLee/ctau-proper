{ lib
, root
, nameMacro ? "xAna_monoZ_preselect.C"
, mainProgram ? "xAna_monoZ_preselect"
, srcRaw ? ./.
}:

root.stdenv.mkDerivation {
  pname = mainProgram;
  version = "0.0.1";

  src = lib.cleanSourceWith {
    filter = (name: type:
      (lib.cleanSourceFilter name type)
      && ! (lib.hasSuffix ".root" name)
    );
    src = srcRaw;
  };

  nativeBuildInputs = [
    root
  ];

  buildInputs = [
  ] ++ root.buildInputs;

  dontUnpack = true;

  buildPhase = ''
    runHook preBuild
    cc -o "${mainProgram}" $(root-config --glibs --cflags) -g -Wall -Wextra -O1 "$src/${nameMacro}"
    runHook postBuild
  '';

  installPhase = ''
    runHook preInstall
    mkdir -p "$out/bin"
    mv "${mainProgram}" "$out/bin/"
    runHook postInstall
  '';

  meta = with lib; {
    description = "The package to build my analysis macro";
    inherit mainProgram;
  };
}
