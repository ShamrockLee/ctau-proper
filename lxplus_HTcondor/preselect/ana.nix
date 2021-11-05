# Nix expression to build the analysis ROOT macro

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
  ];

  dontUnpack = true;

  # Use g++/c++ instead of gcc/cc
  # to prevent the code being falsely identified as C instead of C++
  # See https://discourse.nixos.org/t/how-to-compile-cern-root-macro-analysis-c-code/15695/2
  buildPhase = ''
    runHook preBuild
    c++ -o "${mainProgram}" $(root-config --glibs --cflags) -g -Wall -Wextra -O2 "$src/${nameMacro}" || true
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
