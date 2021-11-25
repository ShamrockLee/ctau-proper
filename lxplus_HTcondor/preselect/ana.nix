# Nix expression to build the analysis ROOT macro

{ lib
, makeWrapper
, root
, nameMacro ? "xAna_monoZ_preselect.C"
, mainProgram ? "xAna_monoZ_preselect"
, srcRaw ? ./.
}:

let
  # The default and unwrapped gcc package
  ccUnwrapped = root.stdenv.cc.cc;
in
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
    makeWrapper
  ];

  buildInputs = [
    root
  ];

  dontUnpack = true;

  # Use g++/c++ instead of gcc/cc
  # to prevent the code being falsely identified as C instead of C++
  # See https://discourse.nixos.org/t/how-to-compile-cern-root-macro-analysis-c-code/15695/2
  buildPhase = ''
    runHook preBuild
    c++ -o "${mainProgram}.o" $(root-config --glibs --cflags) -g -O2 "$src/${nameMacro}"
    runHook postBuild
  '';

  installPhase = ''
    runHook preInstall
    mkdir -p "$out/bin"
    mv "${mainProgram}.o" "$out/bin/"
    makeWrapper "$out/bin/${mainProgram}.o" "$out/bin/${mainProgram}" --prefix LD_LIBRARY_PATH : ${ lib.makeLibraryPath [ root ccUnwrapped ] }
    runHook postInstall
  '';

  meta = with lib; {
    description = "The package to build my analysis macro";
    inherit mainProgram;
  };
}
