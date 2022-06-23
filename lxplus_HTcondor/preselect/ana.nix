# Nix expression to build the analysis ROOT macro

{ lib
, makeWrapper
, root
  # The analysis macro
, macroName
  # Name of the executable
  # Default to the base file name of the main analysis macro without the .C suffix
 , mainProgram ? (builtins.elemAt (builtins.match "(.*/)?([^.]+)[.].+" macroName) 1)
  # List of source files
, srcs
  # Don't try to fetch binary substitute for this drv
, allowSubstitutes ? false
, version ? "0.1.0"
}:

let
  # The default and unwrapped gcc package
  ccUnwrapped = root.stdenv.cc.cc;
in
root.stdenv.mkDerivation {
  pname = mainProgram;
  inherit version srcs;
  inherit allowSubstitutes;

  sourceRoot = "source";

  unpackPhase = ''
    runHook preUnpack
    mkdir -p "$sourceRoot"
    for _src in $srcs; do
      cp "$_src" "$sourceRoot/$(stripHash "$_src")"
    done
    chmod -R u+w -- "$sourceRoot"
    echo "source root is $sourceRoot"
    runHook postUnpack
  '';

  nativeBuildInputs = [
    root
    makeWrapper
  ];

  buildInputs = [
    root
  ];

  # Use g++/c++ instead of gcc/cc
  # to prevent the code being falsely identified as C instead of C++
  # See https://discourse.nixos.org/t/how-to-compile-cern-root-macro-analysis-c-code/15695/2
  buildPhase = ''
    runHook preBuild
    c++ -o "${mainProgram}.o" $(root-config --cflags) -g -O2 "${macroName}" $(root-config --glibs)
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
    description = "My analysis macro";
    inherit mainProgram;
  };
}
