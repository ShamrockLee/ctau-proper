# Nix expression to build the analysis ROOT macro

{ lib
, makeWrapper
, root
  # Sub-path to the analysis macro
, subPathMacro
  # Name of the executable
  # Default to the base file name of the main analysis macro without the .C suffix
, mainProgram ? (builtins.elemAt (builtins.match "(.*/)?([^.]+)[.].+" subPathMacro) 1)
  # The source base of the analysis program
, srcBase ? ./.
  # List of sub-paths to source files or directories relative to srcBase
  # If omitted or specified to `null`, all possible files will be included
, subPaths ? null
  # Don't try to fetch binary substitute for this drv
, allowSubstitutes ? false
, version ? "0.1.0"
}:

let
  # The default and unwrapped gcc package
  ccUnwrapped = root.stdenv.cc.cc;
in
root.stdenv.mkDerivation ((
 if builtins.isList subPaths then {
    unpackPhase = ''
      runHook preUnpack
      sourceRoot="''${sourceRoot:-source}"
      mkdir -p "$sourceRoot"
      for _subPath in ${toString (lib.escapeShellArgs subPaths)}; do
        cp -ar "${srcBase}/$_subPath" "$sourceRoot/$_subPath"
      done
      chmod -R +w "$sourceRoot"
      runHook postUnpack
    '';
  } else {
    src = lib.cleanSourceWith {
      filter = (name: type:
        (lib.cleanSourceFilter name type)
        && ! (lib.hasSuffix ".root" name)
        && ! (lib.hasSuffix ".sif" name)
        && ! (builtins.match "(.*/)?result.*" name != null)
      );
      src = srcBase;
    };
    unpackPhase = ''
      runHook preUnpack
      sourceRoot="''${sourceRoot:-source}"
      mkdir -p "$sourceRoot"
      cp -ar "$src/." "$sourceRoot/"
      chmod -R +w "$sourceRoot"
      runHook postUnpack
    '';
  }
) // {
  pname = mainProgram;
  inherit version;

  inherit allowSubstitutes;


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
    c++ -o "${mainProgram}.o" $(root-config --cflags) -g -O2 "${subPathMacro}" $(root-config --glibs)
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
})
