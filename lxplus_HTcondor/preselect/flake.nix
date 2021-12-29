{
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-21.05";
  inputs.flake-utils.url = "github:numtide/flake-utils";
  inputs.root-source.url = "github:root-project/root/master";
  inputs.root-source.flake = false;
  outputs = inputs@{ self, nixpkgs, flake-utils, root-source, ... }: flake-utils.lib.eachDefaultSystem (system:
    let
      lib = nixpkgs.lib.${system};
      switchFlag = patternToRemove: flagToAdd: flags: (
        let
          flags' = lib.filter (flag: builtins.match patternToRemove flag == null) flags;
        in
        if lib.elem flagToAdd flags' then flags'
        else flags' ++ [ flagToAdd ]
      );
      pkgs = import nixpkgs {
        inherit system;
        overlays = [ (final: prev:
          {
            root = (final.callPackage dependency-builders/root (with final; {
              python = python3;
              inherit (darwin.apple_sdk.frameworks) Cocoa CoreSymbolication OpenGL;
            })).overrideAttrs (oldAttrs: {
              src = root-source;
              cmakeFlags = builtins.foldl'
                (flags: pair:
                  switchFlag (lib.first pair) (lib.last pair) flags
                ) [
                  [ "-Dimt=OFF" "-Dimt=ON" ]
                  [ "-Dssl=OFF" "-Dssl=ON" ]
                  # # The dependencies are not packaged yet.
                  # [ "-Dgfal=OFF" "-Dgfal=ON" ]
                  # [ "-Dxrootd=OFF" "-Dxrootd=ON" ]
                  [ "-DCMAKE_BUILD_TYPE=.*" "-DCMAKE_BUILD_TYPE=RelWithDebInfo" ]
              ];
              buildInputs = oldAttrs.buildInputs ++ (with final; [
                tbb # for implicit multithreading
                openssl # for ssl support
              ]);
            });
          }) ];
      };
      inherit (pkgs) root;
      devShell = pkgs.mkShell {
        buildInputs = (with pkgs; [
          root
          gcc
        ]);
        nativeBuildInputs = (with pkgs;[
          root
          gcc
          gnumake
          cmake
          gdb
          gmock
          gtest
          gawk
          gitAndTools.gitFull
        ]);
      };
      packagesSub = {
        inherit (pkgs) root gcc gnumake cmake gdb gmock gtest;
        inherit (pkgs) gawk;
        inherit (pkgs.gitAndTools) git gitFull;
      };
      run = pkgs.writeShellScriptBin "run" ''
        export PATH="${ with packagesSub; pkgs.lib.makeBinPath [ root gcc gnumake cmake gawk gitFull ]}:$PATH"
        export LD_LIBRARY_PATH="${ pkgs.lib.makeLibraryPath (with pkgs; [ gcc ]) ++ [ root ] }:$LD_LIBRARY_PATH"
        if test -n "${devShell.shellHook}"; then
          . "${devShell.shellHook}";
        fi
        exec "$@"
      '';
      ana = pkgs.callPackage ./ana.nix { inherit (packagesSub) root; };
    in
    {
      legacyPackages = pkgs;
      inherit devShell;
      defaultPackage = run;
      packages = packagesSub // {
        srcRaw = self;
        inherit run ana;
      };
    });
}
