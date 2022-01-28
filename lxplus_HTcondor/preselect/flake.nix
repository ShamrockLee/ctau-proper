{
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-21.11";
  inputs.flake-utils.url = "github:numtide/flake-utils";
  inputs.root-source.url = "github:root-project/root/master";
  inputs.root-source.flake = false;
  inputs.nix-portable-flake.url = "github:DavHau/nix-portable";
  # Use the same nixpkgs as other packages
  inputs.nix-portable-flake.inputs.nixpkgs.follows = "nixpkgs";
  inputs.nix-portable-flake.inputs.flake-utils.follows = "flake-utils";

  outputs = inputs@{ self, nixpkgs, flake-utils, root-source, nix-portable-flake, ... }: flake-utils.lib.eachDefaultSystem (system:
    let
      lib = nixpkgs.lib;
      switchFlag = patternToRemove: flagToAdd: flags: (
        let
          flags' = lib.filter (flag: builtins.match patternToRemove flag == null) flags;
        in
        if lib.elem flagToAdd flags' then flags'
        else flags' ++ [ flagToAdd ]
      );
      pkgs = import nixpkgs {
        inherit system;
        overlays = [
          (final: prev:
            {
              root = (final.callPackage ./dependency-builders/root (with final; {
                python = python3;
                inherit (darwin.apple_sdk.frameworks) Cocoa CoreSymbolication OpenGL;
              })).overrideAttrs (oldAttrs: {
                src = root-source;
                cmakeFlags = builtins.foldl'
                  (flags: pair:
                    switchFlag (lib.head pair) (lib.last pair) flags
                  )
                  oldAttrs.cmakeFlags
                  [
                    [ "-Dimt=ON" "-Dimt=OFF" ]
                    [ "-Dssl=OFF" "-Dssl=ON" ]
                    # # The dependencies are not packaged yet.
                    # [ "-Dgfal=OFF" "-Dgfal=ON" ]
                    # [ "-Dxrootd=OFF" "-Dxrootd=ON" ]
                    [ "-DCMAKE_BUILD_TYPE=.*" "-DCMAKE_BUILD_TYPE=RelWithDebInfo" ]
                  ];
                buildInputs = (lib.subtractLists (with final; [
                  tbb # for static linking
                ]) (oldAttrs.buildInputs or [ ])) ++ (with final; [
                  # tbb # for implicit multithreading
                  openssl # for ssl support
                ]);
              });
            })
        ];
      };
      inherit (pkgs) root;
      devShell = pkgs.mkShell {
        buildInputs = (with pkgs; [
          root
          gcc
        ]);
        nativeBuildInputs = (with pkgs; [
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
      # `nix-shell -c` replacement with LD_LIBRARY_PATH and ./thisroot.sh soucing
      # Use as `nix run .#run -- something to run` or `nix run .# -- something to run`
      run = pkgs.writeShellScriptBin "run" ''
        export PATH="${ lib.makeBinPath (lib.unique (devShell.nativeBuildInputs ++ devShell.buildInputs ++ [ compile-with-root ])) }:$PATH"
        export LD_LIBRARY_PATH="${ pkgs.lib.makeLibraryPath devShell.nativeBuildInputs }:$LD_LIBRARY_PATH"
        if test -n "${devShell.shellHook}"; then
          . "${devShell.shellHook}";
        fi
        . "${packagesSub.root}/bin/thisroot.sh"
        exec "$@"
      '';
      # `nix run .#run -- g++ $(root-config ROOTCONFIG_ARGS) CC_ARGS` equivalence
      # to make life easier in VSCode's launch.json
      # Use as `nix run .# -- nix run .#compile-with-root -- ROOTCONFIG_ARGS -- CC_ARGS`
      # `''${}` is the escaped form of `${}` in a ''...'' string
      compile-with-root = pkgs.writeShellScriptBin "compile-with-root" ''
        declare -a ROOTCONFIG_ARGS=()
        declare -a CC_ARGS=()
        isRootConfigArg=true
        for arg in "$@"; do
          if $isRootConfigArg && [ "$arg" == "--" ]; then
            isRootConfigArg=false
          else
            if $isRootConfigArg; then
              ROOTCONFIG_ARGS+=("$arg")
            else
              CC_ARGS+=("$arg")
            fi
          fi
        done
        echo g++ $(root-config "''${ROOTCONFIG_ARGS[@]}") "''${CC_ARGS[@]}" >2
        g++ $(root-config "''${ROOTCONFIG_ARGS[@]}") "''${CC_ARGS[@]}"
      '';
      ana = pkgs.callPackage ./ana.nix { inherit (packagesSub) root; };
    in
    {
      legacyPackages = pkgs;
      inherit devShell;
      defaultPackage = run;
      packages = packagesSub // {
        srcRaw = self;
        inherit run ana compile-with-root;
      } // lib.optionalAttrs (system == "x86_64-linux") (lib.mapAttrs (
        # Use the same pkgs as other packages
        name: value: value.override { inherit pkgs; }
      ) {
        inherit (nix-portable-flake.packages.${system}) nix-portable nix-portable-aarch64-linux;
      });
    });
}
