{
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-21.11";
  inputs.flake-utils.url = "github:numtide/flake-utils";
  inputs.root-source.url = "github:root-project/root/master";
  inputs.root-source.flake = false;
  inputs.nix-portable-flake.url = "github:DavHau/nix-portable";
  inputs.nix-for-nix-portable.url = "nix/2.5.1";
  inputs.nix-for-nix-portable.inputs.flake-utils.follows = "flake-utils";
  inputs.nix-for-nix-portable.inputs.nixpkgs.follows = "nixpkgs";
  # We don't need to benchmark Nix, but we try to keep the flake.lock unchanged
  inputs.nix-for-nix-portable.inputs.nixpkgs-regression.follows = "nixpkgs";
  # Use the same nixpkgs as other packages
  inputs.nix-portable-flake.inputs.nixpkgs.follows = "nixpkgs";
  inputs.nix-portable-flake.inputs.flake-utils.follows = "flake-utils";
  inputs.nix-portable-flake.inputs.nix.follows = "nix-for-nix-portable";
  # Apptainer is the new name chosen by the Singularity community
  # A symlink to the name singularity at $out/bin is preserved
  inputs.apptainer-source.url = "github:ShamrockLee/apptainer/noroot";
  # This input is to follow the changes at
  # https://github.com/apptainer/apptainer/pull/250
  # inputs.apptainer-source.url = "github:ShamrockLee/apptainer/noroot";
  inputs.apptainer-source.flake = false;
  # inputs.nix-prefetch-flake.url = "github:msteen/nix-prefetch";
  inputs.nix-prefetch-flake.url = "github:ShamrockLee/nix-prefetch/experimental-features";
  inputs.nix-prefetch-flake.inputs.nixpkgs.follows = "nixpkgs";
  inputs.nix-prefetch-flake.inputs.flake-utils.follows = "flake-utils";

  outputs = inputs@{ self
  , nixpkgs
  , flake-utils
  , root-source
  , nix-portable-flake
  , apptainer-source
  , nix-prefetch-flake
  , ... }: flake-utils.lib.eachDefaultSystem (system:
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
                version = root-source.lastModifiedDate;
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
            }
          )
          (final: prev:
            {
              inherit (
                final.callPackage ./dependency-builders/singularity/packages.nix { }
              ) singularity-legacy apptainer singularity-ce;
              singularity = final.apptainer;
              singularity-tools = final.callPackage ./dependency-builders/singularity-tools { };
            }
          )
          (final: prev:
            {
              apptainer = (prev.apptainer.override {
                pname = "apptainer-unstable";
                version = apptainer-source.lastModifiedDate;
                src = apptainer-source;
                vendorSha256 = "sha256-u6iJScCAtLfHunUQyWYz1+xtDMZ51F7i+jnWcfn0KTw=";
              }).overrideAttrs (oldAttrs: {
                postPatch = (oldAttrs.postPatch or "") + ''
                  echo "${apptainer-source.lastModifiedDate}" > VERSION
                '';
              });
            }
          )
        ];
      };
      nix-prefetch = nix-prefetch-flake.packages.${system}.nix-prefetch.override { nix = pkgs.nixFlakes; };
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
          # apptainer
          nix-prefetch
        ]);
      };
      packagesSub = {
        inherit (pkgs) root gcc gnumake cmake gdb gmock gtest;
        inherit (pkgs) gawk;
        inherit (pkgs.gitAndTools) git gitFull;
        # inherit (pkgs) apptainer;
        inherit nix-prefetch;
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
      ana = pkgs.callPackage ./ana.nix (rec {
        subPathMacro = "xAna_monoZ_preselect.C";
        subPaths = [ subPathMacro "skeeto_optparse.h" ];
      });
      mkSingularityBuildscript = package: pkgs.callPackage ./make-singularity-buildscript.nix {
        inherit package;
        singularity = pkgs.apptainer;
      };
      ana-singularity-image = pkgs.singularity-tools.buildImage' {
        name = ana.pname;
        contents = [ ana pkgs.parallel ];
        executableFlags = [ "--verbose" ];
        diskSize = 4096;
        memSize = 2048;
      };
      ana-singularity-buildscript = (pkgs.singularity-tools.buildImageFromDef {
        name = ana.pname;
        contents = [ ana pkgs.parallel ];
        executableFlags = [ "--verbose" ];
        buildImageFlags = [ "--unprivileged" ];
      }).buildscriptPackage;
      apptainer-prefetch-vendorsha256 = (pkgs.writeShellScriptBin "apptainer-prefetch-vendorsha256" ''
        NIX_PATH="nixpkgs=${nixpkgs}" "${nix-prefetch}/bin/nix-prefetch" \
          --extra-experimental-features "nix-command flakes" \
          "let preselect-flake = builtins.getFlake \"git+file://$( (cd ../..; pwd) )?dir=lxplus_HTcondor/preselect\"; preselect-pkgs = preselect-flake.legacyPackages.${system}; in { sha256 }: (preselect-pkgs.apptainer.override { vendorSha256 = sha256; }).go-modules"
      '').overrideAttrs (oldAttrs: {
        meta = (oldAttrs.meta or { }) // {
          mainProgram = "apptainer-prefetch-vendorsha256";
        };
      });
    in
    {
      legacyPackages = pkgs;
      inherit devShell;
      defaultPackage = run;
      packages = packagesSub // {
        srcRaw = self;
        inherit run ana compile-with-root ana-singularity-buildscript ana-singularity-image apptainer-prefetch-vendorsha256;
      } // lib.optionalAttrs (system == "x86_64-linux") (lib.mapAttrs (
        # Use the same pkgs as other packages
        name: value: value.override { inherit pkgs; }
      ) {
        inherit (nix-portable-flake.packages.${system}) nix-portable nix-portable-aarch64-linux;
      });
    });
}
