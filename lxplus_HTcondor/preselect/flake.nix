{
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-21.05";
  inputs.nixpkgs-root.url = "github:ShamrockLee/nixpkgs/root-6-25";
  inputs.flake-utils.url = "github:numtide/flake-utils";
  inputs.root-source.url = "github:root-project/root/master";
  inputs.root-source.flake = false;
  outputs = inputs@{self, nixpkgs, nixpkgs-root, flake-utils, root-source, ...}: flake-utils.lib.eachDefaultSystem (system: let
    pkgs = nixpkgs.legacyPackages.${system};
    pkgs-root = import nixpkgs-root {
      inherit system;
      overlays = [
        (final: prev: {
          root = prev.root.overrideAttrs (oldAttrs: {
            version = "2021-09-01";
            src = root-source;
            cmakeFlags = (map (oldFlag:
              if oldFlag == "-Dimt=OFF" then "-Dimt=ON"
              else if oldFlag == "-Dssl=OFF" then "-Dssl=ON"
              # else if oldFlag == "-Dgfal=OFF" then "-Dgfal=ON"
              # else if oldFlag == "-Dxrootd=OFF" then "-Dxrootd=ON"
              else oldFlag
            ) oldAttrs.cmakeFlags) ++ [
              "-DCMAKE_BUILD_TYPE=RelWithDebInfo"
            ];
            buildInputs = oldAttrs.buildInputs ++ (with pkgs-root; [
              tbb # for implicit multithreading
              openssl # for ssl support
            ]);
          });
        })
      ];
    };
    devShell = pkgs.mkShell {
      buildInputs = (with pkgs-root; [
        root
      ]);
      nativeBuildInputs = (with pkgs-root; [
        gcc
        gnumake
        cmake
        gdb
        gmock
        gtest
      ]) ++ (with pkgs;[
        gawk
        gitAndTools.gitFull
      ]);
    };
    packagesSub = {
      inherit (pkgs-root) root gcc gnumake cmake gdb gmock gtest;
      inherit (pkgs) gawk;
      inherit (pkgs.gitAndTools) git gitFull;
    };
    run = pkgs.writeShellScriptBin "run" ''
      export PATH="${ with packagesSub; pkgs.lib.makeBinPath [ root gcc gnumake cmake gawk gitFull ]}:$PATH"
      if test -n "${devShell.shellHook}"; then
        . "${devShell.shellHook}";
      fi
      exec "$@"
    '';
    ana = pkgs.callPackage ./ana.nix { inherit (packagesSub) root; };
  in{
    legacyPackages = pkgs;
    legacyPackages-root = pkgs-root;
    inherit devShell;
    defaultPackage = run;
    packages = packagesSub // {
      srcRaw = self;
      inherit run ana;
    };
  });
}
