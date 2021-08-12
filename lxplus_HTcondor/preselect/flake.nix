{
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-21.05";
  inputs.nixpkgs-root.url = "github:ShamrockLee/nixpkgs/root-6-25";
  inputs.flake-utils.url = "github:numtide/flake-utils";
  outputs = inputs@{nixpkgs, nixpkgs-root, flake-utils, ...}: flake-utils.lib.eachDefaultSystem (system: let
    pkgs = nixpkgs.legacyPackages.${system};
    pkgs-root = nixpkgs-root.legacyPackages.${system};
  in{
    legacyPackages = pkgs;
    legacyPackages-root = pkgs-root;
    packages = {
      inherit (pkgs-root) root gcc gnumake cmake;
      inherit (pkgs) gawk;
    };
    defaultPackage = pkgs-root.root;
    devShell = pkgs.mkShell {
      buildInputs = (with pkgs-root; [
        root
      ]);
      nativeBuildInputs = (with pkgs-root; [
        gcc
        gnumake
        cmake
      ]);
    };
  });
}
