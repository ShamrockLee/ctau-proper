{ pkgs ? import <nixpkgs> {}
, mkShell ? pkgs.mkShell
, clang ? pkgs.clang
}:
mkShell {
  buildInputs = [ clang ];
}
