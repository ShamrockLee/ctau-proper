{ pkgs ? import <nixpkgs> {}
, mkShell ? pkgs.mkShell
, gdb ? pkgs.gdb
}:
mkShell {
  buildInputs = [ gdb ];
}
