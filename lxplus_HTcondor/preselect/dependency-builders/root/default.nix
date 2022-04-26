{ stdenv
, lib
, fetchurl
, fetchpatch
, makeWrapper
, cmake
, git
, ftgl
, gl2ps
, glew
, gnugrep
, gnused
, gsl
, lapack
, libX11
, libXpm
, libXft
, libXext
, libGLU
, libGL
, libxml2
, llvm_9
, lz4
, xz
, man
, openblas
, pcre
, nlohmann_json
, pkg-config
, python
, xxHash
, zlib
, zstd
, libAfterImage
, giflib
, libjpeg
, libtiff
, libpng
, procps
, tbb
, which
, Cocoa
, CoreSymbolication
, OpenGL
, version ? "6.24.06"
, noSplash ? false
, sw-vers-patch ? ./sw_vers.patch
}:

stdenv.mkDerivation rec {
  pname = "root";
  inherit version;

  src = fetchurl {
    url = "https://root.cern.ch/download/root_v${version}.source.tar.gz";
    sha256 = "sha256-kH9p9LrKHk8w7rSXlZjKdZm2qoA8oEboDiW2u6oO9SI=";
  };

  nativeBuildInputs = [ makeWrapper cmake pkg-config git ];
  buildInputs = [
    ftgl
    gl2ps
    glew
    pcre
    zlib
    zstd
    lapack
    libxml2
    llvm_9
    lz4
    xz
    gsl
    openblas
    xxHash
    libAfterImage
    giflib
    libjpeg
    libtiff
    libpng
    nlohmann_json
    python.pkgs.numpy
    tbb
  ]
  ++ lib.optionals (!stdenv.isDarwin) [ libX11 libXpm libXft libXext libGLU libGL ]
  ++ lib.optionals (stdenv.isDarwin) [ Cocoa CoreSymbolication OpenGL ]
  ;

  patches =
  lib.optional (sw-vers-patch != null) sw-vers-patch
  # Fix builtin_llvm=OFF support
  ++ lib.optional (lib.versionOlder version "6.25.2") (fetchpatch {
    url = "https://github.com/root-project/root/commit/0cddef5d3562a89fe254e0036bb7d5ca8a5d34d2.diff";
    excludes = [ "interpreter/cling/tools/plugins/clad/CMakeLists.txt" ];
    sha256 = "sha256-VxWUbxRHB3O6tERFQdbGI7ypDAZD3sjSi+PYfu1OAbM=";
  })
  ;

  # Fix build against vanilla LLVM 9
  postPatch = ''
    sed \
      -e '/#include "llvm.*RTDyldObjectLinkingLayer.h"/i#define private protected' \
      -e '/#include "llvm.*RTDyldObjectLinkingLayer.h"/a#undef private' \
      -i interpreter/cling/lib/Interpreter/IncrementalJIT.h
  '';

  preConfigure = ''
    rm -rf builtins/*
    substituteInPlace cmake/modules/SearchInstalledSoftware.cmake \
      --replace 'set(lcgpackages ' '#set(lcgpackages '

    # Don't require textutil on macOS
    : > cmake/modules/RootCPack.cmake

    # Hardcode path to fix use with cmake
    sed -i cmake/scripts/ROOTConfig.cmake.in \
      -e 'iset(nlohmann_json_DIR "${nlohmann_json}/lib/cmake/nlohmann_json/")'

    patchShebangs build/unix/
  '' + lib.optionalString noSplash ''
    substituteInPlace rootx/src/rootx.cxx --replace "gNoLogo = false" "gNoLogo = true"
  '' + lib.optionalString stdenv.isDarwin ''
    # Eliminate impure reference to /System/Library/PrivateFrameworks
    substituteInPlace core/CMakeLists.txt \
      --replace "-F/System/Library/PrivateFrameworks" ""
  '';

  cmakeFlags = [
    "-Drpath=ON"
    "-DCMAKE_INSTALL_BINDIR=bin"
    "-DCMAKE_INSTALL_LIBDIR=lib"
    "-DCMAKE_INSTALL_INCLUDEDIR=include"
    "-Dbuiltin_llvm=OFF"
    "-Dbuiltin_nlohmannjson=OFF"
    "-Dbuiltin_openui5=OFF"
    "-Dalien=OFF"
    "-Dbonjour=OFF"
    "-Dcastor=OFF"
    "-Dchirp=OFF"
    "-Dclad=OFF"
    "-Ddavix=OFF"
    "-Ddcache=OFF"
    "-Dfail-on-missing=ON"
    "-Dfftw3=OFF"
    "-Dfitsio=OFF"
    "-Dfortran=OFF"
    "-Dimt=ON"
    "-Dgfal=OFF"
    "-Dgviz=OFF"
    "-Dhdfs=OFF"
    "-Dhttp=ON"
    "-Dkrb5=OFF"
    "-Dldap=OFF"
    "-Dmonalisa=OFF"
    "-Dmysql=OFF"
    "-Dodbc=OFF"
    "-Dopengl=ON"
    "-Doracle=OFF"
    "-Dpgsql=OFF"
    "-Dpythia6=OFF"
    "-Dpythia8=OFF"
    "-Drfio=OFF"
    "-Droot7=OFF"
    "-Dsqlite=OFF"
    "-Dssl=OFF"
    "-Dtmva=ON"
    "-Dvdt=OFF"
    "-Dwebgui=OFF"
    "-Dxml=ON"
    "-Dxrootd=OFF"
  ]
  ++ lib.optional (stdenv.cc.libc != null) "-DC_INCLUDE_DIRS=${lib.getDev stdenv.cc.libc}/include"
  ++ lib.optionals stdenv.isDarwin [
    "-DOPENGL_INCLUDE_DIR=${OpenGL}/Library/Frameworks"
    "-DCMAKE_DISABLE_FIND_PACKAGE_Python2=TRUE"

    # fatal error: module map file '/nix/store/<hash>-Libsystem-osx-10.12.6/include/module.modulemap' not found
    # fatal error: could not build module '_Builtin_intrinsics'
    "-Druntime_cxxmodules=OFF"
  ];

  postInstall = ''
    for prog in rootbrowse rootcp rooteventselector rootls rootmkdir rootmv rootprint rootrm rootslimtree; do
      wrapProgram "$out/bin/$prog" \
        --set PYTHONPATH "$out/lib" \
        --set ${lib.optionalString stdenv.isDarwin "DY"}LD_LIBRARY_PATH "$out/lib"
    done
  '';

  postFixup = ''
    # Patch thisroot.{sh,csh,fish} setxrd.sh
    patchCommandInplace(){
      local file="$1"
      shift
      local sedCommand=""
      local arg
      while [ "$#" -gt 0 ]; do
        local exeName="$1"
        local exePath="$2"
        shift 2
        sedCommand="$sedCommand s|^$exeName\$|$exePath|g;"
        sedCommand="$sedCommand s|^$exeName\\([^A-Za-z_-]\\)|$exePath\\1|g;"
        sedCommand="$sedCommand s|\\([^A-Za-z_-]\\)$exeName\$|\\1$exePath|g;"
        sedCommand="$sedCommand s|\\([^A-Za-z_-]\\)$exeName\\([^A-Za-z_-]\\)|\\1$exePath\\2|g;"
      done
      sed -i "''${sedCommand:1}" "$file"
    }
    patchWhichInplace(){
      for file in "$@"; do
        patchCommandInplace "$file" \
          "which -p" which-p \
          which "${lib.makeBinPath which}/which" \
          which-p "which -p"
      done
    }
    patchCommandInplace "$out/bin/thisroot.sh" \
      "grep" "${lib.makeBinPath gnugrep}/grep" \
      "manpath" "${lib.makeBinPath man}/manpath" \
      "man" "${lib.makeBinPath man}/man" \
      "ps" "${lib.makeBinPath procps}/ps" \
      "sed" "${lib.makeBinPath gnused}/sed"
    patchCommandInplace "$out/bin/thisroot.csh" \
      "grep" "${lib.makeBinPath gnugrep}/grep" \
      "manpath" "${lib.makeBinPath man}/manpath" \
      "man" "${lib.makeBinPath man}/man" \
      "sed" "${lib.makeBinPath gnused}/sed"
    patchCommandInplace "$out/bin/thisroot.fish" \
      "manpath" "${lib.makeBinPath man}/manpath" \
      "man" "${lib.makeBinPath man}/man"
    patchCommandInplace "$out/bin/setxrd.sh" \
      "manpath" "${lib.makeBinPath man}/manpath" \
      "man" "${lib.makeBinPath man}/man" \
      "sed" "${lib.makeBinPath gnused}/sed"
    patchWhichInplace "$out"/bin/{thisroot.sh,thisroot.csh,thisroot.fish,setxrd.sh}

    # Make ldd and sed available to the ROOT executable
    wrapProgram "$out/bin/root" --prefix PATH : "${lib.makeBinPath [ stdenv.cc.libc gnused ]}"
  '';

  setupHook = ./setup-hook.sh;

  meta = with lib; {
    homepage = "https://root.cern.ch/";
    description = "A data analysis framework";
    platforms = platforms.unix;
    maintainers = [ maintainers.veprbl ];
    license = licenses.lgpl21;
    mainProgram = "root";
  };
}
