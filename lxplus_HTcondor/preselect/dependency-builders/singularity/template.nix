{ lib
, util-linux
, gpgme
, openssl
, libuuid
, bash
, coreutils
, which
, makeWrapper
, cryptsetup
, squashfsTools
, nvidia-docker # For nvidia-container-cli
, glibc # For ldconfig needed by nvidia-container-cli
, go
, buildGoModule
, pname
, version
, src
, projectName # "apptainer" or "singularity"
  # These cannot be correctly override with `overrideAttrs`
, vendorSha256 ? null
, deleteVendor ? false
, proxyVendor ? false
, doCheck ? true
, enableNvidiaContainerCli ? true
}:

with lib;
let
  propagatedBuildInputs = [ coreutils squashfsTools bash ];
in
buildGoModule {
  inherit pname version src;

  # `buildGoModule` reads the `go.mod` file,
  # so no need to specify `goPackagePath`

  # All three projects provide vendored source tarball
  # so we set vendorSha256 = null;
  # One who use unvendored source can override this
  # with `override`
  inherit vendorSha256 deleteVendor proxyVendor;

  # apptainer and singularity-ce both want to reference Go
  # to compile the extensions when building container images
  allowGoReference = true;

  passthru = {
    inherit projectName;
  };

  buildInputs = [ gpgme openssl libuuid ];
  nativeBuildInputs = [ util-linux which makeWrapper cryptsetup ];
  inherit propagatedBuildInputs;

  postPatch = ''
    substituteInPlace internal/pkg/build/files/copy.go \
      --replace /bin/cp ${coreutils}/bin/cp
  '';

  # Using `buildGoModule` instead of `buildGoPackage`,
  # We no longer need to `cd go/src/${goPackagePath}`
  postConfigure = ''
    patchShebangs .
    for subPath in cmd/internal/cli/actions.go internal/pkg/util/env/env.go e2e/env/env.go; do
      if [ -f "$subPath" ]; then
        echo "Patching [Dd]efaultPath in \"$subPath\""
        sed -i 's|^\([ \t]*[Dd]efaultPath[ \t]*[:]\?=[ \t]*\)"[^"]*"|\1"${lib.makeBinPath propagatedBuildInputs}"|' "$subPath"
      fi
    done

    ./mconfig -V ${version} -p $out --localstatedir=/var

    # Don't install SUID binaries
    sed -i 's/-m 4755/-m 755/g' builddir/Makefile
  '';

  buildPhase = ''
    runHook preBuild
    make -C builddir
    runHook postBuild
  '';

  installPhase = ''
    runHook preInstall
    make -C builddir install LOCALSTATEDIR=$out/var
    chmod 755 $out/libexec/${projectName}/bin/starter-suid

    # Explicitly configure paths in the config file
    sed -i 's|^# mksquashfs path =.*$|mksquashfs path = ${lib.makeBinPath [squashfsTools]}/mksquashfs|' $out/etc/${projectName}/${projectName}.conf
    sed -i 's|^# unsquashfs path =.*$|unsquashfs path = ${lib.makeBinPath [squashfsTools]}/unsquashfs|' $out/etc/${projectName}/${projectName}.conf
    sed -i 's|^# cryptsetup path =.*$|cryptsetup path = ${lib.makeBinPath [cryptsetup]}/cryptsetup|' $out/etc/${projectName}/${projectName}.conf
    sed -i 's|^# go path =.*$|go path = ${lib.makeBinPath [go]}/go|' $out/etc/${projectName}/${projectName}.conf
    ${optionalString enableNvidiaContainerCli ''
      sed -i 's|^# nvidia-container-cli path =.*$|nvidia-container-cli path = ${lib.makeBinPath [nvidia-docker]}/nvidia-container-cli|' $out/etc/${projectName}/${projectName}.conf
      sed -i 's|^# ldconfig path =.*$|ldconfig path = ${lib.makeBinPath [glibc.bin]}/ldconfig|' $out/etc/${projectName}/${projectName}.conf
    ''}

    runHook postInstall
  '';

  meta = with lib; {
    description = "Application containers for linux";
    longDescription = ''
      Singularity (the upstream) renamed themselves to Apptainer
      to distinguish themselves from a fork made by Sylabs Inc.. See

      https://sylabs.io/2021/05/singularity-community-edition
      https://apptainer.org/news/community-announcement-20211130

      singularity-legacy is from the original repo before the renaming.

      apptainer is from the new repo after the renaming.

      singularity-ce is from the fork of Sylabs Inc..

      The GitHub repo hpcng/singularity now redirects to apptainer/apptainer
      and https://singularity.lbl.gov now redirects to https://apptainer.org/

      Historical versions of the homepage (prior to the redirection)
      are preserved on the Wayback Machine of Internet Archive.
    '';
    license = licenses.bsd3;
    platforms = platforms.linux;
    maintainers = [ maintainers.jbedo ];
    mainProgram = projectName;
  };
}
