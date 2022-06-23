{ lib
, callPackage
, fetchFromGitHub
}:

{
  apptainer = callPackage ./generic.nix rec {
    pname = "apptainer";
    version = "1.0.2";
    projectName = "apptainer";

    src = fetchFromGitHub {
      owner = "apptainer";
      repo = "apptainer";
      rev = "v${version}";
      sha256 = "3yN2bDpgeOtCzFFhtPLf43wwzDNYxPCFnbgcL/1Mweg=";
    };

    # Update by running
    # nix-prefetch -E "{ sha256 }: ((import ./. { }).apptainer.override { vendorSha256 = sha256; }).go-modules"
    # at the root directory of the Nixpkgs repository
    vendorSha256 = "N0vbQToRCASKvNW/1bSLUXG2uwwg+0FO4lT3kG3ab04=";

    # This would soon be false
    # https://github.com/apptainer/apptainer/pull/495
    defaultToSuid = true;

    extraDescription = " (previously known as Singularity)";
    extraMeta.homepage = "https://apptainer.org";
  };

  singularity = callPackage ./generic.nix rec {
    pname = "singularity-ce";
    version = "3.9.8";
    projectName = "singularity";

    src = fetchFromGitHub {
      owner = "sylabs";
      repo = "singularity";
      rev = "v${version}";
      sha256 = "qjM+2JN0KJ8WrzNCCqeK6dyOoOnVMDin/BLzgjyWy50=";
    };

    vendorSha256 = "LkmfDytRr/IeRsI6B9ypSCGAvAqBjK9VFDmICYJ6eco=";

    defaultToSuid = true;

    extraDescription = " (Sylabs Inc's fork of Singularity)";
    extraMeta.homepage = "https://sylabs.io/";
  };
}
