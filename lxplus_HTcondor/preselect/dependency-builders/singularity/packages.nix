{ lib
, callPackage
, fetchurl
}:

{
  singularity-legacy = (callPackage ./template.nix rec {
    pname = "singularity-legacy";
    version = "3.8.5";
    projectName = "singularity";

    src = fetchurl {
      url = "https://github.com/hpcng/singularity/releases/download/v${version}/singularity-${version}.tar.gz";
      sha256 = "sha256-f/94tcB7XU0IJpvSZ6xemUOQ+TMyHlTv1rfIZoMVPOQ=";
    };
  }).overrideAttrs (oldAttrs: {
    meta = oldAttrs.meta // (with lib; {
      description = oldAttrs.meta.description + " (leagacy version before renaming)";
      homepage = "https://singularity.lbl.gov";
    });
  });

  apptainer = (callPackage ./template.nix rec {
    pname = "apptainer";
    version = "1.0.0-rc.1";
    projectName = "apptainer";

    src = fetchurl {
      url = "https://github.com/apptainer/apptainer/releases/download/v${version}/apptainer-${version}.tar.gz";
      sha256 = "1clb9vlw4mv6vd3yxyfvngrzj469qmfpwxih4vzx2q89l2wf4313";
    };
  }).overrideAttrs (oldAttrs: {
    meta = oldAttrs.meta // {
      description = oldAttrs.meta.description + " (the new, renamed Singularity)";
      homepage = "https://apptainer.org";
    };
  });

  singularity-ce = (callPackage ./template.nix rec {
    pname = "singularity-ce";
    version = "3.9.5";
    projectName = "singularity";

    src = fetchurl {
      url = "https://github.com/sylabs/singularity/releases/download/v${version}/singularity-ce-${version}.tar.gz";
      sha256 = "1sb0d9jbj66jd6z4garvb6w8g7azwdhp5zl9jslsy2hxqbm10bvz";
    };
  }).overrideAttrs (oldAttrs: {
    meta = oldAttrs.meta // {
      description = oldAttrs.meta.description + " (Sylabs's fork)";
      homepage = "http://www.sylabs.io";
    };
  });
}
