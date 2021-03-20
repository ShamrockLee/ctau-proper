#!/usr/bin/env sh

nix-shell -p bindfs --run "bindfs /run/media/shamrock/data-ext4/ROOT_data/preselect/dataTest_DYJets $(dirname $0)/dataTest_DYJets"
nix-shell -p bindfs --run "bindfs /run/media/shamrock/data-ext4/ROOT_data/preselect/dataTest_TT $(dirname $0)/dataTest_TT"
