#!/usr/bin/env bash
# Reproduces examples/anim_nyc_2013.mp4 (+ _poster.png, _frames.png): the July 2013 New York event (G34542) of the catalogue.
# DATA: folder of the unpacked Zenodo dataset (doi:10.5281/zenodo.23966834). ERA5 is read from the public ARCO-ERA5 store.
set -euo pipefail
DATA=${1:-na025_heatwave_objects_v1}
cd "$(dirname "$0")"
python ../scripts/anim_event.py --data "$DATA" \
    --point -74.01 40.71 --date 2013-07-18 --name NYC \
    --t2m arco --start 2013-07-10 --end 2013-07-31 --box -97 -60 30 54 \
    --out anim_nyc_2013.mp4
