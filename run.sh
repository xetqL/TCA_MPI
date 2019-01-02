#!/bin/bash
make
mpirun -np $1 bin/CA_road $2 out.jpg $3
ffmpeg -framerate 16 -pattern_type glob -i '*.jpg'   -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4
rm *.jpg

