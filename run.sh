#!/bin/bash
make
mpirun --oversubscribe -np $1 bin/CA_road out.jpg $2 $3 $4
python create_images.py $3 $4
ffmpeg -framerate 16 -pattern_type glob -i '*py_out.jpg'   -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4
rm -I *.jpg
rm -I *_waiting_time.txt