#!/bin/sh


PATH_PVPYTHON="/work/utox/users/helleboid_work_utox/software/ParaView-5.9.0-MPI-Linux-Python3.8-64bit/bin/pvpython --mesa"
PARAVIEW_SCRIPT_PATH="/work/utox/users/helleboid_work_utox/EmpiricalPseudopotential/python/paraview_animation_steady.py"
BSUB_COMMAND="bsub -R \"select[ (rh70) && cpuf >= 220] rusage[mem=20960 ]\" -n 1 -q reg"

NB_BAND=15
NB_FRAMES=450
OUT_DIR="steady_anim"


for Indexband in $(seq 0 $NB_BAND); do bsub -q reg $PATH_PVPYTHON $PARAVIEW_SCRIPT_PATH -f ggg -b ${Indexband} -F ${NB_FRAMES} -o ${OUT_DIR}; done




#avi to mp4
# find . -type f -name "*.avi" -exec ffmpeg -i {} -c:v copy -c:a copy {}.mp4 \;

# Place all the videos filename into a file
# find . -type f -name "*.avi.mp4" -exec echo "file" {} >>  FILEMP4.txt   \;

# Sort the files
# sort FILEMP4.txt -o SFILEMP4.txt

# Concatenate all the videos
# ffmpeg -f concat -i SFILEMP4.txt -c copy compilation_iso.mp4
