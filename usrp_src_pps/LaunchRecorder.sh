#!/bin/bash
gain_val=$1

for ((i=0; i<$(nproc); i++)); do sudo cpufreq-set -c $i -r -g performance; done

echo "Starting USRP"

trc_src="internal"
out_dir="/home/starnav"
OUTFILE="$out_dir/$trc_src.raw"

config_dir="$PWD/config.ini"

#change the absolute path for the creation of the executable
make
sudo chmod +x logger
./logger --config_file="$config_dir" --file="$OUTFILE"  --gain=$gain_val --args="recv_frame_size=16360,num_recv_frames=128"

# make clean

echo "Finished recording"

# Set all SPU cores to power saver
for ((i=0; i<$(nproc); i++)); do sudo cpufreq-set -c $i -r -g powersave; done
