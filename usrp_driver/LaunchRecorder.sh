#!/bin/bash

for ((i=0; i<$(nproc); i++)); do sudo cpufreq-set -c $i -r -g performance; done

echo "Starting USRP"

sudo sysctl -w net.core.wmem_max=2500000

trc_src="internal"
out_dir="/home/kevin"
OUTFILE="$out_dir/$trc_src.raw"

config_dir="$PWD/config.ini"

#change the absolute path for the creation of the executable
make
sudo chmod +x logger
./logger --config_file="$config_dir" --file="$OUTFILE"

# make clean

echo "Finished recording"

# Set all SPU cores to power saver
for ((i=0; i<$(nproc); i++)); do sudo cpufreq-set -c $i -r -g powersave; done
