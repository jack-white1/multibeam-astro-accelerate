#!/bin/sh

# USAGE: ./run_AA.sh

# AA Script USAGE: $AA_dir/scripts/astro-accelerate.sh [ASTRO_ACCELERATE_REPOSITORY_PATH] [OUTPUT_DIR] [ASTRO_ACCELERATE_CONFIGURATION_FILE]

AA_dir=$(pwd)
$AA_dir/scripts/astro-accelerate.sh $AA_dir $AA/AA_Test_output $AA_dir/input_files/GMRT.SPS.10m_newAA.txt