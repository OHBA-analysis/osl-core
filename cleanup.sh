#!/bin/bash

HERE=$(pwd)
THIS=${BASH_SOURCE}
OSLCORE=$(dirname "$THIS")
TARGET=${1:-$OSLCORE}

cd "$TARGET"

# remove temporary files
echo "Removing temporary files..."
find . -type f \
    -name '.DS_Store*' \
    -o -name '._*' \
    -o -name '*~' \
    -delete

# 644 for Matlab files
echo "Setting 644 permissions for Matlab + text files..."
find . -type f \
    -name '*.m' \
    -o -name '*.md' \
    -o -name '*.txt' \
    -o -name '*.json' \
    -o -name '*.conf' \
    -exec chmod 644 {} +

find . -type f \
    -name '*.mat' \
    -o -name '*.dat' \
    -o -name '*.nii*' \
    -exec chmod 644 {} +

# 755 for shell scripts
echo "Setting 755 permissions for shell scripts..."
find . -type f -name '*.sh' -exec chmod 755 {} +

# all done
cd "$HERE"
echo "Done!"
