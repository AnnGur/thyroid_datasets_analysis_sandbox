#!/bin/bash

# Base URL
BASE_URL="https://download.cncb.ac.cn/gsa-human/HRA001684"

# Read md5sum.txt line by line
while read -r md5 filepath; do
    # Remove the initial /HRA001684 from filepath as it's already in BASE_URL
    filepath_trimmed=$(echo "$filepath" | sed 's/^\/HRA001684//')
    
    # Extract directory path
    dir=$(dirname "$filepath")
    # Create directory structure
    mkdir -p ".$dir"
    
    # Construct full URL
    file_url="$BASE_URL$filepath_trimmed"
    
    echo "Downloading: $file_url"
    # Download file to corresponding directory
    wget -P ".$dir" "$file_url"
    
    # Only verify MD5 if download was successful
    if [ -f ".$filepath" ]; then
        calculated_md5=$(md5sum ".$filepath" | cut -d' ' -f1)
        if [ "$calculated_md5" = "$md5" ]; then
            echo "MD5 checksum verified for $filepath"
        else
            echo "WARNING: MD5 checksum mismatch for $filepath"
        fi
    else
        echo "ERROR: Download failed for $filepath"
    fi
done < md5sum.txt