#!/bin/bash

# Check if a file has been provided as an argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <sam_file>"
    exit 1
fi

# Retrieve the SAM file path
sam_file="$1"

# Check if the file exists and is not empty
if [ ! -f "$sam_file" ]; then
    echo "Error: '$sam_file' is not a valid file."
    exit 1
fi

if [ ! -s "$sam_file" ]; then
    echo "Error: The file '$sam_file' is empty."
    exit 1
fi

# Check if the file contains a SAM header
if ! grep -q "^@" "$sam_file"; then
    echo "Error: The file '$sam_file' does not appear to contain a SAM header."
    exit 1
fi

# If all checks pass
echo "The file '$sam_file' is valid and ready for analysis."
nt
echo "Le fichier '$fichier_sam' est un fichier SAM valide."
exit 0
