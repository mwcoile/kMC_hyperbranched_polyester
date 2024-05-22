#!/bin/bash

# Loop through all files in the current directory starting with "p."
for file_path in p.*; do
    # Check if the file exists
    if [ ! -e "$file_path" ]; then
        echo "File not found: $file_path"
        continue
    fi

    # Get the last line of the file
    last_line=$(tail -n 1 "$file_path")

    # Check if the last line starts with "total"
    if [[ "$last_line" == total* ]]; then
        echo "File '$file_path': GOOD: The last line starts with 'total'."
    else
        echo "File '$file_path': BAD: The last line does not start with 'total'."
    fi
done
