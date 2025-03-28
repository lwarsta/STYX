#!/bin/bash

# Output file
OUTPUT_FILE="full_source.txt"

# Ensure the output file is empty
> "$OUTPUT_FILE"

# Declare an associative array to store unique base names
declare -A baseNames

# Process header files: add base names (without .h extension)
for file in include/*.h; do
    if [[ -f "$file" ]]; then
        base=$(basename "$file" .h)
        baseNames["$base"]=1
    fi
done

# Process source files: add base names (without .cpp extension)
for file in src/*.cpp; do
    if [[ -f "$file" ]]; then
        base=$(basename "$file" .cpp)
        baseNames["$base"]=1
    fi
done

# Iterate over all collected base names
for base in "${!baseNames[@]}"; do
    echo "=============================" >> "$OUTPUT_FILE"
    echo "Base: $base" >> "$OUTPUT_FILE"
    echo "=============================" >> "$OUTPUT_FILE"
    
    # Check for header file
    header_file="include/$base.h"
    if [[ -f "$header_file" ]]; then
        echo "Header: $header_file" >> "$OUTPUT_FILE"
        cat "$header_file" >> "$OUTPUT_FILE"
    else
        echo "Header: None" >> "$OUTPUT_FILE"
    fi
    echo -e "\n-----------------------------\n" >> "$OUTPUT_FILE"
    
    # Check for source file
    source_file="src/$base.cpp"
    if [[ -f "$source_file" ]]; then
        echo "Source: $source_file" >> "$OUTPUT_FILE"
        cat "$source_file" >> "$OUTPUT_FILE"
    else
        echo "Source: None" >> "$OUTPUT_FILE"
    fi
    
    echo -e "\n\n" >> "$OUTPUT_FILE"
done

echo "Concatenation complete. Output saved to $OUTPUT_FILE."
