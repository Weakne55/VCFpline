#!/bin/bash

input_file="gencode.v45.annotation.gtf"
output_file="gencode.v45.annotation.NCnames.gtf"

replace=(NC_000001.11 NC_000002.12 NC_000003.12 NC_000004.12 NC_000005.10 NC_000006.12 NC_000007.14 NC_000008.11 NC_000009.12 NC_000010.11 NC_000011.10 NC_000012.12 NC_000013.11 NC_000014.9 NC_000015.10 NC_000016.10 NC_000017.11 NC_000018.10 NC_000019.10 NC_000020.11 NC_000021.9 NC_000022.11 NC_000023.11 NC_000024.10 NC_012920.1)

temp_file=$(mktemp)

cp "$input_file" "$temp_file"
for i in {1..22}; do
	    sed -i "s/^chr${i}\b/${replace[$((i))]}/" "$temp_file"
    done

sed -i "s/^chrX\b/${replace[23]}/" "$temp_file"
sed -i "s/^chrY\b/${replace[24]}/" "$temp_file"
sed -i "s/^chrM\b/${replace[25]}/" "$temp_file"

mv "$temp_file" "$output_file"

echo "All done, results in $output_file"
