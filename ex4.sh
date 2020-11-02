#!/bin/bash

for ((i = 1; i < $#; i++ )); do
	printf '%s\n' "Arg $i: ${!i}"
	if [ ! -f "${!i}" ]; then
		echo "${!i} is not a regular file"
		exit 1
	fi
done

eval "$(conda shell.bash hook)"
conda activate bio

if [ ! -d "bwa_index" ]; then
	mkdir bwa_index

	echo "[!!] BWA indexing... [!!]"
	bwa index -p bwa_index/index "$1" 
else
	echo "[!!] Indexing already found on disk, going straight to aligning [!!]"
fi

 
echo "[!!] BWA aligning... [!!]"
bwa mem bwa_index/index "$2" "$3" > "$4"

echo "[!!] Finished! Check results in $4 [!!]"