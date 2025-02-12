#!/bin/bash

while read -r line 
do 
        ./sratoolkit.3.0.0-centos_linux64/bin/fasterq-dump -O ./out -e 10 -3 ./normal_samples/${line} 
        pigz /media/user/New\ Volume/NAFLD/out/*_1.fastq 
        pigz /media/user/New\ Volume/NAFLD/out/*_2.fastq

done <  sample_names.txt
