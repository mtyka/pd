#!/bin/bash

if [[ $# -lt 1 ]]; then
   echo "Syntax: $0 <benchname>"
   exit 0
fi

tag=$1

## create
cd _$tag
for i in $(ls -d bench_*); do 
 echo Running $i '-->' log.$i
 cd $i
 python run.py &> ../log.$i 
 cd ..
done


