#!/bin/bash

if [[ $# -lt 1 ]]; then
   echo "Syntax $0 <benchname> [ <alternative path to pd.py and _pd.so>]"
   exit 0
fi

tag=$1
binpath="/work/mtyka/trunk/bin"

if [[ $# -gt 1 ]]; then
  binpath=$2
fi

## create
rm _$tag -rf
mkdir _$tag
cd bench
benches=$(ls -d bench_*)
cd ..
for i in $benches; do cp -r bench/$i _$tag; done
for i in $benches; do cp "$binpath"/_pd.so   "$binpath"/pd.py _$tag/$i/; done

echo "A benchmark set has been created in _$tag"
echo "Type 'sh runbench.sh $tag' to run the benchmarks"
