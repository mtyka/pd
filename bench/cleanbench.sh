#!/bin/bash 

rm `find bench/bench_* -name "_pd.so"` `find bench/bench_* -name "pd.py"` `find bench/bench_* -name "inout.pdb"` `find bench/bench_* -name "*.pyc"` -f
