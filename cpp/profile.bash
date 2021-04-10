#!/bin/bash
g++ -std=gnu++2a -g -lprofiler -o $1 $1.cpp -L/usr/local/Cellar/gperftools/2.9.1/lib -I/usr/local/Cellar/gperftools/2.9.1/include/gperftools
env CPUPROFILE=$1.prof CPUPROFILE_FREQUENCY=10000 ./$1
pprof --pdf $1 $1.prof > $1.pdf
#echo "Profiling results: $1.pdf"
