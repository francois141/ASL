gcc benchmark.c -O3 -mfma -fno-tree-vectorize -o benchmark 
./benchmark
rm benchmark