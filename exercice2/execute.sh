gcc main.c -O0 -o main_no_optimisation
gcc main.c -O3 -fno-tree-vectorize -o main_no_vector
gcc main.c -O3 -ffast-math -march=native -o main_vector

echo "RUN NO OPTIMISATION"

./main_no_optimisation 100
./main_no_optimisation 200
./main_no_optimisation 300
./main_no_optimisation 400
./main_no_optimisation 500
./main_no_optimisation 600
./main_no_optimisation 700
./main_no_optimisation 800
./main_no_optimisation 900
./main_no_optimisation 1000
./main_no_optimisation 1100
./main_no_optimisation 1200
./main_no_optimisation 1300
./main_no_optimisation 1400
./main_no_optimisation 1500

echo "RUN OPTIMISATION"

./main_no_vector 100
./main_no_vector 200
./main_no_vector 300
./main_no_vector 400
./main_no_vector 500
./main_no_vector 600
./main_no_vector 700
./main_no_vector 800
./main_no_vector 900
./main_no_vector 1000
./main_no_vector 1100
./main_no_vector 1200
./main_no_vector 1300
./main_no_vector 1400
./main_no_vector 1500


echo "RUN VECTOR"

./main_vector 100
./main_vector 200
./main_vector 300
./main_vector 400
./main_vector 500
./main_vector 600
./main_vector 700
./main_vector 800
./main_vector 900
./main_vector 1000
./main_vector 1100
./main_vector 1200
./main_vector 1300
./main_vector 1400
./main_vector 1500

rm main_no_optimisation
rm main_no_vector
rm main_vector
