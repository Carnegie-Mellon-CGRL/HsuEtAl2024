This code is for the optimization of tissue engineered vascular graft microstructure and material behavior. The directory GNR includes the code for simulating the evolving geometry and composition of an ovine TEVG based on the data/parameters in Blum et al 2022 and Latorre et al 2023. The SMF directory includes a surrogate management framework optimization algorithm for constrained optimization. Export includes wrapper scripts to run optimization simulations and plot includes useful scripts for visualizing results.

Both the GNR and SMF codes are compiled using the Makefile in their respective directories. Running the command 'make' should generate an executable file 'gnr' and 'krig', respectively. Both codes are written in C++ and compiled with the GNU compilers. GNR also includes dependencies on boost and gsl. On Linux, one can run:
sudo apt-get install build-essential
sudo apt-get install libgsl-dev
sudo apt-get install libboost-all-dev

These can also be installed through Brew on macOS or using WSL on Windows. A simple code to test the installation of gsl is available in the GNR directory, 'boost_test.cpp'. It can be compiled and run with:
gcc -c boost_test.c
gcc boost_test.o -lgsl -lgslcblas -lm
./a.out

After compilation, the executables should be moved to the export folder. An optimization can be run from 'run_script_gnr.m', which will setup 10 directories and run 10 instances of the optimization to account for stochasticity in the initialization.
