echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
echo        This script is depreceated.
echo  Run make from the project directory instead.
echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd src
/usr/local/opt/llvm/bin/clang++ -o ../build/plasMD -O3 -std=c++11 -fopenmp -stdlib=libc++ -L/usr/local/opt/llvm/lib -I/usr/local/opt/llvm/include -I/Users/alaric-mba/Programming/include main.cpp Constant.cpp DecayRates.cpp EigenSolver.cpp Grid.cpp HartreeFock.cpp Input.cpp IntegrateRateEquation.cpp Numerics.cpp Potential.cpp RadialWF.cpp ComputeRateParam.cpp Plasma.cpp Wigner/wignerSymbols.cpp
