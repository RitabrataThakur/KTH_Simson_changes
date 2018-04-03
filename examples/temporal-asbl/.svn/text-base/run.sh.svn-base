#!/bin/sh

# Script to run the spatial Falkner-Skan-Cooke example with cross flow vortices
echo '####################################################################'
echo 'Running the temporal asymptotic suction boundary layer'
echo 'This will take a while...'
echo '####################################################################'

export OMP_NUM_THREADS='1'

# Generate initial velocity field bl0.u
echo 'Generating initial velocity field...'
cp bls.localized.i bls.i
../../bls/bls > out2.dat

# Run simulation
echo 'Running simulation...'
../../bla/bla > out3.dat

echo '####################################################################'
echo 'Finished'
echo '####################################################################'

