#!/bin/bash

repetitions=8
outputfile=mass_flux_error
echo -n "# Timestep ala isentropic hydrostatic projected-density" > $outputfile
#for repetitions in 1 2 4 8 16 32 64; do
for maximum_timestep in 15778800e7 7889400e7 3944700e7 1972350e7 986175e7 4930875e6 24654375e5; do
  echo '' >> $outputfile
  echo -n $maximum_timestep >> $outputfile
  for formulation in ala isentropic hydrostatic projected-density; do #projected-density
    data_folder=output-lateral-pipe-transient-repetitions-${repetitions}-${formulation}-${maximum_timestep}
    cat $data_folder/statistics | tail -n 1 | gawk '{printf " %g",sqrt(($21-0.02050001028)*($21-0.02050001028))/0.02050001028}' >> $outputfile
  done
done
