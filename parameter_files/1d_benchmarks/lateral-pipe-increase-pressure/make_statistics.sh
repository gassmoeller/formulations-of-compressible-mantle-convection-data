#!/bin/bash

repetitions=(2 4 8 16 32 64 128)
maximum_timestep=(15778800e7 7889400e7 3944700e7 1972350e7 986175e7 4930875e6 24654375e5)
outputfile=mass_flux_error
echo -n "# Timestep ala isentropic hydrostatic projected-density projected-density-full-pressure" > $outputfile

for index in 0 1 2 3 4 5 6; do #4 8 16 32 64
  echo '' >> $outputfile
  echo -n ${maximum_timestep[$index]} >> $outputfile
  for formulation in ala isentropic hydrostatic projected-density projected-density-full-pressure; do #projected-density
    data_folder=output-lateral-pipe-transient-repetitions-${repetitions[$index]}-${formulation}-${maximum_timestep[$index]}
    echo $data_folder $outputfile
    #expected_mass_flux=4.184501057159e-04
    expected_mass_flux=0.0011382473816094456
    cat $data_folder/statistics | tail -n 1 | gawk "{printf \" %g\",sqrt((\$21+$expected_mass_flux)*(\$21+$expected_mass_flux))/$expected_mass_flux}" >> $outputfile
  done
done



