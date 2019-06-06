#!/bin/bash

# Run all combinations of models in this folder

repetitions=8

#for repetitions in 1 2 4 8 16 32 64; do #4 8 16 32 64
  echo "subsection Geometry model" > repetitions_${repetitions}.prm
  echo "subsection Box" >> repetitions_${repetitions}.prm
  echo "set X repetitions = $repetitions" >> repetitions_${repetitions}.prm
  echo "end" >> repetitions_${repetitions}.prm
  echo "end" >> repetitions_${repetitions}.prm

for maximum_timestep in 31557600e7 15778800e7 7889400e7 3944700e7 1972350e7 986175e7 4930875e6 24654375e5; do #31557600e7 15778800e7 7889400e7 3944700e7 1972350e7 986175e7
  echo "set Maximum time step = $maximum_timestep" > timestep.prm
  for formulation in ala isentropic hydrostatic projected-density; do #ala isothermal hydrostatic projected-density
    output_folder=output-lateral-pipe-transient-repetitions-${repetitions}-${formulation}-${maximum_timestep}
    echo "set Output directory = ${output_folder}" > output.prm
    cat lateral-pipe-transient.prm ${formulation}.prm repetitions_${repetitions}.prm output.prm timestep.prm | ../../../plugins/aspect --
  done
done

rm output.prm
rm repetitions_*.prm

bash make_statistics.sh
