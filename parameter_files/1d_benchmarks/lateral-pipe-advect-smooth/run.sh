#!/bin/bash

# Run all combinations of models in this folder

for repetitions in 2 4 8 16 32 64 128; do #1 2 4 8 16 32 64 128
  echo "subsection Geometry model" > repetitions.prm
  echo "subsection Box" >> repetitions.prm
  echo "set X repetitions = $repetitions" >> repetitions.prm
  echo "end" >> repetitions.prm
  echo "end" >> repetitions.prm

  for temperature in sub-adiabatic; do #sub-adiabatic adiabatic super-adiabatic
    for formulation in ala isentropic hydrostatic projected-density; do #ala isothermal hydrostatic projected-density
      output_folder=output-lateral-pipe-advect-jump-repetitions-${repetitions}-${formulation}-${temperature}
      echo "set Output directory = ${output_folder}" > output.prm
      cat lateral-pipe.prm ${temperature}.prm ${formulation}.prm repetitions.prm output.prm | ../../../plugins/aspect --
    done
  done
done

rm output.prm
rm repetitions.prm

bash make_statistics.sh
