#!/bin/bash

for temperature in sub-adiabatic adiabatic; do
	outputfile=mass_flux_error_${temperature}
	echo -n "# Resolution ala isentropic hydrostatic projected-density" > $outputfile
	  for repetitions in 2 4 8 16 32 64 128; do
		  echo '' >> $outputfile
		  echo -n $repetitions >> $outputfile
  for formulation in ala isentropic hydrostatic projected-density; do #projected-density
    data_folder=output-vertical-pipe-repetitions-${repetitions}-${formulation}-${temperature}
    cat $data_folder/statistics | tail -n 1 | gawk '{printf " %g",sqrt(($22+$23)*($22+$23))/sqrt($23*$23)}' >> $outputfile
    done
  done
  echo '' >> $outputfile
done


