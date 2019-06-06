#!/bin/bash

for temperature in sub-adiabatic; do
	outputfile=mass_flux_error
	echo -n "# Resolution ala isentropic hydrostatic projected-density" > $outputfile
	  for repetitions in 2 4 8 16 32 64 128; do
		  echo '' >> $outputfile
		  echo -n $repetitions >> $outputfile
  for formulation in ala isentropic hydrostatic projected-density; do #projected-density
    data_folder=output-lateral-pipe-advect-jump-repetitions-${repetitions}-${formulation}-${temperature}
    echo $data_folder
    cat $data_folder/statistics | tail -n 1 | gawk '{printf " %g",sqrt(($25-71410462.65337557)*($25-71410462.65337557))/(71410462.65337557)}' >> $outputfile
    done
  done
done
