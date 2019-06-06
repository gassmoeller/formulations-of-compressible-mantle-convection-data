#!/bin/bash

for temperature in sub-adiabatic adiabatic; do
	outputfile=mass_flux_error_${temperature}
	echo -n "# Resolution ala isentropic hydrostatic projected-density" > $outputfile
	  for repetitions in 2 4 8 16 32 64 128; do
		  echo '' >> $outputfile
		  echo -n $repetitions >> $outputfile
  for formulation in ala isentropic hydrostatic projected-density; do #projected-density
    data_folder=output-lateral-pipe-repetitions-${repetitions}-${formulation}-${temperature}
    cat $data_folder/statistics | tail -n 2 | head -n 1 | gawk '{printf " %g",sqrt(($20+$21)*($20+$21))/$21}' >> $outputfile
    done
  done
done
