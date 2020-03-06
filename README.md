# Formulation of Compressible Mantle Convection: Data

This repository accompanies the paper
```
Formulation of Compressible Mantle Convection
Rene Gassmoeller, Juliane Dannberg, Wolfgang Bangerth, Timo Heister, Robert Myhill
Geophysical Journal International, Volume 221, Issue 2, May 2020, Pages 1264–1280.
```
and contains data files and source code to reproduce the computations of the paper available at https://doi.org/10.1093/gji/ggaa078.

Contents
--------

- Parameter files for ASPECT to reproduce the benchmarks and application models
- Source files necessary for the computations
- Postprocessing and plotting scripts to recreate the figures.

see below for the instruction on how to run these cases.

Instructions
------------

The first step for reproducing the results is installing ASPECT, or using a docker image or virtual machine. Instructions are provided in the [ASPECT manual](http://www.math.clemson.edu/%7Eheister/manual.pdf).

Next the plugins in the folder plugins need to be compiled as shared libraries as described [here](http://www.math.clemson.edu/~heister/manual.pdf#sec%3Awrite-plugin).

Now, all models in the parameter_files folder can be run. The benchmarks in parameter_files/1d_benchmarks come with shell scripts called run.sh that start all necessary model variations, and with gnuplot scripts (in the subdirectory figures) that plot the results.

The computations of the paper were produced with the following [ASPECT version](https://github.com/gassmoeller/aspect/releases/tag/compressible-formulations-submission). Official ASPECT releases starting with release 2.2.0 will also support the contained features (but might require minor changes to the parameter files).
