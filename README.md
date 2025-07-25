# Catenary-Identifiablity

### Written by Brett Hajdaj and Garrett Rhoads

## What is this program?

This program is a way to calculate various useful attributes of catenary models including: identifiability, rank, and singular locus. It accomplishes this through a combination of C++ libraries such as GiNaC and CLN as well as making calls to a wolfram mathematica kernel.

## Potential issues

In the code there are a two spots where you may have to change file paths to run the code, you also must have the desktop version of mathematica and have it be activated by a valid account. 

## Compilation

Working from a MacBook, I was not able to get CMake working properly though you are fully encouraged to give it a shot, I was able to get it working on linux so ymmv. I have instead written an extremely basic bash file `make.sh` that can be run with `./make.sh` that compiles the program for my machine, this does include some file paths that will again, most likely have to be changed to fit your needs.

## Use

You can run single catenary models by just running `./catmodel` or if you would like to run all possible catenary models up to a specific graph size and write the data to a .csv file, you can do that by running it with the `-s` command line argument.