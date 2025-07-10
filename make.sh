# @breif Makes the project bc im done writing out a fuck ass terminal prompt and
#        cmake is a fucking WHORE
# @author Garrett Rhoads
# @file make.sh
# @date 07/10/2025

cpp_filename="catmodel.cc"
exac_filename="catmodel"

echo Making $exac_filename exacutable from $cpp_filename
g++ -o $exac_filename $cpp_filename -lginac -lcln
echo $exac_filename exacutable compiled