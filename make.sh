# @breif Makes the project bc im done writing out a fuck ass terminal prompt and
#        cmake is a fucking WHORE
# @author Garrett Rhoads
# @file make.sh
# @date 07/10/2025

cpp_filename="catmodel.cc"
exac_filename="catmodel"

echo Making $exac_filename exacutable from $cpp_filename
g++ -I/Applications/Wolfram.app/Contents/SystemFiles/Links/WSTP/DeveloperKit/MacOSX-x86-64/CompilerAdditions \
    -L/Applications/Wolfram.app/Contents/SystemFiles/Links/WSTP/DeveloperKit/MacOSX-x86-64/CompilerAdditions \
    -lWSTPi4 -framework Foundation \
    $cpp_filename -o $exac_filename \
    -lcln -lginac 
echo $exac_filename exacutable compiled