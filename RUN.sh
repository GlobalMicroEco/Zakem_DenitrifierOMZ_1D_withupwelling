#!/bin/sh
# Simple bash script to run and organize output

#/Applications/Julia-1.9.app/Contents/Resources/julia/bin/julia run_model.jl 
/Users/dccoy/.juliaup/bin/julia run_model.jl
mv out*nc out/.
