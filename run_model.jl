#############################################################################################################
# Run the model!
#############################################################################################################

# Load packages 
using DelimitedFiles
using NCDatasets
using DataFrames
using Printf
using Dates
using SparseArrays, LinearAlgebra

# Import settings
include("settings.jl")

# Calculate dependencies
include("dependencies.jl")

# Import microbial traits
include("inputs/traits.jl")

# Import grazing and mortality parameters 
include("inputs/grazingandmortality.jl")

# Initialize output struct 
include("initialize.jl")

# Import functions
include("functions.jl")

# Run the model and save output
p, b, z, n, d, o = run_npzdb(paras1)

#############################################################################################################
# Functions below 
#############################################################################################################


