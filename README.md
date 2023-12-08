# MicrOMZ 1D ecosystem model 
1D model of an ocean water column coupled to a NPBD model capable of resolving the metabolisms responsible for nitrogen cycling in oxygen minimum zones

## Table of Contents
- [Updates](#updates)
- [Getting started](#getting-started)
- [Code structure](#code-structure)

Requires the Julia programming language

## Updates
* 12/07/23: Danny's first commit of MicrOMZ
* 08/10/23: need to add w to param struct (only wd there now)
* 03/01/23: EJZtinkering beings
* 02/02/23: version to send to Xin Sun
* 01/19/23: add DIN and N-cycling types
* 11/18/22: cont
* 09/14/22: add oxygen, also take out things Xin won't use (ex: crossfeeding)
* 09/14/22: modified from Liang Xu's "NPZDBmodel": //github.com/xl0418/NPZDBmodel/blob/main/run_NPZDBmodel.jl

## Getting started
#### Install necessary packages 
    Call 'getPkgs.jl' to install necessary Julia packages
    You should only need to call this once
        > include("getPkgs.jl")
#### Update the settings in the model 
    Edit 'settings.jl' to toggle user-specific settings (see below)
#### Run the model
    After toggling settings, startup Julia then run the following command
        > include("run_model.jl")

## Code structure 
#### settings.jl 
    Main file where user can toggle things like:
        input/output file names
        runtime, output frequency, and timestep length
        1-D grid settings like total depth, and height of each grid cell
        ecosystem settings, like the number of phytoplankton, chemoautotrophs, and heterotrophs
        how the model represents the physical environment, like light attenuation, temperature, etc.
        how the model represents air-sea gas exchange

#### dependencies.jl  
    Script that will process 'settings.jl' to set up the model run

#### functions.jl  
    Script where timestepping and ecosystem model live

#### inputs/traits.jl 
    File where microbial parameters (e.g., traits) live

#### inputs/grazingandmortality.jl 
    File where grazing and mortality parameters live

#### inputs/w2000_10m.txt
    Text file containing vertical velocity

#### out/
    Folder to store output files, which can be used to 'restart' a new run

#### plots/
    Folder containing post-run plotting scripts

## Support
Contact Daniel McCoy, Xin Sun, or Emily Zakem at Carnegie Science
