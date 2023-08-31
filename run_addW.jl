#8/10/23: need to add w to param struct (only wd there now)
#3/1/23: EJZtinkering beings
#2/2/23: version to send to Xin Sun
#1/19/23: add DIN and N-cycling types
#11/18/22: cont
#9/14/22: add oxygen, also take out things Xin won't use (ex: crossfeeding) 
#9/14/22: modified from Liang Xu's "NPZDBmodel": //github.com/xl0418/NPZDBmodel/blob/main/run_NPZDBmodel.jl

using NCDatasets
using DataFrames

######################################################################################
# save file, date is au
fsave = "out_addW" #fsave = "missing" #use this instead to not save if on, then output isn't saved

# input file for IC
#folder = ""
#finput = "out_addW_20230830_2032.nc" #10 yrs total run (kappa = 5e-5)
#OR
finput = "newIC" #if this is on, then don't load a file for input and specify IC below

# Things vary from run to run
tt = 365*3#0 # days -- run time (days)
nrec = 100 # number of timepoints to record (e.g., if nrec=1, only save the last time point)

## Xin: some choices
#if use different K_n values, use 1, if want same K_n, use anything that is not 1
diffks_n = 0 

#if let OM1 be DOM, use 1, if other than 1, then two POMs.
DOM1 = 0
#if let OM1 be DOM and 90%detritus be DOM1, use 1, if other than 1, then two POMs, 50%, 50%.
#DOM1_90perc = 1 

# set the water column depth and number of 0-D making up the 1-D model
H = 2000  # the total depth of the water column (meters)
dz = 10  # height of each box (meters)

np = 1 # number of phytoplankton
nb = 10  # number of bac (and archaea): o, n1 to n6, aoa, nob, and anx 
# !The two numbers below must add up to nb.
nchemoautos = 3
nhets = 7 # aerobic hets + denitrifiers
if nhets + nchemoautos != nb
    println("PROBLEM: nhets + nchemoautos != nb")
    return
end
#@test nhets + nchemoautos == nb
nd = 1  # number of organic matter pools, if use 2
nz = 1 # number of zooplankton: one for p and one for b for now (not sufficient since b are different sizes)
nn = 5 # NH4, NO2, NO3, N2O, N2 #number of inorganic matter pools
# !order of the nn types matters!! match this with code below!

######################################################################################
#set up grid
ngrid = Int(H/dz) # number of grids or boxes = total height of the water column/height of each box
zc = [dz/2 : dz : H - dz/2;] # centered depth array (; makes it vertical), start with dz/2, interval=dz, end with H - dz/2. (e.g., a surface 10 m box's center depth is 5m)
zf = [0 : dz : H;] #face depth array; 1 number longer than zc ( = total horizontal edges of the boxes)

######################################################################################
## Physical Environment 

# Vertical mixing
mlz=20 # !m mixed layer lengthscale 
# (this is the only parameter to worry about in Vertical mixing, change this number will change how well mixed or stratified the water column is)
kappazmin=5e-5 #(small number = slower mixing) !m2/s --min mixing coefficient -value for most of the deep ocean (higher at top and bottom)
kappazmax=1e-2 # !m2/s --max mixing coefficient -value at top of mixed layer (as well as bottom boundary mixed layer)
kappa_z=(kappazmax .* exp.(-zf/mlz) .+ kappazmin .+ kappazmax .* exp.((zf .- H) / 100.)) .* 3600 .* 24 
kappa_z[1]=0
kappa_z[end]=0
#we also has strong mixing at the bottom, to avoid the OM piled up at the deepest box which will mess up the model. 

#Vertical velocity
using DelimitedFiles
w2D = readdlm("w2000_10m.txt"); #vertical velocity from EJZ's 2D FLOW model
w = w2D[90,:]; #peak upwelling velocity
#Oxygenated: no vertical velocity
w[:] .= 0.
#using Plots
#plot(w[1:end-1],-zc)


# Light (Irradiance I) 
euz=25 # !m -- euphotic zone lengthscale
Light_top = 700   ## the light hits the surface ocean, W/m2 avg incoming PAR = (1400/2)  Light_avg*(cos(t*dt*2*3.1416)+1) for light daily cycle
Light = Light_top .* exp.( -zc ./ euz)

#Oxygen
#air-sea exchange
o2sat=212.1 #mmol/m3 from calc_oxsat(25+273,35) in matlab. WOCE clim-avg surf T at 10S, E. Pac.
#Kgast=3e-5*86400 #m/d: air-sea gas transfer coefficient kw
#from S&G Table 3.3.2:
Kgast = 10/100*24 #m/d from 10 cm/hr
#mlboxes=100/dz #! 100 m is mixed, dz length of the box, this is the # of box constitues mix layer; discrete n of boxes in the mixed layer, close to 100m total sum
#koverh=Kgast/100/mlboxes #gas transfer coefficient for each of the n boxes comprising the ml. why the 100 here?
koverh = Kgast/dz #transfer just into first box!
#eqmask(:)=0D0
#eqmask(3:mlboxes+2)=1D0; !mask for air-sea equilibration
#oxygen relaxation 
o2_deep = 100. #mmol/m3 # assign O2 to the deep water
t_o2relax = 0.01 # unit (1/day), = 0 to turn off (to supply O2 to all boxes via lateral fluxes)

#air-sea exchange for N2 and N2O


###########################################################################################################
#Temperature Function 

#Temperature: A simple profile (fit to N. Pacific oligo gyre data from Zakem et al 2018 Nat Comm)
Temp = 12 .*exp.(-zc./ 150) .+ 12 .*exp.(-zc ./ 500) .+ 2 # temperature profile

# Temperature Modification to Metabolic Rates
TempCoeffArr = 0.8
TempAeArr = -4000
TemprefArr = 293.15   
Tkel = 273.15
TempFun = TempCoeffArr .* exp.(TempAeArr .*(1 ./ (Temp .+ Tkel) .- 1 ./ TemprefArr)) # this modify all rates based on temperature at depth (e.g. 0.5 * rates)


######################################################################################
## Microbial parameters ##############################################################

include("input_traits.jl")

include("input_grazingandmortality.jl")

include("input_OMdetail.jl") #Organic Matter ("d") detail: sinking speeds and probability of generation


###########################################################################################################
#Initial conditions (IC)

if finput == "newIC" 
    pIC = Array{Float64, 2}(undef, ngrid, np) # the initial biomass condition for the np p for all boxes at time 0
    pIC = fill!(pIC, 0.1)
    bIC = Array{Float64, 2}(undef, ngrid, nb) # the initial biomass condition for the nb b for all boxes at time 0
    bIC = fill!(bIC, 0.1/nb)
    #bIC[:,2] .= 0 #turn off type 2 
    zIC = Array{Float64, 2}(undef, ngrid, nz) # the initial biomass condition for the nz z for all boxes at time 0
    zIC = fill!(zIC, 0.1)
    nIC = Array{Float64, 2}(undef, ngrid, nn) # the initial concentration condition for the nn n for all boxes at time 0
    nIC = fill!(nIC, 0.001)
    nIC[:,3] .= 30 #fill NO3 only, the order of nutrients; the  . matches different dimension: means sign the entire vector to a number
    dIC = Array{Float64, 2}(undef, ngrid, nd) # the initial condition for the nd d for all boxes at time 0
    dIC = fill!(dIC, 0.001)
    oIC = Array{Float64, 2}(undef, ngrid, 1) # the initial condition for the nd d for all boxes at time 0
    oIC = fill!(oIC, 10) #oxygenated
else
    println("Loading IC from: ",finput) # input from an output of a previous run.
    ds = NCDataset(folder*finput)
    pIC = ds["p"][:,:,end] #just end point (the last time point of the output of a previous run)
    bIC = ds["b"][:,:,end] 
    zIC = ds["z"][:,:,end] 
    nIC = ds["n"][:,:,end] 
    dIC = ds["d"][:,:,end] 
    #dIConepool = ds["d"][:,:,end] #to input run with just one DIC 
    oIC = ds["o"][:,end]
    oIC = reshape(oIC, length(oIC), 1) #vector into matrix to work with the model
end

# assign the same number to the two OM pool from one OM pool input
#dIC = Array{Float64, 2}(undef, ngrid, nd)
#dIC .= dIConepool #fills all columns with the one pool

#then reset any desired IC (make sure initial O2 is always 1 everywhere even if we use the output O2 from a prevous run):
#oIC .= 10. # uM, set the O2 in each box the same.
#bIC .= 0.01

#reset N2 to see if it is still accumulating, or if it's just residual from spin-up
#nIC[:,5] .= 0.001 # µM set N2 low as initial condition to make sure N2 is not accumulate due to initial spin up.
#nIC[:,4] .= 0.010 # µM set N2O low as initial condition to make sure N2 is not accumulate due to initial spin up.

#！reset denitrifiers to 1e-10
#bIC[:,2:nhets] .= 1e-10
###########################################################################################################
#load model including paras struct:

include("model_addW.jl")

paras1 = paras(tt, nrec, H, dz, np, nb, nhets, nz, nn, nd, pIC, bIC, zIC, nIC, dIC, oIC,
                umax_p, K_pn, m_lp, m_qp, CM, y_d, y_n, y_o, ftoNO3_anx, 
                vmax_d, K_d, vmax_n, K_n, pcoef, m_lb, m_qb,
                g_max, K_g, γ, m_lz, m_qz, prob_generate_d, kappa_z, w, wd, fsave, Light, TempFun,
                K_I, ngrid, GrM,
                e_o, koverh, o2sat, t_o2relax, o2_deep)

###########################################################################################################
#run model and save into fsave_*date*.nc:
#run_npzdb(paras1)
p, b, z, n, d, o = run_npzdb(paras1) # this is calling the function

