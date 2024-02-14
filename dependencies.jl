#############################################################################################################
# Script to calculate dependent fields (i.e., 1D grid based on settings in settings.jl)
#############################################################################################################

###################################################################
# Process 1D grid settings
###################################################################
# Get number of grid cells
ngrid = Int(H/dz) 
# Get grid centers
zc = [dz/2 : dz : H - dz/2;] 
# Get grid faces
zf = [0 : dz : H;] 

###################################################################
# Get number of bacteria 
###################################################################
nb = nchemoautos + nhets

###################################################################
# Get diffusion profile ('law of the wall' means this goes to 0 at model surface and bottom)
###################################################################
kappa_z = (kappazmax .* exp.(-zf/mlz) .+ kappazmin .+ kappazmax .* exp.((zf .- H) / 100.)) .* 3600 .* 24 
kappa_z[1] = 0
kappa_z[end] = 0

###################################################################
# Get light attenuation profile 
###################################################################
Light = Light_top .* exp.( -zc ./ euz)

###################################################################
# Get gas exchange at surface grid cell 
###################################################################
koverh = Kgast/dz 

###################################################################
# Get temperature profile 
###################################################################
# A simple profile (fit to N. Pacific oligo gyre data from Zakem et al 2018 Nat Comm)
Temp = temp_surface .*exp.(-zc./ thermocline_shallow) .+ temp_surface .*exp.(-zc ./ thermocline_deep) .+ temp_min 

###################################################################
# Get temperature modification to metabolic rates 
###################################################################
TempFun = TempCoeffArr .* exp.(TempAeArr .*(1 ./ (Temp .+ Tkel) .- 1 ./ TemprefArr)) 

###################################################################
# Get organic matter distribution(?) 
###################################################################
# distribution of m to d (assign dead cell to the number of OM pools I have. For now, equal assign to the 2 pools)
prob_generate_d = zeros(nd) #need to make it a vector 
pdfx = ones(nd) #uniform distribution
prob_generate_d[:] .= pdfx/sum(pdfx)
# UNUSED FOR NOW
#prob_generate_d = prob_generate_d / sum(prob_generate_d)
#prob_generate_d = transpose(prob_generate_d)
#Lognormal distribution for regularly spaced vmax
#x = log.(vmax_i)
#μ = 0 #log mean
#σ = 1 #log stdev
#pdfx = 1 ./(sqrt(2*pi*σ^2)).*exp.(-(x .- μ).^2 ./ (2*σ^2) )

###################################################################
# Get sinking speed of ocean (ws) and POM (wd) in m/d
###################################################################
ws = sink_speed*ones(nd) 
# Xin_let 90% of OM be DOM, 10% be POM
if DOM1  == 1 #DOM1_90perc
    ws[1] = 0 #, the 1st OM not sinking (Dissolve)
elseif DOM1 == 0
    #prob_generate_d[1] = 0.9
    #prob_generate_d[2] = 0.1
end
wd         = transpose(repeat(ws, 1, ngrid + 1)) + repeat(w, 1, nd)  
wd[1,:]   .= 0   # no flux boundary at surface (nothing to sink to the first box)
wd[end,:] .= 0   # no flux boundary at  bottom (accumulates D, allow nothing to sink out of domain)

###################################################################
# Get initial conditions (IC)
###################################################################
if finput == "newIC" 
    pIC = Array{Float64, 2}(undef, ngrid, np) # the initial biomass condition for the np p for all boxes at time 0
    pIC = fill!(pIC, 0.1)
    bIC = Array{Float64, 2}(undef, ngrid, nb) # the initial biomass condition for the nb b for all boxes at time 0
    bIC = fill!(bIC, 0.1/nb)
    #bIC[:,2] .= 0 #turn off type 2 
    zIC = Array{Float64, 2}(undef, ngrid, nz) # the initial biomass condition for the nz z for all boxes at time 0
    zIC = fill!(zIC, 0.1)
    nIC = Array{Float64, 2}(undef, ngrid, nn) # the initial concentration condition for the nn n for all boxes at time 0
    #nIC = fill!(nIC, 0.001)
    nIC[:,1] = nh4_obs # set to Bianchi N tracer obs 
    nIC[:,2] = no2_obs # set to Bianchi N tracer obs
    nIC[:,3] = no3_obs # set to WOA-18 
    nIC[:,4] = n2o_obs # set to Bianchi N tracer obs (convert)
    nIC[:,5] .= 0
    dIC = Array{Float64, 2}(undef, ngrid, nd) # the initial condition for the nd d for all boxes at time 0
    dIC = fill!(dIC, 0.001)
    oIC = Array{Float64, 2}(undef, ngrid, 1) # the initial condition for the nd d for all boxes at time 0
    oIC = o2_obs # set to WOA-18
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
