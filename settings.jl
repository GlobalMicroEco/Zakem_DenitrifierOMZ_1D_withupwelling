
###################################################################
# Set input and output filenames
###################################################################
# Set output filename
# Set to "missing" if output isn't to be saved
folder = "out/"
fsave  = "out_addW"                  
# Set 'restart' file
# Set to "newIC" if fresh restart is requested
finput = "out_addW_20230907_1051.nc" 

###################################################################
# Runtime settings (+ number of timesteps to be output) 
###################################################################
# Set run time (days)
tt = 365*10 
# Set number of timesteps to record
# NOTE: Use 1 to only save the last timestep
nrec = 100 
# Set timestep (per day)
dt = 1e-3 

###################################################################
# 1-D grid settings 
###################################################################
# Set total depth (meters) of 1D water column
H = 2000  
# Set dz (meters) of each cell
dz = 10   

###################################################################
# Ecosystem settings 
###################################################################
# Set number of phytoplankton (np) and zooplankton (nz)
np = 1           
nz = 1           
# Set number of chemoautotrophs and heterotrophs
nchemoautos = 3  
nhets       = 7        
# Switch for variable (1) or constant (not 1) half-saturation concentrations (K_n)
diffks_n = 1  
# Temperature Modification to Metabolic Rates
TempCoeffArr = 0.8
TempAeArr    = -4000
TemprefArr   = 293.15   
Tkel         = 273.15

##################################################################
# Settings for organic matter 
##################################################################
# Set number of organic matter pools
nd = 1  
# Set sinking speed 
sink_speed = 8
# if let OM1 be DOM, use 1, if other than 1, then two POMs.
DOM1 = 0
# UNUSED FOR NOW
#   if let OM1 be DOM and 90%detritus be DOM1, use 1, if other than 1, then two POMs, 50%, 50%.
#   DOM1_90perc = 1 

#################################################################
# Settings for DIN 
#################################################################
# Set number of DIN pools (default = 5, corresponding to [NH4, NO2, NO3, N2O, N2]
nn = 5 

#################################################################
# Settings for physical environment 
#################################################################
# Set temperature (degC) profile settings
temp_surface        = 12   # max @ surface
temp_min            = 2    # min @ depth
thermocline_shallow = 150  # defines thermocline depth
thermocline_deep    = 500  # defines bottom of thermocline
# Set mixed layer length scale (meters) 
mlz = 20  
# Set deep ocean diffusion coefficient (m2/s)
# NOTE: This should be small number, since diffusion is weak in @ depth due to weak vertical gradients
kappazmin = 5e-5 
# Set surface ocean diffusion coefficient (m2/s)
# NOTE: This should be greater than kappazmin due to strong vertical gradients near surface
# NOTE: This is also applied @ bottom boundary layer
kappazmax = 1e-2 
# Set vertical velocity
w2D = readdlm("inputs/w2000_10m.txt"); # vertical velocity from EJZ's 2D FLOW model
w   = w2D[90,:];                # peak upwelling velocity
# Set euphotic zone lengthscale (meters) for light attenuation  
euz = 25 
# Set surface irradiance (W/mw)
# NOTE: If incoming PAR = (1400/2), then Light_avg*(cos(t*dt*2*3.1416)+1) for light daily cycle
Light_top = 700   

#################################################################
# Settings for oxygen (air-sea gas exchange, deep ocean values) 
#################################################################
# Set oxygen saturation at surface (mmol/m3) assuming T = 25C, S = 35 PSU (~WOCE values in ETSP)
o2sat = 212.1 
# Set gas transfer coefficient (kw, in meters/day)
# NOTE: See Sarmiento & Gruber (2006) Table 3.3.2
Kgast = 10/100*24 #m/d from 10 cm/hr
# Set deep ocean oxygen (mmol/m3)
o2_deep = 10
# Set timescale (1/day) for oxygen relaxation
# NOTE: Set to 0 to turn off
t_o2relax = 0.01


