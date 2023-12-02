#1/23/22
#input_traits.jl for input into run_addDIN.jl

######################################################################################
# Phytoplankton growth (one type of phytoplankton--Proch)
umax_p = fill(0.515, np) # max growth is from Proch.
K_pn = fill(0.003, np) # half saturation constant for nitrogen uptake
K_I = ones(np) * 10 #light half-saturation constant
e_o = 150/16 #production of O2 (Excretion). mol O2/mol N uptake

######################################################################################
# Heterotrophic and chemoautotrophic microbial growth

# Organic matter Consumption Matrix (CM): i X j, bacteria i feed on d j 
#UPTAKE of OM
#CM = zeros(nd,nb)
#CM = [i == j for i = 1:nd, j=1:nb] #diagonals for specialists -- if nb > nd, the rest of the nb aren't assigned to any d (chemoautos!)
#for one pool of d to start:
CM = zeros(nd, nb)
CM[:,1:nhets] .= 1

#input yield values for organic matter (nd x 1:nhets) and both yields and excretions for din (two matrices, each nn x nb)

#initialize yields
y_d = zeros(nd,nb) #organic matter
y_n = zeros(nn,nb) #inorganic nitrogen: NH4, NO2, NO3, N2O, N2
y_o = zeros(nb) #oxygen

#Organic matter yields for hets (Xin updated 09/2023, P=0.16)
y_d[:,1] .= 0.218 #0.2 #aerobe--keep it low, 0.2 or 0.3

y_d[:,2] .= 0.150 #y_d[:,1]*penalty #NO3 to NO2 
y_d[:,3] .= 0.148 #y_d[:,1]*penalty #NO3 to N2O 
y_d[:,4] .= 0.138 #y_d[:,1]*penalty #NO3 to N2 
y_d[:,5] .= 0.201 #y_d[:,1]*penalty #NO2 to N2O
y_d[:,6] .= 0.197 #y_d[:,1]*penalty #NO2 to N2
y_d[:,7] .= 0.298 #y_d[:,1]*penalty #N2O to N2

#Anaerobic heterotroph yields on DIN (Xin updated 09/2023, P=0.16) 
# N2O and N2 are in the unit of N
y_n[3,2] = 0.011 #0.039 #NO3 to NO2 
y_n[3,3] = 0.023 #0.078 #NO3 to N2O 
y_n[3,4] = 0.026 #0.098 #NO3 to N2 
y_n[2,5] = 0.016 #0.039 #NO2 to N2O 
y_n[2,6] = 0.024 #0.059 #NO2 to N2
y_n[4,7] = 0.013 #0.020 #N2O to N2

#AOA and NOB yields (Xin updated 09/2023):
y_n[1, nhets + 1] = 0.0245 # AOA: #mol B/mol NH4 (Bayer et al. 2022; Zakem et al. 2022)
y_n[2, nhets + 2] = 0.0126 # NOB: mol B/mol NO2 (Bayer et al. 2022)

#Anammox yields and excretion factor (Xin updated 09/2023)
y_n[1, nhets + 3] = 1. ./ 75 #mol B/mol NH4 (Lotti et al. 2014 Water Research)
y_n[2, nhets + 3] = 1. ./ 89 #mol B/mol NO2 (Lotti et al. 2014 Water Research)
ftoNO3_anx = 0.08 #percent of excreted DIN that goes to NO3 instead of N2 (Lotti et al. 2014 Water Research)

#Oxygen yields for all aerobes (Xin updated 09/2023)
y_o[1] = 0.035 #0.091 
y_o[nhets + 1] = 0.018 # nhets = 7, # 8 is AOA
y_o[nhets + 2] = 0.022 # nhets = 7, # 9 is NOB

######################################################################################
#UPTAKE KINETICS 

#Organic matter (all the same right now):
vmax_d = zeros(nd,nb)
K_d = zeros(nd,nb)
#slow and fast pools, but all hets eat both
vmax_d[1,1:nhets] .= 1 # set the max OM1 uptake rate for the same for all hets
K_d[1,1:nhets] .= 0.1 # half saturation constant of OM1 the same for all hets

#vmax_d[2,1:nhets] .= 0.25 # same max OM2 uptake rate for all hets
#K_d[2,1:nhets] .= 0.5 # same OM2 half saturation constant

#DIN (just allow all right now):
vmax_n = ones(nn,nb)*50.8 # Xin: we used the same vmax (50.8) in chemostat model
K_n = ones(nn,nb)*0.133 
# nn = 5 # NH4, NO2, NO3, N2O, N2 #number of inorganic matter pools
# nb = 10 
    #1. Aerobic
    #2. NO3 to NO2
    #3. NO3 to N2O
    #4. NO3 to N2 
    #5. NO2 to N2O
    #6. NO2 to N2
    #7. N2O to N2
    #8. AOO #9. NOB #10. anx 
if diffks_n == 1
    K_n[2:3, 2:7] .= 4.0  # NO3 or NO2 uptake by denitrifiers based on 4 – 25 µM NO2 for denitrifiers (Almeida et al. 1995)
    K_n[4, 7] = 0.3*2    # N2O uptake by denitrifier (*2 convert µM N2O to µM N-N2O) based on (Sun et al., 2021) in ODZ k = 0.3 µM N2O, in oxic layer = 1.4~2.8 µM
    K_n[1, 8] = 0.1      # NH4 uptake by AOO based on Martens-Habbena et al. 2009 Nature
    K_n[2, 9] = 0.1      # NO2 uptake by NOB reported from OMZ (Sun et al. 2017) and oligotrophic conditions (Zhang et al. 2020)
    K_n[1:2, 10] .= 0.45  # NH4 and NO2 uptake by anx based on Awata et al. 2013 for Scalindua (K_no2 of 3.0 uM, but this excludes anammox completely in our 0-D)
end

#O2: # O2 is taken up by diffusing int the cell here. I can change it 
pcoef = zeros(nb) #the factor that we multiply [O2] with to get the diffusion-limited uptake rate of O2
#po_coef=985400/2*1D0, & !m3/mol/day - coefficient for diffusive uptake using D=1.5 cm^2/s and diameter = 1um and Litchman's QminN. see October2014_OvsNO3.m with divided by 2 following Stolper 2010
pcoef[:] .= 1000 #above: order 1,000,000 in m3/mol/d, so that x 1e-3 mol/mmol = 1000

## change k_o2 for below when model O2 not as diffusing into cells but M-M
#K_o2_aer = 0.2   # (µM-O2) 10 to 200 nM at extremely low oxygen concentrations (Tiano et al., 2014); 200 nM (low affinity terminal oxidases) (Morris et al., 2013)
#K_o2_aoo = 0.333  # (µM-O2) 333 nM at oxic-anoxic interface (Bristow et al., 2016)
#K_o2_noo = 0.8   # (µM-O2) 778 nM at oxic-anoxic interface or 0.5 nM and 1750 nM breaking into two M-M curves (Bristow et al., 2016)
