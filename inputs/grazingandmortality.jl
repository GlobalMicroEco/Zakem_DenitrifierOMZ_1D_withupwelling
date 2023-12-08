#1/23/22
#input_grazingandmortality.jl for input into run_addDIN.jl

######################################################################################
# Grazing
# grazing rate of z
g_max = ones(nz)
# the half saturation rate of z
K_g = ones(nz)
# the fraction of digestion
γ = ones(nz)*0.5

# Graze matrix
#GrM = ones(nz, np+nb)
# Graze matrix
GrM = ones(Bool, (nz, np+nb))
GrM[:,1:np] .= 0 #no pp for one grazer # turn off grazer for phytoplankton
# ! change: we could turn off grazers for bacteria too to test the grazing conditions.

if nz == 2
    GrM[1,np+1:end] .= 0 #1st z is only for p
    GrM[2,1:np] .= 0 #2nd z is only for b
    #GrM[2,1:np+1] .= 0 #1st b, the PA, is not grazed
end

######################################################################################
# Mortality
m_lp = ones(np) * 0.01  # linear mort: m3/mmol/d 
m_qp = ones(np) * 1.#1e-1 # quadratic mort # give a higher rate because phytoplankton has no grazer in the code above.

m_lb = ones(nb) * 0.01 #linear mort: 1/d
#less for slower-growing chemoautotrophs:
#m_lb[nhets+1:end] .= 0.01
m_qb = ones(nb) * 0.1 #quadratic mort: m3/mmol/d

# ! change--zooplankton mortality, important! because it decided the whole system's magnitude (relative relationship should not be affected)
m_lz = ones(nz) * 0.01 #1e-2
m_qz = ones(nz) * 0.7
