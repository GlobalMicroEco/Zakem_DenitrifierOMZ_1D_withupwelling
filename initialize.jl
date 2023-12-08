# 12/07/23: version updated by Danny McCoy
# 02/02/23: version to send to Xin Sun

# Announce model start
println("Starting DM model version (adapted from EJZ)")

# Set subnormals to zero to avoid slowdown
# Shis way, don't have to cap b and d at a lower limit
set_zero_subnormals(true)

# Make a Composite Type using "struct" to organize parameters (see settings.jl)
# Int64   = integer
# Float64 = floating-point number, e.g., 3.14 
struct paras
    tt::Int64                          # days -- run time
    nrec::Float64                      # the number of timesetps to output
    H::Int64                           # the total depth of the sea
    dz::Int64                          # the height of the box, a parcel of water
    np::Int64                          # the number of the phytoplankton
    nb::Int64                          # the number of the bacteria
    nhets::Int64                       # the number of heterotrophic bacteria 
    nz::Int64                          # the number of the zooplankton
    nn::Int64                          # the number of the inorganic matter pools
    nd::Int64                          # the number of the organic matter pools
    pIC::Array{Float64,2}              # the initial condition for the np p for all boxes at time 0
    bIC::Array{Float64,2}              # the initial condition for the nb b for all boxes at time 0
    zIC::Array{Float64,2}              # the initial condition for the nz z for all boxes at time 0
    nIC::Array{Float64,2}              # the initial condition for the nn n for all boxes at time 0
    dIC::Array{Float64,2}              # the initial condition for the nd d for all boxes at time 0
    oIC::Array{Float64,2}              # the initial condition for the nd d for all boxes at time 0
    umax_p::Array{Float64,1}           # the max growth rate of p 
    K_pn::Array{Float64,1}             # the half saturating rate for each p feeding on n (nutrient)
    m_lp::Array{Float64,1}             # the linear mort of p
    m_qp::Array{Float64,1}             # the quadratic mort of p
    CM::Array{Float64,2}               # Consumption matrix: nd X nb
    y_d::Array{Float64,2}              # the yield rate of bacteria j on d i
    y_n::Array{Float64,2}              # the yield rate of bacteria j on n i
    y_o::Array{Float64,1}              # the yield rate of bacteria j on oxygen
    ftoNO3_anx::Float64                # fraction of anammox excretion that goes to N2 instead of NO3
    vmax_d::Array{Float64,2}           # max uptake rate of bacteria j on d i
    K_d::Array{Float64,2}              # half sat of bacteria j on d i
    vmax_n::Array{Float64,2}           # max uptake rate of bacteria j on n i
    K_n::Array{Float64,2}              # half sat of bacteria j on n i
    pcoef::Array{Float64,1}            # max uptake rate of bacteria j on oxygen
    m_lb::Array{Float64,1}             # the linear mort of b
    m_qb::Array{Float64,1}             # the quadratic mort of b
    g_max::Array{Float64,1}            # the max grazing rate of z
    K_g::Array{Float64,1}              # the half saturation rate of z
    γ::Array{Float64,1}                # the fraction of assimilation for z 
    m_lz::Array{Float64,1}             # the linear mort of z
    m_qz::Array{Float64,1}             # the quadratic mort of z
    prob_generate_d::Array{Float64, 1} # the distribution from the total gain of d to each d pool from mortality
    kappa_z::Array{Float64, 1}         # the vertical eddy diffusivities
    w::Array{Float64, 1}               # Vertical velocity of water
    wd::Array{Float64, 2}              # Sinking rate for POM, INCLUDING vertical velocity already (wd = w + ws) 
    fsave::String                      # saving file name
    Light::Array{Float64, 1}           # Light (irradiance) W/m2
    TempFun::Array{Float64, 1}         # Temperature fun
    K_I::Array{Float64, 1}             # Half-sat rate of light
    number_box::Int64                  # number of boxes at depth
    GrM::Array{Float64, 2}             # Grazing matrix: nz x np+nb
    e_o::Float64
    koverh::Float64
    o2sat::Float64
    t_o2relax::Float64
    o2_deep::Float64
end

paras1 = paras(tt, nrec, H, dz, np, nb, nhets, nz, nn, nd, pIC, bIC, zIC, nIC, dIC, oIC,
            umax_p, K_pn, m_lp, m_qp, CM, y_d, y_n, y_o, ftoNO3_anx, 
            vmax_d, K_d, vmax_n, K_n, pcoef, m_lb, m_qb,
            g_max, K_g, γ, m_lz, m_qz, prob_generate_d, kappa_z, w, wd, fsave, Light, TempFun,
            K_I, ngrid, GrM,
            e_o, koverh, o2sat, t_o2relax, o2_deep)
