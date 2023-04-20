#1/23/22
#input_OMdetail.jl for input into run_addDIN.jl

######################################################################################
# Organic matter detail

# distribution of m to d (assign dead cell to the number of OM pools I have. For now, equal assign to the 2 pools)
prob_generate_d = zeros(nd) #need to make it a vector 
#prob_generate_d = prob_generate_d / sum(prob_generate_d)
#prob_generate_d = transpose(prob_generate_d)
#Lognormal distribution for regularly spaced vmax
#x = log.(vmax_i)
#μ = 0 #log mean
#σ = 1 #log stdev
#pdfx = 1 ./(sqrt(2*pi*σ^2)).*exp.(-(x .- μ).^2 ./ (2*σ^2) )
pdfx = ones(nd) #uniform distribution
prob_generate_d[:] .= pdfx/sum(pdfx)


# Sinking rate for POM (particulate organic matter)
ws = 10*ones(nd) # !m/day --sinking speed of POM (Particulate organic matter) 
# ths is how to define one OM to dissolved:

# Xin_let 90% of OM be DOM, 10% be POM
if DOM1  == 1 #DOM1_90perc
    #prob_generate_d[1] = 0.9
    #prob_generate_d[2] = 0.1
    ws[1] = 0 #, the 1st OM not sinking (Dissolve)
end

#wd=ws .+ w # sinking velocity + water velocity gives vector of w for POM
wd = transpose(repeat(ws, 1, ngrid + 1)) + repeat(w, 1, nd) #ngrid+1 x nd  
wd[1,:].=0 #no flux boundary at surface (for running reason, nothing to sink to the first box)
wd[end,:].=0 # no flux boundary: bottom box accumulates D (allow nothing to sink out of the last box)
