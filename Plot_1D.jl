using NCDatasets
#using CairoMakie
#using Plots, ColorSchemes
using Plots
#gr() #Plots.GRBackend() for gradients in line plots below
#using Printf
#using PlotlyJS, 
using DataFrames

default(show = true)

#read nc data
folder = ""
#folder = "/Users/xinsun/Dropbox/!Carnegie/!Research/!Project_xin_N2O/Emily_1D_Julia/toXin_020223/"
#result_file = "out_addDIN_20230125_1734.nc" #as above but with o2deep = 10. #good! still very small OMZ

folder = "../Sun_DenitrifierOMZ_1D"
result_file = "out_addDIN_20230302.nc" #1 d -- looks promising! got to keep OM
result_file = "out_addDIN_20230302_0736.nc" #1 year -- still spinning up but looks promising! 
result_file = "out_addDIN_20230302_0806.nc" #3 more years run FROM above 1 yr. goes away.
result_file = "out_addDIN_20230302_1356.nc" #3 years total with km = 0.5 for slower OM pool. bad. 
result_file = "out_addDIN_20230302_1850.nc" #1 d for 10 d
result_file = "out_addDIN_20230306.nc" #1 d for 3 yrs -- new SS
#result_file = "out_addDIN_20230306_1543.nc" #1 d for 3 yrs -- new SS
result_file = "out_addDIN_20230313.nc" #kd to 0.1 from 0.5, got rid of mlboxes in params structure
#result_file = "out_addDIN_20230313_1349.nc" #
result_file = "out_addDIN_20230313_1358.nc" #3yrs: higher prod but still no anoxic zone!
##result_file = "out_addDIN_20230313_1527.nc" #3yrs: vmax to 0.5 from 1 for aerobic hets (like 2D model). pile up on bottom, lower pp
#result_file = "out_addDIN_20230313_1529.nc" #3yrs: ws to 3 from 10, and vmax to 0.5 from 1 for aerobic hets (like 2D model). seems good.
result_file = "out_addDIN_20230313_1534.nc" #3yrs: mlz to 20, and ws to 3 from 10, and vmax to 0.5 from 1 for aerobic hets (like 2D model). seems good, not too diff than mlz 25.
##result_file = "out_addDIN_20230313_1557.nc" #3yrs: vmax at 1, mlz to 20, and ws to 3 from 10. NOT good: nitrifiers at surface!! 
#
#result_file = "out_addDIN_20230313_1617.nc" #3yrs: kzmin to 5e-5, vmax at 0.5, mlz to 20, and ws to 3 from 10 (slow model with kzmin and mlz). o2 0 at 500m (instead of 600)
#result_file = "out_addDIN_20230313_1623.nc" #3yrs: vmax at 1, ws = 10, kzmin to 5e-5, mlz to 20 (fast model with kzmin and mlz). same. o2 0 at 500m
result_file = "out_addDIN_20230313_1650.nc" #3yrs: same as fast model but resetting O2 to 1uM everywhere (is it just the IC for O2 that keeps it high below?). o2 to 300m! oh, no zoo o2 cons! 
#result_file = "out_addDIN_20230313_1656.nc" #3yrs: same as slow model but resetting O2 to 1uM everywhere (is it just the IC for O2 that keeps it high below?). also o2 to 300
result_file = "out_addDIN_20230313_1726.nc" #3yrs: including Zoo O2 consumption! from slow model with O2 reset. o2 still to 300m. 
#
result_file = "out_addDIN_20230313_2122.nc" #3yrs: including Zoo O2 consumption! from fast model with O2 reset. again same as slow model (keep fast since more realistic params? or at least 1 and 10 so "order of mag") 


#Adding upwelling!
folder = "../Sun_DenitrifierOMZ_1D/AddUpwelling/"
result_file = "out_addDIN_20230320_1241.nc" #1 yr to test: looks OK! 
result_file = "out_addDIN_20230320_1256.nc" #3 yrs with kappa to 1e-4 from 5e-5
#revisiting in Aug 2023 (following upload to github)
#result_file = "out_addW_20230830.nc" 
result_file = "out_addW_20230830_1832.nc" #running for 1 year, mlz is bigger (25m), relax only for bottom box. smaller anoxic zone than above
#result_file = "out_addW_20230830_1835.nc" #running for 1 year, mlz is bigger (25m), relax for all boxes (changed model file). NO anoxic zone!
result_file = "out_addW_20230830_1903.nc" #running for 1 year, back to mlz=20, kappa 5e-5, relax only for bottom box. much better! significant denitrifiers!
#result_file = "out_addW_20230830_1904.nc" #running for 1 year, back to mlz=20, kappa 5e-5, relax for all boxes (changed model file). nope -- no anoxic zone! 
#result_file = "out_addW_20230830_1922.nc" #running 3 years FROM the good 1 yr run (out_addW_20230830_1903.nc) to see how things evolve

#Change to working from github repository: Zakem_DenitrifierOMZ_1D_withupwelling (this folder)
folder = ""
result_file = "out_addW_20230830.nc"
result_file = "out_addW_20230830_1922.nc" #3y on top of 1y spinup/IC. Woo! Coming along! No NO2 accum yet, as expected
result_file = "out_addW_20230830_2032.nc" #adding another 6 years for a 10y run total. cool! N2O peak at top and bottom! No NO2 accum.
result_file = "out_addW_20230831_1143.nc" #oxygenated


println(result_file)

ds = NCDataset(folder*result_file)

#just end point
p = ds["p"][:,:,end]
b = ds["b"][:,:,end]
z = ds["z"][:,:,end]
n = ds["n"][:,:,end]
d = ds["d"][:,:,end]
o = ds["o"][:,end]

#all
ball = ds["b"][:]
dall = ds["d"][:]
oall = ds["o"][:]
nall = ds["n"][:]

np = size(p)[2]
nb = size(b)[2]
nz = size(z)[2]
nn = size(n)[2]
nd = size(d)[2]

nhets = 7

#Traits
#vmax_i = ds["vmax_i"][:]
#vmax_ij = ds["vmax_ij"][:]
#K_ij = ds["K_ij"][:]
#aff = vmax_ij./K_ij
#CM = ds["CM"][:]

wd = ds["wd"][:]

H = ds["H"][:]
dz = ds["dz"][:]
zc = [dz/2:dz:(H-dz/2)]


### plot stacked biomasses and nuts:
gr(size = (1000, 600))
## start creating subplots
pPZ = plot(sum(p, dims = 2), -zc, linecolor = "green", label = "Total P",  xlabel = "µM-BiomassN")
plot!(sum(z, dims = 2), -zc, linecolor = "Black", label = "Total Z") 
plot!(b[:, 1], -zc, linecolor = "blue", label = "AerHet")
plot!(sum(b[:,2:7],dims = 2), -zc, linecolor = "grey", label = "Total denit hets")  # order of microbes: o, n1 to n6, aoa, nob, and anx 

pOM = plot(d, -zc, palette = :lightrainbow, label = ["OM1" "OM2"], xlabel = "µM")

pDenit = plot(b[:,2:nhets]*1e3, -zc, palette = :tab10, label = ["NO3->NO2" "NO3->N2O" "NO3->N2" "NO2->N2O" "NO2->N2" "N2O->N2"], xlabel = "nM-BiomassN", legend=:right)

pAerobe = plot(b[:, [8,9,10]]*1e3, -zc, palette = :lightrainbow, label = ["AOO" "NOB" "Anx"], xlabel = "nM-BiomassN", legend=:right) #, xlims = (0, 5)
plot!(sum(b[:,2:7],dims = 2)*1e3, -zc, linecolor = "grey", label = "Total denit hets")  # order of microbes: o, n1 to n6, aoa, nob, and anx 

# DIN order: NH4, NO2, NO3, N2O, N2
p3 = plot(o, -zc, linecolor = "black", label = "O2", legend=:bottomleft, xlabel = "µM") # OR: label = false
plot!(n[:,3],-zc, linecolor = "grey", label = "NO3-")

p4 = plot(n[:,2],-zc, linecolor = "green", label = "NO2-", legend=:bottomleft, xlabel = "µM", xlims=(0,0.1)) # OR: label = false
plot!(n[:,1],-zc, linecolor = "blue", label = "NH4+", legend=:bottomleft, xlabel = "µM") # OR: label = false

p5 = plot(n[:,1]*1000, -zc, linecolor = "blue", label = "NH4+", legend=:bottomleft, xlabel = "nM")
plot!(n[:,5]*1000,-zc, linecolor = "grey", label = "N2")

p6 = plot(n[:,4]*1000, -zc, linecolor = "red", label = "N2O", legend=:bottomleft, xlabel = "nM")

plot(pOM, pPZ, pAerobe, pDenit, p3, p4, p5, p6,
    linewidth = 2,
    layout = grid(2,4),
    ylims=(-1000,0)
    #title = ["Biomass" "D" "O2+N" "N" "N2O"]
)
## save plot, change file name!
#png("Fig_20230209_2309_10yrs_sameK_n")
#png("Fig_20230209_2126_10yrs_diffK_n")
#png("Fig_20230209_2307_1day_sameK_n")
#png("Fig_20230214_10yrs_diffK_n")
#png("Fig_20230214_2024_10yrs_diffK_n_denitICsmall")
#png("Fig_20230218_0900_10yrs_diffK_n_denitICsmall_DOM1")
#png("Fig_20230218_1344_10yrs_diffK_n_denitICsmall_90%DOM1")
#png("Fig_20230218_1623_10yrs_diffK_n_denitICsmall_DOM1_O2relax0.001")
#png("Fig_20230223_1544_10yrs_diffK_n_O2relaxDeepOnly")
#png("Fig_20230223_10yrs_diffK_n_DetaGYom_O2relaxDeepOnly")

#plot rates:
#npp = mu_p .* p
#resp = mu_b .* b .* (1/y - 1)  

# did not work: savefig(Results1Dplot.png)
