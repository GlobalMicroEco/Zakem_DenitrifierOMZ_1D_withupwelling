using NCDatasets
#using CairoMakie
using Plots, ColorSchemes
#gr() #Plots.GRBackend() for gradients in line plots below
#using Printf
#using PlotlyJS, 
using DataFrames
#test
default(show = true)

#read nc data
folder = "/Users/xinsun/Dropbox/!Carnegie/!Research/Xin_Zakem_DenitrifierOMZ_1D_withupwelling/"
result_file = "out_addW_20230908.nc"
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

#plot profiles
#for 1 for everyone:
#plot([sum(p, dims = 2), sum(b, dims = 2), sum(z, dims = 2), sum(d, dims = 2), n], -zc; 
 #    layout = 5, 
  #   linewidth = 2, 
   #  legend = false, 
    # title = ["Total P" "Total B" "Total Z" "Total D" "Total N"]
#)

### plot stacked biomasses and nuts:
gr(size = (1000, 600))
## start creating subplots
pPZ = plot(sum(p, dims = 2), -zc, linecolor = "green", label = "Total P",  xlabel = "µM-BiomassN")
plot!(sum(z, dims = 2), -zc, linecolor = "black", label = "Total Z") 
plot!(sum(b[:,2:7],dims = 2), -zc, linestyle =:dash, linecolor = "grey", label = "Total denit")  # order of microbes: o, n1 to n6, aoa, nob, and anx 
#plot!(sum(b[:,1:nhets], dims = 2), -zc, linecolor = "grey", label = "Total B")  # order of microbes: o, n1 to n6, aoa, nob, and anx 
plot!(b[:, 1], -zc, linestyle =:dash, linecolor = "blue", label = "AerHet")

pOM = plot(d, -zc, palette = :lightrainbow, label = ["OM1" "OM2"], xlabel = "µM")

#pDenit = plot(b[:,2:nhets]*1e10, -zc, palette = :tab10, label = ["NO3->NO2" "NO3->N2O" "NO3->N2" "NO2->N2O" "NO2->N2" "N2O->N2"], xlabel = "(1e-10)µM-BiomassN", legend=:right) #xlims = (0, 2)
pDenit = plot(b[:,2:nhets]*1e3, -zc, palette = :tab10, linestyle = [:solid :dash :dot :solid :dash :solid], label = ["NO3->NO2" "NO3->N2O" "NO3->N2" "NO2->N2O" "NO2->N2" "N2O->N2"], xlabel = "nM-BiomassN", legend=:bottomright)

#pAerobe = plot(b[:, [8,9,10]], -zc, palette = :lightrainbow, label = ["AOO" "NOB" "Anx"], xlabel = "µM-BiomassN", legend=:left)
#plot!(sum(b[:,2:7],dims = 2), -zc, linecolor = "grey", label = "Total denit")  # order of microbes: o, n1 to n6, aoa, nob, and anx 
pAerobe = plot(b[:, [8,9,10]]*1e3, -zc, palette = :lightrainbow, label = ["AOO" "NOB" "Anx"], xlabel = "nM-BiomassN", legend=:bottomright) #, xlims = (0, 5)
# palette = :phase OR :lightrainbow OR :rainbow

# DIN order: NH4, NO2, NO3, N2O, N2
pOxy = plot(o, -zc, linecolor = "black", label = "O2", legend=:bottomright, xlabel = "µM") # OR: label = false


pDIN = plot(n[:,2]*10,-zc, linecolor = "green", label = "NO2-*10", legend=:bottomright, xlabel = "µM", xlims = (-0.5, 30)) # OR: label = false
plot!(n[:,3],-zc, linecolor = "grey", label = "NO3-")
plot!(n[:,1]*1000, -zc, linecolor = "blue", label = "NH4+*1000")

# plot!(n[:,2]./n[:,1],-zc, linecolor = "black", label = "NO2-/NH4+", legend=:bottomleft) 

p5 = plot(n[:,5],-zc, linecolor = "grey", label = "N2", legend=:bottomright, xlabel = "µM") #, xlims = (-1, 10)

pN2O = plot(n[:,4]*1000/2, -zc, linecolor = "red", label = "N2O", legend=:bottomright, xlabel = "nM-N2O", xlims = (-0.5, 100))
ylims=(-1500,0)

plot(pOM, pPZ, pAerobe, pDenit, pOxy, pDIN, p5, pN2O,
    linewidth = 2,
    layout = grid(2,4),
    #ylims=(-1500,0)
    #title = ["Biomass" "D" "O2+N" "N" "N2O"]
)
## save plot, change file name!
basename = replace(result_file, ".nc" => "")
basename = replace(basename, "out_addW_" => "")
# change this custom text for figure name
custom_text = "_NairseaN2Oamm_diffK_60y=run30y_from0907_1051" #"_NairseaN2Oamm_diffK_20y=run10y_from09061504"
figure_name = "Fig_$basename$custom_text"
# save
png("$figure_name.png") 

#plot rates:
#npp = mu_p .* p
#resp = mu_b .* b .* (1/y - 1)  

# did not work: savefig(Results1Dplot.png)
