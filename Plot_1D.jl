using NCDatasets
#using CairoMakie
using Plots, ColorSchemes
#gr() #Plots.GRBackend() for gradients in line plots below
#using Printf
#using PlotlyJS, 
using DataFrames

default(show = true)

#read nc data
folder = ""

#just added upwelling -- first test:
result_file = "out_addDIN_20230320_1256.nc" #3 yrs with kappa to 1e-4 from 5e-5

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
plot!(sum(z, dims = 2), -zc, linecolor = "grey", label = "Total Z") 
plot!(sum(b[:,1:nhets], dims = 2), -zc, linecolor = "black", label = "Total B")  # order of microbes: o, n1 to n6, aoa, nob, and anx 
plot!(b[:, 1], -zc, linecolor = "blue", label = "AerHet")

pOM = plot(d, -zc, palette = :lightrainbow, label = ["OM1" "OM2"], xlabel = "µM")

pDenit = plot(b[:,2:nhets]*1e3, -zc, palette = :tab10, label = ["NO3->NO2" "NO3->N2O" "NO3->N2" "NO2->N2O" "NO2->N2" "N2O->N2"], xlabel = "nM-BiomassN", legend=:right)

pAerobe = plot(b[:, [8,9,10]]*1e3, -zc, palette = :lightrainbow, label = ["AOO" "NOB" "Anx"], xlabel = "nM-BiomassN", legend=:right) #, xlims = (0, 5)
plot!(sum(b[:,2:7],dims = 2)*1e3, -zc, linecolor = "grey", label = "Total denit")  # order of microbes: o, n1 to n6, aoa, nob, and anx 

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
