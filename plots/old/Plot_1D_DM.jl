using NCDatasets
using Plots, ColorSchemes
using DataFrames
default(show = false)

# read nc data
result_file = "/Users/DcCoy/Documents/GitHub/MicrOMZ/out/out_addW_20231207.nc"
println(result_file)
ds = NCDataset(result_file)

# extract final timestep
p = ds["p"][:,:,end]
b = ds["b"][:,:,end]
z = ds["z"][:,:,end]
n = ds["n"][:,:,end]
d = ds["d"][:,:,end]
o = ds["o"][:,end]

# extract totals
ball = ds["b"][:]
dall = ds["d"][:]
oall = ds["o"][:]
nall = ds["n"][:]

# extract size of arrays
np = size(p)[2]
nb = size(b)[2]
nz = size(z)[2]
nn = size(n)[2]
nd = size(d)[2]

# set number of heterotrophs
nhets = 7

# Load grid data
wd  = ds["wd"][:]        # sinking rate
H   = first(ds["H"][:])  # water column depth
dz  = first(ds["dz"][:]) # grid spacing
zc  = collect((dz/2):dz:(H-(dz/2)))

#plot C profiles
l = @layout [grid(1,4)]
p1 = plot(sum(p,dims=2),-zc,
     linewidth = 3, 
     legend = false, 
     ylimits=(-2000,0),
     title = "Final phyto biomass"
    )
p2 = plot(sum(z,dims=2),-zc,
     linewidth = 3, 
     legend = false, 
     ylimits=(-2000,0),
     ytickfontcolor=:white,
     title = "Final zoo biomass"
    )
p3 = plot(sum(d,dims=2),-zc,
     linewidth = 3, 
     legend = false, 
     ylimits=(-2000,0),
     ytickfontcolor=:white,
     title = "Final orgC"
    )
p4 = plot(sum(b,dims=2),-zc,
     linewidth = 3, 
     legend = false, 
     ylimits=(-2000,0),
     ytickfontcolor=:white,
     title = "Final total C"
    )
plot(p1, p2, p3, p4, layout = l)
savefig("Results_final_C.png")

# plot N profiles
# [NH4, NO2, NO3, N2O, N2]
l = @layout [grid(1,5)]
p1 = plot(n[:,1],-zc,
     linewidth = 3, 
     legend = false, 
     ylimits=(-2000,0),
     title = "Final NH4"
    )
p2 = plot(n[:,2],-zc,
     linewidth = 3, 
     legend = false, 
     ylimits=(-2000,0),
     ytickfontcolor=:white,
     title = "Final NO2"
    )
p3 = plot(n[:,3],-zc,
     linewidth = 3, 
     legend = false, 
     ylimits=(-2000,0),
     ytickfontcolor=:white,
     title = "Final NO3"
    )
p4 = plot(n[:,4],-zc,
     linewidth = 3, 
     legend = false, 
     ylimits=(-2000,0),
     ytickfontcolor=:white,
     title = "Final N2O"
    )
p5 = plot(n[:,5],-zc,
     linewidth = 3, 
     legend = false, 
     ylimits=(-2000,0),
     ytickfontcolor=:white,
     title = "Final N2"
    )
plot(p1, p2, p3, p4, p5, layout = l)
savefig("Results_final_N.png")

fail
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
custom_text = "_diffK_run10yfromNewIC" #"_NairseaN2Oamm_diffK_40y=run10y_from0907_1051"
figure_name = "Fig_$basename$custom_text"
# save
png("$figure_name.png") 

#plot rates:
#npp = mu_p .* p
#resp = mu_b .* b .* (1/y - 1)  

