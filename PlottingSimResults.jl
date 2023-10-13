using NCDatasets
using CairoMakie, ColorSchemes
using DataFrames

include("MakiePlotFunctions.jl")

#read nc data
result_file = "out_addW_test_newmodel1_0D_20231013_0902.nc" #"out_addW_20230908.nc"out_addW_test_newmodel1_0D_20231013_0902.nc
println(result_file)

ds = NCDataset(result_file)

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
zall = ds["z"][:]
pall = ds["p"][:]

np = size(p)[2]
nb = size(b)[2]
nz = size(z)[2]
nn = size(n)[2]
nd = size(d)[2]


#plot time series given the depth of the water box 
depth =7
conc = nall[:,[1,2],:]

plot_time(conc, depth, "Tplotat7", "Inorganic matter time series at depth 7")


#plot depth profiles given the time slice
tslice = 155
conc = nall[:,[1],:]

plot_depth(conc, tslice, "Vplotat155", "Inorganic matter vertical profile at time 155")


# plot w against time
# w = ds["w"][:]
# scatter(1:length(w), w, title = "w against time")
