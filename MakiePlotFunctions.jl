using CairoMakie, ColorSchemes 
using Printf

println("Loading plot_depth function")
function plot_depth(conc, tslice, figsave, title)
    nz, nc, nt = size(conc)

    if nc == 1
        cmap = cgrad(:roma)
    else
        cmap = cgrad(:roma, nc, categorical = true)
    end

    if !ismissing(title)
        pc = Figure(resolution = (1000, 800))
        pc1 = Axis(pc[1,1], xlabel = "Concentration", ylabel = "Depth", title = title)

        labels = []
        for ll in 1:nc
            push!(labels, string(ll))
        end 

        for j in 1:nc
            scatterlines!(pc1, conc[:,j,tslice], 0 .- [1:1:nz;], color = cmap[j], linewidth=1, label = labels[j])
        end
    else 
        pc = Figure(resolution = (1000, 800))
        pc1 = Axis(pc[1,1], xlabel = "Concentration", ylabel = "Depth")
        for j in 1:nc
            scatterlines!(pc1, conc[:,j,tslice], 0 .- [1:1:nz;], color = cmap[j], linewidth=1)
        end
    end

    yticks = collect(0:-100:-400)
    yticklabels = [ @sprintf("%5.1f", x * 5) for x in yticks ]

    pc[1,2] = Legend(pc, pc1, "Lines", framevisible = false)

    if !ismissing(figsave)
        fs = string(figsave, ".png")
        # save the figure
        save(fs, pc)
    end

    return pc

end



println("Loading plot_time function")
function plot_time(conc, depth, figsave, title)
    nz, nc, nt = size(conc)

    if nc == 1
        cmap = cgrad(:roma)
    else
        cmap = cgrad(:roma, nc, categorical = true)
    end
    
    if !ismissing(title)
        pc = Figure(resolution = (1000, 800))
        pc1 = Axis(pc[1,1], xlabel = "Time", ylabel = "Concentration", title = title)

        labels = []
        for ll in 1:nc
            push!(labels, string(ll))
        end 

        for j in 1:nc
            scatterlines!(pc1, 1:1:nt, conc[depth,j,:], color = cmap[j], linewidth=2, label = labels[j])
        end
    else 
        pc = Figure(resolution = (1000, 800))
        pc1 = Axis(pc[1,1], xlabel = "Time", ylabel = "Concentration", title = title)
        for j in 1:nc
            scatterlines!(pc1, 1:1:nt, conc[depth,j,:], color = cmap[j], linewidth=2)
        end
    end

    pc[1,2] = Legend(pc, pc1, "Lines", framevisible = false)

    if !ismissing(figsave)
        fs = string(figsave, ".png")
        # save the figure
        save(fs, pc)
    end

    return pc

end


