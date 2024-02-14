###################################################################
# Main model timestepping code 
###################################################################
#
# Imports parameter values
# Initializes output matrixes
# Integrated model forward in time
#
println("Loading function: run_npzdb")
function run_npzdb(params)
    # Import parameter values
    println("Importing parameters")
    @printf("np = %5.0f, nb = %5.0f, nz = %5.0f, nn = %5.0f, nd = %5.0f, T = %5.0f \n", 
             params.np,  params.nb,  params.nz,  params.nn,  params.nd,  params.tt)
    
    # Initialize log file
    open("out/jlog.txt","w") do f
        write(f,@sprintf("np = %5.0f, nb = %5.0f, nz = %5.0f, nn = %5.0f, nd = %5.0f, T = %5.0f \n", 
                          params.np,  params.nb,  params.nz,  params.nn,  params.nd,  params.tt))
    end
    tst = now()

    # Start output file
    if params.fsave == "missing"
        # Do nothing
    else
        # Starting time and saving location
        fsaven = string(params.fsave,"_", Dates.format(tst, "yyyymmdd"), ".nc")
        if isfile(fsaven) #if i do more than 1 in 1 day
            fsaven = string(params.fsave,"_",Dates.format(tst, "yyyymmdd_HHMM"), ".nc")
        end
        println("Will be saved as: ", fsaven)
    end

    # Announce starting time
    println("Starting time: ", tst)
    println("Thread id: ", Threads.threadid() )

    # Prep model for integration
    Cs = sparse(params.CM) # sparse matrix (most 0 in the matrix)
    (II, JJ, _) = findnz(Cs) #this gives the list of all "on" indices in the matrix
    println("CM[1:10,1:10]")
    for i=1:min(params.nd,10)
        println(params.CM[i,1:min(params.nb,10)])
    end

    # Get number of timesteps and output frequency
    nt = Int(params.tt/dt)
    trec = nt/params.nrec 

    # Initialize matrices to track tendencies over time
    # This allows tracking the dyanmics of p, b, z, n, d
    # 3D arrays, where i = numDepths, j = num of P/B/Z/N variants, and k = num timesteps)
    nrec1 = Int(params.nrec + 1)
    track_p = Array{Float64, 3}(undef, params.number_box, params.np, nrec1)  
    track_b = Array{Float64, 3}(undef, params.number_box, params.nb, nrec1) 
    track_z = Array{Float64, 3}(undef, params.number_box, params.nz, nrec1)   
    track_n = Array{Float64, 3}(undef, params.number_box, params.nn, nrec1)  
    track_d = Array{Float64, 3}(undef, params.number_box, params.nd, nrec1)  
    track_o = Array{Float64, 3}(undef, params.number_box, 1, nrec1)          
    track_time = Array{Float64,1}(undef, nrec1)   

    # Set initial conditions (hence why nrec1 = ncrec + 1)
    track_p[:,:,1] .= params.pIC
    track_b[:,:,1] .= params.bIC
    track_z[:,:,1] .= params.zIC
    track_n[:,:,1] .= params.nIC
    track_d[:,:,1] .= params.dIC
    track_o[:,:,1] .= params.oIC
    track_time[1] = 0

    # Initialize the variables, change with each time step
    track_p_temp = copy(params.pIC) 
    track_b_temp = copy(params.bIC)
    track_z_temp = copy(params.zIC) 
    track_n_temp = copy(params.nIC) 
    track_d_temp = copy(params.dIC) 
    track_o_temp = copy(params.oIC) 

    # Begin timestepping
    println("Start timestepping")
    @time for t = 1:nt

       # Runge-Kutta 4th order, a way of solving differential equations, numerical solver
       # First calculation
       track_dpdt1, track_dbdt1, track_dzdt1, track_dndt1, track_dddt1, track_dodt1 = one_step_ecof(
            track_p_temp, track_b_temp, track_z_temp, track_n_temp, track_d_temp, track_o_temp, params, II, JJ)

       # Collect tendency over timestep
       track_p1 = track_p_temp .+ dt/2 .* track_dpdt1
       track_b1 = track_b_temp .+ dt/2 .* track_dbdt1
       track_z1 = track_z_temp .+ dt/2 .* track_dzdt1
       track_n1 = track_n_temp .+ dt/2 .* track_dndt1
       track_d1 = track_d_temp .+ dt/2 .* track_dddt1
       track_o1 = track_o_temp .+ dt/2 .* track_dodt1

       # Second calculation
       track_dpdt2, track_dbdt2, track_dzdt2, track_dndt2, track_dddt2, track_dodt2 = one_step_ecof(
           track_p1, track_b1, track_z1, track_n1, track_d1, track_o1, params, II, JJ)

       # Collect tendency over timestep
       track_p2 = track_p_temp .+ dt/2 .* track_dpdt2
       track_b2 = track_b_temp .+ dt/2 .* track_dbdt2
       track_z2 = track_z_temp .+ dt/2 .* track_dzdt2
       track_n2 = track_n_temp .+ dt/2 .* track_dndt2
       track_d2 = track_d_temp .+ dt/2 .* track_dddt2
       track_o2 = track_o_temp .+ dt/2 .* track_dodt2

       # Third calculation
       track_dpdt3, track_dbdt3, track_dzdt3, track_dndt3, track_dddt3, track_dodt3 = one_step_ecof(
           track_p2, track_b2, track_z2, track_n2, track_d2, track_o2, params, II, JJ)

       # Collect tendency over timestep
       track_p3 = track_p_temp .+ dt .* track_dpdt3
       track_b3 = track_b_temp .+ dt .* track_dbdt3
       track_z3 = track_z_temp .+ dt .* track_dzdt3
       track_n3 = track_n_temp .+ dt .* track_dndt3
       track_d3 = track_d_temp .+ dt .* track_dddt3
       track_o3 = track_o_temp .+ dt .* track_dodt3

       # Fourth calculation
       track_dpdt4, track_dbdt4, track_dzdt4, track_dndt4, track_dddt4, track_dodt4 = one_step_ecof(
           track_p3, track_b3, track_z3, track_n3, track_d3, track_o3, params, II, JJ)

       # Apply Runge-Kutta formula to the 4 calculations
       track_p_temp .+= (track_dpdt1 .+ 2 .* track_dpdt2 .+ 2 .* track_dpdt3 .+ track_dpdt4) .* (dt / 6)
       track_b_temp .+= (track_dbdt1 .+ 2 .* track_dbdt2 .+ 2 .* track_dbdt3 .+ track_dbdt4) .* (dt / 6)
       track_z_temp .+= (track_dzdt1 .+ 2 .* track_dzdt2 .+ 2 .* track_dzdt3 .+ track_dzdt4) .* (dt / 6)
       track_n_temp .+= (track_dndt1 .+ 2 .* track_dndt2 .+ 2 .* track_dndt3 .+ track_dndt4) .* (dt / 6)
       track_d_temp .+= (track_dddt1 .+ 2 .* track_dddt2 .+ 2 .* track_dddt3 .+ track_dddt4) .* (dt / 6)
       track_o_temp .+= (track_dodt1 .+ 2 .* track_dodt2 .+ 2 .* track_dodt3 .+ track_dodt4) .* (dt / 6)
       
       # If timestep belongs to number of output records, record the values
       if mod(t, trec)==0
            j = Int(t/trec + 1)
            t_id = t.*dt
            track_p[:,:,j] .= track_p_temp
            track_b[:,:,j] .= track_b_temp 
            track_z[:,:,j] .= track_z_temp 
            track_n[:,:,j] .= track_n_temp 
            track_d[:,:,j] .= track_d_temp
            track_o[:,:,j] .= track_o_temp
            track_time[j] = t_id 

            # Announce recording of results (then record it in log file)
            @printf("Day %7.1f out of %5.0f = %4.0f%% done at %s \n", 
                t_id, params.tt, t_id/params.tt*100, now())
            open("out/jlog.txt","a") do f
                write(f,@sprintf("Day %7.1f out of %5.0f = %4.0f%% done at %s \n", t_id, params.tt, t_id/params.tt*100, now()))
            end

            # Ensure conservation of total N. This simple summing works only because the height of all boxes is the same.
            println("Total N: ", sum(track_p_temp) + sum(track_b_temp) + sum(track_n_temp) + sum(track_d_temp) + sum(track_z_temp))
        end 
    end #end timestepping
    tfn = now() #finish time
    
    # Write files--save after a loop for now. 
    if params.fsave == "missing"
    else
        savetoNC(fsaven, track_p, track_b, track_z, track_n, track_d, track_o,  tst, tfn, params)
    end
    return track_p_temp, track_b_temp, track_z_temp, track_n_temp, track_d_temp, track_o_temp
end

#################################################################################
# Microbial ecosystem model 
#################################################################################
#
# inputs =  p, b, z, n, d 
#
function one_step_ecof(p, b, z, n, d, o, params, II, JJ) # II and JJ record locations of non zero in the sparse matrix (OM * het)

    # Transport, all are for the diffusion term of the function.
    dpdt = DIFF(p, kappa_z, dz) .- UPWIND(p, w, dz)  
    dbdt = DIFF(b, kappa_z, dz) .- UPWIND(b, w, dz)
    dzdt = DIFF(z, kappa_z, dz) .- UPWIND(z, w, dz)
    dndt = DIFF(n, kappa_z, dz) .- UPWIND(n, w, dz)
    dddt = DIFF(d, kappa_z, dz) .- UPWIND(d, wd, dz) # upwind is for sinking term
    dodt = DIFF(o, kappa_z, dz) .- UPWIND(o, w, dz)

    d_gain_total = zeros(params.number_box) # collect OM and then distribute into diff OM pools

    #(below is biological term = growth, nutrient uptake, loss):
    #phytoplankton
    for j = 1:np

        #growth and loss:
        #should do a more complicated N limitation, but OK for now:
        totn = sum(n[:,1:3], dims = 2)
        uptake = params.TempFun .* params.umax_p[j] .* min.(totn ./ (totn .+ params.K_pn[j]), Light ./ (Light .+ params.K_I[j])) .* p[:,j]
        #println("n = ",n[1]," and mu_p = ",uptake[1]./p[1,j])
        mort = (params.m_lp[j] .+ params.m_qp[j] .* p[:,j]) .* p[:,j] # mortality

        #track: (S = S + x can be writen as S += x)
        dndt[:,1] += -uptake .* n[:,1] ./ totn
        dndt[:,2] += -uptake .* n[:,2] ./ totn
        dndt[:,3] += -uptake .* n[:,3] ./ totn
        dodt += uptake*params.e_o # oxygen excretion, e_o in input_traits (use Redfield ratio for now), related to how much O2 is produced per phyto growth, 
        dpdt[:,j] += uptake - mort # uptake * 1/y = growth of biomass, but y = 1 for phyto, no reminerilization.
        d_gain_total += mort #addition to total mortality

    end

    #heterotrophic bacteria growth:
    #@inbounds for n = axes(II, 1)
    for k = axes(II, 1)
        #k from 1 to n of interactions,
        #II[k] values are from 1 to nd (IF there is an interaction at i=nd!) # (OM)
        #JJ[k] values are from 1 to nb # (heterotrophs)
        
        #set up as:
        #1. Aerobic
        #2. NO3 to NO2
        #3. NO3 to N2O
        #4. NO3 to N2 
        #5. NO2 to N2O
        #6. NO2 to N2
        #7. N2O to N2

        #first aerobic heterotrophs
        #if params.y_o[JJ[k]] > 0 #aerobic hets
        if JJ[k] == 1 #index by j in 1:nb # this find the first one: aerobic heterotrophs
       
            mu = params.TempFun .* 
                 min.(params.vmax_d[II[k],JJ[k]] .* d[:,II[k]] ./ (d[:,II[k]] .+ params.K_d[II[k],JJ[k]]) .* params.y_d[II[k],JJ[k]], 
                      params.pcoef[JJ[k]] .* o .* params.y_o[JJ[k]]) #growth rate limited by either d or oxygen
            dodt += -mu .* b[:,JJ[k]] ./ params.y_o[JJ[k]]

        #now anaerobic heterotrophs
        elseif JJ[k] > 1 #all anaerobes (need this separate from aerobes for DIN matrix below)
            
            if JJ[k] == 2 #NO3 to NO2 
                ndinin = 3 #index for DIN required (for respiration)
                ndinex = 2 #index for DIN excreted (from respiration: see NH4 excretion below)
            elseif JJ[k] == 3 #NO3 to N2O 
                ndinin = 3 #index for DIN required (for respiration)
                ndinex = 4 #index for DIN excreted (from respiration: see NH4 excretion below)
            elseif JJ[k] == 4 #NO3 to N2 
                ndinin = 3 #index for DIN required (for respiration)
                ndinex = 5 #index for DIN excreted (from respiration: see NH4 excretion below)
            elseif JJ[k] == 5 #NO2 to N2O
                ndinin = 2 #index for DIN required (for respiration)
                ndinex = 4 #index for DIN excreted (from respiration: see NH4 excretion below)
            elseif JJ[k] == 6 #NO2 to N2
                ndinin = 2 #index for DIN required (for respiration)
                ndinex = 5 #index for DIN excreted (from respiration: see NH4 excretion below)
            elseif JJ[k] == 7 #N2O to N2
                ndinin = 4 #index for DIN required (for respiration)
                ndinex = 5 #index for DIN excreted (from respiration: see NH4 excretion below)
            end
       
            #calculate mu and DIN respiration input and output for anaerobic heterotrophs using ndinin and ndinex:
            mu = params.TempFun .* 
                 min.(params.vmax_d[II[k],JJ[k]] .* d[:,II[k]] ./ (d[:,II[k]] .+ params.K_d[II[k],JJ[k]]) .* params.y_d[II[k],JJ[k]], #this line is for OM
                      params.vmax_n[ndinin,JJ[k]] .* n[:,ndinin] ./ (n[:,ndinin] .+ params.K_n[ndinin,JJ[k]]) .* params.y_n[ndinin,JJ[k]]) # this line is for DIN uptake rate
            dndt[:,ndinin] += -mu .* b[:,JJ[k]] ./ params.y_n[ndinin,JJ[k]] # uptake of DIN
            dndt[:,ndinex] += mu .* b[:,JJ[k]] ./ params.y_n[ndinin,JJ[k]] # excretion of DIN

        end # end of anaerobic heterotrophs 

        #similar terms for all heterotrophs: growth, uptake of d, and excretion of NH4, using specifically-calculated mu:
        dbdt[:,JJ[k]] += mu .* b[:,JJ[k]] #biomass synthesis/growth
        dddt[:,II[k]] += -mu .* b[:,JJ[k]] ./ params.y_d[II[k],JJ[k]] #uptake of organic matter
        dndt[:,1] += mu .* b[:,JJ[k]] .* (1. ./ params.y_d[II[k],JJ[k]] - 1.) #excretion of NH4 

    end # end of all heterotrophs.

    #AOO
    j = params.nhets + 1 #which b?
    k = 1 #which n?
    mu = params.TempFun .* min.(params.vmax_n[k,j] .* n[:,k] ./ (n[:,k] .+ params.K_n[k,j]) .* params.y_n[k,j], pcoef[j] .* o .* params.y_o[j])
    dbdt[:,j] += b[:,j] .* mu  
    dndt[:,k] += -b[:,j] .* mu ./ params.y_n[k,j]
    dndt[:,k+1] += b[:,j] .* mu .* (1. / params.y_n[k,j] - 1.)
    dodt += -b[:,j] .* mu ./ params.y_o[j]
    # Xin add(09/2023): N2O production from ammonia 
    dndt[:,4] += (0.002 ./max.(0.0685, dodt) .+ 0.0008) .* max.(0, b[:,j] .* mu .* (1. / params.y_n[k,j] - 1.)) #(Ji et al., 2018)
 

    #NOO
    j = params.nhets + 2 #which b?
    k = 2 #which n?
    mu = params.TempFun .* min.(params.vmax_n[k,j] .* n[:,k] ./ (n[:,k] .+ params.K_n[k,j]) .* params.y_n[k,j], pcoef[j] .* o .* params.y_o[j])
    dbdt[:,j] += b[:,j] .* mu  
    dndt[:,k] += -b[:,j] .* mu ./ params.y_n[k,j]
    dndt[:,k+1] += b[:,j] .* mu ./ params.y_n[k,j]
    dndt[:,1] += -b[:,j] .* mu #uses NH4 or simple org N for synthesis (negligible, so just have as NH4) 
    dodt += -b[:,j] .* mu ./ params.y_o[j]
    
    #Anammox
    j = params.nhets + 3 #which b?
    k1 = 1 #which n?
    k2 = 2 
    mu = params.TempFun .* min.(params.vmax_n[k1,j] .* n[:,k1] ./ (n[:,k1] .+ params.K_n[k1,j]) .* params.y_n[k1,j], 
                                params.vmax_n[k2,j] .* n[:,k2] ./ (n[:,k2] .+ params.K_n[k2,j]) .* params.y_n[k2,j])
    dbdt[:,j] += b[:,j] .* mu  
    dndt[:,k1] += -b[:,j] .* mu ./ params.y_n[k1,j] #uptake of NH4
    dndt[:,k2] += -b[:,j] .* mu ./ params.y_n[k2,j] #uptake of NO2
    dndt[:,3] += b[:,j] .* mu .* (1. ./params.y_n[k2,j] + 1. ./params.y_n[k2,j] - 1.)*params.ftoNO3_anx #excretion of NO3
    dndt[:,5] += b[:,j] .* mu .* (1. ./params.y_n[k2,j] + 1. ./params.y_n[k2,j] - 1.)*(1. - params.ftoNO3_anx) #excretion of N2

    for j = 1:nz
        #Grazing JUST for bacteria:
        #if sum(GrM[j,np+1:end]) > 0 #bacteria grazers only here.
        #    prey = sum(GrM[j,np+1:end]'.*b, dims = 2) 
        #    g = params.g_max[j] .* prey ./ (prey .+ params.K_g[j])
        #    dzdt[:,j] += params.γ[j] .* g .* z[:,j]
        #    dndt[:,1] += (1 - params.γ[j]) .* g .* z[:,j]
        #    dbdt +=  -g .* z[:,j] .* GrM[j,np+1:end]' .* b ./ prey
        #end
        preyp = sum(GrM[j,1:np]'.*p, dims = 2) #from pp 
        prey = preyp .+ sum(GrM[j,np+1:end]'.*b, dims = 2) #from b
        g = params.g_max[j] .* prey ./ (prey .+ params.K_g[j])
        dzdt[:,j] += params.γ[j] .* g .* z[:,j]
        dndt[:,1] += (1 - params.γ[j]) .* g .* z[:,j]
        dodt[:,1] += -(1 - params.γ[j]) .* g .* z[:,j] * params.e_o
        dpdt +=  -g .* z[:,j] .* GrM[j,1:np]' .* p ./ prey
        dbdt +=  -g .* z[:,j] .* GrM[j,np+1:end]' .* b ./ prey
    end

    #bacterial mortality and addition to total OM
    mort = (transpose(params.m_lb) .+ transpose(params.m_qb) .* b) .* b
    dbdt += -mort
    d_gain_total += sum(mort, dims = 2)
    
    #zooplankton mortality and addition to total OM
    mort = (transpose(params.m_lz) .+ transpose(params.m_qz) .* z) .* z # mort overwrite the above mort, now it is diff, it is for zoo
    dzdt += -mort
    d_gain_total += sum(mort, dims = 2)
    
    #split accumulated OM into nd pools according to probability of generation, now it is 50%, 50%
    dddt += d_gain_total .* transpose(prob_generate_d)

    #oxygen
    dodt[1] += params.koverh.*(params.o2sat .- o[1]) #air-sea flux
    if params.relax_o2 == 1
        dodt += params.t_o2relax.*(params.o2_obs .- o) # O2 relaxation via lateral flux, entire profile
    end
    #dodt[1:params.mlboxes] += params.koverh.*(params.o2sat .- o[1:params.mlboxes]) # ! (may need to change with Emily) oxygen input from air-sea
    #dodt += params.t_o2relax.*(params.o2_deep .- o) #O2 relaxation via lateral flux, not only to the bottom box, but to all boxes
    #dodt[end] += params.t_o2relax.*(params.o2_deep .- o[end]) #O2 relaxation via lateral flux, only to the bottom box.
    
    #N2O and N2 air-sea flux
    dndt[1,4] += params.koverh.*(0.0122 .- dndt[1,4])#N2O in the 1st box (i.e., surface)
    dndt[1,5] += params.koverh.*(385.7*2 .- dndt[1,5])#N2 in the 1st box (i.e., surface)

    #no3
    if params.relax_no3 == 1
        dndt[:,3] += params.t_no3relax.*(params.no3_obs .- dndt[:,3]) # NO3 relaxation via lateral flux, entire profile
    end
    
    return dpdt, dbdt, dzdt, dndt, dddt, dodt
end

#################################################################################
# Diffusion function 
#################################################################################
function DIFF(c, kappa_z, dz)
    c = vcat(transpose(zeros(size(c)[2])), c, transpose(zeros(size(c)[2])))
    kappa_z_len = size(kappa_z)[1]
    f_up = kappa_z[1:kappa_z_len-1] .* (c[2:kappa_z_len,:] .- c[1:kappa_z_len-1,:]) ./dz 
    f_down = kappa_z[2:kappa_z_len] .* (c[3:kappa_z_len+1,:] .- c[2:kappa_z_len,:]) ./dz
    diff = (f_down .- f_up) ./ dz
    return diff
end

#################################################################################
# Upwind function 
#################################################################################
function UPWIND(c, wd, dz) #written for wd, but also works for 1D w!
    c = vcat(zeros(2,size(c)[2]), c, zeros(2, size(c)[2]))
    positive_wd = (wd + abs.(wd))./2. 
    negative_wd = (wd - abs.(wd))./2.
    Fu = positive_wd[2:end,:].*c[3:end-2,:] + negative_wd[2:end,:].*c[4:end-1,:]
    Fd = positive_wd[1:end-1,:].*c[2:end-3,:] + negative_wd[1:end-1,:].*c[3:end-2,:]
    adv = (Fu .- Fd)./dz 
    return adv 
end

#################################################################################
# Netcdf writing function 
#################################################################################
function savetoNC(fsaven, p, b, z, n, d, o, tst, tfn, params)

    println("Saving to: ",fsaven)
    
    #if isfile(fsaven)
    #    f = NCDataset(fsaven, "a") #a for annotate -- this COPIES OVER the previous
    #else
    f = NCDataset(fsaven, "c") #c for create
    #end

    total_tt = Int(params.tt+1) #bc i added time 0

    # define the dim of p, b, z, n, d
    defDim(f,"np",params.np)
    defDim(f,"nb",params.nb)
    defDim(f,"nz",params.nz)
    defDim(f,"nn",params.nn)
    defDim(f,"nd",params.nd)

    # define the dim of the depth
    defDim(f,"ndepth",params.number_box)
    defDim(f,"ndepth1",params.number_box+1)

    # define the dim of the time length
    #defDim(f,"nrec",total_tt)
    nrec1 = Int(params.nrec+1) #bc i added time 0
    defDim(f,"nrec",nrec1)
   
    # info
    f.attrib["title"] = "Community consumption pbznd model i/o"
    f.attrib["Start time"] = string(tst)
    f.attrib["End time"] = string(tfn)
    f.attrib["Run time"] = string(tfn - tst) 

    # simulated results
    w = defVar(f,"p",Float64,("ndepth" ,"np","nrec"))
    w[:,:,:] = p
    w.attrib["units"] = "mmol/m3 C biomass"

    w = defVar(f,"b",Float64,("ndepth" ,"nb","nrec"))
    w[:,:,:] = b
    w.attrib["units"] = "mmol/m3 C biomass"

    w = defVar(f,"z",Float64,("ndepth" ,"nz","nrec"))
    w[:,:,:] = z
    w.attrib["units"] = "mmol/m3 C biomass"
    
    w = defVar(f,"n",Float64,("ndepth" ,"nn","nrec"))
    w[:,:,:] = n
    w.attrib["units"] = "mmol/m3 C OM"

    w = defVar(f,"d",Float64,("ndepth" ,"nd","nrec"))
    w[:,:,:] = d
    w.attrib["units"] = "mmol/m3 C OM"
    
    w = defVar(f,"o",Float64,("ndepth" ,"nrec"))
    w[:,:] = o
    w.attrib["units"] = "mmol/m3 O2"
    
    # w = defVar(f,"v",Float64,("nd","nb"))
    # w[:,:,:] = v
    # w.attrib["units"] = "per d; specific uptake matrix with PA impact"
    
    # w = defVar(f,"uptake",Float64,("nd","nb"))
    # w[:,:,:] = uptake
    # w.attrib["units"] = "mmol/m3 C per d; uptake matrix"
    
    w = defVar(f,"pIC",Float64,("ndepth","np"))
    w[:,:] = params.pIC
    w.attrib["units"] = "mmol/m3 C biomass"

    w = defVar(f,"bIC",Float64,("ndepth","nb"))
    w[:,:] = params.bIC
    w.attrib["units"] = "mmol/m3 C biomass"

    w = defVar(f,"zIC",Float64,("ndepth","nz"))
    w[:,:] = params.zIC
    w.attrib["units"] = "mmol/m3 C biomass"

    w = defVar(f,"nIC",Float64,("ndepth","nn"))
    w[:,:] = params.nIC
    w.attrib["units"] = "mmol/m3 C OM"

    w = defVar(f,"dIC",Float64,("ndepth","nd"))
    w[:,:] = params.dIC
    w.attrib["units"] = "mmol/m3 C OM"
    
    w = defVar(f,"oIC",Float64,("ndepth","nd"))
    w[:,:] = params.dIC
    w.attrib["units"] = "mmol/m3 O2"

    # w = defVar(f, "Ctot", Float64, ("nrec",))
    # w[:] = Ctot
    # w.attrib["units"] = "mmol/m3 total C"
    
    # w = defVar(f,"Cw",Float64,("nd",))
    # w[:,:] = params.Cw 
    # w.attrib["units"] = "Ind C supply weight: probability"
    
    # w = defVar(f,"Cwall",Float64,("nd","nrec"))
    # w[:,:] = Cwall 
    # w.attrib["units"] = "Ind C supply weight: over time"
    
    # w = defVar(f,"PA",Float64,("nb",))
    # w[:,:] = params.PA 
    # w.attrib["units"] = "Ind C supply weight: over time"
    
    # w = defVar(f,"PAall",Float64,("nb","nrec"))
    # w[:,:] = PAall 
    # w.attrib["units"] = "Ind C supply weight: over time"

    # w = defVar(f,"pen",Float64,("nb",))
    # w[:,:] = params.pen
    # w.attrib["units"] = "penalty"
    
    #w = defVar(f, "timet", Float64, ("nrec",))
    #w[:] = timet
    #w.attrib["units"] = "days"
    
    # w = defVar(f, "Pt", Float64, ())
    # w[:] = params.Pt
    # w.attrib["units"] = "Total C supply"

    w = defVar(f, "H", Int, ())
    w[:] = params.H
    w.attrib["units"] = "m; total height"

    w = defVar(f, "dz", Int, ())
    w[:] = params.dz
    w.attrib["units"] = "m; box height"

    w = defVar(f, "umax_p", Float64, ("np",))
    w[:] = params.umax_p
    w.attrib["units"] = "m3/mmol/d; max growth rate of p"

    w = defVar(f, "K_pn", Float64, ("np",))
    w[:] = params.K_pn
    w.attrib["units"] = "m3/mmol; half-sat rate of p"

#    w = defVar(f, "m_lp", Float64, ("ndepth" ,"np"))
#    w[:] = params.m_lp
#    w.attrib["units"] = "m3/mmol; linear death rate of p"

#    w = defVar(f, "m_qp", Float64, ("ndepth" ,"np"))
#    w[:] = params.m_qp
#    w.attrib["units"] = "m3/mmol; quadratic death rate of p"

    w = defVar(f,"CM",Float64,("nd","nb"))
    w[:,:] = params.CM
    w.attrib["units"] = "Community consumption Matrix"

#    # w = defVar(f, "y_i", Float64, ("nd", "nb"))
#    # w[:] = params.y_i
#    # w.attrib["units"] = "mol B/mol C; yield"

    w = defVar(f, "y_d", Float64, ("nd","nb"))
    w[:,:] = params.y_d
    w.attrib["units"] = "per d; max yield rate"

#    #w = defVar(f, "vmax_i", Float64, ("nd",))
#    #w[:,:] = params.vmax_i
#    #w.attrib["units"] = "per d; max uptake rate"
    
    w = defVar(f, "vmax_d", Float64, ("nd", "nb"))
    w[:,:] = params.vmax_d
    w.attrib["units"] = "per d; max uptake rate"
    
    w = defVar(f, "K_d", Float64, ("nd", "nb"))
    w[:] = params.K_d
    w.attrib["units"] = "mmol/m3; half-sat"
    
#    w = defVar(f, "m_lb", Float64, ("ndepth" ,"nb"))
#    w[:] = params.m_lb
#    w.attrib["units"] = "m3/mmol; linear death rate of b"
#
#    w = defVar(f, "m_qb", Float64, ("ndepth" ,"nb"))
#    w[:] = params.m_qb
#    w.attrib["units"] = "m3/mmol; quadratic death rate of b"
#    
#    w = defVar(f, "mu_mz", Float64, ("ndepth" ,"nz"))
#    w[:] = params.mu_mz
#    w.attrib["units"] = "m3/mmol; grazing rate of z"
    
    w = defVar(f, "K_g", Float64, ("nz",))
    w[:] = params.K_g
    w.attrib["units"] = "m3/mmol; half-sat rate of z"
    
    w = defVar(f, "γ", Float64, ("nz",))
    w[:] = params.γ
    w.attrib["units"] = "fraction of digestion of z"
    
#    w = defVar(f, "m_lz", Float64, ("ndepth" ,"nz"))
#    w[:] = params.m_lz
#    w.attrib["units"] = "m3/mmol; linear death rate of z"
#
#    w = defVar(f, "m_qz", Float64, ("ndepth" ,"nz"))
#    w[:] = params.m_qz
#    w.attrib["units"] = "m3/mmol; quadratic death rate of z"
    
    w = defVar(f, "prob_generate_d", Float64, ("nd",))
    w[:] = params.prob_generate_d
    w.attrib["units"] = "distribution from m to each d"
    
    w = defVar(f, "kappa_z", Float64, ("ndepth1",))
    w[:] = params.kappa_z
    w.attrib["units"] = "vertical water velocity"
    
    w = defVar(f, "wd", Float64, ("ndepth1","nd"))
    w[:] = params.wd
    w.attrib["units"] = "sinking rate"

    close(f)
end

