using Base.Threads
using Distributed
using CSV, DataFrames
using Distributions
using Random
using Plots
using Distances

#initialize structures
begin
    struct Cell
        x::Float64
        y::Float64
        r::Float64
    end
    
    struct Track
        x::Float64
        y::Float64
        Rk::Float64
    end
    
    struct Ion
        ion::String
        E::Float64
        A::Int64
        LET::Float64
        rho::Float64
    end
    
    struct Irrad
        dose::Float64
        kR::Float64
        doserate::Float64
    end
    
    struct GSM2
        r::Float64
        a::Float64
        b::Float64
        rd::Float64
    end
    
    struct NTGSM2
        r::Float64
        a::Float64
        b::Float64
        rd::Float64
    
        rM::Float64
        aM::Float64
        bM::Float64
        rdM::Float64
    
        pa::Float64
        qa::Float64
    end
end

include("utilities_spatial.jl")

global R = r_nucleus = 8.0;
global cell = Cell(0.,0., R)

global ion = Ion("4He", 56.0, 1, 4.5, 1.0)
global irrad = Irrad(1.0, 0.8, 0.18)

global (Rc, Rp, Rk) = ATRadius(ion, irrad)
println("Rc=",Rc,"\nRp=",Rp,"\nRk=",Rk)

x, y = GenerateHit(cell, Rk)
track = Track(x, y, Rk)

global DoseRate_h = irrad.doserate*3600
global F = irrad.dose/(1.602*10^(-9)*ion.LET)
global Npar = round(Int, F*(pi*(R+Rk)^2*10^(-8)))
# if no simulation of time
#Np = round(Int, Dose/(1.602*10^(-9)*LET/(pi*(R+Rk)^2*10^(-8))));
#Np = Npar;

#zF = 1.602*10^(-9)*LET/(pi*(R+Rk)^2*10^(-8));
global zF = irrad.dose/Npar
global D = DoseRate_h/zF
global T = irrad.dose/(zF*D)*3600
#######################


#Simulation dose
#Nd is the dimension fo the space, cell is assumed to be a cylinder. The geometry can be more complicated but for now it is fine
Nd = 3;

#X and Y are matrices with the cartesian coordinates of the damages
X = Array{Float64}(undef, 0, Nd);
Y = Array{Float64}(undef, 0, Nd);


#run code to simulate the spatial distribution of the damage for all the tracks
#the function calculate_damage is not optimized but for the moment is only slighlty slower than the code that calculates the dose 
#e.g. the dose deposited by 1046 tracks is computed in 0.56 seconds
#the dose deposited by 1109 tracks and the corresponding spatial distribution of the damages is computed in 0.62 seconds
@time begin
    Np=rand(Poisson(Npar))
    cell = Cell(0.,0., R)
    global DOSE_tot = 0. ;
    global GYR_tot  = 0. ;
    
    println("Number of particles = $Np")

    for i in 1:Np
        x, y = GenerateHit(cell, Rk);
        track = Track(x, y, Rk);
        dose, area, Gyr = distribute_dose(cell,track);
        DOSE_tot += dose
        GYR_tot += Gyr
    end
    println("Total imparted dose = $GYR_tot")
end

@time begin
    Np=rand(Poisson(Npar))
    cell = Cell(0.,0., R)
    global DOSE_tot = 0. ;
    global GYR_tot  = 0. ;

    println("Number of particles = $Np")

    @everywhere DOSE_tot = Array{Float64}(undef, 0);
    @everywhere GYR_tot = Array{Float64}(undef, 0);

    Threads.@threads for i in 1:Np
        global x, y = GenerateHit(cell, Rk);
        global track = Track(x, y, Rk);
        @everywhere dose, area, Gyr = distribute_dose(cell,track);
        push!(DOSE_tot,dose)
        push!(GYR_tot,Gyr)
    end
    sum(DOSE_tot)
    dose_cell = sum(GYR_tot)

    println("Total imparted dose = $dose_cell")

end


#############This works fine

global Np = rand(Poisson(Npar))
global cell = Cell(0.,0., R)
global DOSE_tot = 0. ;
global GYR_tot  = 0. ;

println("Number of particles = $Np")
X = Array{Float64}(undef, 0, Nd);
Y = Array{Float64}(undef, 0, Nd);
@time begin
    for i in 1:Np
        #println("ii = $ii on thread $(Threads.threadid())")
        x, y = GenerateHit(cell, Rk);
        track = Track(x,y,Rk)
        integral, theta, Gyr, radius = distribute_dose_vector(ion, cell, track);
    
        X_, Y_ = calculate_damage(ion, cell, integral, theta, Gyr, radius);

        dist = sqrt.(X_[:,1].*X_[:,1] .+ X_[:,2].*X_[:,2])
        if size(dist[dist .> 8],1) != 0
            println("Error")
            return x, y 
        end

        X = vcat(X,X_)
        Y = vcat(Y,Y_)
        #DOSE_tot+=dose
        GYR_tot+=Gyr
    end
    println("Total imparted dose = $GYR_tot")
end

dist = sqrt.(X[:,1].*X[:,1] .+ X[:,2].*X[:,2]);
if size(dist[dist .> 8],1) != 0
    println("Error")
else
    println("Everything looks good")
end



#from here, if I run the code alone it works, if I run the function it doesn't!

cell = Cell(0.,0., R)
gsm2 = GSM2(4.,0.1,0.1,0.8)

Np = rand(Poisson(Npar))
GYR_tot  = 0. ;

X = Array{Float64}(undef, 0, Nd);
Y = Array{Float64}(undef, 0, Nd);

function simulate_GSM2(ion::Ion, cell::Cell, gsm2::GSM2, Np::Int64, Nd::Int64, Rk::Float64, GYR_tot::Float64, X::Matrix{Float64}, Y::Matrix{Float64})

    for i in 1:Np
        #println("ii = $ii on thread $(Threads.threadid())")
        x, y = GenerateHit(cell, Rk);
        track = Track(x,y,Rk)
        integral, theta, Gyr, radius = distribute_dose_vector(ion, cell, track);
    
        X_, Y_ = calculate_damage(ion, cell, integral, theta, Gyr, radius);

        dist = sqrt.(X_[:,1].*X_[:,1] .+ X_[:,2].*X_[:,2])
        if size(dist[dist .> 8],1) != 0
            println("Error")
            return x, y 
        end

        X = vcat(X,X_)
        Y = vcat(Y,Y_)
        #DOSE_tot+=dose
        GYR_tot += Gyr
    end
    println("Total imparted dose = $GYR_tot")

    survP = spatial_GSM2_fast(X, Y, gsm2)

    return survP

end

Np = rand(Poisson(Npar))
GYR_tot  = 0. ;

X = Array{Float64}(undef, 0, Nd);
Y = Array{Float64}(undef, 0, Nd);

simulate_GSM2(ion, cell, gsm2, Npar, Np, Rk, GYR_tot, X, Y)

@time begin
    survP = simulate_GSM2(ion, cell, gsm2, Npar, Np, Rk, GYR_tot, X, Y)
end


sim = 10^3;
global survP = Array{Float64};
Threads.@threads for i in 1:sim
    survP[i] = simulate_GSM2(ion, cell, gsm2, Npar, Np, Rk)
end
