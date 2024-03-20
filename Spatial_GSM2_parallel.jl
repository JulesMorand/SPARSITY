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
        z::Float64
		r::Float64### r_nucleus
        R::Float64### R_cell 
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
global R_cell=10.;
global cell = Cell(0.0, 0.0, 0.0, r_nucleus, R_cell)
global ion = Ion("4He", 56.0, 1, 4.5, 1.0)
global irrad = Irrad(1.0, 0.8, 0.18)

global (Rc, Rp, Rk) = ATRadius(ion, irrad)
println("Rc=", Rc, "\nRp=", Rp, "\nRk=", Rk)

x, y = GenerateHit(cell,Rk)
track = Track(x, y, Rk)

global DoseRate_h = irrad.doserate * 3600
global F = irrad.dose / (1.602 * 10^(-9) * ion.LET)
global Npar = round(Int, F * (pi * (R + Rk)^2 * 10^(-8)))
# if no simulation of time
#Np = round(Int, Dose/(1.602*10^(-9)*LET/(pi*(R+Rk)^2*10^(-8))));
#Np = Npar;

#zF = 1.602*10^(-9)*LET/(pi*(R+Rk)^2*10^(-8));
global zF = irrad.dose / Npar
global D = DoseRate_h / zF
global T = irrad.dose / (zF * D) * 3600
#######################


#Simulation dose
#Nd is the dimension fo the space, cell is assumed to be a cylinder. The geometry can be more complicated but for now it is fine
Nd = 3;

#X and Y are matrices with the cartesian coordinates of the damages
X = Array{Float64}(undef, 0, Nd);
Y = Array{Float64}(undef, 0, Nd);

#############This works fine
global survP = Array{Float64};
global Np       = rand(Poisson(Npar))
global DOSE_tot = 0.0
global GYR_tot  = 0.0
println("Number of particles = $Np")
global X = Array{Float64}(undef, 0, Nd)
global Y = Array{Float64}(undef, 0, Nd)

#from here, if I run the code alone it works, if I run the function it doesn't!

global gsm2 = GSM2(4.0, 0.1, 0.1, 0.8)


global Np = rand(Poisson(Npar))
global GYR_tot = 0.0;

@everywhere global X = Array{Float64}(undef, 0, Nd);
@everywhere global Y = Array{Float64}(undef, 0, Nd);

simulate_GSM2(ion, cell, gsm2, Np, Rk, GYR_tot, X, Y)

global sim = 10;
#@everywhere survP = zeros(sim);
@time begin
	@everywhere global survP = Array{Float64}(undef, 0);
	Threads.@threads for i_ in 1:sim
		println("ii = $i_ on thread $(Threads.threadid())")
		local Np = rand(Poisson(Npar))
		local GYR_tot = 0.0;

		local X_ = Array{Float64}(undef, 0, Nd); 
		local Y_ = Array{Float64}(undef, 0, Nd);

		@everywhere survPi_, X_, Y_ = simulate_GSM2(ion, cell, gsm2, Np, Rk, GYR_tot, X, Y);

		push!(survP, survPi_);
		#print(X_)
		#push!(Y,Y_)
		#survP[i_,:] = survPi;
	end
end
