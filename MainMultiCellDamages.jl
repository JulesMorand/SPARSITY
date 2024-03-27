using Base.Threads
using Distributed
using CSV, DataFrames
using Distributions
using Random
using Plots
using Distances
using WebIO,PlotlyJS
#Nd is the dimension fo the space, the nucleus of the cell is assumed to be a cylinder. Thde cell a sphere around the centerThe geometry can be more complicated but for now it is fine
Nd = 3;
#initialize structures
begin
	mutable struct Cell
		x::Float64
		y::Float64
        z::Float64
		r::Float64### r_nucleus    
        R::Float64### R_cell
		Dam_X::Array{Float64}
		Dam_Y::Array{Float64}
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
include("./FunctionsCreationCells.jl") 
#include("./FunctionsVisualisation.jl")

###Parameters 
###Dimension of the hit box
global X_box=500.;  #1000µmm side size of the square box
###Size Cell nucleus
global R = r_nucleus = r_nucl= 8.0;; #8µm
###Size Cell
global R_cell=30.; #30µm
#############CHOOSE squared or triangluar Lattice####
#N,nodes_positions=generate_cells_positions_squaredlattice(X_box,R_cell)
N,nodes_positions = generate_cells_positions_triangularlattice(X_box,R_cell) 
println(N)
#############Create Cells Array######################
arrayOfCell = Creation_ArrayOfCell(N, nodes_positions, r_nucl, R_cell);
##########################################################
global ion = Ion("4He", 56.0, 1, 4.5, 1.0)
global irrad = Irrad(1.0, 0.8, 0.18)
global (Rc, Rp, Rk) = ATRadius(ion, irrad)
println("Rc=", Rc, "\nRp=", Rp, "\nRk=", Rk)

global DoseRate_h = irrad.doserate * 3600
global F = irrad.dose / (1.602 * 10^(-9) * ion.LET)
#########################################

global R_beam=100.; #Radius of the beam
global x_cb,y_cb=125.,125.; #Coordiante of the center of the beam
#### dose for circle beam
global Npar = round(Int, F * (pi * (R_beam)^2 * 10^(-8)))
#### dose for squarte box
#global Npar = round(Int, F * ((X_box^2) * 10^(-8)))

# if no simulation of time
#Np = round(Int, Dose/(1.602*10^(-9)*LET/(pi*(R+Rk)^2*10^(-8))));
#Np = Npar;

#zF = 1.602*10^(-9)*LET/(pi*(R+Rk)^2*10^(-8));
global zF = irrad.dose / Npar
global D = DoseRate_h / zF
global T = irrad.dose / (zF * D) * 3600
# ############Test function###
@time begin 
    Np = rand(Poisson(Npar))
    #Np=10;
    #local DOSE_tot = 0. ;
    global GYR_tot = zeros(length(arrayOfCell));
    println(Np)
	global X = Array{Float64}(undef, 0, Nd);
 	global Y = Array{Float64}(undef, 0, Nd);
	#global TrackArray = Array{Float64}(undef, 0, Nd);

	#@everywhere global cell_origin = Cell(0.0, 0.0, 0.0, r_nucleus, R_cell);
	@time Threads.@threads for i in 1:Np
		#println("ii = $i on thread $(Threads.threadid())")
        #local x, y = GenerateHit_BOX(X_box);
        local x, y = GenerateHit_Circle(x_cb,y_cb,R_beam )
		local track = Track(x, y, Rk);
		#global TrackArray=vcat(TrackArray,[x y 0.])
		#println(x," ",y)
		#track=Track(30,30,Rk);
        for j in 1:length(arrayOfCell) 

            cell = arrayOfCell[j];

			local integral = Array{Float64}(undef,0)
			local theta = Array{Float64}(undef,0)
			local radius = Array{Float64}(undef,0)	
			local Gyr = 0
            #println(cell)
			if (cell.x - x)^2 + (cell.y - y)^2 < (cell.r + Rk)^2
                
            	integral, theta, Gyr, radius= distribute_dose_vector(ion,cell,track);
				#println(integral)
            	local X_, Y_ = calculate_damage(ion, cell, track, integral, theta, Gyr, radius);
                dist = sqrt.((X_[:, 1] .- cell.x).^2 .+ (X_[:, 2] .-cell.y).^2)
				if size(dist[dist .> cell.r], 1) != 0
					println("Error")
					println(j)
			#return X_i
				end
				global X = vcat(X,X_)
            	global Y = vcat(Y,Y_)
				arrayOfCell[j].Dam_X=vcat(arrayOfCell[j].Dam_X,X_)
				arrayOfCell[j].Dam_Y=vcat(arrayOfCell[j].Dam_Y,Y_)
			end
			#DOSE_tot+=dose
            GYR_tot[j] += Gyr
			
			# if length(Y_cell)!=0
			# 	println(arrayOfCell[j].Dam_Y)
			# end
        end
    end
# println(DOSE_tot)
println(GYR_tot)
end 
######Plots Damages###### 
let
	include("./FunctionsVisualisation.jl")
	plt = Plot_Lattice_Cells(arrayOfCell,X_box)
	Xnc,Ync,Znc=CylinderShape(x_cb,y_cb,0.,R_beam,R_beam);
	Plots.surface!(plt,
			Xnc,Ync,Znc, 
			opacity=0.7, 
			#color=cgrad(:matter, N, categorical = true)[i], 
			color= :yellow,
			legend=false)
	display(plt)
	Plots.savefig(plt,"PlotDamages.png")
end#dfX = DataFrame(TrackArray,:auto)
#Plots.scatter!(plt,dfX.x1,dfX.x2,dfX.x3, mode="dots",markersize=0.5 , color= :yellow )
