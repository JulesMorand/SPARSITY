using Base.Threads
using Distributed
using CSV, DataFrames
using Distributions
using Random
using Plots
using Distances
using WebIO,PlotlyJS
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
include("./FunctionsCreationCells.jl") 
include("./FunctionsVisualisation.jl")

###Parameters 
###Dimension of the hit box
global X_box=250.;  #1000µmm side size of the square box
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
#Plot_Lattice_Cells(arrayOfCell)
##########################################################
global ion = Ion("4He", 56.0, 1, 4.5, 1.0)
global irrad = Irrad(1.0, 0.8, 0.18)
global (Rc, Rp, Rk) = ATRadius(ion, irrad)
println("Rc=", Rc, "\nRp=", Rp, "\nRk=", Rk)
#global x_center,y_center,z_center=X_box/2.,X_box/2.,R_cell
#global targetCell=Cell(x_center,y_center,z_center,r_nucl,r_nucl)
# x, y = GenerateHit(X_box)
# track = Track(x, y, Rk)
global DoseRate_h = irrad.doserate * 3600
global F = irrad.dose / (1.602 * 10^(-9) * ion.LET)
#global Npar = round(Int, F * (pi * (R + Rk)^2 * 10^(-8)))
global Npar = round(Int, F * ((X_box^2) * 10^(-8)))

# if no simulation of time
#Np = round(Int, Dose/(1.602*10^(-9)*LET/(pi*(R+Rk)^2*10^(-8))));
#Np = Npar;

#zF = 1.602*10^(-9)*LET/(pi*(R+Rk)^2*10^(-8));
global zF = irrad.dose / Npar
global D = DoseRate_h / zF
global T = irrad.dose / (zF * D) * 3600
# ############Test function###
#Nd is the dimension fo the space, the nucleus of the cell is assumed to be a cylinder. Thde cell a sphere around the centerThe geometry can be more complicated but for now it is fine
Nd = 3;
@time begin
    Np = rand(Poisson(Npar))
    #Np=10;
    #local DOSE_tot = 0. ;
    @everywhere global GYR_tot = zeros(length(arrayOfCell));
    println(Np)
	@everywhere global X = Array{Float64}(undef, 0, Nd);
    @everywhere global Y = Array{Float64}(undef, 0, Nd);
    #@everywhere global cell_origin = Cell(0.0, 0.0, 0.0, r_nucleus, R_cell);
	println(X)
	for i in 1:Np
        local x, y = GenerateHit_BOX(X_box);
        local track = Track(x, y, Rk);
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
                                 
               # track_x_tr = track.x - cell.x;
                #track_y_tr = track.y - cell.y;

            	integral, theta, Gyr, radius= distribute_dose_vector(ion,cell,track);
				#println(integral)
            	local X_, Y_ = calculate_damage(ion, cell, track, integral, theta, Gyr, radius);

                # X_[:, 1] .+= cell.x;
                # X_[:, 2] .+= cell.y;
                # Y_[:, 1] .+= cell.x;
                # Y_[:, 2] .+= cell.y;   

				dist = sqrt.((X_[:, 1] .- cell.x).^2 .+ (X_[:, 2] .-cell.y).^2)
				if size(dist[dist .> cell.r], 1) != 0
					println("Error")
					println(j)
			#return X_i
				end
				#X_=hcat(X_,repeat([j],size(X_,1)))
				#Y_=hcat(Y_,repeat([j],size(Y_,1)))
				#print(size(X_),"\n",size(Y_))
				global X = vcat(X,X_)
            	global Y = vcat(Y,Y_)
			end
			#DOSE_tot+=dose
            GYR_tot[j] += Gyr
        end
    end
# println(DOSE_tot)
println(GYR_tot)
end 

df = DataFrame(X,:auto)
println(df) 
plt = Plot_Lattice_Cells(arrayOfCell)
Plots.scatter!(plt,df.x1,df.x2,df.x3,mode="markers",markersize=1)
# plotlyjs()
display(plt)