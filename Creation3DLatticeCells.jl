####Dimension of the hit box
X_box=250  #1000µmm side size of the square box
r_nucl=10. #10µm
R_cell=30. #


using Plots,WebIO,PlotlyJS,Distributions
plotlyjs()
function SphereShape(x,y,z,r)# (the cell)
    N = 50
    u = range(0, stop=2π, length=N)
    v = range(0, stop=π, length=N)
    X = r.*cos.(u) .* sin.(v)'.+x
    Y = r.*sin.(u) .* sin.(v)'.+y
    Z = repeat(r.*cos.(v)',outer=[N, 1]).+z
   
    return X,Y,Z
end

function CylinderShape(x,y,z,r,h)# (the nucleus)
    r_cyl = r
    h_cyl  = h
    m_cyl , n_cyl  =50, 50
    u  = range(0, 2pi, length=n_cyl )
    v = range(z-h_cyl/2.,z+h_cyl/2., length=m_cyl )
    us = ones(m_cyl).*transpose(u)
    vs = v.*transpose(ones(n_cyl))
    #Surface parameterization
    X  = r_cyl*cos.(us).+x
    Y  = r_cyl*sin.(us).+y
    Z  = vs
    return X,Y,Z
end
    
mutable struct Cell3D
    x::Float64
    y::Float64
    z::Float64
    r_nucl::Float64
    R_cell::Float64
end


function generate_cells_positions_squaredlattice(X_box::Int,R_cell::Float64)
    N_CellsSide=convert(Int64,floor(X_box/(2*R_cell)))
    nodes_positions = Vector{Tuple{Float64, Float64}}()
    for i in 1:N_CellsSide
        for j in 1:N_CellsSide
            push!(nodes_positions, (R_cell+(i-1)*2*R_cell,R_cell+(j-1)*2*R_cell))
                
        end
    end
    local N=N_CellsSide^2
    return N,nodes_positions
end
function generate_cells_positions_triangularlattice(X_box::Int,R_cell::Float64)
    N_CellsSide=convert(Int64,floor(X_box/(2*R_cell)))
    N_CellsSide2=convert(Int64,floor((X_box)/(R_cell*sqrt(3))))
    println("\n Nx=", N_CellsSide,"\n Ny=", N_CellsSide2)
    nodes_positions = Vector{Tuple{Float64, Float64}}()
    for i in 1:N_CellsSide
        for j in 1:N_CellsSide2
            if rem(j,2)==1
                push!(nodes_positions, (R_cell+(i-1)*2*R_cell,R_cell+(j-1)*R_cell*sqrt(3)))
            else
                push!(nodes_positions, ((i)*2*R_cell,R_cell+(j-1)*R_cell*sqrt(3)))
            end
        end
    end
    local N=N_CellsSide*N_CellsSide2
    return N,N_CellsSide2,nodes_positions
end

#############CHOOSE squared or triangluar Lattice####
########To chech!! height betwwen layer###
#N,nodes_positions=generate_cells_positions_squaredlattice(X_box,R_cell)
N,N_CellsSide2,nodes_positions=generate_cells_positions_triangularlattice(X_box,R_cell)
N_CellsSide3=convert(Int64,floor((X_box)/(R_cell*sqrt(3/2))))
println("N_CellsSide3=",N_CellsSide3)
Ntot=N*N_CellsSide3;
println("N=",N)
println("Ntot=",Ntot)
#############Create Cells Array######
arrayOfCells = Array{Cell3D, 1}(undef, Ntot)
println(length(arrayOfCells)) 
for k in 1:N_CellsSide3
    for i in 1:N
        local x,y=nodes_positions[i]
        if rem(k,2)==1
            cell=Cell3D(x,y,R_cell+sqrt(3/2)*R_cell*(k-1),r_nucl,R_cell)
        else
            cell=Cell3D(x+R_cell,y+R_cell*sqrt(3)/2,sqrt(3/2)*R_cell*(k),r_nucl,R_cell)
        end
        arrayOfCells[i+(k-1)*N]=cell
       # println(i+(k-1)*N)
    end
end
#######Plot the cell lattice############
#Initialise the plot in ploting the first cell
let
    x,y,z,r,R=arrayOfCells[1].x,arrayOfCells[1].y,arrayOfCells[1].z,arrayOfCells[1].r_nucl,arrayOfCells[1].R_cell
    X,Y,Z=SphereShape(x,y,z,R)
    plt=Plots.surface(
        X, Y, Z, 
        size=(700,700),
        opacity=0.6, 
        color=cgrad(:matter, Ntot, categorical = true)[1],
        legend=false,
        xlims=(0,X_box+R_cell+R_cell*sqrt(3/2)),
        ylims=(0,X_box+R_cell+R_cell*sqrt(3/2)),
        zlims=(0,X_box+R_cell+R_cell*sqrt(3/2)))
    Xnc,Ync,Znc=CylinderShape(x,y,z,r,8.)
        Plots.surface!(
        Xnc,Ync,Znc, 
        opacity=1, 
        color=cgrad(:matter, Ntot, categorical = true)[1],
        legend=false,
        aspect_ratio=1.0)

    for i in 2:Ntot
        println(i)
        x,y,z,r,R=arrayOfCells[i].x,arrayOfCells[i].y,arrayOfCells[i].z,arrayOfCells[i].r_nucl,arrayOfCells[i].R_cell
        X,Y,Z=SphereShape(x,y,z,R) 
        Plots.surface!(
        X, Y, Z, 
        opacity=0.6, 
        color=cgrad(:matter, Ntot, categorical = true)[i], legend=false)
        Xnc,Ync,Znc=CylinderShape(x,y,z,r,8.)
        Plots.surface!(
        Xnc,Ync,Znc, 
        opacity=1, 
        color=cgrad(:matter, Ntot, categorical = true)[i], legend=false)

    end
    display(plt)
end