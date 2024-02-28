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

####Dimension of the hit box
X_box=500  #1000µmm side size of the square box
r_nucl=10. #10µm
R_cell=30.  #

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
    println(N_CellsSide," ", N_CellsSide2)
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
    return N,nodes_positions
end

#############CHOOSE squared or triangluar Lattice####

#N,nodes_positions=generate_cells_positions_squaredlattice(X_box,R_cell)
N,nodes_positions=generate_cells_positions_triangularlattice(X_box,R_cell)

#############Create Cells Array######
arrayOfCells = Array{Cell3D, 1}(undef, N) 
for i in 1:N
    x,y=nodes_positions[i]
    cell=Cell3D(x,y,R_cell,r_nucl,R_cell)
    arrayOfCells[i]=cell
end

#######Plot the cell lattice############
#Initialise the plot in ploting the first cell
x,y,z,r,R=arrayOfCells[1].x,arrayOfCells[1].y,arrayOfCells[1].z,arrayOfCells[1].r_nucl,arrayOfCells[1].R_cell
X,Y,Z=SphereShape(x,y,z,R)
plt=Plots.surface(
    X, Y, Z, 
    size=(700,700),
    opacity=0.4, 
    color=cgrad(:matter, N, categorical = true)[1],
    legend=false,
    xlims=(0,X_box+R_cell),
    ylims=(0,X_box+R_cell),
    zlims=(0,X_box+R_cell))
Xnc,Ync,Znc=CylinderShape(x,y,z,r,8.)
    Plots.surface!(
    Xnc,Ync,Znc, 
    opacity=1, 
    color=cgrad(:matter, N, categorical = true)[1], legend=false)

for i in 2:length(arrayOfCells)
    x,y,z,r,R=arrayOfCells[i].x,arrayOfCells[i].y,arrayOfCells[i].z,arrayOfCells[i].r_nucl,arrayOfCells[i].R_cell
    X,Y,Z=SphereShape(x,y,z,R)
    Plots.surface!(
    X, Y, Z, 
    opacity=0.4, 
    color=cgrad(:matter, N, categorical = true)[i], legend=false)
    Xnc,Ync,Znc=CylinderShape(x,y,z,r,8.)
    Plots.surface!(
    Xnc,Ync,Znc, 
    opacity=1, 
    color=cgrad(:matter, N, categorical = true)[i], legend=false)

end
display(plt)