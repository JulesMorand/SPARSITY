using Distributions
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
#############Create Cells Array################
function Creation_ArrayOfCell(N,nodes_positions)
    arrayOfCells = Array{Cell3D, 1}(undef, N) 
        for i in 1:N
            x,y=nodes_positions[i]
            cell=Cell3D(x,y,R_cell,r_nucl,R_cell)
            arrayOfCells[i]=cell
        end
    return arrayOfCells
end