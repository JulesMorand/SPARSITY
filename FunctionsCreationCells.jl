######################################################################
function generate_cells_positions_squaredlattice(X_box::Float64,R_cell::Float64)
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
######################################################################
function generate_cells_positions_triangularlattice(X_box::Float64,R_cell::Float64)
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
function Creation_ArrayOfCell(N::Int64, nodes_positions::Vector{Tuple{Float64, Float64}}, r__nucl::Float64, R__cell::Float64)
    arrayOfCells = Array{Cell, 1}(undef, N) 
        for i in 1:N
            x,y=nodes_positions[i]
            ######z=R_cell for this layer of cell
            cell=Cell(x,y,R__cell,r__nucl,R__cell,Array{Float64}(undef, 0, Nd),Array{Float64}(undef, 0, Nd))
            arrayOfCells[i]=cell
        end
    return arrayOfCells
end