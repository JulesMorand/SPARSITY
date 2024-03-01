using Base.Threads
include("./CreationLatticeCells.jl") 


N,nodes_positions=generate_cells_positions_triangularlattice(X_box,R_cell)
println(N)