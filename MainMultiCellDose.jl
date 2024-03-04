####Parameters 
###Dimension of the hit box
X_box=500  #1000µmm side size of the square box
###Size Cell nucleus
r_nucl=10. #10µm
###Size Cell
R_cell=30. #
using Base.Threads
include("./CreationLatticeCells.jl") 
include("./FunctionVisualisation.jl")
#############CHOOSE squared or triangluar Lattice####
#N,nodes_positions=generate_cells_positions_squaredlattice(X_box,R_cell)
N,nodes_positions=generate_cells_positions_triangularlattice(X_box,R_cell) 
println(N)
#############Create Cells Array######################
arrayOfCell=Creation_ArrayOfCell(N,nodes_positions)
Plot_Lattice_Cells(arrayOfCell)
