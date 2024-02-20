using Plots,WebIO,PlotlyJS,Distributions
plotlyjs()
function SphereShape(x,y,z,r)# (the cell)
    N = 100
    u = range(0, stop=2π, length=N)
    v = range(0, stop=π, length=N)
    X = r.*cos.(u) .* sin.(v)'.+x
    Y = r.*sin.(u) .* sin.(v)'.+y
    Z = repeat(r.*cos.(v)',outer=[N, 1]).+z
   
    return X,Y,Z
end
mutable struct Cell3D
    x::Float64
    y::Float64
    z::Float64
    r::Float64
end
function GenerateCells(r_min,r_max,R_tumor)
    radius = (R_tumor-r_max)*sqrt((rand(Uniform(0,1))));
    theta = 2*π*rand(Uniform(0,1));
    phi=π*rand(Uniform(0,1))
    X = radius.*cos.(theta) .* sin.(phi)'
    Y = radius.*sin.(theta) .* sin.(phi)'
    Z = radius.*cos.(phi)'
    R=r_min+rand(Uniform(0,1))*(r_max-r_min)
    return Cell3D(X,Y,Z,R)
end   
let 
    X,Y,Z=SphereShape(0,0,0,30.)
    plt2=Plots.surface(
        X, Y, Z, 
        size=(600,600),
        opacity=0.3, 
        cbar=:none, 
        legend=false)
    N = 15 
    # create an uninitialized 1D array of size N
    arrayOfCells = Array{Cell3D, 1}(undef, N) 
    for i in 1:N
        arrayOfCells[i]=GenerateCells(7.5,8.5,30.)
    end
   

    for i in 1:length(arrayOfCells)
        x,y,z,r=arrayOfCells[i].x,arrayOfCells[i].y,arrayOfCells[i].z,arrayOfCells[i].r
        X,Y,Z=SphereShape(x,y,z,r)
        Plots.surface!(
        X, Y, Z, 
        size=(600,600),
        opacity=1, 
        color=cgrad(:matter, N, categorical = true)[i], 
        legend=false)
        
    end
display(plt2)
end
