
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
#############Plot the cell fro the arrayOfCells ############
function Plot_Lattice_Cells(arrayOfCells)
    #Initialise the plot in ploting the first cell
    local x,y,z,r,R=arrayOfCells[1].x,arrayOfCells[1].y,arrayOfCells[1].z,arrayOfCells[1].r_nucl,arrayOfCells[1].R_cell
    local X,Y,Z=SphereShape(x,y,z,R)
    plt=Plots.surface(
    X, Y, Z, 
    size=(700,700),
    opacity=0.4, 
    color=cgrad(:matter, N, categorical = true)[1],
    legend=false,
    xlims=(0,X_box+R_cell),
    ylims=(0,X_box+R_cell),
    zlims=(0,X_box+R_cell))
    local Xnc,Ync,Znc=CylinderShape(x,y,z,r,8.)
    Plots.surface!(
    Xnc,Ync,Znc, 
    opacity=1, 
    color=cgrad(:matter, N, categorical = true)[1], legend=false)

    for i in 2:length(arrayOfCells)
        local x,y,z,r,R=arrayOfCells[i].x,arrayOfCells[i].y,arrayOfCells[i].z,arrayOfCells[i].r_nucl,arrayOfCells[i].R_cell
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
    return plt
end
 