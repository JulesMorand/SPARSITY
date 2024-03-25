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
@everywhere function Plot_Lattice_Cells(arrayOfCells,X_box)
    plotlyjs()
    #Initialise the plot in ploting the first cell
    local x,y,z,r,R,Dam_X,Dam_Y=arrayOfCells[1].x,arrayOfCells[1].y,arrayOfCells[1].z,arrayOfCells[1].r,arrayOfCells[1].R,arrayOfCells[1].Dam_X,arrayOfCells[1].Dam_Y
    local dfX = DataFrame(Dam_X,:auto)
    local dfY = DataFrame(Dam_Y,:auto)

    #local X,Y,Z=SphereShape(x,y,z,R)
    #println(length(Dam_Y)==0)
    if length(Dam_Y)== 0
        local X,Y,Z=SphereShape(x,y,z,R)
        plt=Plots.surface(
        X, Y, Z,
        opacity=0.3,  
        #color=cgrad(:matter, N, categorical = true)[i], 
        color= :red,
        legend=false)
        local Xnc,Ync,Znc=CylinderShape(x,y,z,r,r)
        Plots.surface!(
        Xnc,Ync,Znc, 
        opacity=0.9, 
        #color=cgrad(:matter, N, categorical = true)[i], 
        color= :red,
        legend=false,
        xlim = (0, X_box+30),
        ylim = (0, X_box+30),
        zlim = (0, X_box+30))
    else
        local X,Y,Z=SphereShape(x,y,z,R)
        plt=Plots.surface(
        X, Y, Z, 
        opacity=0.3, 
        #color=cgrad(:matter, N, categorical = true)[i], 
        color= :black,
        legend=false)
        local Xnc,Ync,Znc=CylinderShape(x,y,z,r,r)
        Plots.surface!(
        Xnc,Ync,Znc, 
        opacity=0.9, 
        #color=cgrad(:matter, N, categorical = true)[i], 
        color= :black,
        legend=false,
        xlim = (0, X_box+30),
        ylim = (0, X_box+30),
        zlim = (0, X_box+30))
    end
    Plots.scatter!(plt,dfX.x1,dfX.x2,dfX.x3,mode="markers",markersize=0.5 , color= :purple ) 
    Plots.scatter!(plt,dfY.x1,dfY.x2,dfY.x3,mode="markers",markersize=3 , color= :purple ) 
    for i in 2:length(arrayOfCells)
        local x,y,z,r,R,Dam_X,Dam_Y=arrayOfCells[i].x,arrayOfCells[i].y,arrayOfCells[i].z,arrayOfCells[i].r,arrayOfCells[i].R,arrayOfCells[i].Dam_X,arrayOfCells[i].Dam_Y
        local dfX = DataFrame(Dam_X,:auto)
        local dfY = DataFrame(Dam_Y,:auto)
        if length(Dam_Y)==0
            local X,Y,Z=SphereShape(x,y,z,R)
            Plots.surface!(
            X, Y, Z, 
            opacity=0.3,  
            #color=cgrad(:matter, N, categorical = true)[i], 
            color= :red,
            legend=false)
            local Xnc,Ync,Znc=CylinderShape(x,y,z,r,r)
            Plots.surface!(
            Xnc,Ync,Znc, 
            opacity=0.9, 
            #color=cgrad(:matter, N, categorical = true)[i], 
            color= :red,
            legend=false)
        else
            local X,Y,Z=SphereShape(x,y,z,R)
            Plots.surface!(
            X, Y, Z, 
            opacity=0.3, 
            #color=cgrad(:matter, N, categorical = true)[i], 
            color= :black,
            legend=false)
            local Xnc,Ync,Znc=CylinderShape(x,y,z,r,r)
            Plots.surface!(
            Xnc,Ync,Znc, 
            opacity=0.9, 
            #color=cgrad(:matter, N, categorical = true)[i], 
            color= :black,
            legend=false)
        end
        Plots.scatter!(plt,dfX.x1,dfX.x2,dfX.x3,mode="markers",markersize=0.5 , color= :green ) 
        Plots.scatter!(plt,dfY.x1,dfY.x2,dfY.x3,mode="markers",markersize=3 , color= :purple ) 
    end
    plot!(size=(800,800))
    return plt
end
 