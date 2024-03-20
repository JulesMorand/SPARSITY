@everywhere function ATRadius(ion::Ion, irrad::Irrad)

    Rc = 0.01;

    Rp = 0.05*((ion.E/ion.A)^(1.7));

    if irrad.kR < 1.0
        Rk = Rc*exp((irrad.kR*(1+2*log(Rp/Rc))-1)/2);
    else
        Rk = Rp;
    end

    return Rc, Rp, Rk

end
#Generate uniform Hit in a box from 0 to X_box in x and y
@everywhere function GenerateHit_BOX(X_box::Float64)
    x0 = rand(Uniform(0,1))*X_box;
    y0 = rand(Uniform(0,1))*X_box;
    return x0, y0
end
#Generate uniform Hit in a circle center in 0 0 from with radius cell.r+Rk
@everywhere function GenerateHit(cell::Cell, Rk::Float64)
    radius = (cell.r+Rk)*sqrt((rand(Uniform(0,1))));
    theta = 2*pi*rand(Uniform(0,1));
    x0 = radius*cos(theta);
    y0 = radius*sin(theta);
    return x0, y0
end
#######################
@everywhere function GetRadialLinearDose(r::Float64, ion::Ion)
    #LET normalized to Rk ???
    LETk = ion.LET*0.1602;####JM: why thi is usefull?
    D_arc=0.
    if r <= Rc
        D_arc = (1/(pi*Rc^2))*(LETk/(1*(1+2*log(Rk/Rc))));
    elseif r <= Rk
        D_arc = (1/(pi*r*r))*(LETk/(1*(1+2*log(Rk/Rc))));
    end
    return D_arc
end

####Original functions

@everywhere function arc_intersection_angle(r::Float64, b::Float64, r_nucleus::Float64)
    if b < r_nucleus
        if r <= r_nucleus - b
            return 2 * π
        elseif r < b + r_nucleus
            arg = b/(2*r)  + r/(2*b)  - r_nucleus * r_nucleus / (2 * b * r)
            arg = min(arg, 1.0) #not sure if usefull?
            arg = max(arg, -1.0)
            return 2 * acos(arg)
        end
    else #(b>=r_nucleus)
        if r <= b - r_nucleus
            return 0.0
        elseif r < b + r_nucleus
            arg =b/(2*r)  + r/(2*b)  - r_nucleus * r_nucleus / (2 * b * r)
            arg = min(arg, 1.0) #not sure if usefull?
            arg = max(arg, -1.0)
            return 2 * acos(arg)
        end
    end
    return 0.0
end
##############
@everywhere function integrate_weighted_radial_track( rMin::Float64, rMax::Float64, b::Float64, r_nucleus::Float64, step::Int64)
    r1, r2, log_r2, log_rMin, log_rMax, log_step = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    f1, f2, f, arc_weight1, arc_weight2 = 0.0, 0.0, 0.0, 0.0, 0.0
    integral = 0.0
    if rMin > 0
        log_rMin = log10(rMin)
    else
        log_rMin = -5.
    end
    log_rMax = log10(rMax)
    nSteps = step
    log_step = (log_rMax - log_rMin) / nSteps
    #println(log_step)
    if nSteps < 3
        log_step = (log_rMax - log_rMin) / 3
        nSteps = 3
    end
    area = 0.0
    arc_weight2 = arc_intersection_angle(rMin, b, r_nucleus)
    f2 = GetRadialLinearDose(rMin, ion) * rMin * arc_weight2
    r2 = rMin
    for i in 1:nSteps - 1
        log_r2 = log_rMin + log_step * (i + 1)
        f1 = f2
        r1 = r2
        arc_weight1 = arc_weight2
        #println(r2 ,"\n")
        r2 = 10^log_r2
        arc_weight2 = arc_intersection_angle(r2, b, r_nucleus)
        f2 = GetRadialLinearDose(r2, ion) * r2 * arc_weight2
        f = (r2 - r1) * (f1 / 2.0 + f2 / 2.0)
        integral += f
        area += (r2 - r1) * (arc_weight1 * r1 / 2.0 + arc_weight2 * r2 / 2.0)
    end
    ####Pass 1<-2
    f1 = f2
    r1 = r2
    arc_weight1 = arc_weight2
    #### Fill 2
    r2 = rMax
    arc_weight2 = arc_intersection_angle(r2, b, r_nucleus)
    f2 = GetRadialLinearDose(r2, ion) * r2 * arc_weight2
    f = (r2 - r1) * (f1 / 2.0 + f2 / 2.0)
    integral += f
    area += (r2 - r1) * (arc_weight1 * r1 / 2.0 + arc_weight2 * r2 / 2.0)
    return area, integral, integral/area
end

@everywhere function distribute_dose(cell::Cell, track::Track)
    x_track, y_track = track.x,track.y
    x_track = (x_track - cell.x) #* 1e3  # mm -> um ??
    y_track = (y_track - cell.y) #* 1e3  # mm -> um
    b = sqrt(x_track^2 + y_track^2)

    rMax = min(Rk, b + cell.r)

    area1 = area2 = area3 = 0.0

    if b <= cell.r
        #nucleus.inNucleusCount += 1
        rMin = 0.0

        if b + track.Rk < cell.r
            r_intersection = track.Rk
        else
            r_intersection = cell.r - b
        end

        area1 = π * r_intersection^2
        area,integral,Gyr=integrate_weighted_radial_track(0., r_intersection, b, cell.r,1000);
        dose=integral
        # dose = track.getRadialIntegral(0.0, r_intersection) * area1

        if rMax > r_intersection
            area,integral,Gyr=integrate_weighted_radial_track(r_intersection, rMax, b, cell.r,1000);
            area2 = area #track, r_intersection, rMax, b, 0.01
            dose += integral
        end

        if rMax == track.Rk
            if track.Rk > cell.r - b
               theta1 = acos((b/(2*rMax) + rMax/(2*b) - (cell.r^2) / (2 * b * rMax)))
               theta2 = acos((b/(2*cell.r) - (rMax*rMax)/(2*b*cell.r) + (cell.r) / (2 * b)))
               area3 = π * cell.r^2 - (theta1 * track.Rk^2 + theta2 * cell.r^2 - track.Rk * b * sin(theta1))
            else
                area3 = π * (cell.r^2 - r_intersection^2)
            end
        end

        #dose /= area1 + area2 + area3
        Gyr=dose/(area1+area2+ area3)
       # println("area1: ",area1," area2: ",area2," area3: ",area3,"\n1 sumareas:",area1+area2+area3," aera cell:",cell.r*cell.r*π, "\n",area1+area2+area3-cell.r*cell.r*π)

        #nucleus.totalNucleusDose += dose
       # push!(nucleus.doses, dose)
        #push!(nucleus.times, track.getTime())

    elseif b <= cell.r + track.Rk
        #nucleus.intersectionCount += 1
        rMin = b - cell.r
        area, integral, Gyr = integrate_weighted_radial_track(rMin, rMax, b, cell.r,1000)
        dose = integral
        area2=area
        if rMax == track.Rk
            theta1 = acos((b/(2*rMax) + rMax/(2*b) - (cell.r^2) / (2 * b * rMax)))
            theta2 = acos((b/(2*cell.r) - (rMax*rMax)/(2*b*cell.r) + (cell.r) / (2 * b)))
            area3 = π * cell.r^2 - (theta1 * track.Rk^2 + theta2 * cell.r^2 - track.Rk * b * sin(theta1))
        end

        Gyr = dose/(area2 + area3)

        #nucleus.totalNucleusDose += dose
        #return dose, area1+area2+area3, Gyr
       # push!(cell.doses, dose)
        #push!(nucleus.times, track.getTime())
    end
    return dose, area1+area2+area3, Gyr
end



#Added function as NameFunction_vector that return vectors instead of numbers as before. 
#They are used to caluclate the spatal distribution of damages that is proportional to the dose

@everywhere function arc_intersection_angle_vector(r::Float64, b::Float64, r_nucleus::Float64)
    if b < r_nucleus
        if r <= r_nucleus - b
            return 2 * π, 0
        elseif r < b + r_nucleus
            arg = b/(2*r)  + r/(2*b)  - r_nucleus * r_nucleus / (2 * b * r)
            arg = min(arg, 1.0) #not sure if usefull?
            arg = max(arg, -1.0)
            return 2 * acos(arg), arg*r_nucleus
        end
    else #(b>=r_nucleus)
        if r <= b - r_nucleus
            return 0.0, 0.0
        elseif r < b + r_nucleus
            arg =b/(2*r)  + r/(2*b)  - r_nucleus * r_nucleus / (2 * b * r)
            arg = min(arg, 1.0) #not sure if usefull?
            arg = max(arg, -1.0)
            return 2 * acos(arg), arg*r_nucleus
        end
    end
    return 0.0, 0.0
end
##############

@everywhere function integrate_weighted_radial_track_vector(ion::Ion, rMin::Float64, rMax::Float64, b::Float64, r_nucleus::Float64, nSteps::Int64)
    local r1, r2, log_r2, log_rMin, log_rMax, log_step = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    local f1, f2, f, arc_weight1, arc_weight2 = 0.0, 0.0, 0.0, 0.0, 0.0
    local integral = Array{Float64}(undef,0)
    local theta = Array{Float64}(undef,0)
    local radius = Array{Float64}(undef,0)
    if rMin > 0
        log_rMin = log10(rMin)
    else
        log_rMin = -5.
    end
    log_rMax = log10(rMax)
    #nSteps = step
    log_step = (log_rMax - log_rMin) / nSteps
    #println(log_step)
    if nSteps < 3
        log_step = (log_rMax - log_rMin) / 3
        nSteps = 3
    end
    area = 0.0
    arc_weight2 = arc_intersection_angle(rMin, b, r_nucleus)
    #push!(theta,arc_weight2/2)
    f2 = GetRadialLinearDose(rMin, ion) * rMin * arc_weight2
    r2 = rMin
    for i in 1:nSteps - 1
        log_r2 = log_rMin + log_step * (i + 1)
        f1 = f2
        r1 = r2
        arc_weight1 = arc_weight2
        #println(r2 ,"\n")
        r2 = 10^log_r2
        arc_weight2 = arc_intersection_angle(r2, b, r_nucleus)
        push!(theta, arc_weight1)
        f2 = GetRadialLinearDose(r2, ion) * r2 * arc_weight2
        f = (r2 - r1) * (f1 / 2.0 + f2 / 2.0)
        push!(integral,f)
        area += (r2 - r1) * (arc_weight1 * r1 / 2.0 + arc_weight2 * r2 / 2.0)
        push!(radius,r1)
    end
    ####Pass 1<-2
    f1 = f2
    r1 = r2
    arc_weight1 = arc_weight2
    #### Fill 2
    r2 = rMax
    arc_weight2= arc_intersection_angle(r2, b, r_nucleus)
    push!(theta, arc_weight1)
    push!(theta, arc_weight2)
    f2 = GetRadialLinearDose(r2, ion) * r2 * arc_weight2
    f = (r2 - r1) * (f1 / 2.0 + f2 / 2.0)
    push!(integral,f)
    area += (r2 - r1) * (arc_weight1 * r1 / 2.0 + arc_weight2 * r2 / 2.0)
    push!(radius,r1)
    push!(radius,r2)
    return integral, theta, area, radius
end

#####################################

@everywhere function distribute_dose_vector(ion::Ion, cell::Cell, track::Track)
    local x_track, y_track = track.x,track.y
    local x_track = (x_track - cell.x) #* 1e3  # mm -> um ??
    local y_track = (y_track - cell.y) #* 1e3  # mm -> um
    local b = sqrt(x_track^2 + y_track^2)

    local rMax = min(Rk, b + cell.r)

    local area1 = area2 = area3 = 0.0

    if b <= cell.r
        #nucleus.inNucleusCount += 1
        local rMin = 0.0

        if b + track.Rk < cell.r
            local r_intersection = track.Rk
        else
            local r_intersection = cell.r - b
        end

        local area1 = π * r_intersection^2
        local integral1, theta11, area1, radius1 = integrate_weighted_radial_track_vector(ion, 0., r_intersection, b, cell.r, 1000);
        local dose = sum(integral1)
        # dose = track.getRadialIntegral(0.0, r_intersection) * area1

        if rMax > r_intersection
            local integral2, theta21, area, radius2 = integrate_weighted_radial_track_vector(ion, r_intersection, rMax, b, cell.r,1000);
            local area2 = area #track, r_intersection, rMax, b, 0.01
            local dose += sum(integral2)
            local integral = [integral1; integral2];
            local theta = [theta11; theta21];
            local radius = [radius1; radius2];
        else
            local integral = integral1;
            local theta = theta11;
            local radius = radius1;
        end

        if rMax == track.Rk
            if track.Rk > cell.r - b
                local theta1 = acos((b/(2*rMax) + rMax/(2*b) - (cell.r^2) / (2 * b * rMax)))
                local theta2 = acos((b/(2*cell.r) - (rMax*rMax)/(2*b*cell.r) + (cell.r) / (2 * b)))
                local area3 = π * cell.r^2 - (theta1 * track.Rk^2 + theta2 * cell.r^2 - track.Rk * b * sin(theta1))
            else
                local area3 = π * (cell.r^2 - r_intersection^2)
            end
        end

        #dose /= area1 + area2 + area3
        local Gyr = dose/(area1+area2+ area3)
       # println("area1: ",area1," area2: ",area2," area3: ",area3,"\n1 sumareas:",area1+area2+area3," aera cell:",cell.r*cell.r*π, "\n",area1+area2+area3-cell.r*cell.r*π)

        #nucleus.totalNucleusDose += dose
       # push!(nucleus.doses, dose)
        #push!(nucleus.times, track.getTime())

    elseif b <= cell.r + track.Rk
        #nucleus.intersectionCount += 1
        local rMin = b - cell.r
        local integral, theta, area, radius = integrate_weighted_radial_track_vector(ion, rMin, rMax, b, cell.r, 1000)
        local dose = sum(integral)
        local area2=area
        if rMax == track.Rk
            local theta1 = acos((b/(2*rMax) + rMax/(2*b) - (cell.r^2) / (2 * b * rMax)))
            local theta2 = acos((b/(2*cell.r) - (rMax*rMax)/(2*b*cell.r) + (cell.r) / (2 * b)))
            local area3 = π * cell.r^2 - (theta1 * track.Rk^2 + theta2 * cell.r^2 - track.Rk * b * sin(theta1))
        end

        local Gyr = dose/(area2 + area3)

        #nucleus.totalNucleusDose += dose
        #return dose, area1+area2+area3, Gyr
       # push!(cell.doses, dose)
        #push!(nucleus.times, track.getTime())
    end

    local integral[integral .< 0] .= 0
    local theta = minimum([theta[1:end-1]./2 theta[2:end]./2], dims = 2)

    return integral, theta, Gyr, radius
end

#function to calculate the average yield of damage per unit Gy
@everywhere function calculate_kappa(ion::Ion)

    if ion.ion == "12C"
        p1 = 6.8;
        p2 = 0.156;
        p3 = 0.9214;
        p4 = 0.005245;
        p5 = 1.395;
    elseif ion.ion == "4He"
        p1 = 6.8;
        p2 = 0.1471;
        p3 = 1.038;
        p4 = 0.006239;
        p5 = 1.582;
    elseif ion.ion == "3He"
        p1 = 6.8;
        p2 = 0.1471;
        p3 = 1.038;
        p4 = 0.006239;
        p5 = 1.582;
    elseif ion.ion == "1H"
        p1 = 6.8;
        p2 = 0.1773;
        p3 = 0.9314;
        p4 = 0;
        p5 = 1;
    elseif ion.ion == "2H"
        p1 = 6.8;
        p2 = 0.1773;
        p3 = 0.9314;
        p4 = 0;
        p5 = 1;
    elseif ion.ion == "16O"
        p1 = 6.8;
        p2 = 0.1749;
        p3 = 0.8722;
        p4 = 0.004987;
        p5 = 1.347;
    else
        println("Unknown ion specie")
        return -1
    end

    yield = (p1 + (p2*ion.LET)^p3)/(1+ (p4*ion.LET)^p5);

    return yield
end

#function to calculate the spatial position of the damages
#for each track calculate the dose on the cell and ditribute the damages around the cell nucleus

#@everywhere function calculate_damage(ion::Ion, integral::Vector{Float64}, theta::Matrix{Float64}, Gyr::Float64)

@everywhere function calculate_damage(ion_::Ion, cell_::Cell, integral_, theta_, Gyr_, radius_, x_, y_)
    
    #theta__ = [theta_[1:end-1]./2 theta_[2:end]./2]
    #theta_ = minimum(theta__, dims = 2)
    local b = sqrt(x_*x_ + y_*y_)
    
    local kappa_DSB = 9*calculate_kappa(ion_);
    local lambda_DSB = kappa_DSB*10^-3;
    
    local x0d = rand(Poisson(kappa_DSB*Gyr_))
    local y0d = rand(Poisson(lambda_DSB*Gyr_))
    
    if (x0d == 0) & (y0d == 0)
        local X_CD = Array{Float64}(undef, 0, Nd);
        local Y_CD = Array{Float64}(undef, 0, Nd);

        return X_CD, Y_CD
    end
    
    if x0d > 0
        local X_CD = Array{Float64}(undef, 0, Nd);
        #local X_CD = zeros(x0d, Nd);

        local radius__xP = rand(Categorical((integral_/sum(integral_))), x0d)
        local radius__x = (radius_[radius__xP .+ 1] - radius_[radius__xP]).*rand(Uniform(0,1),x0d) .+ radius_[radius__xP];
        if (x_ >= 0)
            local theta__x = 3*π/2 .- acos.(y_/b) .+ theta_[radius__xP].*rand(Uniform(0,1),x0d).*[-1 ,1][rand(Bernoulli(),x0d) .+ 1];
        elseif (x_ < 0)
            local theta__x = 3*π/2 .+ acos.(y_/b) .+ theta_[radius__xP].*rand(Uniform(0,1),x0d).*[-1 ,1][rand(Bernoulli(),x0d) .+ 1];
        end
        local Xx = radius__x.*cos.(theta__x) .+ x_;
        local Xy = radius__x.*sin.(theta__x) .+ y_;
        for i in 1:x0d
            X_CD = vcat(X_CD,reshape([Xx[i], Xy[i], cell_.R-cell_.r/2+cell_.r*rand(Uniform(0,1),1)[1]], 1, :));
        end
    else
        local X_CD = Array{Float64}(undef, 0, Nd);
    end
    
    if y0d > 0
        local Y_CD = Array{Float64}(undef, 0, Nd);
        #local Y_CD = zeros(y0d, Nd);

        local radius__yP = rand(Categorical((integral_/sum(integral_))), y0d)
        local radius__y = (radius_[radius__yP .+ 1] - radius_[radius__yP]).*rand(Uniform(0,1),y0d) .+ radius_[radius__yP];
        if (x_ >= 0) 
            local theta__y = 3*π/2 .- acos.(y_/b) .+ theta_[radius__yP].*rand(Uniform(0,1),y0d).*[-1 ,1][rand(Bernoulli(),y0d) .+ 1];
        elseif (x_ < 0) 
            local theta__y = 3*π/2 .+ acos.(y_/b) .+ theta_[radius__yP].*rand(Uniform(0,1),y0d).*[-1 ,1][rand(Bernoulli(),y0d) .+ 1];  
        end
        local Yx = radius__y.*cos.(theta__y) .+ x_;
        local Yy = radius__y.*sin.(theta__y) .+ y_;
        for i in 1:x0d
            #global Y_CD = vcat(Y_CD,reshape([Yx[i], Yy[i], cell_.r*rand(Uniform(0,1),1)[1]], 1, :));
            Y_CD = vcat(Y_CD,reshape([Yx[i], Yy[i], cell_.R-cell_.r/2+cell_.r*rand(Uniform(0,1),1)[1]], 1, :));
        end
    else
        local Y_CD = Array{Float64}(undef, 0, Nd);
    end
    
    return X_CD, Y_CD
end

#@everywhere function calculate_damage_NT(ion::Ion, integral::Float64, theta::Float64, Gyr::Float64, NT::Bool)

@everywhere function calculate_damage_NT(ion::Ion, cell::Cell, integral, theta, Gyr, NT::Bool)
    
    X = Array{Float64}(undef, 0, Nd);
    Y = Array{Float64}(undef, 0, Nd);

    theta_ = [theta[1:end-1]./2 theta[2:end]./2]
    theta = minimum(theta_, dims = 2)

    b = sqrt(x*x + y*y)
    
    if NT == false
        kappa_DSB = 9*calculate_kappa(ion);
        lambda_DSB = kappa_DSB*10^-3;
    else
        kappa_DSB = 0.5*9*calculate_kappa(ion);
        lambda_DSB = kappa_DSB*10^-3;
    end
 

    x0d = rand(Poisson(kappa_DSB*Gyr))
    y0d = rand(Poisson(lambda_DSB*Gyr))

    if (x0d == 0) & (y0d == 0)
        return X, Y
    end

    if x0d > 0
        radius_xP = rand(Categorical((integral/sum(integral))), x0d)
        radius_x = (radius[radius_xP .+ 1] - radius[radius_xP]).*rand(Uniform(0,1),x0d) .+ radius[radius_xP];
        if (x >= 0)
            theta_x = 3*π/2 .- acos.(y/b) .+ theta[radius_xP].*rand(Uniform(0,1),x0d).*[-1 ,1][rand(Bernoulli(),x0d) .+ 1];
        elseif (x < 0)
            theta_x = 3*π/2 .+ acos.(y/b) .+ theta[radius_xP].*rand(Uniform(0,1),x0d).*[-1 ,1][rand(Bernoulli(),x0d) .+ 1];
        end
        Xx = radius_x.*cos.(theta_x) .+ x;
        Xy = radius_x.*sin.(theta_x) .+ y;
        for i in 1:x0d
            X = vcat(X,reshape([Xx[i], Xy[i], cell.r*rand(Uniform(0,1),1)[1]], 1, :));
        end
    end
    if y0d > 0
        radius_yP = rand(Categorical((integral/sum(integral))), y0d)
        radius_y = (radius[radius_yP .+ 1] - radius[radius_yP]).*rand(Uniform(0,1),y0d) .+ radius[radius_yP];
        if (x >= 0) 
            theta_y = 3*π/2 .- acos.(y/b) .+ theta[radius_yP].*rand(Uniform(0,1),y0d).*[-1 ,1][rand(Bernoulli(),y0d) .+ 1];
        elseif (x < 0) 
            theta_y = 3*π/2 .+ acos.(y/b) .+ theta[radius_yP].*rand(Uniform(0,1),y0d).*[-1 ,1][rand(Bernoulli(),y0d) .+ 1];  
        end
        Yx = radius_y.*cos.(theta_y) .+ x;
        Yy = radius_y.*sin.(theta_y) .+ y;
        for i in 1:x0d
            Y = vcat(Y,reshape([Yx[i], Yy[i], cell.r*rand(Uniform(0,1),1)[1]], 1, :));
        end
    end

    return X, Y
end

@everywhere function spatial_GSM2(X::Matrix{Float64}, Y::Matrix{Float64}, gsm2::GSM2)
    while size(X)[1] != 0
    
        init_size = size(X)[1];

        Tr = [0,0,0];
        Rr = [rand(Uniform(0,1)), rand(Uniform(0,1)), rand(Uniform(0,1))];
        Sr = log.(1 ./ Rr);
        time_all = Array{Float64}(undef, 0);
        time = 0.0;

        aX = gsm2.a*size(X)[1];
        dist = pairwise(Euclidean(), transpose(X));
        dist[dist .> gsm2.rd] .= 0;
        bXM = gsm2.b.*dist;
        bX = sum(bXM)/2;
        rX = gsm2.r*size(X)[1];

        tau = (Sr .- Tr)./[rX, aX, bX];
        taus = minimum(tau);

        Tr += taus.*[rX, aX, bX];
        time += taus;
        append!(time_all,time);

        if argmin(tau) == 1
            println("r")

            X = X[1:end .!= rand(DiscreteUniform(1,size(X)[1])), :];
        elseif argmin(tau) == 2
            println("a")

            rem = rand(DiscreteUniform(1,size(X)[1]));
            Y = vcat(Y,reshape(X[rem,:], 1, :));
            X = X[1:end .!= rem, :];

        elseif argmin(tau) == 3
            println("b")

            prob_matrix = bXM./sum(bXM);
            selected_index = rand(Categorical(vec(prob_matrix)), 1)[1];
            row = (selected_index ÷ size(bXM)[1]) + 1;
            cols = ifelse(selected_index % size(bXM)[1] == 0, size(bXM)[1], selected_index % size(bXM)[1]);

            Y = vcat(Y,reshape((X[cols,:]+X[row,:])/2, 1, :));
            X = X[setdiff(1:end, (cols, row)), :];
        end
        Sr .+= log.(1/rand(Uniform(0,1)));


        if init_size == size(X)[1]
            println("Error")
            return -1
        end
        if size(Y)[1] > 0
            println("Dead")
            return 0
        end
    end

    return 1
end

@everywhere function spatial_GSM2_fast(X::Matrix{Float64}, Y::Matrix{Float64}, gsm2::GSM2)
    
    dist = pairwise(Euclidean(), transpose(X));
    dist[dist .> gsm2.rd] .= 0;
    p = 1;
    if size(Y)[1] > 0
        return 0
    else
        for i in 1:size(X)[1]
            p *= gsm2.r/(gsm2.r+gsm2.a+gsm2.b*sum(dist[i,:]))
        end
    end

    return p
end

@everywhere function spatial_NTGSM2(X::Matrix{Float64}, Y::Matrix{Float64}, N::Matrix{Float64}, M::Matrix{Float64}, ntgsm2::NTGSM2)   
    while maximum([size(X)[1],size(N)[1]]) != 0
    
        init_size = size(X)[1];

        Tr = [0,0,0,0,0,0];
        Rr = [rand(Uniform(0,1)), rand(Uniform(0,1)), rand(Uniform(0,1)), rand(Uniform(0,1)), rand(Uniform(0,1)), rand(Uniform(0,1))];
        Sr = log.(1 ./ Rr);
        time_all = Array{Float64}(undef, 0);
        time = 0.0;

        aX = ntgsm2.a*size(X)[1];
        dist = pairwise(Euclidean(), transpose(X));
        dist[dist .> ntgsm2.rd] .= 0;
        bXM = ntgsm2.b.*dist;
        bX = sum(bXM)/2;
        rX = ntgsm2.r*size(X)[1];

        aN = ntgsm2.aM*size(N)[1];
        distN = pairwise(Euclidean(), transpose(N));
        distN[distN .> ntgsm2.rdM] .= 0;
        bXMN = bM.*distN;
        bN = sum(bXMN)/2;
        rN = ntgsm2.rM*size(N)[1];

        tau = (Sr .- Tr)./[rX, aX, bX, rN, aN, bN];
        taus = minimum(tau);

        Tr += taus.*[rX, aX, bX, rN, aN, bN];
        time += taus;
        append!(time_all,time);

        if argmin(tau) == 1
            println("r")

            X = X[1:end .!= rand(DiscreteUniform(1,size(X)[1])), :];

        elseif argmin(tau) == 2  
            println("a")

            psim = rand(Uniform(0,1));

            if psim < ntgsm2.pa
                X = X[1:end .!= rand(DiscreteUniform(1,size(X)[1])), :];
            else
                rem = rand(DiscreteUniform(1,size(X)[1]));
                N = vcat(N,reshape(X[rem,:], 1, :));
                X = X[1:end .!= rem, :];
            end

        elseif argmin(tau) == 3 
            println("b")

            prob_matrix = bXM./sum(bXM);
            selected_index = rand(Categorical(vec(prob_matrix)), 1)[1];
            row = (selected_index ÷ size(bXM)[1]) + 1;
            cols = ifelse(selected_index % size(bXM)[1] == 0, size(bXM)[1], selected_index % size(bXM)[1]);

            psim = rand(Uniform(0,1));

            if psim < ntgsm2.qa
                Y = vcat(Y,reshape((X[cols,:]+X[row,:])/2, 1, :));
                X = X[setdiff(1:end, (cols, row)), :];
            else
                N = vcat(N,reshape((X[cols,:]+X[row,:])/2, 1, :));
                X = X[setdiff(1:end, (cols, row)), :];
            end

        elseif argmin(tau) == 4 
            println("rM")

            N = N[1:end .!= rand(DiscreteUniform(1,size(N)[1])), :];

        elseif argmin(tau) == 5 
            println("aM")

            rem = rand(DiscreteUniform(1,size(N)[1]));
            M = vcat(M,reshape(N[rem,:], 1, :));
            N = N[1:end .!= rem, :];

        elseif argmin(tau) == 6 
            println("bM")

            prob_matrix = bXMN./sum(bXMN);
            selected_index = rand(Categorical(vec(prob_matrix)), 1)[1];
            row = (selected_index ÷ size(bXMN)[1]) + 1;
            cols = ifelse(selected_index % size(bXMN)[1] == 0, size(bXMN)[1], selected_index % size(bXMN)[1]);

            psim = rand(Uniform(0,1));

            M = vcat(M,reshape((N[cols,:]+N[row,:])/2, 1, :));
            N = N[setdiff(1:end, (cols, row)), :];
        end
        Sr .+= log.(1/rand(Uniform(0,1)));

        #if init_size == size(X)[1]
        #    println("Error")
        #    break
        #end

        if (size(Y)[1] > 0)
            println("Dead")
            return 0
            break
        end
    end

    if (size(Y)[1] == 0) & (size(M)[1] > 0 )
        return 1
    elseif (size(Y)[1] == 0) & (size(M)[1] == 0 )
        return 2
    end
end

@everywhere function spatial_NTGSM2_fast(X::Matrix{Float64}, Y::Matrix{Float64}, N::Matrix{Float64}, M::Matrix{Float64}, ntgsm2::NTGSM2)   
    
    if ntgsm2.pa != 0
        println("Error pa")
        return -1
    elseif ntgsm2.qa != 0
        println("Error qa")
        return -1
    end

    dist = pairwise(Euclidean(), transpose(X));
    dist[dist .> ntgsm2.rd] .= 0;
    p = 1;
    if size(Y)[1] > 0
        return 0
    else
        for i in 1:size(X)[1]
            p *= ntgsm2.r/(ntgsm2.r+ntgsm2.a+ntgsm2.b*sum(dist[i,:]))
        end
    end

    dist = pairwise(Euclidean(), transpose(N));
    dist[dist .> ntgsm2.rdM] .= 0;
    pM = 1;
    if size(M)[1] > 0
        return 0
    else
        for i in 1:size(N)[1]
            pM *= ntgsm2.rM/(ntgsm2.rM+ntgsm2.aM+ntgsm2.bM*sum(dist[i,:]))
        end
    end

    return p*(1-pM)
end

@everywhere function simulate_GSM2(ion_::Ion, cell_::Cell, gsm2_::GSM2, Np_::Int64, Rk_::Float64, GY_R_tot_::Float64, X_::Matrix{Float64}, Y_::Matrix{Float64})
	for i in 1:Np_
		#println("ii = $ii on thread $(Threads.threadid())")
		local x, y = GenerateHit(cell_, Rk_)
		local track = Track(x, y, Rk_)
		local integral, theta, Gyr, radius = distribute_dose_vector(ion_, cell_, track)

		local X_i, Y_i = calculate_damage(ion_, cell_, integral, theta, Gyr, radius, x, y)

		local dist = sqrt.((X_i[:, 1]).^2 .+ (X_i[:, 2] ).^2)
		if size(dist[dist.>cell_.r], 1) != 0
			println("Error damage", X_i)
			return X_i
		end
		X_ = vcat(X_, X_i)
		Y_ = vcat(Y_, Y_i)
		#DOSE_tot+=dose
		GY_R_tot_ += Gyr
	end
	#println("Total imparted dose = $GY_R_tot_")

	local survP = spatial_GSM2_fast(X_, Y_, gsm2_)

	return survP, X_, Y_
end
