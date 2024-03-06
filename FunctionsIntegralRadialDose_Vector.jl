
#Added function as NameFunction_vector that return vectors instead of numbers as before. 
#They are used to caluclate the spatal distribution of damages that is proportional to the dose
@everywhere function arc_intersection_angle_vector(r, b, r_nucleus)
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
@everywhere function integrate_weighted_radial_track_vector( rMin::Float64, rMax::Float64, b::Float64, r_nucleus::Float64, nSteps::Int64)
    r1, r2, log_r2, log_rMin, log_rMax, log_step = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    f1, f2, f, arc_weight1, arc_weight2 = 0.0, 0.0, 0.0, 0.0, 0.0
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
    f2 = GetRadialLinearDose(rMin) * rMin * arc_weight2
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
        f2 = GetRadialLinearDose(r2) * r2 * arc_weight2
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
    f2 = GetRadialLinearDose(r2) * r2 * arc_weight2
    f = (r2 - r1) * (f1 / 2.0 + f2 / 2.0)
    push!(integral,f)
    area += (r2 - r1) * (arc_weight1 * r1 / 2.0 + arc_weight2 * r2 / 2.0)
    push!(radius,r1)
    push!(radius,r2)
    return integral, theta, area, radius
end

#####################################

@everywhere function distribute_dose_vector(cell::Cell, track::Track)
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
        integral1, theta11, area1, radius1 = integrate_weighted_radial_track_vector(0., r_intersection, b, cell.r, 1000);
        dose = sum(integral1)
        # dose = track.getRadialIntegral(0.0, r_intersection) * area1

        if rMax > r_intersection
            integral2, theta21, area, radius2 = integrate_weighted_radial_track_vector(r_intersection, rMax, b, cell.r,1000);
            area2 = area #track, r_intersection, rMax, b, 0.01
            dose += sum(integral2)
            integral = [integral1; integral2];
            theta = [theta11; theta21];
            radius = [radius1; radius2];
        else
            integral = integral1;
            theta = theta11;
            radius = radius1;
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
        integral, theta, area, radius = integrate_weighted_radial_track_vector(rMin, rMax, b, cell.r, 1000)
        dose = sum(integral)
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
    return integral, theta, Gyr, radius
end
