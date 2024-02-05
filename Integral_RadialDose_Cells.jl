using Distributions
###########
R = 8.0; #radius of the cicle where uniform hit is generated
#Particle
#amorphous track parameters
E = 56.0; #C 149MeV %p 80MeV %p 18.6MeV %C 280.0Mev #He 145.74-56MeV
A = 1;
LET = 4.5; #C 149MeV 19.56 keV/um, p 80MeV 0.86(0.6) keV/um, p 18.6MeV  2.76 keV/um %C 280.0Mev 13.12 keV/um #He 145.74-56MeV 4.5
rho = 1;

par = string("He", "_", string(E));

#Time and number of particle
Dose = 12.0;
kR = 0.8;

function ATRadius(E, A, kR)

    Rc = 0.01;

    Rp = 0.05*((E/A)^(1.7));

    if kR < 1.0
        Rk = Rc*exp((kR*(1+2*log(Rp/Rc))-1)/2);
    else
        Rk = Rp;
    end

    return Rc, Rp, Rk

end

(Rc, Rp, Rk) = ATRadius(E, A, kR);
println("Rc=",Rc,"\nRp=",Rp,"\nRk=",Rk)

DoseRate = 0.18; #10
DoseRate_h = DoseRate*3600;

F = Dose/(1.602*10^(-9)*LET);
Npar = round(Int, F*(pi*(R+Rk)^2*10^(-8)));
# if no simulation of time
#Np = round(Int, Dose/(1.602*10^(-9)*LET/(pi*(R+Rk)^2*10^(-8))));
#Np = Npar;

#zF = 1.602*10^(-9)*LET/(pi*(R+Rk)^2*10^(-8));
zF = Dose/Npar;
D = DoseRate_h/zF;
T = Dose/(zF*D)*3600;
#######################
function GenerateHit(R, Rk)

    radius = (R+Rk)*sqrt((rand(Uniform(0,1))));
    theta = 2*pi*rand(Uniform(0,1));

    x0 = radius*cos(theta);
    y0 = radius*sin(theta);

    return x0, y0

end
#######################
function GetRadialLinearDose(r)
    #LET normalized to Rk ???
    LETk = LET*0.1602;####JM: why thi is usefull?
    D_arc=0.
    if r <= Rc
        D_arc = (1/(pi*Rc^2))*(LETk/(1*(1+2*log(Rk/Rc))));
    elseif r <= Rk
        D_arc = (1/(pi*r*r))*(LETk/(1*(1+2*log(Rk/Rc))));
    end
    return D_arc
end

##############


#Compute the angle theta in rad. of the intersection arc betweentwo circles:
#lenght of the arc=theta*r, cos(Theta)=x/r=arg
#x is here:
#https://mathworld.wolfram.com/Circle-CircleIntersection.html
function arc_intersection_angle(r, b, r_nucleus)
    if b < r_nucleus
        if r <= r_nucleus - b
            return 2 * π
        elseif r < b + r_nucleus
            arg = (b*b  + r*r  - r_nucleus * r_nucleus) / (2 * b * r)
            return 2 * acos(arg)
        end
    else #(b>=r_nucleus)
        if r <= b - r_nucleus
            return 0.0
        elseif r < b + r_nucleus
            arg =(b*b  + r*r  - r_nucleus * r_nucleus) / (2 * b * r)
            #arg = min(arg, 1.0)not sure if usefull?
            return 2 * acos(arg)
        end
    end
    return 0.0
end

#########
r_nucleus=5.;
function integrate_weighted_radial_track( rMin::Float64, rMax::Float64, b::Float64, step::Int64)
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
            f2 = GetRadialLinearDose(r2) * r2 * arc_weight2
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
        f2 = GetRadialLinearDose(r2) * r2 * arc_weight2
        f = (r2 - r1) * (f1 / 2.0 + f2 / 2.0)
        integral += f
        area += (r2 - r1) * (arc_weight1 * r1 / 2.0 + arc_weight2 * r2 / 2.0)

        return area, integral, area/integral
    end