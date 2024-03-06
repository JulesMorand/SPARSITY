using Distributions

###########
R = 8.0; #radius of the beam cicle where uniform hit is generated
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
    #LET normalized to Rk -- the cut off
    LETk = LET*0.1602;####JM: why thi is usefull?
    D_arc=0.
    if r <= Rc
        D_arc = (1/(pi*Rc^2))*(LETk/(1*(1+2*log(Rk/Rc))));
    elseif r <= Rk
        D_arc = (1/(pi*r*r))*(LETk/(1*(1+2*log(Rk/Rc))));
    end
    return D_arc
end