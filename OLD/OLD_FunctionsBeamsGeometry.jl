using Distributions

######amorphous track radius
function ATRadius(E, A, kR)
    local Rc = 0.01;
    local Rp = 0.05*((E/A)^(1.7));
    if kR < 1.0
        local Rk = Rc*exp((kR*(1+2*log(Rp/Rc))-1)/2);
    else
        local Rk = Rp;
    end
    return Rc, Rp, Rk
end
#######################
#amorphous track parameters
global E = 56.0; #C 149MeV %p 80MeV %p 18.6MeV %C 280.0Mev #He 145.74-56MeV
global A = 1;
global kR = 0.8;
global (Rc, Rp, Rk) = ATRadius(E, A, kR);
#println("Rc=",Rc,"\nRp=",Rp,"\nRk=",Rk)
function GenerateHit(R, Rk)
    radius = (R+Rk)*sqrt((rand(Uniform(0,1))));
    theta = 2*pi*rand(Uniform(0,1));
    x0 = radius*cos(theta);
    y0 = radius*sin(theta);
    return x0, y0
end
#######################

function GetRadialLinearDose(r,LET)

    #LET normalized to Rk -- the cut off
    local LETk = LET*0.1602;
    D_arc=0.
    if r <= Rc
        D_arc = (1/(pi*Rc^2))*(LETk/(1*(1+2*log(Rk/Rc))));
    elseif r <= Rk
        D_arc = (1/(pi*r*r))*(LETk/(1*(1+2*log(Rk/Rc))));
    end
    return D_arc
end