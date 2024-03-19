function calculate_kappa(ion,LET)

    if ion == "12C"
        p1 = 6.8;
        p2 = 0.156;
        p3 = 0.9214;
        p4 = 0.005245;
        p5 = 1.395;
    elseif ion == "4He"
        p1 = 6.8;
        p2 = 0.1471;
        p3 = 1.038;
        p4 = 0.006239;
        p5 = 1.582;
    elseif ion == "3He"
        p1 = 6.8;
        p2 = 0.1471;
        p3 = 1.038;
        p4 = 0.006239;
        p5 = 1.582;
    elseif ion == "1H"
        p1 = 6.8;
        p2 = 0.1773;
        p3 = 0.9314;
        p4 = 0;
        p5 = 1;
    elseif ion == "2H"
        p1 = 6.8;
        p2 = 0.1773;
        p3 = 0.9314;
        p4 = 0;
        p5 = 1;
    elseif ion == "16O"
        p1 = 6.8;
        p2 = 0.1749;
        p3 = 0.8722;
        p4 = 0.004987;
        p5 = 1.347;
    end
    yield = (p1 + (p2*LET)^p3)/(1+ (p4*LET)^p5);

    return yield
end

#function to calculate the spatial position of the damages
#for each track calculate the dose on the cell and ditribute the damages around the cell nucleus
function calculate_damage(ion, LET,cell, integral,radius, theta, Gyr)
    
    X = Array{Float64}(undef, 0, Nd);
    Y = Array{Float64}(undef, 0, Nd);

    theta_ = [theta[1:end-1]./2 theta[2:end]./2]
    theta = minimum(theta_, dims = 2)

    b = sqrt(x*x + y*y)
    
    kappa_DSB = 9*calculate_kappa(ion,LET);
    lambda_DSB = kappa_DSB*10^-3;

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
            X = vcat(X,reshape([Xx[i], Xy[i], cell.r_nucl*rand(Uniform(0,1),1)[1]], 1, :));
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
            Y = vcat(Y,reshape([Yx[i], Yy[i], cell.r_nucl*rand(Uniform(0,1),1)[1]], 1, :));
        end
    end

    return X, Y
end