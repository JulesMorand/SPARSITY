double Nucleus_Integral_t::ArcIntersectionWeight(double r,
                                                 double b)
{
    double arg;
    
    if (b < r_nucleus) {
        if (r <= r_nucleus - b) {
            return 2 * M_PI;
        } else if (r < b + r_nucleus) {
            arg = b/(2*r) + r/(2*b) - r_nucleus*r_nucleus / (2*b*r);
            return 2 * acos( arg );
        }
    } else {
        if (r <= b - r_nucleus) {
            return 0.;
        } else if (r < b + r_nucleus) {
            arg = b/(2*r) + r/(2*b) - r_nucleus*r_nucleus / (2*b*r);
            if (arg > 1.)
                arg = 1.;
            return 2 * acos( arg );
        }
    }
    return 0.;
}