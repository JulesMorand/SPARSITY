using Base.MathConstants

module Survival

    export Nucleus_Integral_t

    struct Nucleus_Integral_t
        cellLine::CellLine
        x_nucleus::Float64
        y_nucleus::Float64
        r_nucleus::Float64
        totalNucleusDose::Float64
        inNucleusCount::Int
        intersectionCount::Int
        doses::Vector{Float64}
        times::Vector{Float64}
    end

    function Nucleus_Integral_t(cellLine::CellLine, xPosition::Float64, yPosition::Float64)
        r_nucleus = cellLine.getNucleusRadius()
        new(cellLine, xPosition, yPosition, r_nucleus, 0.0, 0, 0, Float64[], Float64[])
    end

    function Nucleus_Integral_t(nn::Nucleus_Integral_t, cellLine::CellLine)
        totalNucleusDose = nn.totalNucleusDose
        inNucleusCount = nn.inNucleusCount
        intersectionCount = nn.intersectionCount
        doses = copy(nn.doses)
        times = copy(nn.times)
        Nucleus_Integral_t(cellLine, nn.x_nucleus, nn.y_nucleus, nn.r_nucleus, totalNucleusDose, inNucleusCount, intersectionCount, doses, times)
    end

    function addBackgroundDose!(nucleus::Nucleus_Integral_t, dose::Float64, t::Float64)
        nucleus.totalNucleusDose += dose
        push!(nucleus.doses, dose)
        push!(nucleus.times, t)
    end

    function arc_intersection_weight(r::Float64, b::Float64, r_nucleus::Float64)
        if b < r_nucleus
            if r <= r_nucleus - b
                return 2 * π
            elseif r < b + r_nucleus
                arg = b / (2 * r) + r / (2 * b) - r_nucleus * r_nucleus / (2 * b * r)
                return 2 * acos(arg)
            end
        else
            if r <= b - r_nucleus
                return 0.0
            elseif r < b + r_nucleus
                arg = b / (2 * r) + r / (2 * b) - r_nucleus * r_nucleus / (2 * b * r)
                arg = min(arg, 1.0)
                return 2 * acos(arg)
            end
        end
        return 0.0
    end

    function clean_nucleus!(nucleus::Nucleus_Integral_t)
        nucleus.totalNucleusDose = 0.0
        nucleus.inNucleusCount = 0
        nucleus.intersectionCount = 0
        empty!(nucleus.times)
        empty!(nucleus.doses)
    end

    function clone_nucleus(nn::Nucleus_Integral_t, cellLine::CellLine)
        Nucleus_Integral_t(cellLine, nn.x_nucleus, nn.y_nucleus)
    end

    function distribute_dose!(nucleus::Nucleus_Integral_t, track::Track)
        x_track, y_track = track.getPosition()
        x_track = (x_track - nucleus.x_nucleus) * 1e3  # mm -> um
        y_track = (y_track - nucleus.y_nucleus) * 1e3  # mm -> um
        b = sqrt(x_track^2 + y_track^2)

        rMax = min(track.getRadius(), b + nucleus.r_nucleus)

        area1 = area2 = area3 = 0.0

        if b <= nucleus.r_nucleus
            nucleus.inNucleusCount += 1
            rMin = 0.0

            if b + track.getRadius() < nucleus.r_nucleus
                r_intersection = track.getRadius()
            else
                r_intersection = nucleus.r_nucleus - b
            end

            area1 = π * r_intersection^2
            dose = track.getRadialIntegral(0.0, r_intersection) * area1

            if rMax > r_intersection
                area2 = integrate_weighted_radial_track(track, r_intersection, rMax, b, 0.01)
                dose += area2
            end

            if rMax == track.getRadius()
                if track.getRadius() > nucleus.r_nucleus - b
                    theta1 = acos((b^2 + track.getRadius()^2 - nucleus.r_nucleus^2) / (2 * b * track.getRadius()))
                    theta2 = acos((b^2 - track.getRadius()^2 + nucleus.r_nucleus^2) / (2 * b * nucleus.r_nucleus))
                    area3 = π * nucleus.r_nucleus^2 - (theta1 * track.getRadius()^2 + theta2 * nucleus.r_nucleus^2 - track.getRadius() * b * sin(theta1))
                else
                    area3 = π * (nucleus.r_nucleus^2 - r_intersection^2)
                end
            end

            dose /= area1 + area2 + area3

            nucleus.totalNucleusDose += dose
            push!(nucleus.doses, dose)
            push!(nucleus.times, track.getTime())
        elseif b <= nucleus.r_nucleus + track.getRadius()
            nucleus.intersectionCount += 1
            rMin = b - nucleus.r_nucleus
            area2 = integrate_weighted_radial_track(track, rMin, rMax, b, 0.01)
            dose = area2

            if rMax == track.getRadius()
                theta1 = acos((b^2 + track.getRadius()^2 - nucleus.r_nucleus^2) / (2 * b * track.getRadius()))
                theta2 = acos((b^2 - track.getRadius()^2 + nucleus.r_nucleus^2) / (2 * b * nucleus.r_nucleus))
                area3 = π * nucleus.r_nucleus^2 - (theta1 * track.getRadius()^2 + theta2 * nucleus.r_nucleus^2 - track.getRadius() * b * sin(theta1))
            end

            dose /= area2 + area3

            nucleus.totalNucleusDose += dose
            push!(nucleus.doses, dose)
            push!(nucleus.times, track.getTime())
        end
    end

    function distribute_dose!(nucleus::Nucleus_Integral_t, tracks::Tracks)
        for i in 1:length(tracks)
            distribute_dose!(nucleus, tracks[i])
        end
    end

    function get_cell_type(nucleus::Nucleus_Integral_t)
        return nucleus.cellLine.getCellType()
    end

    function get_dose_and_survival(nucleus::Nucleus_Integral_t, dose::Float64, doseUncertainty::Float64, survival::Float64, survivalUncertainty::Float64)
        dose = nucleus.totalNucleusDose
        sum_lethal = nucleus.cellLine.getLogSurvival_X(nucleus.doses, nucleus.times)
        survival = exp(-sum_lethal)
        doseUncertainty = -1
        survivalUncertainty = -1
    end

    function get_dose_and_lethals(nucleus::Nucleus_Integral_t, dose::Float64, doseUncertainty::Float64, lethals::Float64, lethalsUncertainty::Float64)
        dose = nucleus.totalNucleusDose
        lethals = nucleus.cellLine.getLogSurvival_X(nucleus.doses, nucleus.times)
        doseUncertainty = -1
        lethalsUncertainty = -1
    end

    function get_position(nucleus::Nucleus_Integral_t, returnX::Float64, returnY::Float64)
        returnX = nucleus.x_nucleus
        returnY = nucleus.y_nucleus
    end

    function integrate_weighted_radial_track(track::Track, rMin::Float64, rMax::Float64, b::Float64, step::Float64)
        r1, r2, log_r2, log_rMin, log_rMax, log_step = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        f1, f2, f, arc_weight1, arc_weight2 = 0.0, 0.0, 0.0, 0.0, 0.0
        integral = 0.0

        if rMin > 0
            log_rMin = log10(rMin)
        else
            log_rMin = -5
        end

        log_rMax = log10(rMax)
        log_step = step

        nSteps = ceil((log_rMax - log_rMin) / log_step)

        if nSteps < 3
            log_step = (log_rMax - log_rMin) / 3
            nSteps = 3
        end

        area = 0.0

        arc_weight2 = arc_intersection_weight(rMin, b, nucleus.r_nucleus)
        f2 = track.getLocalDose(rMin) * rMin * arc_weight2
        r2 = rMin

        for i in 1:nSteps - 1
            log_r2 = log_rMin + log_step * (i + 1)
            f1 = f2
            r1 = r2
            arc_weight1 = arc_weight2

            r2 = 10^log_r2
            arc_weight2 = arc_intersection_weight(r2, b, nucleus.r_nucleus)
            f2 = track.getLocalDose(r2) * r2 * arc_weight2
            f = (r2 - r1) * (f1 / 2.0 + f2 / 2.0)
            integral += f
            area += (r2 - r1) * (arc_weight1 * r1 / 2.0 + arc_weight2 * r2 / 2.0)
        end

        f1 = f2
        r1 = r2
        arc_weight1 = arc_weight2

        r2 = rMax
        arc_weight2 = arc_intersection_weight(r2, b, nucleus.r_nucleus)
        f2 = track.getLocalDose(r2) * r2 * arc_weight2
        f = (r2 - r1) * (f1 / 2.0 + f2 / 2.0)
        integral += f
        area += (r2 - r1) * (arc_weight1 * r1 / 2.0 + arc_weight2 * r2 / 2.0)

        return integral / area
    end

end  # module Survival
