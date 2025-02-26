function battery_capacity_fade(soc_data::Vector{Float64}, resolution::Int, dt_hr_adjust::Float64)
    temperature = 25  # mean temperature (C)
    timestep_s = resolution * 60  # timestep in sec

    # Calendar aging
    d_cal = calendar_ageing(soc_data, temperature, timestep_s)

    # Analyze input data of battery's state of charge
    ext, exttime = sig2ext(soc_data) # Calculate turning points
    dod, soc, c_rate = rainflow_top(ext, exttime, dt_hr_adjust) # calculate converted parameters
    cycles = rainflow(ext, exttime) # calculate cycles using rainflow algorithm
    if isempty(cycles)
        d_tot = sum(d_cal)
        q_bat = 1 - nonlinear_dynamic(d_tot)
        return q_bat
    end

    # Cycle aging
    cycle_start_time = sort(cycles[4,:])
    d_cyc = Vector{Float64}(undef, size(cycles,2))
    for i in axes(cycles,2)
        time = cycle_start_time[i]
        index = findfirst(isequal(time), exttime)
        soc_stress = soc_effect(soc[index])
        dod_stress = dod_effect(dod[index])
        c_rate_stress = c_rate_effect(c_rate[index])
        d_cyc[i] = soc_stress * dod_stress * c_rate_stress * cycles[3,i]
    end

    # Total degradation
    d_cal_tot = sum(d_cal)
    d_cyc_tot = sum(d_cyc)
    d_tot = d_cyc_tot + d_cal_tot

    # Capacity degredation percentage
    q_bat = 1 - nonlinear_dynamic(d_tot)

    return q_bat
end

function nonlinear_dynamic(d::Float64)
    p_sei = 0.0296
    r_sei = 150.24
    soh = p_sei * exp(-r_sei * d) + (1 - p_sei) * exp(-d)
    return soh
end

function calendar_ageing(soc::Vector{Float64}, t_avg::Int, time_interval_sec::Int)
    k_soc = 0.265
    k_cal = 6.396e-10
    soc_ref = 0.5
    k_t = 0.0693
    t_ref = 25

    d_c = k_cal * exp.(-((soc .- soc_ref) / k_soc) .^ 2 .* (soc_ref ./ soc)) .* time_interval_sec .* 
           exp.(k_t * (t_avg .- t_ref) .* (273 + t_ref) ./ (273 + t_avg))
    
    return d_c
end

function soc_effect(soc::Float64)
    k_soc = 0.2685
    soc_ref = 0.5
    d_s = exp(-((soc - soc_ref) / k_soc) ^ 2 * (soc_ref / soc))
    return d_s
end

# DoD effect on cycle ageing
# Input: 
#   dod - average cycle DoD level
#   threshold - The limit between 0 and which a linear degradation is applied
#   e - degradation at theoretical 0 DoD in linearized degradation
# Output:
#   d_d - cycle degradation damage of each DoD
function dod_effect(dod::Float64)
    # Coefficients
    a1 = 1.108e7
    a2 = -0.7778
    a3 = 2.334e5
    a4 = -0.08211

    d_d = 1 / (a1 * exp(a2 * dod) + a3 * exp(a4 * dod))

    return d_d
end

function c_rate_effect(c_rate::Float64)
    k_c1 = 0.5681  
    k_c2 = 0.5655
    d_c = k_c1 * exp(k_c2 * c_rate)
    return d_c
end

function rainflow_top(ext::Vector{Float64}, exttime::Vector{Float64}, dt_hr_adjust::Float64)
    tot_len = length(ext) - 1
    dod = Vector{Float64}(undef, tot_len)
    soc = Vector{Float64}(undef, tot_len)
    c_rate = Vector{Float64}(undef, tot_len)

    for i in 1:tot_len
        dod[i]  = abs(ext[i + 1] - ext[i])
        # State of Charge
        soc[i] = ext[i] + ((ext[i + 1] - ext[i]) / 2)
        # Crate
        duration_hr = (exttime[i + 1] - exttime[i]) * dt_hr_adjust
        c_rate[i] = dod[i] / duration_hr
    end

    return dod, soc, c_rate
end

# Outputs
#   1: Cycles amplitude, 
#   2: Cycles mean value, 
#   3: Number of cycles (0.5 or 1.0),
#   4: Begining time (when input includes dt or extt data),
#   5: Cycle period (when input includes dt or extt data),
function rainflow(array_ext::Vector{Float64}, array_t::Vector{Float64})
    tot_num = length(array_ext)
    if tot_num != length(array_t)
        throw(ArgumentError("RAINFLOW: Time Array size error."))
    end
    a = zeros(Float64, tot_num)
    t = zeros(Float64, tot_num)
    outputs = Vector{Float64}[]

    j = 0
    cNr = 1
    for index in 1:tot_num
        j += 1
        a[j] = array_ext[index]
        t[j] = array_t[index]

        while j >= 3 && abs(a[j-1] - a[j-2]) <= abs(a[j] - a[j-1])
            ampl = abs((a[j-1] - a[j-2]) / 2)

            if j == 3
                mean = (a[1] + a[2]) / 2
                period = (t[2] - t[1]) * 2
                atime = t[1]
                a[1] = a[2]
                a[2] = a[3]
                t[1] = t[2]
                t[2] = t[3]
                j = 2

                if ampl > 0
                    push!(outputs, [ampl, mean, 0.50, atime, period])
                end
            else
                mean = (a[j-1] + a[j-2]) / 2
                period = (t[j-1] - t[j-2]) * 2
                atime = t[j-2]
                a[j-2] = a[j]
                t[j-2] = t[j]
                j -= 2

                if ampl > 0
                    push!(outputs, [ampl, mean, 1.00, atime, period])
                    cNr += 1
                end
            end
        end
    end

    for index in 1:j-1
        ampl = abs(a[index] - a[index + 1]) / 2
        mean = (a[index] + a[index + 1]) / 2
        period = (t[index + 1] - t[index]) * 2
        atime = t[index]

        if ampl > 0
            push!(outputs, [ampl, mean, 0.50, atime, period])
        end
    end
    resize!(outputs, tot_num - cNr)

    if isempty(outputs)
        return outputs
    end

    output_arr = reduce(hcat, outputs)

    return output_arr
end

function sig2ext(sig::Vector{Float64}; dt::Float64=1.0)
    # Finding local extrema
    w1 = diff(sig)
    w = [true; (w1[1:end-1] .* w1[2:end]) .<= 0; true]
    ext = sig[w]
    exttime = findall(w) * dt

    # Remove triple values
    w1 = diff(ext)
    w = .!([false; (w1[1:end-1] .== 0) .& (w1[2:end] .== 0); false])
    ext = ext[w]
    exttime = exttime[w]

    # Remove double values and adjust times
    w = .!([false; ext[1:end-1] .== ext[2:end]])
    ext = ext[w]
    w1 = (exttime[2:end] - exttime[1:end-1]) ./ 2
    exttime = [exttime[1:end-1] + w1 .* .!w[2:end]; exttime[end]]
    exttime = exttime[w]

    # Final check for extrema
    if length(ext) > 2
        w1 = diff(ext)
        w = [true; w1[1:end-1] .* w1[2:end] .< 0; true]
        ext = ext[w]
        exttime = exttime[w]
    end

    return ext, exttime
end

