# PARAMETERS
# Function to create dictionary to hold all global parameters
function create_global_parameters()
    params = Dict{String, Real}()
    
    # Choose daily optimization resolution in minutes
    params["Resolution"] = 15
    # Slack bus reference value [pu]
    params["V_REF_PU"] = 1.0
    # Nominal network voltage in the lv side [V]
    params["V_NOM"] = 416.0
    # Angle for linear approx. of quadratic thermal constraints [rad]
    params["PHI2"] = pi/4 
    # Voltage Limits
    params["Vmin"] = 0.9
    params["Vmax"] = 1.1
    params["Vmin2"] = params["Vmin"]^2
    params["Vmax2"] = params["Vmax"]^2
    # Thermal limits
    params["ld_max_%"] = 100.0

    # PU Values
    params["V_base"] = params["V_NOM"] #V
    params["S_base"] = 10000.0 #VA  
    params["Z_base"] = params["V_base"]^2 / params["S_base"] #Ohm   

    # for converting energy to the right time scale. It is used for the battery State of Energy.
    # Power is in units MW. So it needs to be multiplied by 0.25 to calculate the volume of energy in a 15-minute period.
    params["Dt_hr_adjust"] = params["Resolution"] / 60.0

    # total timesteps of the daily optimization problem 
    params["timesteps"] = div(60*24, params["Resolution"]) # 96 for 15min resolution

    # Power Factor limit for the control of reactive power of the pv inverters
    pf_lim = 0.93
    params["max_phi"] = acos(pf_lim)

    # lines thermal limits constraints
    params["a_phi2"] = (1-cos(params["PHI2"]))/sin(params["PHI2"])
    params["b_phi2"] = sin(params["PHI2"])

    return params
end