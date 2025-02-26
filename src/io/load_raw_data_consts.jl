const _cluster_columns = [("year", Int), ("id", Int), ("type", String), ("f_value", Float64)]

const _source_columns = [("voltage_kv", Float64), ("v_setpoint_pu", Float64)]

const _trafo_columns = [
    ("name", String), ("bus1", Int), ("bus2", Int),
    ("kv_pri", Float64), ("kv_sec", Float64),
    ("mva", Float64), ("x_100", Float64), ("r_100", Float64)
]

const _bus_columns = [("id", Int), ("name", String)]

const _line_columns = [
    ("from_bus", Int), ("to_bus", Int), ("length", Float64), # in meters
    ("linecode", String)
]

const _linecode_columns = [
    ("name", String), ("type", String),
    ("r1", Float64), ("x1", Float64), ("units", String), ("ampacity_a", Int), 
    ("installation_cost_per_m", Float64), ("capital_cost_per_m", Float64)
]

const _load_columns = [
    ("bus", Int), ("id", Int), ("profile_name", String), 
]
const _load_profile_columns = [
    ("date", Int), ("time", Int), ("consumption_kW", Float64)
]

const _ev_columns = [
    ("bus", Int), ("id", Int), ("profile_name", String), ("max_ac_charging_kW", Float64)
]
const _ev_profile_columns = [
    ("date", Int), ("time", Int), ("charging_power_kW", Float64)
]

const _pv_columns = [
    ("bus", Int), ("id", Int), ("profile_name", String), ("pv_size", Float64)
]
const _pv_profile_columns = [
    ("date", Int), ("time", Int), ("p_w", Float64)
]

const _batterycode_columns = [
    ("code", Int), ("name", String),
    ("c_rate", Float64), ("soc_min", Float64), ("soc_max", Float64), ("efficiency", Float64),
    ("costcapacity_eur_per_kwh", Float64), ("costpower_eur_per_kw", Float64), ("costmaint_eur_per_year_per_kWh", Float64),
    ("lifetime_years", Int), ("cycles", Int), ("end_life_deg_percent", Float64)
]

const _tap_changer_columns = [
    ("code", Int), ("name", String), ("reg", Float64), ("ns", Float64), ("max_n", Float64),
    ("cinv", Float64), ("cmain", Float64)
]

