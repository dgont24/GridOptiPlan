import os
import numpy as np
import pandas as pd
import pandapower as pp
from pandapower.control import ConstControl
from pandapower.timeseries import DFData
from pandapower.timeseries import OutputWriter
from pandapower.timeseries.run_time_series import run_timeseries

######################################################################################
# Time-series Simulation
######################################################################################

# Data source
## Data taken from csv data provided in IEEE benchmark grid

def create_load_data_source(aggregate_df, data_path, year):
    profiles = pd.DataFrame()
    for _, row in aggregate_df.iterrows():
        # add yealy load profile
        load_file_path = os.path.join(data_path,"load_data","15min",row['load_profile_name'])
        load_profile = pd.read_csv(load_file_path).consumption_kW.values * 1e-3  # read in mw

        # add yearly ev profile
        ev_profile = [0] * len(load_profile)
        if pd.notna(row['ev_profile_name']) and row['ev_integration_year'] <= year:
            ev_file_path = os.path.join(data_path,"ev_data","15min",row['ev_profile_name'])
            ev_profile = pd.read_csv(ev_file_path).charging_power_kW.values * 1e-3  # read in mw 

        profiles[row['load_name'] + '_P'] = np.add(load_profile, ev_profile)
        
        # add reactive power column
        profiles[row['load_name'] + '_Q'] = 0 

    ds = DFData(profiles)
    return profiles, ds

def create_gen_data_source(gens_df, data_path, year, timesteps):
# carefull data is in watt
    profiles = pd.DataFrame()
    for _, row in gens_df.iterrows():
        if row['integration_year'] <= year:
            profiles['PV' + str(row['id']) + '_P'] = pd.read_csv(os.path.join(data_path,"pv_data","15min",row['profile_name'])).p_w.values * 1e-6  # read in mw
        else:
            profiles['PV' + str(row['id']) + '_P'] = [0] * timesteps   
                 
    profiles['PV_Q'] = 0

    ds = DFData(profiles)
    return profiles, ds

# Time series controller
## P, Q values entered using P and power factor(cos_phi) data

def create_load_controllers(net, aggregate_df, ds):
    ConstControl(net, element='load', variable='p_mw', element_index=aggregate_df.index,
                 data_source=ds, profile_name=aggregate_df.load_name+'_P')
    ConstControl(net, element='load', variable='q_mvar', element_index=aggregate_df.index,
                 data_source=ds, profile_name=aggregate_df.load_name+'_Q')
    return net

def create_gen_controllers(net, gens_df, ds):
    ConstControl(net, element='sgen', variable='p_mw', element_index=gens_df.index,
                 data_source=ds, profile_name='PV'+gens_df['id'].astype(str)+'_P')
    ConstControl(net, element='sgen', variable='q_mvar', element_index=gens_df.index,
                 data_source=ds, profile_name='PV_Q')
    return net

# Format output
## - Bus voltage magnitudes and angles for each time step
## - Aggregated Sum of P & Q values for each time step

def create_output_writer(net, time_steps, output_dir):
    ow = OutputWriter(net, time_steps, output_path=output_dir, output_file_type=".json")
    ow.log_variable('res_trafo', 'p_lv_mw', index=net.trafo.index)
    ow.log_variable('res_trafo', 'q_lv_mvar', index=net.trafo.index)
    ow.log_variable('res_bus', 'vm_pu')
    ow.log_variable('res_bus', 'va_degree')
    ow.log_variable('res_line', 'loading_percent')
    return ow

def timeseries_calculation(net, n_timesteps, year, output_dir, aggregate_df, gens_df, data_path):
    # create data source
    load_profiles, ds_load = create_load_data_source(aggregate_df, data_path, year)
    sgen_profiles, ds_sgen = create_gen_data_source(gens_df, data_path, year, n_timesteps)

    # create controllers (to control P values of the load and the sgen)
    create_load_controllers(net, aggregate_df, ds_load)
    create_gen_controllers(net, gens_df, ds_sgen)

    # time steps to be calculated. Could also be a list with non-consecutive time steps
    time_steps = range(0, n_timesteps)

    # 4. the output writer with the desired results to be stored to files.
    ow = create_output_writer(net, time_steps, output_dir=output_dir)

    # 5. the main time series function
    run_timeseries(net, time_steps)


# load results
def load_time_series_results(output_dir, vnom):
    pp_vm = os.path.join(output_dir,'res_bus','vm_pu.json')
    df_pp_vm = pd.read_json(pp_vm) * vnom
    pp_ld = os.path.join(output_dir,'res_line','loading_percent.json')
    df_pp_ld = pd.read_json(pp_ld)

    # This is required since json makes keys as string type, the index order is like 1, 10, 100 ,...
    df_pp_vm.index = df_pp_vm.index.astype(np.int64)
    df_pp_vm = df_pp_vm.sort_index()
    df_pp_ld.index = df_pp_ld.index.astype(np.int64)
    df_pp_ld = df_pp_ld.sort_index()

    # Convert the DataFrame to a numpy array
    bus_voltage_data_array = df_pp_vm.values
    line_loading_data_array = df_pp_ld.values

    return bus_voltage_data_array, line_loading_data_array

######################################################################################
# Clustering
######################################################################################

# function to find CRITICAL BUSES for overvoltages or undervoltages
## KEEP ONLY BUSES ABOVE Q3 (75% quartile) of voltage deviation
def find_critical_voltage_buses(max_dv_per_bus):

    # Calculate the 75th percentile (Q3)
    q3 = np.percentile(max_dv_per_bus, 75)

    # Find indices of elements above Q3
    buses_above_q3 = {i: value for i, value in enumerate(max_dv_per_bus) if value > q3}

    # Find median of values above q3
    # q4_median = np.median(list(buses_above_q3.values()))

    return buses_above_q3

def find_median_cb(buses_q4):
    values = list(buses_q4.values())
    median_value = np.median(values)

    # Find the key corresponding to the median value
    key_with_median = next(key for key, value in buses_q4.items() if value == median_value)

    return key_with_median

def find_reference_voltage_buses(bus_voltage_array, vnom):
    # Calculate voltage descriptive statistics

    max_dv_ov_per_bus = [] # the maximum positive voltage deviation that appears in each bus throughout the year
    max_dv_uv_per_bus = [] # the maximum negative voltage deviation that appears in each bus throughout the year

    max_ov_all_buses = [] # list of the maximum voltage that appears in each bus daily
    max_uv_all_buses = [] # list of the minimum voltage that appears in each bus daily

    var_v_all_buses = [] # list of the voltage variance in each bus daily

    for bus_index in range(117):
        # Extract voltage data for the current bus
        bus_voltages = bus_voltage_array[:, bus_index].reshape((-1, 96))

        # Calculate maximum positive overvoltage for each day
        voltage_deviations = bus_voltages - vnom

        # ov - daily
        max_overvoltages = np.maximum(0, np.max(voltage_deviations, axis=1))
        max_ov_all_buses.append(max_overvoltages)
        # append max yearly overvoltages for each bus
        max_dv_ov_per_bus.append(np.max(max_overvoltages))

        # uv
        max_undervoltages = np.maximum(0, - np.min(voltage_deviations, axis=1))
        max_uv_all_buses.append(max_undervoltages)
        # append max yearly undervoltages for each bus
        max_dv_uv_per_bus.append(np.max(max_undervoltages))

        # voltage sigma^2 - daily - for each bus
        var_v_per_day = np.var(bus_voltages, axis=1)
        var_v_all_buses.append(var_v_per_day)

    # Convert the lists to NumPy arrays
    max_ov_all_buses = np.array(max_ov_all_buses)   # (117,365)
    max_uv_all_buses = np.array(max_uv_all_buses)   # (117,365)
    max_dv_ov_per_bus = np.array(max_dv_ov_per_bus) # (117,)
    max_dv_uv_per_bus = np.array(max_dv_uv_per_bus) # (117,)
    var_v_all_buses = np.array(var_v_all_buses)     # (117,365)

    # Find critical buses for ov and uv
    buses_q4_ov = find_critical_voltage_buses(max_dv_ov_per_bus)
    buses_q4_uv = find_critical_voltage_buses(max_dv_uv_per_bus)

    # Merge critical buses
    # Step 1: Extract buses
    keys_dict1 = set(buses_q4_ov.keys())
    keys_dict2 = set(buses_q4_uv.keys())

    # Step 2: Combine keys
    critical_buses = keys_dict1.union(keys_dict2)
    critical_buses = sorted(critical_buses)

    ref_bus_ov = find_median_cb(buses_q4_ov)
    ref_bus_uv = find_median_cb(buses_q4_uv)

    return ref_bus_ov, ref_bus_uv, max_ov_all_buses, max_uv_all_buses

def calculate_ilo(loading_values, threshold=100):
    overloadings = np.maximum(0, loading_values - threshold)
    ilo = np.sum(overloadings, axis=1)
    return ilo # (365,)

def calculate_loading_data(line_loading_array):

    # Calculate loading descriptive statistics
    max_ld_all_lines = [] # list of the maximum overloading that appears in each line daily
    max_ld_per_line_yearly = [] # the maximum overloading that appears in each line throughout the year

    ilo_each_line_daily = [] # list of the integrated line overloading for each line daily

    for line_index in range(line_loading_array.shape[1]):
        # Extract voltage data for the current bus
        line_loadings = line_loading_array[:, line_index].reshape((-1, 96))

        # overloading - daily
        max_daily_loadings = np.max(line_loadings, axis=1)
        max_ld_all_lines.append(max_daily_loadings)
        # append max yearly overvoltages for each bus
        max_ld_per_line_yearly.append(np.max(max_daily_loadings))

        # integrated overloading - daily - for each bus
        line_daily_ilo = calculate_ilo(line_loadings)
        ilo_each_line_daily.append(line_daily_ilo)

    # Convert the lists to NumPy arrays
    max_ld_all_lines = np.array(max_ld_all_lines)               # (115,365)
    max_ld_per_line_yearly = np.array(max_ld_per_line_yearly)   # (115,)
    ilo_each_line_daily = np.array(ilo_each_line_daily)         # (115,365)
    ilo_all_lines_daily = np.sum(ilo_each_line_daily, axis=0)   # (365,)

    # Find days that are overloaded

    days_with_overloading_ind = np.where(ilo_all_lines_daily > 0)[0]
    days_with_overloading = days_with_overloading_ind + 1

    return ilo_all_lines_daily, days_with_overloading_ind, days_with_overloading


# Create a list with the types that each day of the year belongs to
def create_day_category_list(day_dict):
    # Find the maximum day value
    max_day = max(max(days, default=0) for days in day_dict.values())

    # Initialize a list of 'unknown' for all days
    new_day_categories = ['unknown'] * (max_day + 1)  # +1 because list indexing starts at 0

    # Assign categories to corresponding day indices
    for category, days in day_dict.items():
        for day in days:
            new_day_categories[day] = category

    return new_day_categories[1:]  # slice to ignore index 0 since there's no day 0


def create_clusters(bus_voltage_array, line_loading_array, day_categories_names, vnom, soft_volt_limit, upper_volt_limit):

    ref_bus_ov, ref_bus_uv, max_ov_all_buses, max_uv_all_buses = find_reference_voltage_buses(bus_voltage_array, vnom)
    ilo_all_lines_daily, days_with_overloading_ind, days_with_overloading = calculate_loading_data(line_loading_array)

    # CATEGORIZE DAYS FOR REFERENCE BUSES - BASED ON THE VOLTAGES

    ovs = max_ov_all_buses[ref_bus_ov, :]  # daily max overvoltages for ov ref bus
    uvs = max_uv_all_buses[ref_bus_uv, :]  # daily max undervoltages for uv ref bus
    day_categories = {name: [] for name in day_categories_names}

    for day in range(len(ovs)):
        # Assign each day to the respective category
        if  0 <= ovs[day] < soft_volt_limit  and 0 <= uvs[day] < soft_volt_limit and day not in days_with_overloading_ind:
            day_categories["normal"].append(day + 1)
        elif ((soft_volt_limit  <= ovs[day] < upper_volt_limit and 0                <= uvs[day] < upper_volt_limit and day not in days_with_overloading_ind) or 
            (0                <= ovs[day] < soft_volt_limit  and soft_volt_limit  <= uvs[day] < upper_volt_limit and day not in days_with_overloading_ind)):
            day_categories["limited-ov/uv"].append(day + 1)
        elif (upper_volt_limit  <= ovs[day]                    and 0                <= uvs[day] < upper_volt_limit and day in days_with_overloading_ind):
            day_categories["ov + ol"].append(day + 1)
        elif (upper_volt_limit  <= ovs[day]                    and 0                <= uvs[day] < upper_volt_limit and day not in days_with_overloading_ind):
            day_categories["ov"].append(day + 1)
        elif (0                 <= ovs[day] < upper_volt_limit and upper_volt_limit <= uvs[day]                    and day in days_with_overloading_ind):
            day_categories["uv + ol"].append(day + 1)
        elif (0                 <= ovs[day] < upper_volt_limit and upper_volt_limit <= uvs[day]                    and day not in days_with_overloading_ind):
            day_categories["uv"].append(day + 1)
        elif (upper_volt_limit  <= ovs[day]                    and upper_volt_limit <= uvs[day]                    and day in days_with_overloading_ind):
            day_categories["ov + uv + ol"].append(day + 1)
        elif (upper_volt_limit  <= ovs[day]                    and upper_volt_limit <= uvs[day]                    and day not in days_with_overloading_ind):
            day_categories["ov + uv"].append(day + 1)
        else:
            day_categories["rest"].append(day + 1)

    return day_categories, ovs, uvs, ilo_all_lines_daily
    
# Function to normalize an array
def normalize_array(array):
    """
    Normalize an array based on the formula:
    (value - min) / (max - min)
    """
    min_val = np.min(array)
    max_val = np.max(array)
    # Prevent division by zero if max and min are the same
    if max_val - min_val == 0:
        return array
    normalized_array = (array - min_val) / (max_val - min_val)
    return normalized_array

# Function to calculate a weighted average of normalized values
def calculate_weighted_average(ovs, uvs, ilo):
    """
    Calculate a weighted average of normalized values from three arrays.
    The weights are equally distributed.
    """
    n_ovs = normalize_array(ovs)
    n_uvs = normalize_array(uvs)
    n_ilo = normalize_array(ilo)
    return (n_ovs + n_uvs + n_ilo) / 3

# Choose the representative days. 
# This part is different from the v1 as instead from choosing the days based on variance, 
# a normalization function is used with equivalent weight factor for all three dimensions (ov, uv and ol)
def get_representative_days(ovs, uvs, ilo_all_lines_daily, day_categories):
    combined_values = calculate_weighted_average(ovs, uvs, ilo_all_lines_daily)

    # Choose the representative days for each cluster
    results = {}
    for category, days in day_categories.items():
        # Skip certain categories
        if category in ['normal', 'limited-ov/uv']:
            continue

        if not days:
            results[category] = []  # Assign empty list for empty categories
            continue

        # Calculate and sort normalized values for each day
        normalized_values = [combined_values[day - 1] for day in days]
        days_sorted_by_value = sorted(zip(normalized_values, days), key=lambda x: x[0])

        # Determine representative days based on the number of days
        if len(days) == 1:
            indices = [0]
        elif len(days) == 2:
            indices = [0, 1]
        elif len(days) == 3:
            indices = [0, 1, 2]        
        elif len(days) <= 60:
            # Pick indices that divide the list into roughly equal parts (min, median, max)
            indices = [0, len(days_sorted_by_value) // 2, -1]
        elif len(days) <= 90:
            # Pick indices for min, 1/3, 2/3, max
            indices = [0, len(days_sorted_by_value) // 3, 2 * len(days_sorted_by_value) // 3, -1]
        else:
            # For more than 90 days, pick min, 25%, 50%, 75%, max        
            indices = [0, len(days_sorted_by_value) // 4, len(days_sorted_by_value) // 2, 3 * len(days_sorted_by_value) // 4, -1]

        # Extract the days corresponding to the selected indices
        results[category] = [days_sorted_by_value[i][1] for i in indices]

    return combined_values, results

def update_cluster_dataframe(df, year, combined_values, day_categories):
    year_id = [year] * len(combined_values)
    days_id = list(range(1, len(combined_values) + 1))
    day_cluster_types = create_day_category_list(day_categories)
    
    # Create a DataFrame for the new data
    new_data = pd.DataFrame(list(zip(year_id, days_id, day_cluster_types, combined_values)),
                            columns=['year', 'id', 'type', 'f_value'])
    
    # Append the new data to the existing DataFrame
    updated_df = pd.concat([df, new_data], ignore_index=True)
    
    return updated_df    