import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cmocean.cm as cmo

#=========================================== The following functions are for use in the final project for MS-263======================================================

def strength_location_divisions(list_of_data, wind_threshold, bay_threshold, max_pressure):
     '''
    Categorize dissipation profile data based on wind strength and geographic location.

    This function separates profiles into four groups depending on whether wind forcing 
    is above or below a given threshold (strong vs. weak) and whether the profile lies 
    inside or outside a specified longitudinal threshold (in-bay vs. out-of-bay). 
    It returns a DataFrame for each group, containing filtered dissipation rates, 
    buoyancy frequency (N²), pressure, wind speed, time, and geographic coordinates.

    Parameters:
    -----------
    list_of_data : list
        A list of xarray Datasets, each containing microstructure profile data with 
        variables such as 'dissipation', 'pressure', 'N_squared', 'wind', 
        'start_time', 'profile_lat', and 'profile_long'.
    
    wind_threshold : float
        Wind speed threshold (e.g., in m/s) used to distinguish between strong and weak wind conditions.

    bay_threshold : float
        Longitude threshold to divide profiles into "in-bay" (east of threshold) 
        and "out-of-bay" (west of threshold) groups.

    max_pressure : float
        Maximum pressure up to which profiles will be included; deeper data is excluded.

    Returns:
    --------
    df_strong_in : pandas.DataFrame
        Profiles with strong wind and located in-bay.

    df_strong_out : pandas.DataFrame
        Profiles with strong wind and located out-of-bay.

    df_weak_in : pandas.DataFrame
        Profiles with weak wind and located in-bay.

    df_weak_out : pandas.DataFrame
        Profiles with weak wind and located out-of-bay.
    '''

    diss_strong_in, press_strong_in, time_strong_in, wind_strong_in, lat_strong_in, long_strong_in, N2_strong_in= [], [], [], [], [], [], []
    diss_strong_out, press_strong_out, time_strong_out, wind_strong_out, lat_strong_out, long_strong_out, N2_strong_out = [], [], [], [], [], [], []
    diss_weak_in, press_weak_in, time_weak_in, wind_weak_in, lat_weak_in, long_weak_in, N2_weak_in = [], [], [], [], [], [], []
    diss_weak_out, press_weak_out, time_weak_out, wind_weak_out, lat_weak_out, long_weak_out, N2_weak_out = [], [], [], [], [], [], []
    

    # Loop over all datasets
    for ds in list_of_data:
        num_profiles = len(ds['dissipation'].values)
        for i in range(num_profiles):
        # Loop over each profile in the dissipation data (across profiles)    
            long_i = ds['profile_long'].values[i, 0]
            lat_i = ds['profile_lat'].values[i, 0]
            diss_i = ds['dissipation'].values[i, :]
            press_i = ds['pressure'].values[i,:] 
            N2_i = ds['N_squared'].values[i, :]
            
            epoch_time = ds['start_time'].values[i, 0]
            
            if 0 < ds['wind'].values[i,0] < 100:    
                wind_i = ds['wind'].values[i, 0]
            else:
                print(f"Skipping bad wind value: {wind_i} at index {i}")
                continue  # Skip to the next profile
    
            # Only proceed if the timestamp is within a reasonable time range
            if 1e8 < epoch_time < 2e9:  
                utc_time = pd.to_datetime(epoch_time, unit='s', utc=True)
                local_time = utc_time.tz_convert('America/Los_Angeles')
                time_i = local_time
            else:
                print(f"Skipping bad timestamp: {epoch_time} at index {i}")
                continue  # Skip to the next profile
    
            mask = press_i < max_pressure
            diss_i = diss_i[mask]
            press_i = press_i[mask] 
            N2_i = N2_i[mask]
    
    
            if wind_i >= wind_threshold and long_i >= bay_threshold:
                diss_strong_in.append(diss_i.tolist()) 
                press_strong_in.append(press_i.tolist())
                N2_strong_in.append(N2_i.tolist())
                time_strong_in.append(time_i)
                wind_strong_in.append(wind_i)
                lat_strong_in.append(lat_i)
                long_strong_in.append(long_i)
            elif wind_i >= wind_threshold and long_i < bay_threshold:
                diss_strong_out.append(diss_i.tolist())
                press_strong_out.append(press_i.tolist())
                N2_strong_out.append(N2_i.tolist())
                time_strong_out.append(time_i)
                wind_strong_out.append(wind_i)
                lat_strong_out.append(lat_i)
                long_strong_out.append(long_i)
            elif wind_i < wind_threshold and long_i >= bay_threshold:
                diss_weak_in.append(diss_i.tolist())
                press_weak_in.append(press_i.tolist())
                N2_weak_in.append(N2_i.tolist())
                time_weak_in.append(time_i)
                wind_weak_in.append(wind_i)
                lat_weak_in.append(lat_i)
                long_weak_in.append(long_i)
            elif wind_i < wind_threshold and long_i < bay_threshold:
                diss_weak_out.append(diss_i.tolist())
                press_weak_out.append(press_i.tolist())
                N2_weak_out.append(N2_i.tolist())
                time_weak_out.append(time_i)
                wind_weak_out.append(wind_i)
                lat_weak_out.append(lat_i)
                long_weak_out.append(long_i)
               
    df_strong_in = pd.DataFrame({'time': time_strong_in,'wind': wind_strong_in, 'N_squared': N2_strong_in,
        'dissipation': diss_strong_in,'pressure': press_strong_in, 'latitude': lat_strong_in, 'longitude': long_strong_in})
    
    df_strong_out = pd.DataFrame({'time': time_strong_out,'wind': wind_strong_out, 'N_squared': N2_strong_out,
        'dissipation': diss_strong_out,'pressure': press_strong_out, 'latitude': lat_strong_out, 'longitude': long_strong_out})
    
    df_weak_in = pd.DataFrame({'time': time_weak_in,'wind': wind_weak_in, 'N_squared': N2_weak_in,
        'dissipation': diss_weak_in,'pressure': press_weak_in, 'latitude': lat_weak_in, 'longitude': long_weak_in})
    
    df_weak_out = pd.DataFrame({'time': time_weak_out,'wind': wind_weak_out, 'N_squared': N2_weak_out,
        'dissipation': diss_weak_out,'pressure': press_weak_out, 'latitude': lat_weak_out, 'longitude': long_weak_out}) 

    return df_strong_in, df_strong_out, df_weak_in, df_weak_out




def location_divisions(df_strong_in, df_strong_out, df_weak_in, df_weak_out):
    """
    Concatenate DataFrames into 'in bay' and 'out of bay' categories.

    Parameters:
    -----------
    df_strong_in : DataFrame
        DataFrame for strong wind conditions inside the bay.
    df_strong_out : DataFrame
        DataFrame for strong wind conditions outside the bay.
    df_weak_in : DataFrame
        DataFrame for weak wind conditions inside the bay.
    df_weak_out : DataFrame
        DataFrame for weak wind conditions outside the bay.

    Returns:
    --------
    df_in : DataFrame
        Concatenated DataFrame of all 'inside bay' profiles.
    df_out : DataFrame
        Concatenated DataFrame of all 'outside bay' profiles.
    """
    frames_in = [df_strong_in, df_weak_in]
    frames_out = [df_strong_out, df_weak_out]

    df_in = pd.concat(frames_in)
    df_out = pd.concat(frames_out)

    return df_in, df_out
    

def time_divisions(df_strong_in,df_strong_out, df_weak_in, df_weak_out):
    
    """
    Divides input DataFrames into morning and afternoon subsets based on local time.

    Parameters:
    -----------
    df_strong_in : DataFrame
        Strong wind profiles inside the bay.
    df_strong_out : DataFrame
        Strong wind profiles outside the bay.
    df_weak_in : DataFrame
        Weak wind profiles inside the bay.
    df_weak_out : DataFrame
        Weak wind profiles outside the bay.

    Returns:
    --------
    df_in_morning : DataFrame
        Inside-bay profiles collected before noon.
    df_out_morning : DataFrame
        Outside-bay profiles collected before noon.
    df_in_afternoon : DataFrame
        Inside-bay profiles collected at or after noon.
    df_out_afternoon : DataFrame
        Outside-bay profiles collected at or after noon.
    """
    
    frames_in = [df_strong_in, df_weak_in]
    frames_out = [df_strong_out, df_weak_out]
    df_in = pd.concat(frames_in)
    df_out = pd.concat(frames_out)
    
    idx_morning_out = df_out['time'].dt.hour < 12
    idx_afternoon_out = df_out['time'].dt.hour >= 12
    idx_morning_in = df_in['time'].dt.hour < 12
    idx_afternoon_in = df_in['time'].dt.hour >= 12
    
    df_in_morning = df_in[idx_morning_in]
    df_out_morning = df_out[idx_morning_out]
    df_in_afternoon = df_in[idx_afternoon_in]
    df_out_afternoon = df_out[idx_afternoon_out]

    return df_in_morning, df_out_morning, df_in_afternoon, df_out_afternoon


def dissipation_profile_means(df_dict):
    '''
    Compute the mean value for diisipation
    for each profile in each group of DataFrames.

    Parameters:
    -----------
    df_dict : dict
        Dictionary of DataFrames grouped by division conditions.
        Example: {'strong_in': df_strong_in, 'weak_out': df_weak_out, ...}

    Returns:
    --------
    mean_dict : dict
        Dictionary of array of profile means.
    '''
    mean_dict = {}

    for label, df in df_dict.items():
        profile_means = []
        for profile in df['dissipation']:
            profile_array = np.array(profile)
            profile_array = np.log10(profile_array)
            mean_val = np.nanmean(profile_array)
            profile_means.append(mean_val)
        
        # Clean out NaNs and store result
        profile_means = np.array(profile_means)
        profile_means = profile_means[~np.isnan(profile_means)]
        mean_dict[label] = profile_means

    return mean_dict


def mask_flatten(divided_data_frame):
    """
    Flattens profile-based dissipation and stratification data into point-wise lists for plotting or analysis.

    Parameters:
    -----------
    divided_data_frame : pd.DataFrame
        DataFrame with columns ['dissipation', 'pressure', 'longitude', 'N_squared']

    Returns:
    --------
    longitudes : list of float
    depths : list of float
    log_dissipation : list of float (log10 of dissipation)
    N_squareds : list of float
    dissipation : list of float
    """
    df_sorted = divided_data_frame.sort_values('longitude').reset_index(drop=True)

    # Initialize flattened lists
    longitudes = []
    depths = []
    log_dissipation = []
    N_squareds = []
    dissipation = []

    for i, row in df_sorted.iterrows():
        diss = np.array(row['dissipation'])
        depth = np.array(row['pressure'])
        N2 = np.array(row['N_squared'])
        lon = row['longitude']

        # Filter out profiles with mismatched array lengths or invalid longitude
        if len(diss) != len(depth) or len(diss) != len(N2) or np.isnan(lon):
            continue

        # Compute log10 of dissipation
        log_eps = np.log10(diss)

        # Create valid mask
        mask = np.isfinite(log_eps) & np.isfinite(depth) & np.isfinite(N2)

        # Append flattened values
        longitudes.extend([lon] * np.sum(mask))
        depths.extend(depth[mask])
        log_dissipation.extend(log_eps[mask])
        N_squareds.extend(N2[mask])
        dissipation.extend(diss[mask])

    return longitudes, depths, log_dissipation, N_squareds, dissipation

def turbulence_intensity_parameter(divided_data_frame):

    """
    Computes the turbulence intensity parameter I = ε / (N² * ν)

    Parameters:
    -----------
    divided_data_frame : pd.DataFrame
        DataFrame with 'dissipation', 'pressure', and 'N_squared' columns (each row is a profile)

    Returns:
    --------
    I : np.ndarray
        Flattened array of turbulence intensity values
    """
    
    df_sorted = divided_data_frame.sort_values('longitude').reset_index(drop=True)

    # Initialize flattened lists
    longitudes = []
    depths = []
    log_dissipation = []
    N_squareds = []
    dissipation = []
    viscosity = 10**(-6)

    for i, row in df_sorted.iterrows():
        diss = np.array(row['dissipation'])
        depth = np.array(row['pressure'])
        N2 = np.array(row['N_squared'])
        lon = row['longitude']

        # Filter out profiles with mismatched array lengths or invalid longitude
        if len(diss) != len(depth) or len(diss) != len(N2) or np.isnan(lon):
            continue

        # Compute log10 of dissipation
        log_eps = np.log10(diss)

        # Create valid mask
        mask = np.isfinite(log_eps) & np.isfinite(depth) & np.isfinite(N2)

        # Append flattened values
        longitudes.extend([lon] * np.sum(mask))
        depths.extend(depth[mask])
        log_dissipation.extend(log_eps[mask])
        N_squareds.extend(N2[mask])
        dissipation.extend(diss[mask])
    
    I = np.array(dissipation)/(np.array(N_squareds) *viscosity)
    return I                           



#======================================= Plot functions to make notebook cleaner, axis are hardcoded and not made for general use =======================================

def strength_location_division_plots(divided_data_frame, title=''):
    plt.figure(figsize=(12, 10))
    
    plt.subplot(2,2,1)
    for i in range(len(divided_data_frame['N_squared'])):
        plt.plot((divided_data_frame['N_squared'][i]), divided_data_frame['pressure'][i], linewidth=2, marker='o')
        
    plt.gca().invert_yaxis() 
    plt.xlabel('N^2 [s^-2]')
    plt.ylabel('Pressure (dbar)')
    plt.title(title)
    plt.show()

def time_divisions_plot(time_division_df,title=''):
    for i in range(len(time_division_df['dissipation'])):
        plt.plot(np.log10(np.array(time_division_df['dissipation'].iloc[i])), np.array(time_division_df['pressure'].iloc[i]), linewidth=2, marker='o')
        
    plt.gca().invert_yaxis() 
    plt.xlabel('$\log_{10}\epsilon$ [W/kg]')
    plt.ylabel('Pressure (dbar)')
    plt.title(title)
    plt.show()


def ship_vs_buoy_timeseries(location_dividion_df, buoy_df, title=''):
    df_buoy = buoy_df
    buoy_index = (location_dividion_df['datetime_local'] >= location_dividion_df['time'].min()) & 
    (location_dividion_df['datetime_local'] <= location_dividion_df['time'].max())

    plt.plot(df_buoy['datetime_local'][buoy_index], df_buoy['wndsd'][buoy_index], linewidth=2, color = 'orange', label = 'buoy wind speed')
    plt.scatter( df_in['time'], df_in['wind'], label = 'ship wind speed')
    plt.ylabel('wind [m/s]')
    plt.legend()
    plt.title(title)
    plt.show()    


def slo_map():
    line_C_lat =   [35.1184, 35.1217, 35.1246, 35.1275, 35.1305, 35.1339, 35.1368, 35.1393, 35.1425]
    line_C_lon = line_C_coordinates =   [-120.7436, -120.739, -120.7352, -120.7311, -120.7271, -120.7223, -120.7185, -120.7151, -120.7108]  
    
    line_D_lat =   [35.1125, 35.1186, 35.1215, 35.1245, 35.1276, 35.1337, 35.1367]
    line_D_lon =   [-120.7223 , -120.7142, -120.7101, -120.7059, -120.7019, -120.6937, -120.6897]
    
    
    request = cimgt.GoogleTiles(style = 'satellite')
    
    plt.figure()
    ax = plt.axes(projection=ccrs.Mercator())
    
    #ax.coastlines()
    ax.add_image(request, 12) #lower number bigger pixels, bigger number higher resolution (max 18?)
    
    plt.scatter(line_C_lon, line_C_lat,color = 'gold', linewidth=1, marker='o',transform = ccrs.PlateCarree(), label = 'Line C')
    plt.scatter(line_D_lon, line_D_lat,color = 'yellowgreen', linewidth=1, marker='o',transform = ccrs.PlateCarree(), label = 'Line D')
    plt.plot([-120.755, -120.63], [35.16, 35.031], color= 'orangered',  linestyle='--', linewidth=1, transform = ccrs.PlateCarree() , label = 'Inner/Outer Bay Boundary')
    plt.scatter(-120.7408, 35.1698, color = 'r', marker = '*' , linewidth=2, transform = ccrs.PlateCarree(), label = 'MET Bouy')
    plt.scatter(-120.999, 34.937, color = 'r', marker = '*' , linewidth=2,transform = ccrs.PlateCarree())
    plt.legend()
    plt.title('San Luis Obispo Bay Microstructure Stations')
    
    gl = ax.gridlines(alpha=0.25, crs = ccrs.PlateCarree(), draw_labels = True)
    
    ax.set_extent(extent_slo, crs = ccrs.PlateCarree() )
    plt.show