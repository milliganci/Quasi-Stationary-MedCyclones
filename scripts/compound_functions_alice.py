# -*- coding: utf-8 -*-
"""
Created on Thursday Sept 15 2022

@author: Raphael
"""

import numpy as np
import datetime
import pandas as pd
import xarray as xr
from scipy.signal import convolve2d
from scipy.io import loadmat
import multiprocessing
import compound_params as pr
import os.path
import calendar
from scipy import ndimage


## FUNCTIONS TO READ / SELECT CYCLONE TRACKS

def open_tracks():
    """open_tracks creates two dataframes: df_tracks,
    which contains the storm track information loaded from the original .mat files,
    and df_clust, which matches cyclone IDs to cluster numbers.
    It takes no inputs and depends only upond the tracks and clustering data,
    the path to which is hardcoded in the function"""
    # Tracks dataset
    tracks_mat = loadmat("/scratch2/aportal/data/tracks/Filtered_Tracks.mat")
    tmp = tracks_mat["Filtered_Tracks"]
    # Resetting track numbering to be consistent with clusters
    for ii, track_ID in enumerate(np.unique(tmp[:, 0])):
        tmp[tmp[:, 0] == track_ID, 0] = ii + 1
    df_tracks = pd.DataFrame(tmp)
    df_tracks[0] = df_tracks[0].astype(int)
    df_tracks[3] = df_tracks[3].astype(int)
    df_tracks[4] = df_tracks[4].astype(int)
    df_tracks[5] = df_tracks[5].astype(int)
    df_tracks[6] = df_tracks[6].astype(int)
    # df_tracks.columns = ['Event Number','lon','lat','year','mon','day','hour','MSLP']
    # Clustering table
    ind_mat = loadmat("/scratch2/aportal/data/tracks/ind.mat")
    ind_mat = np.transpose(ind_mat["ind"])
    df_clust = pd.DataFrame(np.argwhere(ind_mat) + 1)
    df_clust.columns = ["Event Number", "Cluster Number"]
    return df_tracks, df_clust


def get_storms_alltime(df_tracks, df_clust, cluster_number, id_storm, year_range):
    """takes into input the two dataframes produced by open_tracks,
    along with the selected cluster and number of tracks.
    Setting in input cluster_number=0 selects all clusters and n_storm=0 selects all storms.
    The function also uses the compound_params.py to evaluate where the data is available,
    and then outputs a dataframe of tracks containing the selected storms
    (in get_storms_minp number of storms),
    for the selected clusters at all times within y_range
    (in get_storms_minp filtering of times valid_hours).
    Notice that

    Parameters:

    Returns:

    """
    if isinstance(cluster_number, list):
        event_ID = np.array(
            df_clust.loc[df_clust["Cluster Number"].isin(cluster_number)][
                "Event Number"
            ]
        )
        df_select = df_tracks.loc[df_tracks[0].isin(event_ID)]
    else:
        if cluster_number != 0:
            event_ID = np.array(
                df_clust.loc[df_clust["Cluster Number"] == cluster_number][
                    "Event Number"
                ]
            )
            df_select = df_tracks.loc[
                df_tracks[0].isin(event_ID)
            ]  # Select all times of all events within a given cluster
        else:
            event_ID = np.array(df_clust["Event Number"])
            df_select = df_tracks
    # Second set of track selection criterion
    # Select only years that are included in the dust dataset
    df_select = df_select.loc[
        (df_select[3] >= year_range[0]) & (df_select[3] <= year_range[1])
    ]
    # Reset index to make track selection easier. The event_ID is retained in df_test[0]
    # Note that tracks IDs increasing incrementally by 1
    df_select = df_select.reset_index(drop=True)
    # If number of storms is specified to a non-zero value, else take all storms
    if isinstance(id_storm, list):
        df_select = df_select.loc[df_select[0].isin(id_storm)]
    else:
        if id_storm != 0:
            df_select = df_select.loc[df_select[0].isin([id_storm])]
    return df_select


def get_storms_sometime(df_tracks, df_clust, cluster_number, id_storm, var_time):
    """
    Takes into input the two dataframes produced by open_tracks,
    along with the selected cluster and number of tracks.
    Setting in input cluster_number=0 selects all clusters and n_storm=0 selects all storms.
    The function also uses the compound_params.py to evaluate where the data is available,
    and then outputs a dataframe of tracks containing the selected storms
    (in get_storms_minp number of storms),
    for the selected clusters at the timesteps within var_time.

    """
    # If cluster_number is specified by list or by a number different from 0, else take all clusters
    if isinstance(cluster_number, list):
        event_ID = np.array(
            df_clust.loc[df_clust["Cluster Number"].isin(cluster_number)][
                "Event Number"
            ]
        )
        df_select = df_tracks.loc[df_tracks[0].isin(event_ID)]
    else:
        if cluster_number != 0:
            event_ID = np.array(
                df_clust.loc[df_clust["Cluster Number"] == cluster_number][
                    "Event Number"
                ]
            )
            df_select = df_tracks.loc[
                df_tracks[0].isin(event_ID)
            ]  # Select all times of all events within a given cluster
        else:
            event_ID = np.array(df_clust["Event Number"])
            df_select = df_tracks

    # If number of storms is specified to a list or a non-zero value, else take all storms
    if isinstance(id_storm, list):
        df_select = df_select.loc[df_select[0].isin(id_storm)]
    else:
        if id_storm != 0:
            df_select = df_select.loc[df_select[0].isin([id_storm])]

    # Select based on var_time
    if var_time != []:
        time_select = make_var_time(df_select).astype("datetime64[h]")
        df_select = df_select[np.isin(time_select, var_time)]

    # Re-indexing
    df_select = df_select.reset_index(drop=True)
    return df_select


def indices_tseries_pmin(df_tracks, df_clust, cluster_number, valid_hours, year_range, n_storm, time_steps):
    """
    Indices of timesteps around the minimum pressure of each storm in the selected cluster.
    Identifies indices 'ts_ind' relative to the dataframe object 'df_select_vh'.
    Non existing timesteps with respect to time of minimum pressure are returned as nan values. 
    """
    # Get all events that fall within the selected cluster. If cluster_number==0: select all events
    if cluster_number != 0:
        event_ID = np.array(df_clust.loc[df_clust['Cluster Number'].isin(cluster_number)]['Event Number'])
        df_select = df_tracks.loc[df_tracks[0].isin(event_ID)] # Select all times of all events within a given cluster
    else: 
        event_ID = np.array(df_clust['Event Number'])
        df_select = df_tracks
    # Select only years that are included in the dataset
    df_select = df_select.loc[(df_select[3] >= year_range[0])&(df_select[3] <= year_range[1])]     
    # If number of storms is specified to a non-zero value, else take all storms
    if n_storm != 0:
        df_select = df_select.loc[0:n_storm-1]
    # Extracting all lines where time of the day is compatible with whatever data we want to process
    df_select_vh = df_select.loc[df_select[6].isin(valid_hours)]
    
    # Extract timeseries around ALL-TIME minimum pressure for each event ii
    # for each event_ID take the pmin index + [time_steps]
    list_ts_ind = []
    for ii in event_ID:    
        tmp = df_select.loc[df_select[0]==ii][7]              # pressure at all times for cyclone ID ii
        tmp_vh = df_select_vh.loc[df_select_vh[0]==ii][7]     # pressure at valid hours for cyclone ID ii
        # Index of valid time closest to pmin
        if any(tmp):
            pmin_ind = tmp.idxmin()                           # index of pmin
            ind_diff = tmp_vh.index - pmin_ind
            pmin_ind_vh = tmp_vh.index[np.absolute(ind_diff) == np.absolute(ind_diff).min()][0] 
            # Tseries indices (set to -1 on non-existing time steps)
            ind_th = pd.Index(pmin_ind_vh + time_steps).values
            ind_ev = df_select_vh.loc[df_select_vh.index.isin(ind_th)].index.values
            ind_th[np.isin(ind_th,ind_ev)==False] = -1
            ind_ev = list(ind_th)
            list_ts_ind.append(ind_ev)
    ts_ind = np.where(np.array(list_ts_ind)==-1,np.nan,np.array(list_ts_ind))
    return df_select_vh, ts_ind


def get_variable_storm_alltime(df_tracks, id_storm, list_var):
    """From df_tracks, get all-time variable of the cyclones selected by id_storm.

    Parameters:
    df_tracks: dataframe stormtracks.
    id_storm: list of id cyclones (1 to 3190).
    list_var: list of variables in columns "dataframe".

    Returns:
    var_tracks: variables in list_var for each cyclone at every time-step.

    """
    if isinstance(id_storm, list):
        df_select = df_tracks.loc[df_tracks[0].isin(id_storm)]
    else:
        if id_storm != 0:
            df_select = df_tracks.loc[df_tracks[0].isin([id_storm])]
    var_select = df_select.loc[:, list_var].to_numpy()
    return var_select



## TIME FUNCTIONS

def make_var_time(df_select):
    """make_var_time takes as inputs the dataframe of the selected storm tracks at time of minimum pressure,
    and produces a datetime array corresponding to the times of minimum pressure."""
    # List of row indices
    row_ind = df_select.index.values
    var_time = []
    for tt in row_ind:
        yy = df_select.loc[tt, 3]
        mm = df_select.loc[tt, 4]
        dd = df_select.loc[tt, 5]
        hh = df_select.loc[tt, 6]
        # datetime format
        var_time.append(datetime.datetime(year=yy, month=mm, day=dd, hour=hh))
    var_time = np.array(var_time, dtype="datetime64")
    return var_time


def timerange_datetime(t_start, t_end, t_int, t_res, list_mon):
    """Build ndarray of datetime64 ndarray, filtering months

    Parameters:
    t_start: t_str first time step %Y-%M-%D(T%h:%m:%s).
    t_end: t_str first excluded time step %Y-%M-%D(T%h:%m:%s).
    t_int: number to convert in timedelta64.
    t_res: resolution of t_int.
    list_mon: list months to be selected. All months are taken if list_mon=[]


    Returns:
    var_time: list of timesteps in datetime.datetime format.

    """
    t_int = np.timedelta64(t_int, t_res)
    var_time = np.arange(t_start, t_end, t_int, dtype="datetime64")
    if list_mon != []:
        months = var_time.astype("datetime64[M]").astype(int) % 12 + 1
        var_time = var_time[np.isin(months, list_mon)]
    return var_time



## FUNCTIONS TO READ / OPEN VARIBLES

def generate_fname(var_name, var_time):
    """Creates a string with the path to the file of the variable requested at the time requested.

    Parameters:
    var_name: name of the variable.
    var_time: datetime.datetime object.

    Returns:
    fname: name of the file (string).

    """
    # Set type var_time
    if type(var_time) != datetime.datetime:
        var_time = var_time.astype(datetime.datetime)

    # Define directory
    if var_name == "precip":
        prefix_var = "/scratch2/aportal/data/precip/precip"
        file_date = var_time.strftime("%Y")
        fname = prefix_var + file_date + ".nc"
    elif (var_name == "precip_6h") or (var_name == "convprecip_6h"):
        prefix_var = "/scratch2/aportal/processed_data/precip/"+var_name[:-3]
        file_date = var_time.strftime("%Y")
        fname = prefix_var + file_date + var_name[-3:] + ".nc"
    elif (var_name == "precip_24h"):
        prefix_var = "/scratch2/aportal/processed_data/precip/"+var_name[:-4]
        file_date = var_time.strftime("%Y")
        fname = prefix_var + file_date + var_name[-4:] + ".nc"
    elif var_name == "convprecip_6h_extr":
        prefix_var = "/scratch2/aportal/processed_data/precip/"+var_name[:-8]
        file_date = var_time.strftime("%Y")
        fname = prefix_var + file_date + "_6h_tpextr_cp08.nc"
    elif var_name == "windgust_6h":
        prefix_var = "/scratch2/aportal/processed_data/windgust_6h/windgust"
        file_date = var_time.strftime("%Y")
        fname = prefix_var + file_date + "_6h.nc"
    elif var_name == "swh_6h":
        prefix_var = "/scratch2/aportal/processed_data/swh_6h/swh"
        file_date = var_time.strftime("%Y")
        fname = prefix_var + file_date + "_6h.nc"
    elif var_name == "wind_mag":
        prefix_var = "/scratch2/aportal/data/wind_mag/wmag"
        file_date = var_time.strftime("%Y")
        fname = prefix_var + file_date + ".nc"
    elif var_name == "swell_wave_height":
        prefix_var = "/scratch2/aportal/data/swell_wave_height/swh_"
        file_date = var_time.strftime("%Y")
        fname = prefix_var + file_date + "_ERA5.nc"
    elif var_name=='pm10_24h':
        fname = '/scratch2/aportal/data/dust/CAMS_dust_movmean_24h_regridded.nc'
    elif var_name == "mslp":
        prefix_var = "/scratch2/aportal/data/mslp/mslp"
        file_date = var_time.strftime("%Y")
        fname = prefix_var + file_date + ".nc"
    elif var_name == "WCB":
        prefix_var = "/scratch2/aportal/data/WCB/"
        file_date = var_time.strftime("%Y%m")
        folder_date = var_time.strftime("%Y%m%d_%H")
        fname = prefix_var + file_date + "/" + "T" + folder_date
    elif var_name == "WCB_bool":
        prefix_var = "/scratch2/aportal/processed_data/WCB_bool/"
        folder_date = var_time.strftime("%Y")
        fname = prefix_var + "wcb_gt800-midtrop_bool_" + folder_date + "_mediterranean.nc"
    elif (var_name == "WCBin_bool") or (var_name == "WCBout_bool"):
        prefix_var = "/scratch2/aportal/processed_data/WCB_bool/"
        folder_date = var_time.strftime("%Y")
        fname = prefix_var + "wcb_inout_gt800_midtrop_bool_" + folder_date + "_mediterranean.nc"
    elif var_name == "DI_bool":
        prefix_var = "/scratch2/aportal/processed_data/DI_bool/"
        folder_date = var_time.strftime("%Y")
        fname = prefix_var + "di_gt700_bool_sm4_" + folder_date + "_mediterranean.nc"
    elif var_name == "cold_front_by_year":
        prefix_var = "/scratch2/aportal/data/front_data/by_year/"
        folder_date = var_time.strftime("%Y")
        fname = (
            prefix_var
            + "era5_fronts_expanded_cold_"
            + folder_date
            + "_mediterranean.nc"
        )
    elif (var_name == "cold_front") or (var_name == "warm_front"):
        fname = (
            "/scratch2/aportal/data/front_data/era5_fronts_expanded_"
            + var_name[0:4]
            + "_1979--2020_mediterranean.nc"
        )
    elif var_name[-5:] == "_cycl":
        if var_name[:2] == "IA":
            prefix_var = (
                "/scratch2/aportal/processed_data/impact-area/" + var_name[:4] + "_bool_"
            )
            folder_date = var_time.strftime("%Y")
        else:
            prefix_var = "/scratch2/aportal/processed_data/impact-area/dynfeats_bool_"
            folder_date = var_time.strftime("%Y")
        fname = prefix_var + folder_date + ".nc"
    elif var_name == "prob_lightning_remap":
        prefix_var = "/scratch2/aportal/processed_data/lightning/"
        folder_date = var_time.strftime("%Y")
        fname = prefix_var + "prob_lightning_"+folder_date+"_remapcon.nc"
    elif (var_name == 'lightning_ATDnet'):
        prefix_var = '/scratch2/aportal/data/lightning/ATDnet_'
        file_date = var_time.strftime('%Y')
        fname = prefix_var+file_date+'_025deg.nc'
    else:
        print("Variable name is not known")
    return fname


def str_to_name(varstr):
    """
    Assign variable to collable name (for dataset structure)
    """
    if varstr == "windgust_6h":
        varname = "fg10"
    elif (varstr == "precip") or (varstr == ("precip_6h")) or (varstr == ("precip_24h")):
        varname = "tp"
    elif (varstr == "convprecip_6h") or (varstr == "convprecip_6h_extr"):
        varname = "cp"
    elif varstr == "wind_mag":
        varname = "wmag"
    elif varstr == "swh_6h":
        varname = "swh"
    elif varstr == "pm10_24h":
        varname = "pm10"
    elif varstr == "mslp":
        varname = "mslp"
    elif varstr == "WCB":
        varname = "WCB"
    elif (varstr == "WCB_bool") or (varstr == "WCBin_bool")  or (varstr == "WCBout_bool"):
        varname = varstr[:-5]
    elif varstr == "DI_bool":
        varname = "DI"
    elif varstr == "DI":
        varname = "N"
    elif varstr == "prob_lightning_remap":
        varname = "prob_lightning"
    elif varstr == "lightning_ATDnet":
        varname = "Flash count"
    elif (
        (varstr == "cold_front_by_year")
        or (varstr == "cold_front")
        or (varstr == "warm_front")
    ):
        varname = "fronts"
    elif varstr[-5:] == "_cycl":
        varname = varstr[:-5]
    else:
        print("Variable is not in list")
    return varname


def open_var(var_name, var_time, fname_prev, ds_prev):
    """Returns the xarray dataset of the selected variable, within which the selected time is available.
    If fn_prev is empty or different from the output of generate_fname, reopens a dataset.

    Parameters:
    var_name: name of the variable.
    var_time: datetime.datetime object.
    fname_prev: name of the previous dataset file opened.
    ds_prev: previous dataset file opened.


    Returns:
    ds_var: xarray dataset.
    fname_var: name of the data file (string).

    """
    fname_var = generate_fname(var_name, var_time)
    if fname_var != fname_prev:
        ds_var = xr.open_dataset(fname_var)
    else:
        ds_var = ds_prev
    # return the xarray dataset and the path to the associated file
    return ds_var, fname_var


def process_var(ds, var_name, list_time, lat_obj, lon_obj, l_seltime):
    """Returns the processed xarray dataset ds_obj, with the var_name information condensed in a unique variable.
    If ds coordinates do not correspond to the given lon, lat coordinates, the dataset is interpolated to fit these.

    Parameters:
    ds: dataset.
    var_name: name of the objective variable.
    list_time: time of interest.
    lon_obj, lat_obj: objective coordinates.
    l_seltime: flag for time slicing.


    Returns:
    ds_obj: objective xarray dataset.

    """
    # Select time range to reduce dimension dataset (if necessary)
    if (var_name != "WCB") and (l_seltime == True):
        ds = ds.sel(time=slice(str(list_time[0])[:10], str(list_time[-1])[:10]))

    # Define and/or modify variable and coordinates (if necessary)
    if var_name == "WCB":
        ds0 = ds["MIDTROP"].squeeze("dimz_MIDTROP.INPUT").rename(var_name)
        ds0[var_name] = ds["GT800"].squeeze("dimz_GT800.INPUT") + ds["MIDTROP"].squeeze(
            "dimz_MIDTROP.INPUT"
        )
        ds = ds0.where(ds0[var_name] > 0.0, 1.0, 0.0)  # dataset of 0, 1
    elif var_name in ["cold_front", "warm_front"]:
        ds = ds.fillna(0.0)  # replace nan with 0
        ds = ds.rename({"longitude": "lon", "latitude": "lat"})
    elif var_name in ["cold_front_by_year", "WCB_bool", "WCBin_bool", "WCBout_bool", "DI_bool"]:
        ds = ds.fillna(0.0)  # replace nan with 0
    elif var_name=='prob_lightning_remap':
        vname = str_to_name(var_name)
        ds = ds.where((ds[vname] > 1E-10) | (ds[vname].isnull()), 0.0)  
    elif var_name=='lightning_ATDnet':
        vname = str_to_name(var_name)
        ds = ds.where((ds[vname] < 1.), 1) # lightning hours (1 or 0)

    # Interpolate (if necessary)
    nlat = len(lat_obj)
    nlon = len(lon_obj)
    if (
        ((len(ds["lat"]) == nlat) and (len(ds["lon"]) == nlon))
        and ((ds["lat"] == lat_obj).all())
        and ((ds["lon"] == lon_obj).all())
    ):
        ds_obj = ds
    else:
        if var_name in [
            "cold_front",
            "warm_front",
            "cold_front_by_year",
            "WCB",
            "WCB_bool",
            "WCBin_bool",
            "WCBout_bool",
            "DI",
            "DI_bool",
            "IA01_cycl",
            "IA02_cycl",
            "IA03_cycl",
            "r_500_cycl",
            "r_1000_cycl",
            "DI_cycl",
            "WCB_cycl",
            "WCBin_cycl",
            "WCBout_cycl",
            "fronts_cycl",
            "convprecip_6h_extr"
        ]:
            interp_met = "nearest"  # preserve boolean
        else:
            interp_met = "linear"
        ds_obj = ds.interp(
            lon=lon_obj,
            lat=lat_obj,
            method=interp_met,
            kwargs={
                "bounds_error": False,
                "fill_value": np.nan,
            },  # nan values outside input grid
        )
    return ds_obj


def select_time(ds, var_name, var_time):
    """Returns the array of the selected variable and time.

    Parameters:
    ds: input dataset.
    var_name: name of the variable in the file.
    var_time: datetime.datetime object.


    Returns:
    var_array: ndarray (lon,lat).

    """
    var_date = var_time.astype(datetime.datetime).strftime(
        "%Y-%m-%dT%H:%M:%S" + ".000000000"
    )
    if var_name == "WCB":
        var_array = ds[str_to_name(var_name)].sel(time=0.0).data
    else:
        var_array = ds[str_to_name(var_name)].sel(time=var_date).data
    return var_array




## FUNCTIONS TO COMPUTE EXTREME - CYCLONE STATISTICS

def dist_extreme_cyclones(pos_extr, time_extr, df_cycl, time_cycl):
    """Returns the distance between the position of the extreme
    and the closest cyclone track in df_cycl at time time_extr.

    Parameters:
    pos_extr: position of the extreme (lat, lon).
    time_extr: time of the extreme, type datetime.datetime.
    df_cycl: cyclone tracks, id at column 0, lon at column 1, lat at column 2.
    time_cycl: cyclone tracks times, type datetime.datetime.

    Returns:
    dist: distance of closest cyclone track in km.

    """
    n_cycl = np.count_nonzero(time_cycl == time_extr)
    dist = np.nan
    if n_cycl == 0:
        pass
    # 1 cyclone track
    elif n_cycl == 1:
        pos_cycl = df_cycl[time_cycl == time_extr].iloc[0, 1:3].iloc[::-1]  # lat, lon
        dist = geodesic(pos_extr, pos_cycl).km
    # more than 1 cyclone track
    elif n_cycl > 1:
        pos_cycl = (
            df_cycl[time_cycl == time_extr].iloc[:, 1:3].iloc[:, ::-1]
        )  # lat, lon
        # find closest cyclone
        for jj in range(n_cycl):
            dist = np.append(dist, geodesic(pos_extr, pos_cycl.iloc[jj]).km)
        dist = np.amin(dist)
    return dist


def count_compounds_impactarea(
    var_list,
    ds_var1,
    ds_var2,
    th_extr,
    l_ia,
    var_ia,
    ds_ia,
    t_step,
    df_cycl,
    time_cycl,
    r_cycl,
    grid_tuples,
    c_extr,
    c_comp,
    c_ia0,
    c_extr_ia0,
    c_comp_ia0,
    c_ial,
    c_extr_ial,
    c_comp_ial,
):
    """
    Returns the counts of extreme and compounds and compounds.and.cyclone impact area over lon and lat grid.

    Parameters:
    var_list: list of variables for computing extremes.
    var_ia: name of impact area variable
    th_extr: threshold for extreme definition.
    ds_var1, ds_var2, ds_ia: dataset of 2 variables and impact area on same lon, lat grid.
    l_ia: label indicating impact area type ('' is nothing, '00' is 1000 km circle, '01' is VERSION 01, see versions_impact-area.txt)
    t_step: current timestep, np.datetime64 format.
    df_cycl: cyclone dataframe.
    time_cycl: np.datetime64 array of cyclone timesteps, corresponding to cyclone instances in df_cycl.
    r_cycl: radius of the impact area '00' around the cyclone center.
    grid_tuples: tuples (lat,lon) for each grid point.
    c_extr, c_comp: bins counting extremes (2) and compounds (1) for each grid point.
    c_ia0, c_extr_ia0, c_comp_ia0: bins counting impact area '00' (1), extremes.and.IA00 (2) and compounds.and.IA00 (1) for each grid point.
    c_ial, c_extr_ial, c_comp_ial: bins counting impact area 'l' (1), extremes.and.IAl (2) and compounds.and.IAl (1) for each grid point.

    Returns:
    c_extr, c_comp
    c_ia0, c_extr_ia0, c_comp_ia0
    c_ial, c_extr_ial, c_comp_ial

    """
    # definition dimensions
    nvar = 2
    nlat = c_extr.shape[1]
    nlon = c_extr.shape[2]

    # select variables, count extremes and compounds
    var_arr = np.nan * np.ones((nvar, nlat, nlon))
    for iv, ds_var in enumerate([ds_var1, ds_var2]):
        var_arr[iv] = select_time(ds_var, var_list[iv], t_step)
        # count extremes
        mask_extr = np.where((var_arr[iv] > th_extr[iv]), 1, 0)
        c_extr[iv] += mask_extr
    # count compounds
    mask_comp = np.where((var_arr[0] > th_extr[0]) & (var_arr[1] > th_extr[1]), 1, 0)
    c_comp += mask_comp

    # count compound extremes and cyclone impact area
    mask_time = time_cycl == t_step
    n_cycl = np.count_nonzero(mask_time)
    if n_cycl == 0:
        pass
    elif l_ia == "00" or l_ia == "01":
        df_cycl = df_cycl[mask_time]
        pos_cycl = df_cycl.iloc[:, 1:3]
        if n_cycl == 1:
            distances = haversine(
                pos_cycl[1].squeeze(),
                pos_cycl[2].squeeze(),
                grid_tuples[:, 0],
                grid_tuples[:, 1],
            )
            mask_cycl = distances.reshape(nlat, nlon) < r_cycl
        else:
            mask_cycl = np.zeros((nlat, nlon)).astype(dtype=bool)
            for ic in range(n_cycl):
                distances = haversine(
                    pos_cycl.iloc[ic][1].squeeze(),
                    pos_cycl.iloc[ic][2].squeeze(),
                    grid_tuples[:, 0],
                    grid_tuples[:, 1],
                )
                mask_cycl = mask_cycl | (distances.reshape(nlat, nlon) < r_cycl)
        mask_cycl = np.where(mask_cycl, 1, 0)
        c_ia0 += mask_cycl
        c_extr_ia0 += np.where((mask_cycl[None, :] * c_extr == 1), 1, 0)
        c_comp_ia0 += np.where((mask_cycl * mask_comp == 1), 1, 0)
        # dynamical features
        if l_ia == "01":
            mask_df = select_time(ds_ia, var_ia, t_step)
            c_ial += np.where((mask_cycl == 1) | (mask_df == 1), 1, 0)
            mask_nocycl = np.where((mask_cycl == 0), 1, 0)
            c_extr_ial += np.where((mask_cycl[None, :] * mask_extr == 1), 1, 0) + \
                          np.where((mask_nocycl[None, :] * mask_df[None, :] * c_extr == 1), 1, 0)
            c_comp_ial += np.where((mask_cycl * mask_comp == 1), 1, 0) + \
                          np.where((mask_nocycl * mask_df * mask_comp == 1), 1, 0)
    else:
        print("Impact area was not selected. c_comp_ia0 and c_comp_ial untouched.")
    return c_extr, c_comp, c_ia0, c_extr_ia0, c_comp_ia0, c_ial, c_extr_ial, c_comp_ial


def count_extr_impactarea(
    var_name,
    ds_var,
    th_extr,
    var_ia,
    ds_ia,
    t_step,
    time_cycl,
    c_extr,
    c_ia,
    c_extr_ia,

):
    """
    Returns the counts of extremes and extremes.and.cyclone impact area over lon and lat grid.
    Differently from count_compounds_impactarea, this function uses the input impact area in ds_ia,
    and does not compute the interception with the circle around the cyclone center.

    Parameters:
    ds_var, ds_ia: dataset of variable and impact area on same lon, lat grid
    
    Returns:
    c_extr: extreme count.
    c_comp: compound extreme count.
    c_comp_ia: compound extreme.and.impact area count.

    """

    # select variables, count extremes
    var_arr = select_time(ds_var, var_name, t_step)
    # count extremes
    mask_extr = np.where((var_arr > th_extr), 1, 0).squeeze()
    c_extr += mask_extr

    # count compound extremes and cyclone impact area
    mask_time = time_cycl == t_step
    n_cycl = np.count_nonzero(mask_time)
    if n_cycl == 0:
        pass
    else:
        mask_ia = select_time(ds_ia, var_ia, t_step)
        mask_ia = np.where(mask_ia, 1, 0)
        c_ia += mask_ia
        c_extr_ia += np.where((mask_ia * mask_extr == 1), 1, 0)
    return c_extr, c_ia, c_extr_ia


def count_extr_comp_impactarea(
    var_list,
    ds_var1,
    ds_var2,
    th_extr,
    var_ia,
    ds_ia,
    t_step,
    time_cycl,
    c_extr,
    c_comp,
    c_ia,
    c_extr_ia,
    c_comp_ia,

):
    """
    Returns the counts of extreme and compounds and compounds.and.cyclone impact area over lon and lat grid.
    Differently from count_compounds_impactarea, this function uses the input impact area in ds_ia,
    and does not compute the interception with thecircle around the cyclone center.

    Parameters:
    ds_var1, ds_var2, ds_ia: dataset of 2 variables and impact area on same lon, lat grid
    
    Returns:
    c_extr: extreme count.
    c_comp: compound extreme count.
    c_comp_ia: compound extreme.and.impact area count.

    """
    # definition dimensions
    nvar = 2
    nlat = c_extr.shape[1]
    nlon = c_extr.shape[2]

    # select variables, count extremes and compounds
    var_arr = np.nan * np.ones((nvar, nlat, nlon))
    for iv, ds_var in enumerate([ds_var1, ds_var2]):
        var_arr[iv] = select_time(ds_var, var_list[iv], t_step)
    # count extremes
    mask_extr = np.where((var_arr > th_extr), 1, 0)
    c_extr += mask_extr
    # count compounds
    mask_comp = np.where((var_arr[0] > th_extr[0]) & (var_arr[1] > th_extr[1]), 1, 0)
    c_comp += mask_comp

    # count compound extremes and cyclone impact area
    mask_time = time_cycl == t_step
    n_cycl = np.count_nonzero(mask_time)
    if n_cycl == 0:
        pass
    else:
        mask_ia = select_time(ds_ia, var_ia, t_step)
        mask_ia = np.where(mask_ia, 1, 0)
        c_ia += mask_ia
        c_extr_ia += np.where((mask_ia[None, :, :] * mask_extr == 1), 1, 0)
        c_comp_ia += np.where((mask_ia * mask_comp == 1), 1, 0)
    return c_extr, c_comp, c_ia, c_extr_ia, c_comp_ia


def count_extr_comp_impactarea_conv(
    var_list,
    ds_var1_bool, # boolean
    ds_var2,
    th_extr2,
    var_ia,
    ds_ia,
    t_step,
    time_cycl,
    c_extr,
    c_comp,
    c_extr_ia,
    c_comp_ia,
):
    """
    Returns the counts of convective extremes and compounds and (extremes.or.compounds).and.cyclone impact area over lon and lat grid.
    Differently from count_compounds_impactarea, this function uses the input impact area in ds_ia,
    and does not compute the interception with thecircle around the cyclone center.

    Parameters:
    ds_var1, ds_var2: dataset of 2 variables (boolean and float) on same lon, lat grid
    t_step, time_cycl: current time step, cyclone time steps.

    Returns:
    c_extr: conv.extreme count.
    c_comp: conv.compound extreme count.
    c_extr_ia: conv.extreme.and.impact area count.
    c_comp_ia: conv.compound.and.impact area count.

    """
    # definition dimensions
    nvar = 2
    nlat = c_extr.shape[0]
    nlon = c_extr.shape[1]

    # select variables, count extremes and compounds
    var_arr2 = select_time(ds_var2, var_list[1], t_step)
    mask_extr2 = np.where((var_arr2 > th_extr2), 1, 0)
    
    # count convective extremes
    mask_cextr = select_time(ds_var1_bool, var_list[0], t_step)
    c_extr += mask_cextr

    # count convective compounds
    mask_ccomp = np.where((mask_extr2==1) & (mask_cextr==1), 1, 0)
    c_comp += mask_ccomp

    # count compound extremes and cyclone impact area
    mask_time = time_cycl == t_step
    n_cycl = np.count_nonzero(mask_time)
    if n_cycl == 0:
        pass
    else:
        mask_ia = select_time(ds_ia, var_ia, t_step)
        mask_ia = np.where(mask_ia, 1, 0)
        c_extr_ia += np.where((mask_ia * mask_cextr == 1), 1, 0)
        c_comp_ia += np.where((mask_ia * mask_ccomp == 1), 1, 0)
    return c_extr, c_comp, c_extr_ia, c_comp_ia


def mask_comp_impactarea(
    var_list,
    ds_var1,
    ds_var2,
    th_extr,
    var_ia,
    ds_ia,
    t_step,
    time_cycl,
    nlat, 
    nlon,
    nvar,
):
    """
    Returns masks of compounds, cyclone impact area and compounds.and.cyclone impact area over lon and lat grid.
    See method of count_extr_comp_impactarea.

    Parameters:
    ds_var1, ds_var2, ds_ia: dataset of 2 variables and impact area on same lon, lat grid

    Returns:
    masks
    """

    # select variables, mask compounds
    var_arr = np.nan * np.ones((nvar, nlat, nlon))
    for iv, ds_var in enumerate([ds_var1, ds_var2]):
        var_arr[iv] = select_time(ds_var, var_list[iv], t_step)
    # count compounds
    mask_comp = np.where((var_arr[0] > th_extr[0]) & (var_arr[1] > th_extr[1]), 1, 0)

    # count compound extremes and cyclone impact area
    mask_time = time_cycl == t_step
    n_cycl = np.count_nonzero(mask_time)
    if n_cycl == 0:
        pass
    else:
        mask_ia = select_time(ds_ia, var_ia, t_step)
        mask_ia = np.where(mask_ia, 1, 0)
        mask_comp_ia = np.where((mask_ia * mask_comp == 1), 1, 0)
    return mask_comp, mask_ia, mask_comp_ia



## FUNCTIONS WITH VARIABLE SPECIFIC OUTPUT

def compound_parameters(var_list, l_conv):
    """
    Assign titles and lower limits of quantile-based extremes to variables used to compute compounds
    """
    if var_list == ['precip_6h','windgust_6h'] and l_conv==False:
        qtls_ll = [.002,10] # lower-limit quantile
        var_title = ['Rain6h','Windgust6h']
    elif var_list == ['precip_6h','windgust_6h'] and l_conv==True:
        var_list[0] = 'convprecip_6h_extr'
        qtls_ll = [.8,10] # lower-limit quantile
        var_title = ['ConvRain6h','Windgust6h']
    elif var_list == ['swh_6h','windgust_6h']:
        qtls_ll = [2,10] # lower-limit quantile
        var_title = ['SWheight6h','Windgust6h']
        l_conv = False  # don't count wave-wind convective events
    elif var_list == ['pm10_24h','windgust_6h']:
        qtls_ll = [10E-9,10] # lower-limit quantile
        var_title = ['Dust24h','Windgust6h']
        l_conv = False  # don't count wave-wind convective events
    return qtls_ll, var_title, l_conv


def make_labels(mon_list, var_name, var_qt, var_llim):
    """
    Returns months and variable labels for name files and figures.

    Parameters:
    mon_list: list of month indices (int from 0 to 11).
    var_name: list of names used for variables.
    var_qt: list of quantiles (for each variable) used to compute the threshold for extreme identification.
    var_llim: list of lower limits (for each variable) of the quantile thresholds used for extreme identification.

    Returns:
    lab_mon: '_' + initial letter of each month
    lab_var: list of variable labels.
    """
    # define labels months
    mon_lab = ""
    if mon_list != []:
        mon_lab += "_"
    for mn in mon_list:
        mon_lab += calendar.month_abbr[mn][:1]
    # define labels variables
    var_lab = []
    for iv, vname in enumerate(var_name):
        vstr = vname.lower() + str(var_qt[iv])[-2:]
        if var_llim[iv] != 0:
            if vname == "Rain6h":
                vstr += "inf" + str(int(var_llim[iv] * 1000)) + "mm"
            elif vname == "ConvRain6h":
                vstr = "rain6h98inf2mmconv80pc"
            elif vname == "SWheight6h":
                vstr += "inf" + str(int(var_llim[iv])) + "m"
            elif vname == "Windgust6h":
                vstr += "inf" + str(int(var_llim[iv])) + "ms"
            elif vname == "Dust24h":
                vstr += "inf" + str(int(var_llim[iv] * 1E9)) + "ug"
            elif vname == "lightning_ATDnet":
                vstr += "gt1"
            elif vname == "LightningProb":
                vstr = vname.lower()
            else:
                print('Ups! Variable not recognised!')
        var_lab.append(vstr)
    return mon_lab, var_lab



## FIDDLY FUNCTIONS

def haversine(lon1, lat1, lon2, lat2):
    """
    Haversine formula to compute distance (in km) between lon, lat points (vectors)
    """
    # convert decimal degrees to radians
    lon1 = np.deg2rad(lon1)
    lon2 = np.deg2rad(lon2)
    lat1 = np.deg2rad(lat1)
    lat2 = np.deg2rad(lat2)

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    c = 2 * np.arcsin(np.sqrt(a))
    r = 6371
    return c * r


def smooth(y, box_pts):
    box = np.ones((box_pts, box_pts)) / box_pts**2
    y_smooth = convolve2d(y, box, mode="same", boundary="symm")
    return y_smooth


def string_combination(list_var, list_bool, list_labs):
    """Return a string describing the combination of Y/N variables."""
    n_comb = len(list_var)
    str_comb = ''
    for ii in range(n_comb):
        if list_bool[ii]==False:
            str_comb += 'NO'
        str_comb += list_labs[ii]+'-'
    str_comb = str_comb[:-1]
    return str_comb


def intercecting_connected_components(cc_in, cc_def, mask):
    """Separate cc_in in connected components, and keep the ones that intersect with the mask.
    Return:
    cc_out: connected components that intersect with the mask.
    """
    cc_out = np.zeros_like(cc_in)
    map_cc, n_cc = ndimage.label(cc_in, structure=cc_def)
    # interception cc and mask
    for i_cc in range(1,n_cc+1):
        mask_cc = (map_cc==i_cc)
        if np.any( mask * mask_cc ):
            cc_out += (mask_cc)
    return cc_out


def ds_ratio(ds1, var1, ds2, var2):
    """
    Returns ratio between var1 in ds1 and var2 in ds2
    """
    ds_ratio = ds1[var1] / ds2[var2]
    ds_ratio = ds_ratio.where(ds_ratio<np.Inf,0)
    return ds_ratio


def update_weight_sum(sum, ii_sum, weight, ii_weight):
    """Update the sum and weight of a variable."""
    if np.isnan(ii_sum)==False:
        weight += ii_weight
        sum += ii_sum * ii_weight
    return sum, weight

