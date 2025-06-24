# Used python packages

import numpy as np
import datetime
import pandas as pd
import xarray as xr
from scipy.io import loadmat
import calendar
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker

## TRACK SELECTION

def open_tracks_flaounas(path_name):
    """open_tracks df_tracks,
    which contains the storm track information loaded from the original .dat file 
    of CL5 tracks provided by Flaounas et al. (2023).
    It takes no inputs and depends only upon the tracks,
    the path to which it is hardcoded in the function"""
    # Tracks dataset
    df_tracks = pd.read_csv(path_name, delim_whitespace=True, header=None)
    # Resetting track numbering to be consistent with clusters
    for icol in [0,3,4,5,6]:
        df_tracks[icol] = df_tracks[icol].astype(int)
    return df_tracks


def get_storms_sometime(df_tracks, id_storm, var_time):
    """
    From the two dataframes produced by open_tracks,
    selects all the time steps of the storms falling in the selected track IDs, and list of timesteps.
    Note that the track availability goes from 1979 to 2020.
    Setting in input id_storm=0 selects all storms.

    Parameters:
    df_tracks, df_clust: dataframes with track obtained from open_tracks_flaounas
    id_storm: storm ID, corresponding to column [0] in df_tracks.
    var_time: array of datetime64 time steps.

    Returns:
    df_select: selection of df_tracks based on id_storm and var_time. The ID count is reset.

    """

    # If number of storms is specified to a list or a non-zero value, else take all storms
    if isinstance(id_storm, list):
        df_select = df_tracks.loc[df_select[0].isin(id_storm)]
    else:
        if id_storm != 0:
            df_select = df_tracks.loc[df_tracks[0].isin([id_storm])]
        else:
            df_select = df_tracks
    # Select based on var_time
    if var_time != []:
        # Call function make_var_time
        time_select = make_var_time(df_select).astype("datetime64[h]")
        df_select = df_select[np.isin(time_select, var_time)]

    # Re-indexing
    df_select = df_select.reset_index(drop=True)
    return df_select


## TIME FUNCTIONS

def make_var_time(df_select):
    """
    Takes as inputs the dataframe of the selected storm tracks,
    and produces a corresponding datetime array var_time.
    """
    # List of row indices
    row_ind = df_select.index.values
    time_select = []
    for tt in row_ind:
        yy = df_select.loc[tt, 'year']
        mm = df_select.loc[tt, 'month']
        dd = df_select.loc[tt, 'day']
        hh = df_select.loc[tt, 'time']
        # datetime format
        time_select.append(datetime.datetime(year=yy, month=mm, day=dd, hour=hh))
    time_select = np.array(time_select, dtype="datetime64")
    return time_select


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


## HAVERSINE formula to compute distance (in km) between lon, lat points (vectors)
    
def haversine(lon1, lat1, lon2, lat2):
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


## FIND_NEAREST_INDEX: Formula
def find_nearest_index(lon_ref, lat_ref, lon_input, lat_input):
    xi = np.nanargmin((lon_ref-lon_input)**2)
    yi  = np.nanargmin((lat_ref-lat_input)**2)
    return (xi, yi)
