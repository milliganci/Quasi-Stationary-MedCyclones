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
from scipy.stats import rankdata

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


## VARIABLE SELECTION

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
    data_dir = "/media/alice/Crucial X9/portal/data_UNIBE/Bern_data/"
    if var_name == "precip":
        prefix_var = data_dir + "data/precip/precip"
        file_date = var_time.strftime("%Y")
        fname = prefix_var + file_date + ".nc"
    elif (var_name == "precip_6h") or (var_name == "convprecip_6h"):
        prefix_var = data_dir + "processed_data/precip/"+var_name[:-3]
        file_date = var_time.strftime("%Y")
        fname = prefix_var + file_date + var_name[-3:] + ".nc"
    elif (var_name == "precip_24h"):
        prefix_var = data_dir + "processed_data/precip/"+var_name[:-4]
        file_date = var_time.strftime("%Y")
        fname = prefix_var + file_date + var_name[-4:] + ".nc"
    elif var_name == "convprecip_6h_extr":
        prefix_var = data_dir + "processed_data/precip/"+var_name[:-8]
        file_date = var_time.strftime("%Y")
        fname = prefix_var + file_date + "_6h_tpextr_cp08.nc"
    else:
        print("Variable name is not known")
    return fname


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


## ASSIGN QUANTILES to Full-Track (FT) and Along-Track (AT) stationarity definitions

def assign_quantiles_and_categories_FT(values, df_original, name):
    # compute quantile ranks (0 to 1)
    ranks = rankdata(values, method='ordinal')
    quantiles = np.round((ranks - 1) / (len(values) - 1), 3)

    # assign categories: 0 = unclassified, 1 = low, 2 = medium, 3 = high
    # quantiles can be changed to suit needs of user
    categories = np.zeros_like(quantiles)
    categories[quantiles <= 0.10] = 1
    categories[(quantiles >= 0.45) & (quantiles <= 0.55)] = 2
    categories[quantiles >= 0.90] = 3
    categories[np.isnan(values)] = np.nan

    return pd.DataFrame({
        'id': df_original['id'].unique()[:len(values)],  # match IDs
        f'{name}_v': np.round(values, 3),
        f'{name}_q': quantiles,
        f'{name}_c': categories
    })


def categorise_distances_AT(df_dist, df_original, tstep_window, prefix):
    # categorisation function, reusable for any dist DataFrame with ['id', 'lon', 'lat', 'dist_sum']
    # calculate percentiles (user may change these thresholds according to preference)
    all_dist = df_dist['dist_sum'].values
    perc_90 = np.percentile(all_dist, 90)
    perc_45 = np.percentile(all_dist, 45)
    perc_55 = np.percentile(all_dist, 55)
    perc_10 = np.percentile(all_dist, 10)

    # full values per ID with padding NaNs for trailing points
    full_values = []
    for ID_unique in np.unique(df_original.id.values):
        vals = np.array([v for v in df_dist.loc[df_dist['id'] == ID_unique, 'dist_sum']])
        vals = np.round(vals, 3)
        vals_padded = np.append(vals, [np.nan]*tstep_window)  # pad for alignment
        full_values.append(vals_padded)
    full_values_concat = np.concatenate(full_values)

    # calculate quantile ranks
    ranks = rankdata(full_values_concat, method='ordinal', nan_policy='omit')
    quantiles = np.where(np.isnan(full_values_concat), np.nan, (ranks - 1) / (len(full_values_concat[~np.isnan(full_values_concat)]) - 1))
    quantiles = np.round(quantiles, 3)

    # categorise based on quantiles
    categories = np.zeros_like(quantiles)
    categories[quantiles <= 0.1] = 1
    categories[(quantiles >= 0.45) & (quantiles <= 0.55)] = 2
    categories[quantiles >= 0.9] = 3
    categories[np.isnan(quantiles)] = np.nan

    # create DataFrame to join back to 'original' df_Medcrossers
    new_cols = pd.DataFrame({
        f'{prefix}_v': full_values_concat,
        f'{prefix}_q': quantiles,
        f'{prefix}_c': categories    
    })
    new_cols.index = df_original.index  # align index

    # return categorised dataframe and the three tracks filtered by percentile groups
    return new_cols