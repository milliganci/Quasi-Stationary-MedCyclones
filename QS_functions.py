# Used python packages

import numpy as np
import datetime
import pandas as pd
import xarray as xr
from scipy.signal import convolve2d
from scipy.io import loadmat
import multiprocessing
import os.path
import calendar
from scipy import ndimage
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker

## TRACK SELECTION

def open_tracks_givon():
    """open_tracks creates two dataframes: df_tracks,
    which contains the storm track information loaded from the original .mat files,
    and df_clust, which matches cyclone IDs to cluster numbers.
    It takes no inputs and depends only upond the tracks and clustering data,
    the path to which is hardcoded in the function"""
    # Tracks dataset
    tracks_mat = loadmat("insert_your_path!!!/ind.mat")
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
    ind_mat = loadmat("insert_your_path!!!/ind.mat")
    ind_mat = np.transpose(ind_mat["ind"])
    df_clust = pd.DataFrame(np.argwhere(ind_mat) + 1)
    df_clust.columns = ["Event Number", "Cluster Number"]
    return df_tracks, df_clust


def open_tracks_flaounas():
    """open_tracks df_tracks,
    which contains the storm track information loaded from the original .dat file 
    of CL5 tracks provided by Flaounas et al. (2023).
    It takes no inputs and depends only upon the tracks,
    the path to which it is hardcoded in the function"""
    # Tracks dataset
    df_tracks = pd.read_csv("insert_your_path!!!/TRACKS_CL5.dat", delim_whitespace=True, header=None)
    # Resetting track numbering to be consistent with clusters
    for icol in [0,3,4,5,6]:
        df_tracks[icol] = df_tracks[icol].astype(int)
    # df_tracks.columns = ['Event Number','lon','lat','year','mon','day','hour','MSLP']
    return df_tracks


def get_storms_alltime(df_tracks, id_storm, year_range):
    """
    From the two dataframes produced by open_tracks,
    selects all the time steps of the storms falling in the selected clusters, track IDs, and range of years.
    Note that the track availability goes from 1979 to 2020.
    Setting in input cluster_number=0 selects all clusters and n_storm=0 selects all storms.

    Parameters:
    df_tracks: dataframe with track obtained from open_tracks_flaounas
    id_storm: storm ID, corresponding to column [0] in df_tracks.
    year_range: list of two years (numbers) as in [year_first, year_last].

    Returns:
    df_select: selection of df_tracks based on cluster_number and id_storm. The ID count is reset.
    """
    # Select only years included in year_range
    df_select = df_tracks.loc[
        (df_tracks[3] >= year_range[0]) & (df_tracks[3] <= year_range[1])
    ]
    
    # Select specific storm numbers if id_storm is specified to a non-zero value.
    # Else take all storms.
    if isinstance(id_storm, list):
        df_select = df_select.loc[df_select[0].isin(id_storm)]
    else:
        if id_storm != 0:
            df_select = df_select.loc[df_select[0].isin([id_storm])]
    
    # Reset index to make tracks IDs increase incrementally by 1.
    df_select = df_select.reset_index(drop=True)
    
    return df_select


def get_storms_sometime(df_tracks, id_storm, var_time):
    """
    From the two dataframes produced by open_tracks,
    selects all the time steps of the storms falling in the selected clusters, track IDs, and list of timesteps.
    Note that the track availability goes from 1979 to 2020.
    Setting in input cluster_number=0 selects all clusters and n_storm=0 selects all storms.

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
            df_select = df_tracks.loc[df_select[0].isin([id_storm])]
    # Select based on var_time
    if var_time != []:
        # Call function make_var_time
        time_select = make_var_time(df_select).astype("datetime64[h]")
        df_select = df_select[np.isin(time_select, var_time)]

    # Re-indexing
    df_select = df_select.reset_index(drop=True)
    return df_select


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


def indices_tseries_pmin(df_tracks, valid_hours, year_range, time_steps):
    """
    Indices of timesteps around the minimum pressure of each storm in the selected cluster.
    Identifies indices 'ts_ind' relative to the dataframe object 'df_select_vh'.
    Non existing timesteps with respect to time of minimum pressure are returned as nan values. 
    """

    # Select only years that are included in the dataset
    df_select = df_tracks.loc[(df_select[3] >= year_range[0]) & (df_select[3] <= year_range[1])]
    # Extracting all lines where time of the day is compatible with whatever data we want to process
    df_select_vh = df_select.loc[df_select[6].isin(valid_hours)]
    
    # Extract timeseries around ALL-TIME minimum pressure for each event ii
    # for each event_ID take the pmin index + [time_steps]
    list_ts_ind = []
    event_ID = np.unique(np.array(df_select[0]))
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
    # Convert -1 to nan
    ts_ind = np.where(np.array(list_ts_ind)==-1,np.nan,np.array(list_ts_ind))

    return df_select_vh, ts_ind


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
        yy = df_select.loc[tt, 3]
        mm = df_select.loc[tt, 4]
        dd = df_select.loc[tt, 5]
        hh = df_select.loc[tt, 6]
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
    var_array = ds[str_to_name(var_name)].sel(time=var_date).data
    return var_array



## FUNCTIONS TO READ / OPEN VARIBLES

def generate_fname(var_name, var_time):   ## NOTE: modify with path to the data
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
        prefix_var = "insert_your_path!!!//precip"
        file_date = var_time.strftime("%Y")
        fname = prefix_var + file_date + ".nc"
    elif (var_name == "precip_6h") or (var_name == "precip_24h"):
        prefix_var = "insert_your_path!!!/precip/"+var_name[:-3]
        file_date = var_time.strftime("%Y")
        fname = prefix_var + file_date + var_name[-3:] + ".nc"
    else:
        print("Variable name is not known")
    return fname


def str_to_name(varstr):
    """
    Assign variable name to collable name (for dataset structure)
    """
    if (varstr == "precip_6h") or (varstr == "precip_24h") or (varstr == ("precip")):
        varname = "tp"
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


## OTHER FUNCTIONS

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


def subplot_Europe(
    fig, nums_plot, lat, lon, field_s, dict_s, field_c, dict_c, ext_plot, title, title_font
):
    """
    Draw subplot of standard Europe map, with shading (and contours)

    Parameters
    ------------
    nums_plot: [nrows,ncols,id_subplot]
    lon, lat: vectors
    field_s,  dict_s: field, dictionary for shading
    field_c, dict_c: field, dictionary for contours
    ext_plot: [lonmin,lonmax,lonmin,lonmax]
    title: title of individual subplot
    """
    ax = fig.add_subplot(
        nums_plot[0], nums_plot[1], nums_plot[2], projection=ccrs.PlateCarree()
    )
    if ext_plot == []:
        ax.set_extent([-30, 55, 10, 70], crs=ccrs.PlateCarree())
    else:
        ax.set_extent(ext_plot, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE)
    ax.set_title(title, fontsize=title_font)

    # shading
    if "colors" in dict_s:
        h_s = ax.contourf(
            lon,
            lat,
            field_s,
            colors=dict_s["colors"],
            levels=dict_s["levels"],
            transform=ccrs.PlateCarree(),
            zorder=0,
            extend=dict_s["extend"],
            spacing=dict_s["spacing"],
        )
    elif "extend" in dict_s:
        h_s = ax.contourf(
            lon,
            lat,
            field_s,
            cmap=dict_s["cmap"],
            levels=dict_s["levels"],
            transform=ccrs.PlateCarree(),
            zorder=0,
            extend=dict_s["extend"],
            spacing=dict_s["spacing"],
        )
    else:
        h_s = ax.contourf(
            lon,
            lat,
            field_s,
            cmap=dict_s["cmap"],
            levels=dict_s["levels"],
            transform=ccrs.PlateCarree(),
            zorder=0,
        )
    # white
    if dict_s["levels"].min() > 0:
        ax.contourf(
            lon,
            lat,
            field_s,
            levels=[0, dict_s["levels"].min()],
            colors=dict_s["cfill"],
            transform=ccrs.PlateCarree(),
        )
    # contours
    h_c = ax.contour(
        lon,
        lat,
        field_c,
        colors=dict_c["color"],
        levels=dict_c["levels"],
        linewidths=dict_c["linewidth"],
        alpha=dict_c["alpha"],
        transform=ccrs.PlateCarree(),
    )
    if dict_c["clabel"] == True:
        ax.clabel(h_c)
    return ax, h_s, h_c


def add_Europe(ax, lat, lon, fl_s, field_s, dict_s, fl_c, field_c, dict_c):
    """
    Adds shading and/or contours to pre-existing axis

    Parameters
    ------------
    ax: axis where to plot
    lon, lat: vectors
    fl_s, field_s, dict_s: if fl_s is True, plots field using dictionary (shading)
    fl_c, field_c, dict_c: if fl_c is True, plots field using dictionary (contours)
    """
    if fl_s == True:
        # shading
        h_s = ax.contourf(
            lon,
            lat,
            field_s,
            cmap=dict_s["cmap"],
            levels=dict_s["levels"],
            transform=ccrs.PlateCarree(),
            zorder=0,
        )
        # white
        ax.contourf(
            lon,
            lat,
            field_s,
            levels=[0, dict_s["levels"].min()],
            colors=dict_s["cfill"],
            transform=ccrs.PlateCarree(),
        )
    else:
        h_s = None

    if fl_c == True:
        # contours
        h_c = ax.contour(
            lon,
            lat,
            field_c,
            colors=dict_c["color"],
            levels=dict_c["levels"],
            linewidths=dict_c["linewidth"],
            alpha=dict_c["alpha"],
            transform=ccrs.PlateCarree(),
        )
        if dict_c["clabel"] == True:
            ax.clabel(h_c)
    else:
        h_c = None

    return ax, h_s, h_c

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

