"""
PAL 9/8/2014
plot time series of p and T
"""

from scipy.io import netcdf as nc
import matplotlib.pyplot as plt
import numpy as np
import os
import datetime
import matplotlib.patches as mpatches
from matplotlib.dates import date2num


lat_cv, lon_cv = 38.237, -121.566
lat_bay, lon_bay = 38.053, -122.395
lat_colusa, lon_colusa = 39.157, -121.971
utc_offset = 8.  # hrs to PST
mintimestep = 48  # start time step to plot (30 min time steps)
maxtimestep = 144  # end time step to plot (30 min time steps)

def get_box_index(f, lat, lon):
    dif_lat = f.variables['XLAT'][:, :]-lat
    dif_lon = f.variables['XLONG'][:, :]-lon
    distance = dif_lat**2 + dif_lon**2
    ixlat, ixlon = np.unravel_index(np.argmin(distance), distance.shape)
    return ixlat, ixlon

def convert_times(f):
    times = []
    for ii in xrange(f.variables['Times'].shape[0]):
        time_arr = f.variables['Times'][ii, :]
        year = int(''.join(time_arr[:4]))
        month = int(''.join(time_arr[5:7]))
        day = int(''.join(time_arr[8:10]))
        hour = int(''.join(time_arr[11:13]))
        minute = int(''.join(time_arr[14:16]))
        second = int(''.join(time_arr[17]))
        times.append(datetime.datetime(year, month, day, hour, minute, second))
    return np.array(times) - datetime.timedelta(hours=utc_offset)

def calc_distance(ij1, ij2, f):
    # calculate distance between ij1 and ij2
    # i is lon/x, j is lat/y
    di = ij1[0]-ij2[0]
    dj = ij1[1]-ij2[1]
    dx = f.variables['DX'][0]
    dy = f.variables['DY'][0]
    distance = ( (di*dx)**2 + (dj*dy)**2 )**0.5
    distance = distance/1000.  # convert m -> km
    return distance

def plot_var(times, data, levels, ax):
    
    # data: time x level
    # levels: list of numeric vertical level heights

    colors = ['b', 'c', 'g', 'y', 'r', 'm', 'k']

    # loop through levels and plot
    handles = []
    labels = []
    for i_lev, lev in enumerate(levels):
        h = ax.plot(times[mintimestep:maxtimestep], data[mintimestep:maxtimestep, i_lev], color=colors[i_lev])
        handles.append(h[0])
        labels.append(str(int(lev))+" m")

    ax.set_xlabel('date in PST')

    return handles, labels

def plot_p(times, p_bay, p_cv, p_col, d_baycv, d_cvcol, levels):
    # p_bay, p_cv, and p_col are time x level

    fig, ax = plt.subplots(nrows=len(levels), ncols=2)

    pgrad_cv_bay = (p_cv - p_bay)/d_baycv
    pgrad_col_cv = (p_col - p_cv)/d_cvcol

    for i_lev, lev in enumerate(levels):

        ax[i_lev, 0].plot(times[mintimestep:maxtimestep], p_bay[mintimestep:maxtimestep, i_lev], 'b', label='bay')
        ax[i_lev, 0].plot(times[mintimestep:maxtimestep], p_cv[mintimestep:maxtimestep, i_lev], 'g', label='delta')
        ax[i_lev, 0].plot(times[mintimestep:maxtimestep], p_col[mintimestep:maxtimestep, i_lev], 'r', label='colusa')
        if i_lev == 0:
            ax[i_lev, 0].legend(loc='lower left', prop={'size':7})
        ax[i_lev, 0].set_ylabel(str(int(lev))+"m\np, hPa")

        ax[i_lev, 1].plot(times[mintimestep:maxtimestep], pgrad_cv_bay[mintimestep:maxtimestep, i_lev], 'c', label='delta-bay')
        ax[i_lev, 1].plot(times[mintimestep:maxtimestep], pgrad_col_cv[mintimestep:maxtimestep, i_lev], 'm', label='colusa-delta')
        if i_lev == 0:
            ax[i_lev, 1].legend(loc='lower left', prop={'size':7})
        ax[i_lev, 1].set_ylabel(str(int(lev))+"m\ngrad p, hPa/km")

    return fig, ax


if __name__ == "__main__":
    
    domains = ["d01", "d02"]
    root_dir = "/scratch2/scratchdirs/plink/WRF/output"

    run_name = "CA-0.2"

    p_plot = False
    t_plot = False
    p_t_plot = True

    plot_dir = os.path.join(root_dir, "timeseries_p_T")
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    for domain in domains:

        print domain

        # open nc file
        f = nc.netcdf_file(os.path.join(root_dir, run_name, "wrfout_interp_"+domain+".nc"))
        times = convert_times(f)

        # get coordinates
        ilat_cv, ilon_cv = get_box_index(f, lat_cv, lon_cv)
        ilat_bay, ilon_bay = get_box_index(f, lat_bay, lon_bay)
        ilat_colusa, ilon_colusa = get_box_index(f, lat_colusa, lon_colusa)

        ### P plot
        p_bay = f.variables["P"][:, :, ilat_bay, ilon_bay]
        p_cv = f.variables["P"][:, :, ilat_cv, ilon_cv]
        p_col = f.variables["P"][:, :, ilat_colusa, ilon_colusa]
        p_ocn = f.variables["P"][:, :, 49, 3]
        d_baycv = calc_distance((ilon_bay, ilat_bay), (ilon_cv, ilat_cv), f)
        d_cvcol = calc_distance((ilon_colusa, ilat_colusa), (ilon_cv, ilat_cv), f)

        if p_plot:
            fig, ax = plot_p(times, p_bay, p_cv, p_col, d_baycv, d_cvcol, f.variables["levels"][:])
            
            # save p fig (make it tall)
            fig.set_size_inches(10, 10)
            fig.suptitle(run_name+" "+domain)
            fig.savefig(os.path.join(plot_dir, "p_"+run_name+"_"+domain+".png"))

            plt.show()
            #1/0
            plt.close(fig)

        ### T plot
        if t_plot:
            fig, ax = plt.subplots(nrows=3, ncols=1)
            
            #ylims_T = 
            handles, labels = plot_var(times, f.variables["T"][:, :, ilat_bay, ilon_bay], f.variables["levels"][:], ax[0])
            h = ax[0].plot(times[mintimestep:maxtimestep], f.variables["T2"][mintimestep:maxtimestep, ilat_bay, ilon_bay]-273., 'b--')
            handles.append(h[0])
            labels.append("2 m")
            h = ax[0].plot(times[mintimestep:maxtimestep], f.variables["TSK"][mintimestep:maxtimestep, ilat_bay, ilon_bay]-273., 'r--')
            handles.append(h[0])
            labels.append("skin")
            ax[0].legend(handles, labels, loc="lower right", prop={"size": 6})
            ax[0].set_ylabel("T bay")
            #ax[0].set_ylim(ylims_T)
            
            handles, labels = plot_var(times, f.variables["T"][:, :, ilat_cv, ilon_cv], f.variables["levels"][:], ax[1])
            ax[1].plot(times[mintimestep:maxtimestep], f.variables["T2"][mintimestep:maxtimestep, ilat_cv, ilon_cv]-273., 'b--')
            ax[1].plot(times[mintimestep:maxtimestep], f.variables["TSK"][mintimestep:maxtimestep, ilat_cv, ilon_cv]-273., 'r--')
            ax[1].set_ylabel("T cv")
            #ax[1].set_ylim(ylims_T)

            handles, labels = plot_var(times, f.variables["T"][:, :, ilat_colusa, ilon_colusa], f.variables["levels"][:], ax[2])
            ax[2].plot(times[mintimestep:maxtimestep], f.variables["T2"][mintimestep:maxtimestep, ilat_colusa, ilon_colusa]-273., 'b--')
            ax[2].plot(times[mintimestep:maxtimestep], f.variables["TSK"][mintimestep:maxtimestep, ilat_colusa, ilon_colusa]-273., 'r--')
            ax[2].set_ylabel("T colusa")
            #ax[2].set_ylim(ylims_T)

            # save T fig
            fig.set_size_inches(8, 9)
            fig.suptitle(run_name+" "+domain)
            fig.savefig(os.path.join(plot_dir, "T_"+run_name+"_"+domain+".png"))

            plt.show()
            #1/0
            plt.close(fig)

        ### P250 and T together
        if p_t_plot:
            fig, ax = plt.subplots(nrows=3, ncols=1)

            h = ax[0].plot(times[mintimestep:maxtimestep], p_bay[mintimestep:maxtimestep, 2]-p_ocn[mintimestep:maxtimestep, 2], \
                            'b', lw=1.5)
            ax[0].invert_yaxis()
            ax0T = ax[0].twinx()
            handles, labels = plot_var(times, f.variables["T"][:, :, ilat_bay, ilon_bay], f.variables["levels"][:], ax0T)
            handles.append(h[0])
            labels.append("P")
            ax0T.legend(handles, labels, loc="lower left", prop={"size": 6})
            ax[0].set_ylabel("P250")
            ax0T.set_ylabel("T")

            ax[1].plot(times[mintimestep:maxtimestep], p_cv[mintimestep:maxtimestep, 2]-p_ocn[mintimestep:maxtimestep, 2], \
                            'b', lw=1.5)
            ax[1].invert_yaxis()
            ax1T = ax[1].twinx()
            handles, labels = plot_var(times, f.variables["T"][:, :, ilat_cv, ilon_cv], f.variables["levels"][:], ax1T)
            ax[1].set_ylabel("P250")
            ax1T.set_ylabel("T")

            ax[2].plot(times[mintimestep:maxtimestep], p_col[mintimestep:maxtimestep, 2]-p_ocn[mintimestep:maxtimestep, 2], \
                            'b', lw=1.5)
            ax[2].invert_yaxis()
            ax2T = ax[2].twinx()
            handles, labels = plot_var(times, f.variables["T"][:, :, ilat_colusa, ilon_colusa], f.variables["levels"][:], ax2T)
            ax[2].set_ylabel("P250")
            ax2T.set_ylabel("T")

            fig.suptitle(run_name+" "+domain)

            plt.show()
            1/0
