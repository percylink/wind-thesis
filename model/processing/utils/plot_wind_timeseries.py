"""
PAL 9/8/2014
plot time series of wind at solano for each test
"""

from scipy.io import netcdf as nc
import matplotlib.pyplot as plt
import numpy as np
import os
import datetime
import matplotlib.patches as mpatches
from matplotlib.dates import date2num


lat_solano, lon_solano = 38.166, -121.817
utc_offset = 8.  # hrs to PST
burn_in = 48  # number of half-hour output timesteps that constitute burn-in period
#coords_solano = {'d01': (49, 51), 'd02': (91, 96), 'd03': (127, 143)}

def get_box_index(f, lat, lon):
    dif_lat = f.variables['XLAT'][:, :]-lat
    dif_lon = f.variables['XLONG'][:, :]-lon
    distance = dif_lat**2 + dif_lon**2
    ixlat, ixlon = np.unravel_index(np.argmin(distance), distance.shape)
    return ixlat, ixlon

def get_winds(f, coords_solano):
    usolano = f.variables['U'][:, :2, coords_solano[1], coords_solano[0]]
    vsolano = f.variables['V'][:, :2, coords_solano[1], coords_solano[0]]
    speed = (usolano**2 + vsolano**2)**0.5  # calculate speed
    return speed

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

def plot_wind(ax0, ax1, case):
    
    # ax0: axes obj for level 0 wind
    # ax1: axes obj for level 1 wind

    colors = ['g', 'r', 'c', 'm', 'y', 'k']

    # open control file
    fc = nc.netcdf_file(os.path.join(root_dir, case["ctrl"], "wrfout_interp_"+domain+".nc"))

    # find solano coords
    ixlat, ixlon = get_box_index(fc, lat_solano, lon_solano)
    coords_solano = (ixlon, ixlat)
    
    # get wind speed at solano and times in PST
    speed_c = get_winds(fc, coords_solano)
    times = convert_times(fc)

    # plot control
    h = ax0.plot(times, speed_c[:,0], 'b', lw=1.5)
    handles = [h[0]]
    labels = [case["ctrl"]]
    ax1.plot(times, speed_c[:,1], 'b', lw=1.5)

    # plot each test
    for i_test, test in enumerate(case["tests"]):
        print test
        # open test file
        ft = nc.netcdf_file(os.path.join(root_dir, test, "wrfout_interp_"+domain+".nc"))
        # get solano wind
        speed_t = get_winds(ft, coords_solano)
        # plot
        h = ax0.plot(times, speed_t[:, 0], colors[i_test])
        handles += h#[0]
        labels += [test]
        ax1.plot(times, speed_t[:, 1], colors[i_test])

    ax0.legend(handles, labels, prop={'size': 8}, loc='lower center', ncol=2)
    ax1.legend(handles, labels, prop={'size': 8}, loc='lower center', ncol=2)

    ax0.set_ylim(2, 14)
    ax1.set_ylim(2, 16)
    ax0.set_xlabel('date in PST')
    ax1.set_xlabel('date in PST')

    rect0 = mpatches.Rectangle([date2num(times[0]), 2], date2num(times[burn_in])-date2num(times[0]), 12, ec="none", fc="gray", alpha=0.5)
    ax0.add_patch(rect0)
    rect1 = mpatches.Rectangle([date2num(times[0]), 2], date2num(times[burn_in])-date2num(times[0]), 14, ec="none", fc="gray", alpha=0.5)
    ax1.add_patch(rect1)

    return handles, labels


if __name__ == "__main__":
    
    domains = ["d01", "d02"]
    root_dir = "/scratch2/scratchdirs/plink/WRF/output"

    case_wetrg = {"ctrl": "CA-0.1", "tests": ["CA-wetCR", "CA-wetCV", "CA-wetSN"], "test_id": "wet_regions"}
    case_dryrg = {"ctrl": "CA-0.25", "tests": ["CA-dryCR", "CA-dryCV", "CA-drySN"], "test_id": "dry_regions"}
    case_CVwet = {"ctrl": "CA-0.2", \
                  "tests": ["CA-CV0.05w", "CA-CV0.1w", "CA-CV0.15w", "CA-CV0.25w", "CA-CV0.3w", "CA-CV0.35w"], \
                  "smois": [0.05, 0.1, 0.15, 0.25, 0.3, 0.35], \
                  "smois_c": 0.2, \
                  "test_id": "CV_0.2"}
    case_CVdry = {"ctrl": "CA-0.1", \
                  "tests": ["CA-CV0.05", "CA-CV0.15", "CA-CV0.2", "CA-CV0.25", "CA-CV0.3", "CA-CV0.35"], \
                  "smois": [0.05, 0.15, 0.2, 0.25, 0.3, 0.35], \
                  "smois_c": 0.1, \
                  "test_id": "CV_0.1"}

    coords_solano = {}

    plot_dir = os.path.join(root_dir, "timeseries_wind")
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    for domain in domains:

        print domain

        # Regional tests
        fig1, ax1 = plt.subplots(nrows=2, ncols=1)  # level 0 (110 m)
        fig2, ax2 = plt.subplots(nrows=2, ncols=1)  # level 1 (150 m)

        # plot dry bkgd on top panels
        handles, labels = plot_wind(ax1[0], ax2[0], case_wetrg)
        ax1[0].set_ylabel('speed (m/s) at 60 m AGL')
        ax2[0].set_ylabel('speed (m/s) at 100 m AGL')
        #ax1[0].legend(handles, labels, prop={'size': 8}, loc='lower center', ncols=2)
        #ax2[0].legend(handles, labels, prop={'size': 8})
        ax1[0].set_title(domain+" dry background wet perturbed regions, Solano wind speed")
        ax2[0].set_title(domain+" dry background wet perturbed regions, Solano wind speed")

        # plot wet bkgd on bottom panels
        handles, labels = plot_wind(ax1[1], ax2[1], case_dryrg)
        ax1[1].set_ylabel('speed (m/s) at 60 m AGL')
        ax2[1].set_ylabel('speed (m/s) at 100 m AGL')
        #ax1[1].legend(handles, labels, prop={'size': 8})
        #ax2[1].legend(handles, labels, prop={'size': 8})
        ax1[1].set_title(domain+" wet background dry perturbed regions, Solano wind speed")
        ax2[1].set_title(domain+" wet background dry perturbed regions, Solano wind speed")

        fig1.set_size_inches(10, 7)
        fig2.set_size_inches(10, 7)
        fig1.subplots_adjust(hspace=0.35)
        fig2.subplots_adjust(hspace=0.35)
        fig1.savefig(os.path.join(plot_dir, "solano_wind_rg_"+domain+"_level0.png"))
        fig2.savefig(os.path.join(plot_dir, "solano_wind_rg_"+domain+"_level1.png"))

        #plt.show()
        #1/0
        plt.close(fig1)
        plt.close(fig2)

        # CV tests
        fig1, ax1 = plt.subplots(nrows=2, ncols=1)  # level 0 (110 m)
        fig2, ax2 = plt.subplots(nrows=2, ncols=1)  # level 1 (150 m)

        # plot wet bkgd on top panels
        handles, labels = plot_wind(ax1[0], ax2[0], case_CVwet)
        ax1[0].set_ylabel('speed (m/s) at 60 m AGL')
        ax2[0].set_ylabel('speed (m/s) at 100 m AGL')
        #ax1[0].legend(handles, labels, prop={'size': 8})
        #ax2[0].legend(handles, labels, prop={'size': 8})
        ax1[0].set_title(domain+" wet background vary CV, Solano wind speed")
        ax2[0].set_title(domain+" wet background vary CV, Solano wind speed")

        # plot dry bkgd on bottom panels
        handles, labels = plot_wind(ax1[1], ax2[1], case_CVdry)
        ax1[1].set_ylabel('speed (m/s) at 60 m AGL')
        ax2[1].set_ylabel('speed (m/s) at 100 m AGL')
        #ax1[1].legend(handles, labels, prop={'size': 8})
        #ax2[1].legend(handles, labels, prop={'size': 8})
        ax1[1].set_title(domain+" dry background vary CV, Solano wind speed")
        ax2[1].set_title(domain+" dry background vary CV, Solano wind speed")

        fig1.set_size_inches(10, 7)
        fig2.set_size_inches(10, 7)
        fig1.subplots_adjust(hspace=0.35)
        fig2.subplots_adjust(hspace=0.35)
        fig1.savefig(os.path.join(plot_dir, "solano_wind_CV_"+domain+"_level0.png"))
        fig2.savefig(os.path.join(plot_dir, "solano_wind_CV_"+domain+"_level1.png"))

        #plt.show()
        #1/0
        plt.close(fig1)
        plt.close(fig2)


# # plot the lower two levels next to each other
# fig3, ax3 = plt.subplots(nrows=3, ncols=1)
# ax3[0].plot(times, usolano)
# ax3[1].plot(times, vsolano)
# ax3[2].plot(times, speed)
# ax3[0].set_ylabel('u')
# ax3[1].set_ylabel('v')
# ax3[2].set_ylabel('speed')
# ax3[0].set_title(run_name+" "+domain)
# fig3.savefig(os.path.join(plot_dir2, "wind_tseries_"+run_name+"_"+domain+".png"))
# plt.close(fig3)



