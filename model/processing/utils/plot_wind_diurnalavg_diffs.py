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
burn_in = 0  # number of half-hour output timesteps that constitute burn-in period
#coords_solano = {'d01': (49, 51), 'd02': (91, 96), 'd03': (127, 143)}

def get_box_index(f, lat, lon):
    dif_lat = f.variables['XLAT'][:, :]-lat
    dif_lon = f.variables['XLONG'][:, :]-lon
    distance = dif_lat**2 + dif_lon**2
    ixlat, ixlon = np.unravel_index(np.argmin(distance), distance.shape)
    return ixlat, ixlon

def get_winds(f, coords_solano):
    mask_levels = (f.variables['levels'][:] == 110.) | (f.variables['levels'][:] == 150.)
    usolano = f.variables['U'][:, mask_levels, coords_solano[1], coords_solano[0]]
    vsolano = f.variables['V'][:, mask_levels, coords_solano[1], coords_solano[0]]
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

def remove_daily_mean(times, speed):
    monthmin = times[0].month
    daymin = times[0].day
    monthmax = times[-1].month
    daymax = times[-1].day
    daysdelta = (datetime.datetime(2009, monthmax, daymax)-datetime.datetime(2009, monthmin, daymin)).days
    sp_rem_dmean = np.zeros(len(speed))+np.nan
    for i in xrange(daysdelta):
        mask = (times >= times[0]+datetime.timedelta(days=i)) & (times < times[0]+datetime.timedelta(days=i+1))
        davg = np.mean(speed[mask])
        sp_rem_dmean[mask] = speed[mask]-davg
    return sp_rem_dmean

def avg_diurnal(times, speed, flag_abs=True):

    hours = np.array([d.hour+float(d.minute)/60. for d in times])
    hrU = list(set(hours))
    hrU.sort()
    sp_av = np.zeros(len(hrU))
    sp_sd = np.zeros(len(hrU))

    if flag_abs:
        speed = remove_daily_mean(times, speed)

    for i, hr in enumerate(hrU):
        mask = hours == hr
        sp_av[i] = np.mean(speed[mask])
        sp_sd[i] = np.std(speed[mask])

    return hrU, sp_av, sp_sd

def plot_wind(ax, level, case):

    # ax has three panels: ax[0] plot sp_av, ax[1] plot diff sp_av, ax[2] plot diff sp_sd
    
    if case["test_id"] in ["CV_0.2", "CV_0.1"]:
        colors = ['m', 'r', 'y', 'g', 'c', 'b']
    else:
        colors = ['b', 'g', 'r', 'c', 'm', 'y']

    # open control file
    fc = nc.netcdf_file(os.path.join(root_dir, case["ctrl"], "wrfout_interp_lower_"+domain+".nc"))

    # find solano coords
    print "get coords"
    ixlat, ixlon = get_box_index(fc, lat_solano, lon_solano)
    coords_solano = (ixlon, ixlat)
    
    # get wind speed at solano and times in PST
    print "get control wind"
    speed_c = get_winds(fc, coords_solano)
    times = convert_times(fc)
    hr_c, sp_av_c, sp_sd_c = avg_diurnal(times, speed_c[:, level], flag_abs=True)
        
    # plot control
    #h = ax[0].plot(hr_c, sp_av_c, 'k', lw=1.5)
    handles = []  #[h[0]]
    labels = []  #[case["ctrl"]]

    # plot each test
    for i_test, test in enumerate(case["tests"]):
        print test
        # open test file
        ft = nc.netcdf_file(os.path.join(root_dir, test, "wrfout_interp_lower_"+domain+".nc"))
        # get solano wind
        speed_t = get_winds(ft, coords_solano)
        # average over diurnal cycle
        # hr_t, sp_av_t, sp_sd_t = avg_diurnal(times, speed_t[:, level], flag_abs=True)
        hr_d, sp_av_d, sp_sd_d = avg_diurnal(times, speed_t[:, level]-speed_c[:, level], flag_abs=False)
        # plot
        # h = ax[0].plot(hr_t, sp_av_t, colors[i_test])
        # ax[0].fill_between(hr_t, sp_av_t-sp_sd_t, sp_av_t+sp_sd_t, color=colors[i_test], alpha=0.4)
        # handles += h#[0]
        # labels += [test]
        # plot diff
        ax.fill_between(hr_d, sp_av_d-sp_sd_d, sp_av_d+sp_sd_d, color=colors[i_test], alpha=0.2)
        h = ax.plot(hr_d, sp_av_d, colors[i_test], lw=2)
        handles += h#[0]
        labels += [test]
        # plot std dev
        # ax[2].plot(hr_d, sp_sd_d, colors[i_test])

    ax.legend(handles, labels, prop={'size': 8}, loc='lower right', ncol=2)
    # ax[0].legend(handles, labels, prop={'size': 8}, loc='lower right', ncol=2)
    #ax1.legend(handles, labels, prop={'size': 8}, loc='lower center', ncol=2)

    ax.set_xticks([6, 12, 18])
    ax.set_xlim(0, 24)
    ax.grid(b=True, axis='x')

    ax.plot([0, 24], [0, 0], "g:")
    # ax[1].plot(ax[1].get_xlim(), [0, 0], "g:")

    ax.set_xlabel('PST hour')
    ax.set_ylabel("avg $\\Delta$ speed, m/s")
    
    #ax0.set_ylim(2, 14)
    #ax1.set_ylim(2, 16)
    # ax[0].set_xlabel('PST hour')
    # ax[1].set_xlabel('PST hour')
    # ax[0].set_ylabel("avg speed, m/s")
    # ax[1].set_ylabel("avg $\\Delta$ speed, m/s")
    # ax[2].set_ylabel("std dev $\\Delta$ speed, m/s")
    #ax[0].xaxis.grid(True)
    #ax[1].xaxis.grid(True)

    #rect0 = mpatches.Rectangle([date2num(times[0]), 2], date2num(times[burn_in])-date2num(times[0]), 12, ec="none", fc="gray", alpha=0.5)
    #ax0.add_patch(rect0)
    #rect1 = mpatches.Rectangle([date2num(times[0]), 2], date2num(times[burn_in])-date2num(times[0]), 14, ec="none", fc="gray", alpha=0.5)
    #ax1.add_patch(rect1)

    return handles, labels

def plot_control_wind(ax, lev, case, subtract_mean=False):

    # open control file
    fc = nc.netcdf_file(os.path.join(root_dir, case["ctrl"], "wrfout_interp_lower_"+domain+".nc"))

    # find solano coords
    print "get coords"
    ixlat, ixlon = get_box_index(fc, lat_solano, lon_solano)
    coords_solano = (ixlon, ixlat)
    
    # get wind speed at solano and times in PST
    print "get control wind"
    speed_c = get_winds(fc, coords_solano)
    times = convert_times(fc)

    monthmin = times[0].month
    daymin = times[0].day
    monthmax = times[-1].month
    daymax = times[-1].day
    daysdelta = (datetime.datetime(2009, monthmax, daymax)-datetime.datetime(2009, monthmin, daymin)).days

    if subtract_mean:
        speed_plot = remove_daily_mean(times, speed_c[:, lev])
    else:
        speed_plot = speed_c[:, lev]
    
    for i in xrange(daysdelta):
        mask = (times >= datetime.datetime(2009, monthmin, daymin)+datetime.timedelta(days=i)) & \
               (times < datetime.datetime(2009, monthmin, daymin)+datetime.timedelta(days=i+1))
        hours = times[mask]-(datetime.datetime(2009, monthmin, daymin)+datetime.timedelta(days=i))
        hours = np.array([h.seconds/3600. for h in hours])
        ax.plot(hours, speed_plot[mask], 'k', alpha=0.5)
        
    ax.set_xticks([6, 12, 18])
    ax.set_xlim(0, 24)
    ax.grid(b=True, axis='x')

    if subtract_mean:
        ax.set_ylabel("Wind speed minus daily mean, m/s")
    else:
        ax.set_ylabel("Wind speed, m/s")

    ax.set_xlim(0, 24)
    ax.set_xlabel("PST hour")

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
    case_CVrgs_dry = {"ctrl": "CA-0.2", \
                  "tests": ["CA-CV0.05w", "CA-CV0.05a", "CA-CV0.05d", "CA-CV0.05n", "CA-CV0.05s"], \
                  "smois_c": 0.2, \
                  "test_id": "CV_rgs_dry"}
    case_CVrgs_wet = {"ctrl": "CA-0.2", \
                  "tests": ["CA-CV0.35w", "CA-CV0.35a", "CA-CV0.35d", "CA-CV0.35n", "CA-CV0.35s"], \
                  "smois_c": 0.2, \
                  "test_id": "CV_rgs_wet"}

    coords_solano = {}

    case_list = [case_wetrg, case_dryrg, case_CVwet]

    plot_dir = os.path.join(root_dir, "timeseries_wind")
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    for domain in domains:

        print domain

        for lev in xrange(2):

            print lev

            for case in case_list:

                print case["test_id"]

                # fig, ax = plt.subplots(nrows=3, ncols=1)
                fig, ax = plt.subplots()

                # plot: 
                handles, labels = plot_wind(ax, lev, case)
                ax.set_title(domain+" "+case["test_id"]+", Solano wind speed")
                # ax[0].set_title(domain+" "+case["test_id"]+", Solano wind speed")

                # fig.set_size_inches(15, 7)
                # fig.set_size_inches(10, 8)
                # fig.subplots_adjust(hspace=0.35)
                fig.savefig(os.path.join(plot_dir, "solano_diurnalwind_"+case["test_id"]+"_"+domain+"_level"+str(lev)+".png"))

                # plt.show()
                # 1/0
                plt.close(fig)

                # plot control wind
                fig, ax = plt.subplots()
                plot_control_wind(ax, lev, case)
                ax.set_title(case["ctrl"]+" "+domain)
                fig.savefig(os.path.join(plot_dir, "solano_controlwind_"+case["ctrl"]+"_"+domain+"_level"+str(lev)+".png"))

                # plt.show()
                # 1/0
                plt.close(fig)

                # plot control wind
                fig, ax = plt.subplots()
                plot_control_wind(ax, lev, case, subtract_mean=True)
                ax.set_title(case["ctrl"]+" "+domain)
                fig.savefig(os.path.join(plot_dir, "solano_controlwind_minusmean_"+case["ctrl"]+"_"+domain+"_level"+str(lev)+".png"))

                # plt.show()
                # 1/0
                plt.close(fig)
 
 
