"""
PAL 9/10/2014
quantify changes in max, min winds and in timing of incr and decr
"""

from scipy.io import netcdf as nc
import matplotlib.pyplot as plt
import numpy as np
import os
import datetime


# global settings
utc_offset = 8  # hours
ix_time_start = 48  # burn-in period: 24 hrs, 30 min timesteps
min_start_hr = 10
min_end_hr = 12
max_start_hr = 20
max_end_hr = 2
lat_solano, lon_solano = 38.166, -121.817
level = 0  # which wind level to use (0=60mAGL, 1=100mAGL)


def get_box_index(f, lat, lon):
    dif_lat = f.variables['XLAT'][:, :]-lat
    dif_lon = f.variables['XLONG'][:, :]-lon
    distance = dif_lat**2 + dif_lon**2
    ixlat, ixlon = np.unravel_index(np.argmin(distance), distance.shape)
    return ixlat, ixlon

def get_winds(f, coords_solano):
    # coords_solano are (lon, lat)
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
    return np.array(times)-datetime.timedelta(hours=utc_offset)

def smooth(timeseries, window=2):
    t2 = np.zeros_like(timeseries)
    for i in xrange(len(timeseries)):
        i_start = max(0, i-window)
        i_end = min(len(timeseries), i+window)
        tmp = timeseries[i_start:i_end]
        mask = np.isfinite(tmp)
        if len(tmp[mask])>0:
            t2[i] = np.mean(tmp[mask])
        else:
            t2[i] = np.nan
    return t2

def calc_max_min(wind, dt):
    # burn-in should already be removed
    # note: assumes all days in the time period are from same year and month

    year = dt[-1].year
    month = dt[-1].month
    day_arr = np.array([d.day for d in dt])
    hr_arr = np.array([d.hour for d in dt])

    mask_finite = np.isfinite(wind)

    days = np.arange(np.min(day_arr), np.max(day_arr)+1)
    mins = np.zeros_like(days)+np.nan
    maxes = np.zeros_like(days)+np.nan
    hour_rampup = np.zeros_like(days)+np.nan
    hour_rampdn = np.zeros_like(days)+np.nan
    for i, d in enumerate(days):
        mask_min = (dt >= datetime.datetime(year, month, d, min_start_hr)) & \
                   (dt <= datetime.datetime(year, month, d, min_end_hr))
        mask_max = (dt >= datetime.datetime(year, month, d, max_start_hr)) & \
                   (dt <= datetime.datetime(year, month, d+1, max_end_hr))
        if len(wind[mask_min & mask_finite]) > 0:
            mins[i] = np.mean(wind[mask_min & mask_finite])
        if len(wind[mask_max & mask_finite]) > 0:
            maxes[i] = np.mean(wind[mask_max & mask_finite])

        # check if this day's min and max could be calculated
        #if dt[0] <= datetime.datetime(year, month, d, min_end_hr):
        if np.isfinite(mins[i]) and np.isfinite(maxes[i]):
            mask_rampup = (dt >= datetime.datetime(year, month, d, min_end_hr)) & \
                          (dt <= datetime.datetime(year, month, d, max_start_hr))
            tmp_smooth = smooth(wind[mask_rampup])
            tmp_dt = dt[mask_rampup]
            ix = np.where(tmp_smooth >= mins[i]+(maxes[i]-mins[i])/2.)
            hour_rampup[i] = tmp_dt[ix[0][0]].hour+tmp_dt[ix[0][0]].minute/60.

        # check if this day's min and previous day's max could be calculated
        #if dt[-1] >= datetime.datetime(year, month, d, min_start_hr):
        if i > 0:
            if np.isfinite(mins[i]) and np.isfinite(maxes[i-1]):
                mask_rampdn = (dt >= datetime.datetime(year, month, d, max_end_hr)) & \
                              (dt <= datetime.datetime(year, month, d, min_start_hr))
                tmp_smooth = smooth(wind[mask_rampdn])
                tmp_dt = dt[mask_rampdn]
                ix = np.where(tmp_smooth <= mins[i]+(maxes[i-1]-mins[i])/2.)
                hour_rampdn[i] = tmp_dt[ix[0][0]].hour+tmp_dt[ix[0][0]].minute/60.

    return days, mins, maxes, hour_rampup, hour_rampdn

def calc_mean_std_diffs(test_vals, ctrl_vals):
    diffs = test_vals-ctrl_vals
    diffs = diffs[np.isfinite(diffs)]
    if len(diffs)>2:
        mu = np.mean(diffs)
        sigma = np.std(diffs)
    else:
        mu = np.nan
        sigma = np.nan
    return mu, sigma

def calc_diff_metrics(d):
    # d should be a dict with these fields: max_c, max_t, min_c, min_v, ru_c, ru_t, rd_c, rd_t
    r = {}  # dictionary for return values
    r['mu_max'], r['sig_max'] = calc_mean_std_diffs(d['max_t'], d['max_c'])
    r['mu_min'], r['sig_min'] = calc_mean_std_diffs(d['min_t'], d['min_c'])
    r['mu_ru'], r['sig_ru'] = calc_mean_std_diffs(d['ru_t'], d['ru_c'])
    r['mu_rd'], r['sig_rd'] = calc_mean_std_diffs(d['rd_t'], d['rd_c'])
    return r

def plot_diffs(fig, ax, case):
    
    if case["test_id"] in ["CV_0.1", "CV_0.2"]:
        x = case["smois"]
    else:
        x = np.arange(len(case["tests"]))
    if case["test_id"] in ["CV_0.1", "wet_regions"]:
        color = 'r'
        symbol = 's'
    elif case["test_id"] in ["CV_0.2", "dry_regions"]:
        color = 'b'
        symbol = 'D'

    # open control file
    fc = nc.netcdf_file(os.path.join(root_dir, case["ctrl"], "wrfout_interp_"+domain+".nc"))

    # find solano coords
    ixlat, ixlon = get_box_index(fc, lat_solano, lon_solano)
    coords_solano = (ixlon, ixlat)
    
    # get wind speed at solano and times in PST
    speed_c = get_winds(fc, coords_solano)
    times = convert_times(fc)

    # get metrics for control
    m_c = {}
    m_c['days_c'], m_c['min_c'], m_c['max_c'], m_c['ru_c'], m_c['rd_c'] = \
            calc_max_min(speed_c[:, level], times)

    # make lists for saving diffs for each test and each metric
    dif_mu_max = []
    dif_mu_min = []
    dif_mu_ru = []
    dif_mu_rd = []
    dif_sd_max = []
    dif_sd_min = []
    dif_sd_ru = []
    dif_sd_rd = []

    # get metrics for each test
    for test in case["tests"]:
        print test
        m_t = m_c.copy()  # copy control data into a temp dict
        # open test file
        ft = nc.netcdf_file(os.path.join(root_dir, test, "wrfout_interp_"+domain+".nc"))
        # get solano wind
        speed_t = get_winds(ft, coords_solano)
        # get metrics for each day for the test
        m_t['days_t'], m_t['min_t'], m_t['max_t'], m_t['ru_t'], m_t['rd_t'] = \
                calc_max_min(speed_t[:, level], times)
        # calculate mean and std dev of the difference in the metrics, and save
        tmp = calc_diff_metrics(m_t)
        dif_mu_max.append(tmp["mu_max"])
        dif_mu_min.append(tmp["mu_min"])
        dif_mu_ru.append(tmp["mu_ru"])
        dif_mu_rd.append(tmp["mu_rd"])
        dif_sd_max.append(tmp["sig_max"])
        dif_sd_min.append(tmp["sig_min"])
        dif_sd_ru.append(tmp["sig_ru"])
        dif_sd_rd.append(tmp["sig_rd"])

    h = ax[0, 0].errorbar(x, dif_mu_max, yerr=dif_sd_max, color=color, marker=symbol, linestyle='None')
    ax[0, 1].errorbar(x, dif_mu_min, yerr=dif_sd_min, color=color, marker=symbol, linestyle='None')
    ax[1, 0].errorbar(x, dif_mu_ru, yerr=dif_sd_ru, color=color, marker=symbol, linestyle='None')
    ax[1, 1].errorbar(x, dif_mu_rd, yerr=dif_sd_rd, color=color, marker=symbol, linestyle='None')

    if case["test_id"] in ["CV_0.1", "CV_0.2"]:
        for a in ax.flatten():
            a.plot(case["smois_c"], 0, 'o', mec=color, mfc='none', mew=1.5)
            a.set_xlim(min(case["smois"])-0.05, max(case["smois"])+0.05)
            a.set_xticks(np.arange(min(case["smois"]), max(case["smois"])+0.05, 0.05))
            a.set_xlabel("CV soil moisture")
    else:
        for a in ax.flatten():
            a.set_xlim(-0.5, 2.5)
            a.set_xticks(np.arange(3))
            a.set_xticklabels(["CR", "CV", "SN"])

    for a in ax.flatten():
        lims = a.get_xlim()
        a.plot(lims, [0, 0], ':', color='gray')

    ax[0, 0].set_ylabel("$\\Delta$ av. max. (m/s)")
    ax[0, 1].set_ylabel("$\\Delta$ av. min. (m/s)")
    ax[1, 0].set_ylabel("$\\Delta$ ramp up (hr)")
    ax[1, 1].set_ylabel("$\\Delta$ ramp down (hr)")

    ax[0, 0].set_ylim(-2, 2)
    ax[0, 1].set_ylim(-2, 2)
    ax[1, 0].set_ylim(-1, 2)
    ax[1, 1].set_ylim(-1, 2)

    return h


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
    

    plot_dir = os.path.join(root_dir, "shifts_wind")
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    for domain in domains:

        print domain
        fig1, ax1 = plt.subplots(nrows=2, ncols=2)
        fig2, ax2 = plt.subplots(nrows=2, ncols=2)

        # wet_regions and dry_regions cases
        hdrybkg = plot_diffs(fig1, ax1, case_wetrg)
        hwetbkg = plot_diffs(fig1, ax1, case_dryrg)
        ax1[0, 0].legend([hdrybkg, hwetbkg], ["dry background", "wet background"], numpoints=1, prop={"size": 8})
        plt.tight_layout()
        fig1.savefig(os.path.join(plot_dir, "shifts_regions_"+domain+".png"))

        # CV_0.2 and CV_0.1 cases
        hwet = plot_diffs(fig2, ax2, case_CVwet)
        hdry = plot_diffs(fig2, ax2, case_CVdry)
        hfakectrl = plt.Line2D([0], [0], marker='o', mec='k', mfc='none', linestyle='')
        ax2[0, 0].legend([hwet, hdry, hfakectrl], ["mtn smois 0.2", "mtn smois 0.1", "control case"], numpoints=1, prop={"size": 8})
        plt.tight_layout()
        fig2.savefig(os.path.join(plot_dir, "shifts_CVsmois_"+domain+".png"))
       
        #plt.show()
        #1/0
