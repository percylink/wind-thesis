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


# global variables
lat_solano, lon_solano = 38.166, -121.817
utc_offset = 8.  # hrs to PST
burn_in = 48  # number of half-hour output timesteps that constitute burn-in period
box_pts = {"bay": {"d01": [(42, 48), (44, 52)], "d02": [(70, 88), (76, 100)]}, 
           "cv": {"d01": [(48, 50), (53, 54)], "d02": [(89, 94), (101, 104)]}}

root_dir = "/scratch2/scratchdirs/plink/WRF/output"
plot_dir = os.path.join(root_dir, "timeseries_forcing")
if not os.path.exists(plot_dir):
    os.mkdir(plot_dir)

# get masks for regions
maskCR = {}
f = nc.netcdf_file(os.path.join(root_dir, "CA-dryCR", "wrfout_d01_2009-07-01_00:00:00"))
maskCR["d01"] = f.variables["SMOIS"][0, 0, :, :] == 0.1
f = nc.netcdf_file(os.path.join(root_dir, "CA-dryCR", "wrfout_d02_2009-07-01_00:00:00"))
maskCR["d02"] = f.variables["SMOIS"][0, 0, :, :] == 0.1
maskCV = {}
f = nc.netcdf_file(os.path.join(root_dir, "CA-dryCV", "wrfout_d01_2009-07-01_00:00:00"))
maskCV["d01"] = f.variables["SMOIS"][0, 0, :, :] == 0.1
f = nc.netcdf_file(os.path.join(root_dir, "CA-dryCV", "wrfout_d02_2009-07-01_00:00:00"))
maskCV["d02"] = f.variables["SMOIS"][0, 0, :, :] == 0.1
maskSN = {}
f = nc.netcdf_file(os.path.join(root_dir, "CA-drySN", "wrfout_d01_2009-07-01_00:00:00"))
maskSN["d01"] = f.variables["SMOIS"][0, 0, :, :] == 0.1
f = nc.netcdf_file(os.path.join(root_dir, "CA-drySN", "wrfout_d02_2009-07-01_00:00:00"))
maskSN["d02"] = f.variables["SMOIS"][0, 0, :, :] == 0.1


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

def get_winds(f, lev):
    u = f.variables['U'][:, lev, :, :]
    v = f.variables['V'][:, lev, :, :]
    speed = (u**2 + v**2)**0.5  # calculate speed
    return speed

def calc_mask_avg(f, data, mask):
    # data is shape time x lat x lon
    # mask is dim lat x lon
    tdim = data.shape[0]
    av = np.zeros(tdim) + np.nan
    for i in xrange(tdim):
        tmp = data[i, :, :]
        tmp = tmp[mask]
        av[i] = np.mean(tmp[np.isfinite(tmp)])
    return av

def calc_box_av(f, data, box_name, domain):
    # data is shape time x lat x lon
    # corner indices are lon, lat
    llcnr, urcnr = box_pts[box_name][domain]
    dslice = data[:, llcnr[1]:urcnr[1], llcnr[0]:urcnr[0]]
    rtrn = np.zeros(data.shape[0])+np.nan
    for i in xrange(len(rtrn)):
        tmp = dslice[i, :, :]
        rtrn[i] = np.mean(tmp[np.isfinite(tmp)])
    return rtrn

def plot_tseries(case):
    
    if case["test_id"] in ["CV_0.2", "CV_0.1"]:
        colors = ['m', 'r', 'y', 'g', 'c', 'b']
    else:
        colors = ['b', 'g', 'r', 'c', 'm', 'y']

    variables = ['SMOIS', 'TSK', 'T2', 'speed_upper']

    for domain in ["d01", "d02"]:

        # open control file
        fc = nc.netcdf_file(os.path.join(root_dir, case["ctrl"], "wrfout_"+domain+"_2009-07-01_00:00:00"))
        times = convert_times(fc)

        # open test files
        fts = {}
        for test in case["tests"]:
            fts[test] = nc.netcdf_file(os.path.join(root_dir, test, "wrfout_"+domain+"_2009-07-01_00:00:00"))

        for v in variables:

            # set up plot
            fig, ax = plt.subplots(nrows=5, ncols=2)
            fig.suptitle(case["test_id"]+" "+domain+" "+v)
            savename = "forcings_"+case["test_id"]+"_"+v+"_"+domain+".png"

            # get data for control
            if v == 'SMOIS':
                data_c = fc.variables[v][:, 0, :, :]
            elif v == 'speed_upper':
                data_c = get_winds(fc, 15)  # level 15 is a little more than 4 km in central valley; ~4-5 km over cr, 5-6 over sn
            else:
                data_c = fc.variables[v][:]

            # get averages for the regions for control
            avCR_c = calc_mask_avg(fc, data_c, maskCR[domain])
            avCV_c = calc_mask_avg(fc, data_c, maskCV[domain])
            avSN_c = calc_mask_avg(fc, data_c, maskSN[domain])
            avbay_c = calc_box_av(fc, data_c, "bay", domain)
            avsac_c = calc_box_av(fc, data_c, "cv", domain)

            # plot control
            h = ax[0, 0].plot(times, avCR_c, 'k', lw=1.5)
            ax[0, 0].set_ylabel("CR")
            handles = [h[0]]
            labels = [case["ctrl"]]
            ax[1, 0].plot(times, avCV_c, 'k', lw=1.5)
            ax[1, 0].set_ylabel("CV")
            ax[2, 0].plot(times, avSN_c, 'k', lw=1.5)
            ax[2, 0].set_ylabel("SN")
            ax[3, 0].plot(times, avbay_c, 'k', lw=1.5)
            ax[3, 0].set_ylabel("bay box")
            ax[4, 0].plot(times, avsac_c, 'k', lw=1.5)
            ax[4, 0].set_ylabel("sac box")
            
            # plot each test
            for i_test, test in enumerate(case["tests"]):

                print test

                # get test file
                ft = fts[test]

                # get data
                if v == 'SMOIS':
                    data_t = ft.variables[v][:, 0, :, :]
                elif v == 'speed_upper':
                    data_t = get_winds(ft, 15)  # level 15 is a little more than 4 km in central valley; ~4-5 km over cr, 5-6 over sn
                else:
                    data_t = ft.variables[v][:]

                # get averages for the regions for control
                avCR_t = calc_mask_avg(ft, data_t, maskCR[domain])
                avCV_t = calc_mask_avg(ft, data_t, maskCV[domain])
                avSN_t = calc_mask_avg(ft, data_t, maskSN[domain])
                avbay_t = calc_box_av(ft, data_t, "bay", domain)
                avsac_t = calc_box_av(ft, data_t, "cv", domain)

                # plot control
                h = ax[0, 0].plot(times, avCR_t, colors[i_test])
                handles += h
                labels.append(test)
                ax[1, 0].plot(times, avCV_t, colors[i_test])
                ax[2, 0].plot(times, avSN_t, colors[i_test])
                ax[3, 0].plot(times, avbay_t, colors[i_test])
                ax[4, 0].plot(times, avsac_t, colors[i_test])

                ax[0, 1].plot(times, avCR_t-avCR_c, colors[i_test])
                ax[1, 1].plot(times, avCV_t-avCV_c, colors[i_test])
                ax[2, 1].plot(times, avSN_t-avSN_c, colors[i_test])
                ax[3, 1].plot(times, avbay_t-avbay_c, colors[i_test])
                ax[4, 1].plot(times, avsac_t-avsac_c, colors[i_test])
                
            ax[0, 0].legend(handles, labels, prop={'size': 6}, loc='lower center', ncol=2)
            ax[0, 0].set_title('abs')
            ax[0, 1].set_title('test - ctrl')
            
            for a in ax.flatten():
                a.set_xlabel('PST')

            fig.set_size_inches(10, 11)
            plt.show()
            1/0
            fig.savefig(os.path.join(plot_dir, savename))
            plt.close(fig)


if __name__ == "__main__":
    
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

    # call plotting function
    plot_tseries(case_wetrg)
    plot_tseries(case_dryrg)
    plot_tseries(case_CVwet)
    plot_tseries(case_CVdry)
    plot_tseries(case_CVrgs_wet)
    plot_tseries(case_CVrgs_dry)



