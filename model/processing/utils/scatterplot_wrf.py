"""
PAL 8/26/2014
Scatter plots of WRF output
"""

import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.io import netcdf as nc


if __name__ == "__main__":

    ctrl_name = "CA-0.24"
    #test_names = ["CA-dryCR", "CA-dryCV", "CA-drySN"]
    #domains = ["d01", "d02", "d03"]
    test_names = ["CA-dryCR"]
    domains = ["d01"]
    sm_test = 0.05  # SMOIS value in the test region
    utc_offset = 8  # offset to subtract to get local time
    root_dir = "/scratch2/scratchdirs/plink/WRF/output"
    mintime = 42
    maxtime = -1
    boxes = ['nCV', 'cCV', 'sCV', 'nCR', 'sCR', 'SN']

    for test in test_names:

        for domain in domains:

            # make plot output directory if it doesn't exist
            plot_dir = os.path.join(root_dir, test, "plots", domain, "scatter_sfc")
            if not os.path.exists(plot_dir):
                os.mkdir(plot_dir)

            # open the two netcdf files
            ft = nc.netcdf_file(os.path.join(root_dir, test, "wrfout_"+domain+"_2009-07-01_00:00:00"), "r")
            fi = nc.netcdf_file(os.path.join(root_dir, test, "wrfout_interp_"+domain+".nc"), "r")
            #fc = nc.netcdf_file(os.path.join(root_dir, ctrl_name, "wrfout_"+domain+"_2009-07-01_00:00:00"), "r")

            # get the arrays of HFX data and T0 (T at the first model level)
            hfx_t = ft.variables["HFX"][mintime:maxtime, :, :]
            hfx_c = fc.variables["HFX"][mintime:maxtime, :, :]
            t0_t = (ft.variables["T"][mintime:maxtime, 0, :, :].T + ft.variables["T00"][mintime:maxtime]).T - 273.  # convert to C
            t0_c = (fc.variables["T"][mintime:maxtime, 0, :, :].T + fc.variables["T00"][mintime:maxtime]).T - 273.
            xlon = ft.variables["XLONG"][mintime:maxtime, :, :]
            xlat = ft.variables["XLAT"][mintime:maxtime, :, :]
            hgt = ft.variables["HGT"][mintime:maxtime, :, :]

            # get hour of day array from times array, and expand to size of HFX array
            hours = []
            times_arr = ft.variables['Times'][mintime:maxtime]
            for ii in xrange(times_arr.shape[0]):
                time_arr = times_arr[ii, :]
                hour = int(''.join(time_arr[11:13]))
                hour = hour-utc_offset
                if hour < 0:
                    hour = hour+24
                hours.append(hour)
            hours = np.array(hours)

            hours_arr = np.tile(np.array([[hours]]).T, (1, hfx_t.shape[1], hfx_t.shape[2]))
            print "hfx shape:", hfx_t.shape
            print "hours_arr shape:", hours_arr.shape

            # get landmask
            mask_land = ft.variables["LANDMASK"][0, :, :] == 1
            mask_land = np.tile(np.array([mask_land]).T, (1, hfx_t.shape[0])).T
            print "mask_land shape:", mask_land.shape

            # get mask of test region, using inital soil moisture
            mask_test = ft.variables["SMOIS"][0, 0, :, :] == sm_test
            mask_test = np.tile(np.array([mask_test]).T, (1, hfx_t.shape[0])).T
            print "mask_test shape:", mask_test.shape

            # flatten the HFX and T0 arrays and the hour of day array
            hfx_t = np.ravel(hfx_t)
            hfx_c = np.ravel(hfx_c)
            t0_t = np.ravel(t0_t)
            t0_c = np.ravel(t0_c)
            hours_arr = np.ravel(hours_arr)
            mask_test = np.ravel(mask_test)
            mask_land = np.ravel(mask_land)
            xlon = np.ravel(xlon)
            xlat = np.ravel(xlat)
            hgt = np.ravel(hgt)

            #xlatmasked = xlat[mask_land & mask_test]
            #xlonmasked = xlon[mask_land & mask_test]
            #maskbox = (xlat == xlatmasked[400]) & (xlon == xlonmasked[400])
            maskbox = (xlat >= 37.8) & (xlat <= 38.1)

            # scatter plot HFXctrl vs HFXtest, and color by hour of day
            f, ax = plt.subplots()
            h = ax.scatter(hfx_c[mask_land & mask_test & maskbox], hfx_t[mask_land & mask_test & maskbox], c=hours_arr[mask_land & mask_test & maskbox])
            #hours_arr[mask_land & mask_test])
            ax.set_xlabel("HFX ctrl")
            ax.set_ylabel("HFX test")
            f.colorbar(h)

            plt.show()
            #1/0

            # calculate dHFX and dT0 (test - control)
            dhfx = hfx_t-hfx_c
            dt0 = t0_t-t0_c

            # scatter plot dHFX vs dT0, and color by hour of day
            f, ax = plt.subplots()
            h = ax.scatter(dhfx[mask_land & mask_test & maskbox], dt0[mask_land & mask_test & maskbox], c=hours_arr[mask_land & mask_test & maskbox])
            #hours_arr[mask_land & mask_test])
            ax.set_xlabel("dHFX, test-ctrl")
            ax.set_ylabel("dT0, test-ctrl")
            f.colorbar(h)

            plt.show()

            # 