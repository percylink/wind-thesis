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
    test_names = ["CA-dryCR", "CA-dryCV", "CA-drySN"]
    domains = ["d01", "d02"]#, "d03"]
    #test_names = ["CA-dryCR"]
    #domains = ["d01"]
    utc_offset = 8  # offset to subtract to get local time
    root_dir = "/scratch2/scratchdirs/plink/WRF/output"
    mintime = 42
    maxtime = -1
    boxes = ['nCV', 'cCV', 'sCV', 'nCR', 'sCR'] #, 'SN']

    for test in [ctrl_name]+test_names:

        print test

        for domain in domains:

            print domain

            # make plot output directory if it doesn't exist
            plot_dir = os.path.join(root_dir, test, "plots", domain, "scatter_sfc")
            if not os.path.exists(plot_dir):
                os.mkdir(plot_dir)

            # open the two netcdf files
            ft = nc.netcdf_file(os.path.join(root_dir, test, "wrfout_"+domain+"_2009-07-01_00:00:00"), "r")
            fi = nc.netcdf_file(os.path.join(root_dir, test, "wrfout_interp_"+domain+".nc"), "r")
            #fc = nc.netcdf_file(os.path.join(root_dir, ctrl_name, "wrfout_"+domain+"_2009-07-01_00:00:00"), "r")

            xlon = ft.variables["XLONG"][0, :, :]
            xlat = ft.variables["XLAT"][0, :, :]
            hgt = ft.variables["HGT"][0, :, :]

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

            hours_arr = np.tile(np.array([[hours]]).T, (1, hgt.shape[0], hgt.shape[1]))
            #print "hours_arr shape:", hours_arr.shape

            # get landmask
            mask_land = ft.variables["LANDMASK"][0, :, :] == 1

            for box in boxes:

                print box

                if box == 'nCV':
                    mask_box = (xlat>39) & (xlon>-122.5) & (hgt<200)
                elif box == 'cCV':
                    mask_box = (xlat>37.7) & (xlat<38.5) & (xlon>-122) & (hgt<200)
                elif box == 'sCV':
                    mask_box = (xlat<37) & (xlon>-120.9) & (hgt<200)
                elif box == 'nCR':
                    mask_box = (xlat>38) & (xlon<-122) & (hgt>200)
                elif box == 'sCR':
                    mask_box = (xlat<37) & (xlon<-120) & (hgt>200)

                hours_box = hours_arr[:, mask_box & mask_land]
                #print "hours_box shape", hours_box.shape
                hfx = ft.variables["HFX"][mintime:maxtime, mask_box & mask_land]
                #print "hfx shape", hfx.shape
                ts = ft.variables["TSK"][mintime:maxtime, mask_box & mask_land]-273.
                #print "ts shape", ts.shape
                t2 = ft.variables["T2"][mintime:maxtime, mask_box & mask_land]-273.
                #print "t2 shape", t2.shape
                tair = {}
                for i_lev, lev in enumerate(fi.variables["levels"][:]):
                    tmp = fi.variables["T"][mintime:maxtime, i_lev, :, :].copy()
                    tmp[tmp>1e36] = np.nan
                    tair[lev] = tmp[:, mask_box & mask_land]
                    #print "tair shape", lev, tair[lev].shape
                    if lev == 250:
                        tmpp = fi.variables["P"][mintime:maxtime, i_lev, :, :].copy()
                        tmpp[tmpp>1e36] = np.nan
                        p250 = tmpp[:, mask_box & mask_land]
                        #print "p250 shape", p250.shape

                # scatter plot HFX vs Ts
                fig, ax = plt.subplots()
                h = ax.scatter(np.ravel(hfx), np.ravel(ts), c=np.ravel(hours_box), edgecolor=None)
                fig.colorbar(h)
                ax.set_xlabel('HFX')
                ax.set_ylabel('Ts')
                ax.set_title(box)
                fig.savefig(os.path.join(plot_dir, "HFX_vs_Ts_"+box+".png"))

                # scatter plot HFX vs T2m
                fig, ax = plt.subplots()
                h = ax.scatter(np.ravel(hfx), np.ravel(t2), c=np.ravel(hours_box), edgecolor=None)
                fig.colorbar(h)
                ax.set_xlabel('HFX')
                ax.set_ylabel('T 2m')
                ax.set_title(box)
                fig.savefig(os.path.join(plot_dir, "HFX_vs_T2_"+box+".png"))

                # scatter plot Ts vs T2m
                fig, ax = plt.subplots()
                h = ax.scatter(np.ravel(ts), np.ravel(t2), c=np.ravel(hours_box), edgecolor=None)
                fig.colorbar(h)
                ax.set_xlabel('Ts')
                ax.set_ylabel('T 2m')
                ax.set_title(box)
                fig.savefig(os.path.join(plot_dir, "Ts_vs_T2_"+box+".png"))

                for lev in fi.variables["levels"][:]:

                    print lev
                    
                    if len(tair[lev][np.isfinite(tair[lev])]):

                        # scatter T2m vs Tair for each level
                        fig, ax = plt.subplots()
                        h = ax.scatter(np.ravel(t2), np.ravel(tair[lev]), c=np.ravel(hours_box), edgecolor=None)
                        fig.colorbar(h)
                        ax.set_xlabel('T 2m')
                        ax.set_ylabel('T '+str(lev)+' m')
                        ax.set_title(box)
                        fig.savefig(os.path.join(plot_dir, "T2_vs_T"+str(lev)+"_"+box+".png"))

                        # scatter Tair for each level vs P250
                        fig, ax = plt.subplots()
                        h = ax.scatter(np.ravel(tair[lev]), np.ravel(p250), c=np.ravel(hours_box), edgecolor=None)
                        fig.colorbar(h)
                        ax.set_xlabel('T '+str(lev)+' m')
                        ax.set_ylabel('P 250 m')
                        ax.set_title(box)
                        fig.savefig(os.path.join(plot_dir, "T"+str(lev)+"_vs_p250_"+box+".png"))

                plt.close('all')

