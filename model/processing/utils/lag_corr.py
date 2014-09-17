"""
PAL 9/11/14
calculate correlation between pairs of variables for a range of lags
"""

import scipy.io.netcdf as nc
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import os


fill_val_thresh = 1e36
lat_solano, lon_solano = 38.166, -121.817
box_pts = {"bay": {"d01": [(42, 48), (44, 52)], "d02": [(70, 88), (76, 100)]},
           "cv_small": {"d01": [(46, 49), (49, 52)], "d02": [(83, 93), (91, 98)]},
           "cv_big": {"d01": [(48, 50), (53, 54)], "d02": [(89, 94), (101, 104)]},
           "cv_sac": {"d01": [(50, 48), (54, 55)], "d02": [(95, 90), (105, 107)]}}
mintime = 42
maxtime = -1


def calc_lag_corr(v1, v2, lags):
    # v1 and v2: dimension time
    # v1 leads v2
    #lags = np.arange(12)
    lag_corr = np.zeros_like(lags)+np.nan
    lag_slope = np.zeros_like(lags)+np.nan
    for i_lag, lag in enumerate(lags):
        if lag > 0:
            v1tmp = v1[:-lag]
        else:
            v1tmp = v1
        v2tmp = v2[lag:]
        mask = np.isfinite(v1tmp) & np.isfinite(v2tmp) & (v1tmp < fill_val_thresh) & (v2tmp < fill_val_thresh)
        if len(v1tmp[mask]) > 5:  # arbitrarily picked this number
            lag_slope[i_lag], intc, lag_corr[i_lag], pval, se = linregress(v1tmp[mask], v2tmp[mask])
    return lag_slope, lag_corr

def box_center(llcnr, urcnr):
    x = int((llcnr[0]+urcnr[0])/2)
    y = int((llcnr[1]+urcnr[1])/2)
    return (x, y)

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

def get_winds(f):
    dif_lat = f.variables['XLAT'][:, :]-lat_solano
    dif_lon = f.variables['XLONG'][:, :]-lon_solano
    distance = dif_lat**2 + dif_lon**2
    ixlat, ixlon = np.unravel_index(np.argmin(distance), distance.shape)
    usolano = f.variables['U'][mintime:maxtime, :2, ixlat, ixlon]
    vsolano = f.variables['V'][mintime:maxtime, :2, ixlat, ixlon]
    speed = (usolano**2 + vsolano**2)**0.5  # calculate speed
    return usolano, vsolano, speed

def avg_press(pressure, box, domain):
    # corner indices are lon, lat
    # pressure is time x lat x lon
    llcnr = box_pts[box][domain][0]
    urcnr = box_pts[box][domain][1]
    pslice = pressure[:, llcnr[1]:urcnr[1], llcnr[0]:urcnr[0]]
    prtrn = np.zeros(pslice.shape[0])+np.nan  # dimension time
    for i in xrange(len(prtrn)):
        ptmp = pslice[i, :, :]
        ptmp = ptmp[np.isfinite(ptmp) & (ptmp < fill_val_thresh)]
        if len(ptmp)>0:
            prtrn[i] = np.mean(ptmp)
    return prtrn



# function to loop through levels, calculate average for the box(es), call calc_lag_corr, plot
if __name__ == "__main__":

    #run_name = "CA-0.2"
    run_list = ["CA-0.2", "CA-0.1", "CA-0.25", "CA-dryCR", "CA-drySN", "CA-CV0.05w", "CA-CV0.35w"]
    domains = ["d01", "d02"]
    lags = np.arange(12)
    root_dir = "/scratch2/scratchdirs/plink/WRF/output"
    timestep = 0.5  # output timestep in hours
    colors = ['b', 'c', 'g', 'y', 'r', 'm', 'k']

    for run_name in run_list:

        print run_name

        figD, axD = plt.subplots(nrows=2, ncols=2)
        plot_dir_D = os.path.join(root_dir, run_name, "plots", "lag_corr")
        if not os.path.exists(plot_dir_D):
            os.mkdir(plot_dir_D)

        for i_dom, domain in enumerate(domains):

            print domain

            fi = nc.netcdf_file(os.path.join(root_dir, run_name, "wrfout_interp_"+domain+".nc"))

            plot_dir = os.path.join(root_dir, run_name, "plots", domain, "lag_corr")
            if not os.path.exists(plot_dir):
                os.mkdir(plot_dir)

            # get u, v at solano
            usolano, vsolano, speed = get_winds(fi)
            fig, ax = plt.subplots(nrows=3, ncols=2)
            handles = []
            labels = []

            for i_lev, level in enumerate(fi.variables["levels"]):

                print level

                # get pressure av timeseries for each of the boxes
                pressure = fi.variables["P"][mintime:maxtime, i_lev, :, :]
                p_bay = avg_press(pressure, "bay", domain)
                p_cvsmall = avg_press(pressure, "cv_small", domain)
                p_cvbig = avg_press(pressure, "cv_big", domain)
                p_cvsac = avg_press(pressure, "cv_sac", domain)
                
                # get locations of box centers and distances
                ij_bay = box_center(*box_pts["bay"][domain])
                ij_cvsmall = box_center(*box_pts["cv_small"][domain])
                ij_cvbig = box_center(*box_pts["cv_big"][domain])
                ij_cvsac = box_center(*box_pts["cv_sac"][domain])
                d_baycvsmall = calc_distance(ij_bay, ij_cvsmall, fi)
                d_baycvbig = calc_distance(ij_bay, ij_cvbig, fi)
                d_baycvsac = calc_distance(ij_bay, ij_cvsac, fi)
                
                # get lag slope and lag corr
                lag_slope_cvsmall_u0, lag_corr_cvsmall_u0 = calc_lag_corr((p_cvsmall-p_bay)/d_baycvsmall, speed[:, 0], lags)
                lag_slope_cvsmall_u1, lag_corr_cvsmall_u1 = calc_lag_corr((p_cvsmall-p_bay)/d_baycvsmall, speed[:, 1], lags)

                lag_slope_cvbig_u0, lag_corr_cvbig_u0 = calc_lag_corr((p_cvbig-p_bay)/d_baycvbig, speed[:, 0], lags)
                lag_slope_cvbig_u1, lag_corr_cvbig_u1 = calc_lag_corr((p_cvbig-p_bay)/d_baycvbig, speed[:, 1], lags)

                lag_slope_cvsac_u0, lag_corr_cvsac_u0 = calc_lag_corr((p_cvsac-p_bay)/d_baycvsac, speed[:, 0], lags)
                lag_slope_cvsac_u1, lag_corr_cvsac_u1 = calc_lag_corr((p_cvsac-p_bay)/d_baycvsac, speed[:, 1], lags)

                # plot (rows: cv small vs cv big, cols: slope and corr)
                h0 = ax[0, 0].plot(lags*timestep, lag_slope_cvsmall_u0, linestyle="-", color=colors[i_lev])
                h1 = ax[0, 0].plot(lags*timestep, lag_slope_cvsmall_u1, linestyle="--", color=colors[i_lev])
                handles += [h0[0]]#, h1[0]]
                labels += ["p "+str(int(level))+" m"] #, u 60 m", "p "+str(int(level))+" m, u 100 m"]

                ax[0, 1].plot(lags*timestep, lag_corr_cvsmall_u0, linestyle="-", color=colors[i_lev])
                ax[0, 1].plot(lags*timestep, lag_corr_cvsmall_u1, linestyle="--", color=colors[i_lev])

                ax[1, 0].plot(lags*timestep, lag_slope_cvbig_u0, linestyle="-", color=colors[i_lev])
                ax[1, 0].plot(lags*timestep, lag_slope_cvbig_u1, linestyle="--", color=colors[i_lev])

                ax[1, 1].plot(lags*timestep, lag_corr_cvbig_u0, linestyle="-", color=colors[i_lev])
                ax[1, 1].plot(lags*timestep, lag_corr_cvbig_u1, linestyle="--", color=colors[i_lev])

                ax[2, 0].plot(lags*timestep, lag_slope_cvsac_u0, linestyle="-", color=colors[i_lev])
                ax[2, 0].plot(lags*timestep, lag_slope_cvsac_u1, linestyle="--", color=colors[i_lev])

                ax[2, 1].plot(lags*timestep, lag_corr_cvsac_u0, linestyle="-", color=colors[i_lev])
                ax[2, 1].plot(lags*timestep, lag_corr_cvsac_u1, linestyle="--", color=colors[i_lev])

                # add to multi-domain plot
                axD[0, i_dom].plot(lags*timestep, lag_slope_cvsac_u0, linestyle="-", color=colors[i_lev])
                axD[0, i_dom].plot(lags*timestep, lag_slope_cvsac_u1, linestyle="--", color=colors[i_lev])

                axD[1, i_dom].plot(lags*timestep, lag_corr_cvsac_u0, linestyle="-", color=colors[i_lev])
                axD[1, i_dom].plot(lags*timestep, lag_corr_cvsac_u1, linestyle="--", color=colors[i_lev])

            # format single-domain figure
            ax[0,0].set_ylabel("u vs pcvsmall-pbay")
            ax[1,0].set_ylabel("u vs pcvbig-pbay")
            ax[2,0].set_ylabel("u vs pcvsac-pbay")
            ax[0,0].set_title("slope (m/s (hPa/km)$^{-1}$)")
            ax[0,1].set_title("r")
            for a in ax.flatten():
                a.set_xlabel("lag (wind behind pressure, hr)")
            
            # make color/pressure legend
            fig.legend(handles, labels, loc='lower center', ncol=3, prop={'size': 10})

            # make dummy linestyle lines for u legend
            solidline = plt.Line2D([0], [0], linestyle='-', color='k')
            dashedline = plt.Line2D([0], [0], linestyle='--', color='k')
            fig.legend([solidline, dashedline], ["speed 60 m AGL", "speed 100 m AGL"], loc='lower right', prop={"size": 10})

            fig.subplots_adjust(bottom=0.17, hspace=0.35)
            fig.suptitle(run_name+" "+domain)
            fig.savefig(os.path.join(plot_dir, "lag_corr_p_u_"+domain+".png"))
            #plt.show()
            #1/0
            plt.close(fig)

        # format multi-domain figure
        axD[0,0].set_ylabel("slope (m/s (hPa/km)$^{-1}$)")
        axD[1,0].set_ylabel("r")
        axD[0,0].set_title(domains[0])
        axD[0,1].set_title(domains[1])
        for a in ax.flatten():
            a.set_xlabel("lag (wind behind pressure, hr)")

        axD[0,0].set_ylim(-9, 0)
        axD[0,1].set_ylim(-9, 0)
        axD[1,0].set_ylim(-0.9, 0)
        axD[1,1].set_ylim(-0.9, 0)
        
        # make color/pressure legend
        figD.legend(handles, labels, loc='lower center', ncol=3, prop={'size': 10})

        # make dummy linestyle lines for u legend
        #solidline = plt.Line2D([0], [0], linestyle='-', color='k')
        #dashedline = plt.Line2D([0], [0], linestyle='--', color='k')
        figD.legend([solidline, dashedline], ["speed 60 m AGL", "speed 100 m AGL"], loc='lower right', prop={"size": 10})

        figD.subplots_adjust(bottom=0.17, hspace=0.35)
        figD.suptitle("pgrad sac-bay vs u "+domain)
        figD.savefig(os.path.join(plot_dir_D, "lag_corr_p_u.png"))

        plt.close(figD)

