"""
PAL 9/5/2014
map of correlation between pressure and wind at solano
"""

import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.io import netcdf as nc
from scipy.stats import linregress, pearsonr
from mpl_toolkits.basemap import Basemap


coords_solano = {'d01': (49, 51), 'd02': (91, 96), 'd03': (127, 143)}
lat_solano, lon_solano = 38.166, -121.817
grad_pts = {'bay': {'d01': (43,49), 'd02': (70,88), 'd03': (61,115)},
            'ocn': {'d01': (35,45), 'd02': (46,76), 'd03': (25,99)},
            'sac': {'d01': (53,53), 'd02': (100,100), 'd03': (151,151)},
            'fth': {'d01': (60,56), 'd02': (121,109), 'd03': (183,168)},
            'red': {'d01': (45,77), 'd02': (84,150), 'd03': (None,None)},
            'bkf': {'d01': (77,18), 'd02': (151,32), 'd03': (None,None)}}
box_pts = {"bay": {"d01": [(42, 48), (44, 52)], "d02": [(70, 88), (76, 100)]}, 
           "cv": {"d01": [(48, 50), (53, 54)], "d02": [(89, 94), (101, 104)]}}
mintime = 42
maxtime = -1
utc_offset = 8  # offset to subtract to get local time
fill_val_thresh = 1e36

def get_hours(f):
    hours = []
    times_arr = f.variables['Times'][mintime:maxtime]
    for ii in xrange(times_arr.shape[0]):
        time_arr = times_arr[ii, :]
        hour = int(''.join(time_arr[11:13]))
        hour = hour-utc_offset
        if hour < 0:
            hour = hour+24
        hours.append(hour)
    hours = np.array(hours)
    return hours

def get_winds(f, domain):
    usolano = f.variables['U'][mintime:maxtime, :2, coords_solano[domain][1], coords_solano[domain][0]]
    vsolano = f.variables['V'][mintime:maxtime, :2, coords_solano[domain][1], coords_solano[domain][0]]
    speed = (usolano**2 + vsolano**2)**0.5  # calculate speed
    return usolano, vsolano, speed

def avg_press(pressure, box, domain):
    # corner indices are lon, lat
    llcnr = box_pts[box][domain][0]
    urcnr = box_pts[box][domain][1]
    pslice = pressure[:, llcnr[1]:urcnr[1], llcnr[0]:urcnr[0]]
    prtrn = np.zeros(pslice.shape[0])+np.nan
    for i in xrange(len(prtrn)):
        ptmp = pslice[i, :, :]
        ptmp = ptmp[np.isfinite(ptmp) & (ptmp < fill_val_thresh)]
        if len(ptmp)>0:
            prtrn[i] = np.mean(ptmp)
    return prtrn

def calc_corr(T, p, lag):
    # T: time x lat x lon
    # p: time

    corr_coefs = np.zeros((T.shape[1], T.shape[2])) + np.nan
    slope_vals = np.zeros((T.shape[1], T.shape[2])) + np.nan
    for i_lat in xrange(T.shape[1]):
        for j_lon in xrange(T.shape[2]):
            Ttmp = T[:, i_lat, j_lon]
            if lag > 0:
                Ttmp = Ttmp[:-lag]
            ptmp = p[lag:]
            mask = np.isfinite(ptmp) & np.isfinite(Ttmp) & (ptmp < fill_val_thresh) & (Ttmp < fill_val_thresh)
            if len(ptmp[mask]) > 5:  # arbitrarily picked this number
                slope, intc, r, pval, se = linregress(ptmp[mask], Ttmp[mask])
                slope_vals[i_lat, j_lon] = slope
                corr_coefs[i_lat, j_lon] = r
    return corr_coefs, slope_vals

def setup_map(f, nrows=0, ncols=0):

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols)
    ax = ax.flatten()
    mlist = []
    for a in ax:
        m = Basemap(width=f.DX*1.2*getattr(f,'WEST-EAST_GRID_DIMENSION'),\
            height=f.DY*1.2*getattr(f,'SOUTH-NORTH_GRID_DIMENSION'),resolution='l',\
            projection='lcc',lat_1=f.TRUELAT1,lat_2=f.TRUELAT2,lat_0=f.CEN_LAT,lon_0=f.CEN_LON, ax=a)
        mlist.append(m)
    return fig, ax, mlist

def decorate_map(m=None, a=None, h=None, fig=None, maxmin=None):

    m.drawparallels(np.arange(maxmin['minlat'], maxmin['maxlat'], 5))
    m.drawmeridians(np.arange(maxmin['minlon'], maxmin['maxlon'], 5))
    m.drawcoastlines()
    m.drawstates()
    fig.colorbar(h, ax=a)

def plot_corr(f, T, p, lag, maxmin):
    # T is time x lat x lon
    # p is time

    fig, ax, mlist = setup_map(f, nrows=1, ncols=2)
    
    corr_coefs, slope = calc_corr(T, p, lag)
    
    for i in xrange(2):
        m = mlist[i]
        a = ax[i]
        if i == 0:
            data = np.ma.masked_array(slope, np.isnan(slope))
            title = "slope"
            cc = False  # flag for corr coef
        elif i == 1:
            data = np.ma.masked_array(corr_coefs, np.isnan(corr_coefs))
            title = "corr. coef.s"
            cc = True  # flag for corr coef

        x,y = m(f.variables['XLONG'][0, :, :], f.variables['XLAT'][0, :, :])
        h = m.pcolormesh(x, y, data)
        if cc:
            h.set_clim(vmin=-0.8, vmax=0.8)
        #else:
        #    h.set_clim(vmin=0, vmax=0.3)
        a.set_title(title)
        decorate_map(m=m, a=a, h=h, fig=fig, maxmin=maxmin)
        xsolano, ysolano = m(lon_solano, lat_solano)
        a.plot(xsolano, ysolano, 'md')

    return fig, ax


if __name__ == "__main__":

    plevel = 2  # 2=250mASL

    #ctrl_name = "CA-0.25"
    #test_names = ["CA-dryCR", "CA-dryCV", "CA-drySN"]
    ctrl_name = "CA-0.2"
    test_names = ["CA-CV0.05w", "CA-CV0.35w"]
    domains = ["d02"]  #["d01", "d02"]
    lags = [0] # [0, 2, 3]
    root_dir = "/scratch2/scratchdirs/plink/WRF/output"
    timestep = 30  # output timestep in minutes

    #for run_name in [ctrl_name]+test_names:
    for domain in domains:

        print domain
        #if run_name == ctrl_name:
        #    control_vals = {}

        for run_name in [ctrl_name]+test_names:

            print run_name
            if run_name == ctrl_name:
                control_vals = {}

            fi = nc.netcdf_file(os.path.join(root_dir, run_name, "wrfout_interp_"+domain+".nc"))
            fr = nc.netcdf_file(os.path.join(root_dir, run_name, "wrfout_"+domain+"_2009-07-01_00:00:00"))

            plot_dir = os.path.join(root_dir, run_name, "plots", domain, "corr_Tair_p")
            if not os.path.exists(plot_dir):
                os.mkdir(plot_dir)

            # get pressure in box
            pressure = fi.variables["P"][mintime:maxtime, plevel, :, :]
            p_box = avg_press(pressure, "cv", domain)
            if run_name == ctrl_name:
                control_vals["f"] = fi
                control_vals["p_box"] = p_box
                control_vals["T"] = {}
                
                # get min/max lat and lon
                maxmin = {}
                maxmin['minlon'] = np.floor(np.min(fi.variables['XLONG'][:,:]/10.))*10.
                maxmin['maxlon'] = np.ceil(np.max(fi.variables['XLONG'][:,:]/10.))*10.
                maxmin['minlat'] = np.floor(np.min(fi.variables['XLAT'][:,:]/10.))*10.
                maxmin['maxlat'] = np.ceil(np.max(fi.variables['XLAT'][:,:]/10.))*10.

            for i_lev, level in enumerate(fi.variables["levels"]):

                print level

                temperature = fi.variables["T"][mintime:maxtime, i_lev, :, :].copy()
                temperature[temperature >= fill_val_thresh] = np.nan
                if run_name == ctrl_name:
                    control_vals["T"][level] = temperature

                for lag in lags:

                    print lag

                    # correlation map pbox vs temperature
                    fig, ax = plot_corr(fr, temperature, p_box, lag, maxmin)
                    fig.suptitle(run_name+" pressure cv box vs T "+str(level)+" m, lag "+str(lag*timestep)+" min")
                    fig.savefig(os.path.join(plot_dir, "corr_pbox_T_lev"+str(round(level))+"_lag"+str(lag)+".png"))

                    # for test cases:
                    if run_name != ctrl_name:

                        # correlation map dpbox vs dtemperature
                        fig, ax = plot_corr(fr, temperature-control_vals["T"][level], p_box-control_vals["p_box"], lag, maxmin)
                        fig.suptitle(run_name+" dpbox vs dT "+str(level)+" m, lag "+str(lag*timestep)+" min")
                        #plt.show()
                        #1/0
                        fig.savefig(os.path.join(plot_dir, "corr_dpbox_dT_lev"+str(round(level))+"_lag"+str(lag)+".png"))

                    plt.close("all")

                
