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
box_pts = {"bay": {"d01": [(42, 48), (44, 52)], "d02": [(70, 88), (76, 100)]}}
mintime = 42
maxtime = -1
utc_offset = 8  # offset to subtract to get local time
fill_val_thresh = 1e36

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

def calc_distance_matrix(ij1, f):
    # calculate distance between ij1 and ij2
    # i is lon/x, j is lat/y
    i_arr = np.arange(0, f.variables["XLAT"].shape[-1])
    j_arr = np.arange(0, f.variables["XLAT"].shape[-2])
    i_grid, j_grid = np.meshgrid(i_arr, j_arr)
    di = i_grid-ij1[0]
    dj = j_grid-ij1[1]
    dx = f.variables['DX'][0]
    dy = f.variables['DY'][0]
    distance = ( (di*dx)**2 + (dj*dy)**2 )**0.5
    distance = distance/1000.  # convert m -> km
    distance[distance == 0] = np.nan
    return distance

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

def avg_press(pressure, domain):
    # corner indices are lon, lat
    llcnr = box_pts["bay"][domain][0]
    urcnr = box_pts["bay"][domain][1]
    pslice = pressure[:, llcnr[1]:urcnr[1], llcnr[0]:urcnr[0]]
    prtrn = np.zeros(pslice.shape[0])+np.nan
    for i in xrange(len(prtrn)):
        ptmp = pslice[i, :, :]
        ptmp = ptmp[np.isfinite(ptmp) & (ptmp < fill_val_thresh)]
        if len(ptmp)>0:
            prtrn[i] = np.mean(ptmp)
    return prtrn

def calc_corr(p, u, lag):
    # pressure: time x lat x lon
    # u: time

    corr_coefs = np.zeros((p.shape[1], p.shape[2])) + np.nan
    slope_vals = np.zeros((p.shape[1], p.shape[2])) + np.nan
    for i_lat in xrange(p.shape[1]):
        for j_lon in xrange(p.shape[2]):
            ptmp = p[:, i_lat, j_lon]
            if lag > 0:
                ptmp = ptmp[:-lag]
            utmp = u[lag:]
            mask = np.isfinite(ptmp) & np.isfinite(utmp) & (ptmp < fill_val_thresh) & (utmp < fill_val_thresh)
            if len(ptmp[mask]) > 5:  # arbitrarily picked this number
                slope, intc, r, pval, se = linregress(ptmp[mask], utmp[mask])
                #corr_coefs[i_lat, j_lon], p_vals[i_lat, j_lon] = pearsonr(ptmp[mask], utmp[mask])
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

def plot_corr(f, p, u, lag, maxmin):

    fig, ax, mlist = setup_map(f, nrows=2, ncols=2)
    
    corr_coefs0, slope0 = calc_corr(p, u[:,0], lag)
    corr_coefs1, slope1 = calc_corr(p, u[:,1], lag)

    for i in xrange(4):
        m = mlist[i]
        a = ax[i]
        if i == 0:
            data = np.ma.masked_array(slope0, np.isnan(slope0))
            title = "slope, u lev 0"
            cc = False  # flag for corr coef
        elif i == 1:
            data = np.ma.masked_array(corr_coefs0, np.isnan(corr_coefs0))
            title = "corr. coef.s, u lev 0"
            cc = True  # flag for corr coef
        elif i == 2:
            data = np.ma.masked_array(slope1, np.isnan(slope1))
            title = "slope, u lev 1"
            cc = False  # flag for corr coef
        elif i == 3:
            data = np.ma.masked_array(corr_coefs1, np.isnan(corr_coefs1))
            title = "corr. coef.s, u lev 1"
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

    #plt.show()
    #1/0

    return fig, ax #, slope1, corr_coefs1


if __name__ == "__main__":

    ctrl_name = "CA-0.25"
    test_names = ["CA-dryCR", "CA-dryCV", "CA-drySN"]
    #ctrl_name = "CA-0.2"
    #test_names = ["CA-CV0.05w", "CA-CV0.35w"]
    domains = ["d02"]  #["d01", "d02"]
    lags = [1] # [0, 2, 3]
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

            plot_dir = os.path.join(root_dir, run_name, "plots", domain, "corr_p_u")
            if not os.path.exists(plot_dir):
                os.mkdir(plot_dir)

            # get hour of day array from times array
            #hours = get_hours(f)

            # get box center, and distance for each point to bay, sac, baybox
            ij_box = box_center(*box_pts["bay"][domain])
            d_bay = calc_distance_matrix(grad_pts["bay"][domain], fi)
            d_sac = calc_distance_matrix(grad_pts["sac"][domain], fi)
            d_box = calc_distance_matrix(ij_box, fi)

            # get u, v at solano
            usolano, vsolano, speed = get_winds(fi, domain)
            if run_name == ctrl_name:
                control_vals["usolano"] = usolano
                control_vals["vsolano"] = vsolano
                control_vals["speed"] = speed
                control_vals["f"] = fi
                control_vals["pressure"] = {}
                control_vals["pgrad_sac"] = {}
                control_vals["pgrad_bay"] = {}
                control_vals["pgrad_baybox"] = {}

                # get min/max lat and lon
                maxmin = {}
                maxmin['minlon'] = np.floor(np.min(fi.variables['XLONG'][:,:]/10.))*10.
                maxmin['maxlon'] = np.ceil(np.max(fi.variables['XLONG'][:,:]/10.))*10.
                maxmin['minlat'] = np.floor(np.min(fi.variables['XLAT'][:,:]/10.))*10.
                maxmin['maxlat'] = np.ceil(np.max(fi.variables['XLAT'][:,:]/10.))*10.

            for i_lev, level in enumerate(fi.variables["levels"]):

                print level

                # get pressure and ocean pressure (if ctrl run) or delta pressure (if test run)
                pressure = fi.variables["P"][mintime:maxtime, i_lev, :, :]
                #ocnpressure = pressure[:, grad_pts["ocn"][domain][1], grad_pts["ocn"][domain][0]]
                baypressure = pressure[:, grad_pts["bay"][domain][1], grad_pts["bay"][domain][0]]
                sacpressure = pressure[:, grad_pts["sac"][domain][1], grad_pts["sac"][domain][0]]
                if run_name == ctrl_name:
                    control_vals["pressure"][level] = pressure

                # get pressure timeseries averaged over bay box
                p_baybox = avg_press(pressure, domain)

                # get grad pressure from sac (if ctrl) or delta grad pressure (if test run)
                #pressure_m_ocn = pressure-ocnpressure.reshape(ocnpressure.shape[0], 1, 1)
                pressure_m_sac = (pressure-sacpressure.reshape(sacpressure.shape[0], 1, 1))/d_sac
                pressure_m_bay = (pressure-baypressure.reshape(baypressure.shape[0], 1, 1))/d_bay
                pressure_m_baybox = (pressure-p_baybox.reshape(p_baybox.shape[0], 1, 1))/d_box
                if run_name == ctrl_name:
                    control_vals["pgrad_sac"][level] = pressure_m_sac
                    control_vals["pgrad_bay"][level] = pressure_m_bay
                    control_vals["pgrad_baybox"][level] = pressure_m_baybox

                for lag in lags:

                    print lag

                    # correlation map u vs pressure
                    #fig, ax = plot_corr(fr, pressure, speed, lag, maxmin)
                    #fig.suptitle(run_name+" u vs pressure "+str(level)+" m, lag "+str(lag*timestep)+" min")
                    #plt.show()
                    #1/0
                    #fig.savefig(os.path.join(plot_dir, "corr_u_p_lev"+str(round(level))+"_lag"+str(lag)+".png"))

                    # # correlation map u vs pressure minus bay pressure
                    # fig, ax = plot_corr(fr, pressure_m_bay, speed, lag, maxmin)
                    # fig.suptitle(run_name+" u vs pressure minus bay "+str(level)+" m, lag "+str(lag*timestep)+" min")
                    # fig.savefig(os.path.join(plot_dir, "corr_u_p-bay_lev"+str(round(level))+"_lag"+str(lag)+".png"))

                    # correlation map u vs pressure minus bay box pressure
                    fig, ax = plot_corr(fr, pressure_m_baybox, speed, lag, maxmin)
                    #1/0
                    fig.suptitle(run_name+" u vs pressure minus bay box "+str(level)+" m, lag "+str(lag*timestep)+" min")
                    fig.savefig(os.path.join(plot_dir, "corr_u_p-baybox_lev"+str(round(level))+"_lag"+str(lag)+".png"))

                    # # correlation map u vs pressure gradient from sac point
                    # fig, ax = plot_corr(fr, pressure_m_sac, speed, lag, maxmin)
                    # fig.suptitle(run_name+" u vs pressure minus sac "+str(level)+" m, lag "+str(lag*timestep)+" min")
                    # fig.savefig(os.path.join(plot_dir, "corr_u_p-sac_lev"+str(round(level))+"_lag"+str(lag)+".png"))

                    # for test cases:
                    if run_name != ctrl_name:

                        # correlation map du vs dpressure
                        fig, ax = plot_corr(fr, pressure-control_vals["pressure"][level], speed-control_vals["speed"], lag, maxmin)
                        fig.suptitle(run_name+" du vs dpressure "+str(level)+" m, lag "+str(lag*timestep)+" min")
                        fig.savefig(os.path.join(plot_dir, "corr_du_dp_lev"+str(round(level))+"_lag"+str(lag)+".png"))

                        # correlation map du vs d(pressure gradient from sac point)
                        fig, ax = plot_corr(fr, pressure_m_sac-control_vals["pgrad_sac"][level], speed-control_vals["speed"], lag, maxmin)
                        fig.suptitle(run_name+" du vs d(pressure minus sac) "+str(level)+" m, lag "+str(lag*timestep)+" min")
                        fig.savefig(os.path.join(plot_dir, "corr_du_d(p-sac)_lev"+str(round(level))+"_lag"+str(lag)+".png"))

                        # correlation map du vs d(pressure gradient from bay point)
                        fig, ax = plot_corr(fr, pressure_m_bay-control_vals["pgrad_bay"][level], speed-control_vals["speed"], lag, maxmin)
                        fig.suptitle(run_name+" du vs d(pressure minus bay) "+str(level)+" m, lag "+str(lag*timestep)+" min")
                        fig.savefig(os.path.join(plot_dir, "corr_du_d(p-bay)_lev"+str(round(level))+"_lag"+str(lag)+".png"))

                        # correlation map du vs d(pressure gradient from bay box)
                        fig, ax = plot_corr(fr, pressure_m_baybox-control_vals["pgrad_baybox"][level], speed-control_vals["speed"], lag, maxmin)
                        fig.suptitle(run_name+" du vs d(pressure minus baybox) "+str(level)+" m, lag "+str(lag*timestep)+" min")
                        fig.savefig(os.path.join(plot_dir, "corr_du_d(p-baybox)_lev"+str(round(level))+"_lag"+str(lag)+".png"))

                    plt.close("all")

                
