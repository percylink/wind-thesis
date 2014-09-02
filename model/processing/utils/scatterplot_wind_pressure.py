"""
PAL 9/1/2014
Scatter plots of WRF output
"""

import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.io import netcdf as nc


coords_solano = {'d01': (49, 51), 'd02': (91, 96), 'd03': (127, 143)}
grad_pts = {'bay': {'d01': (43,49), 'd02': (70,88), 'd03': (61,115)},
            'ocn': {'d01': (35,45), 'd02': (46,76), 'd03': (25,99)},
            'sac': {'d01': (53,53), 'd02': (100,100), 'd03': (151,151)},
            'fth': {'d01': (60,56), 'd02': (121,109), 'd03': (183,168)},
            'red': {'d01': (45,77), 'd02': (84,150), 'd03': (None,None)},
            'bkf': {'d01': (77,18), 'd02': (151,32), 'd03': (None,None)}}
mintime = 42
maxtime = -1


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

def get_p(f, pt, domain):
    ij = grad_pts[pt][domain]
    pressure = f.variables['P'][mintime:maxtime, :, ij[1], ij[0]].copy()
    mask = pressure > 1e36
    pressure[mask] = np.nan
    return pressure

def calc_distance(pt, domain, f):
    # calculate distance between ij_pt and ij_sac
    ij_pt = grad_pts[pt][domain]
    ij_sac = grad_pts['sac'][domain]
    di = ij_pt[0]-ij_sac[0]
    dj = ij_pt[1]-ij_sac[1]
    dx = f.variables['DX'][0]
    dy = f.variables['DY'][0]
    distance = ( (di*dx)**2 + (dj*dy)**2 )**0.5
    distance = distance/1000.  # convert m -> km
    return distance

def plot_p_vs_wind(p, u, v, speed, f, run_name, pt, plot_dir):
    fig, ax = plt.subplots(nrows=p.shape[1], ncols=3)
    for i_lev in xrange(p.shape[1]):
        ax[i_lev, 0].plot(p[:, i_lev], u, '.')
        ax[i_lev, 1].plot(p[:, i_lev], v, '.')
        ax[i_lev, 2].plot(p[:, i_lev], speed, '.')
        ax[i_lev, 0].set_ylabel(str(f.variables['levels'][i_lev])+" m")
    ax[0, 0].set_title('U')
    ax[0, 1].set_title('V')
    ax[0, 2].set_title('speed')
    fig.suptitle(run_name+" x: pressure "+pt+", y: wind solano")
    fig.set_size_inches(6, p.shape[1]*2)
    fig.savefig(os.path.join(plot_dir, "pressure_vs_wind_"+pt+".png"))
    plt.close('all') 


if __name__ == "__main__":

    ctrl_name = "CA-0.24"
    test_names = ["CA-dryCR", "CA-dryCV", "CA-drySN"]
    #domains = ["d01", "d02", "d03"]
    #test_names = ["CA-dryCR"]
    domains = ["d02", "d03"]
    sm_test = 0.05  # SMOIS value in the test region
    utc_offset = 8  # offset to subtract to get local time
    root_dir = "/scratch2/scratchdirs/plink/WRF/output"


    for run_name in [ctrl_name]+test_names:

        print run_name
        if run_name == ctrl_name:
            control_vals = {}

        for domain in domains:

            print domain
            if run_name == ctrl_name:
                control_vals[domain] = {}

            f = nc.netcdf_file(os.path.join(root_dir, run_name, "wrfout_interp_"+domain+".nc"))
            if run_name != ctrl_name:
                fc = nc.netcdf_file(os.path.join(root_dir, ctrl_name, "wrfout_interp_"+domain+".nc"))

            plot_dir = os.path.join(root_dir, run_name, "plots", domain, "scatter_wind_pressure")
            if not os.path.exists(plot_dir):
                os.mkdir(plot_dir)

            # get hour of day array from times array
            hours = get_hours(f)

            # get u, v at solano
            usolano, vsolano, speed = get_winds(f, domain)

            # get p at each gradient point
            pdict = {}
            for pt in grad_pts:

                ij = grad_pts[pt][domain]
                if ij[0] is None:
                    continue

                pdict[pt] = get_p(f, pt, domain)
                
                # for each of the p locations, scatter u, v, speed vs each level of pressure
                plot_p_vs_wind(pdict[pt], usolano, vsolano, speed, f, run_name, pt, plot_dir)

                fig, ax = plt.subplots(nrows=pdict[pt].shape[1], ncols=3)
                for i_lev in xrange(pdict[pt].shape[1]):
                    ax[i_lev, 0].plot(pdict[pt][:, i_lev], usolano, '.')
                    ax[i_lev, 1].plot(pdict[pt][:, i_lev], vsolano, '.')
                    ax[i_lev, 2].plot(pdict[pt][:, i_lev], speed, '.')
                    ax[i_lev, 0].set_ylabel(str(f.variables['levels'][i_lev])+" m")
                ax[0, 0].set_title('U')
                ax[0, 1].set_title('V')
                ax[0, 2].set_title('speed')
                fig.suptitle(run_name+" x: pressure "+pt+", y: wind solano")
                fig.set_size_inches(6, pdict[pt].shape[1]*2)
                fig.savefig(os.path.join(plot_dir, "pressure_vs_wind_"+pt+".png")) 

                # if this is a diff case, scatter diffs in u, v, speed vs diffs in pressure
                if run_name != ctrl_name:
                    pdiff = pdict[pt]-control_vals[domain]['p'][pt]
                    udiff = usolano-control_vals[domain]['u']
                    vdiff = vsolano-control_vals[domain]['v']
                    spdiff = speed-control_vals[domain]['speed']

                    # for each of the p locations, scatter diffs of u, v, speed vs each level of pressure
                    fig, ax = plt.subplots(nrows=pdiff.shape[1], ncols=3)
                    for i_lev in xrange(pdiff.shape[1]):
                        ax[i_lev, 0].plot(pdiff[:, i_lev], udiff, '.')
                        ax[i_lev, 1].plot(pdiff[:, i_lev], vdiff, '.')
                        ax[i_lev, 2].plot(pdiff[:, i_lev], spdiff, '.')
                        ax[i_lev, 0].set_ylabel(str(f.variables['levels'][i_lev])+" m")
                    ax[0, 0].set_title('U diff')
                    ax[0, 1].set_title('V diff')
                    ax[0, 2].set_title('speed diff')
                    fig.suptitle(run_name+" minus ctrl\n x: pressure diff "+pt+", y: wind diff solano")
                    fig.set_size_inches(6, pdiff.shape[1]*2)
                    fig.savefig(os.path.join(plot_dir, "pressure_vs_wind_DIFF_"+pt+".png")) 

                else:
                    control_vals[domain]['p'] = pdict
                    control_vals[domain]['u'] = usolano
                    control_vals[domain]['v'] = vsolano
                    control_vals[domain]['speed'] = speed

            plt.close('all')

            for pt in grad_pts:

                # for each of the gradients, scatter u, v, speed vs each level of gradient
                if (pt != 'sac') and (pt in pdict):
                    
                    distance = calc_distance(pt, domain, f)
                    gradP = (pdict[pt] - pdict['sac']) / distance

                    fig, ax = plt.subplots(nrows=gradP.shape[1], ncols=3)
                    for i_lev in xrange(gradP.shape[1]):
                        ax[i_lev, 0].plot(gradP[:, i_lev], usolano, '.')
                        ax[i_lev, 1].plot(gradP[:, i_lev], vsolano, '.')
                        ax[i_lev, 2].plot(gradP[:, i_lev], speed, '.')
                        ax[i_lev, 0].set_ylabel(str(f.variables['levels'][i_lev])+" m")
                    ax[0, 0].set_title('U')
                    ax[0, 1].set_title('V')
                    ax[0, 2].set_title('speed')
                    fig.suptitle(run_name+" x: press grad "+pt+", y: wind solano")
                    fig.set_size_inches(6, gradP.shape[1]*2)
                    fig.savefig(os.path.join(plot_dir, "pgrad_vs_wind_"+pt+".png")) 

                    # if this is a diff case, scatter u, v, speed vs diffs in gradient
                    if run_name != ctrl_name:
                        gradPc = (control_vals[domain]['p'][pt]-control_vals[domain]['p']['sac'])/distance
                        udiff = usolano-control_vals[domain]['u']
                        vdiff = vsolano-control_vals[domain]['v']
                        spdiff = speed-control_vals[domain]['speed']
                        gradPdiff = gradP-gradPc

                        # for each of the p locations, scatter diffs of u, v, speed vs each level of pressure
                        fig, ax = plt.subplots(nrows=gradP.shape[1], ncols=3)
                        for i_lev in xrange(gradPdiff.shape[1]):
                            ax[i_lev, 0].plot(gradPdiff[:, i_lev], udiff, '.')
                            ax[i_lev, 1].plot(gradPdiff[:, i_lev], vdiff, '.')
                            ax[i_lev, 2].plot(gradPdiff[:, i_lev], spdiff, '.')
                            ax[i_lev, 0].set_ylabel(str(f.variables['levels'][i_lev])+" m")
                        ax[0, 0].set_title('U diff')
                        ax[0, 1].set_title('V diff')
                        ax[0, 2].set_title('speed diff')
                        fig.suptitle(run_name+" minus ctrl\nx: press grad diff "+pt+", y: wind diff solano")
                        fig.set_size_inches(6, gradP.shape[1]*2)
                        fig.savefig(os.path.join(plot_dir, "pgrad_vs_wind_DIFF_"+pt+".png")) 

            plt.close('all')
            

