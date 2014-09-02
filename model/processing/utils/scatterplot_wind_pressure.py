"""
PAL 9/1/2014
Scatter plots of WRF output
"""

import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.io import netcdf as nc
from scipy.stats import linregress


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

def plot_p_vs_wind(p, u, v, speed, f, run_name, pt, plot_dir, diff=False, grad=False, resid=False):

    fig, ax = plt.subplots(nrows=p.shape[1], ncols=3)
    if p.shape[1] > 1:
        for i_lev in xrange(p.shape[1]):
            ax[i_lev, 0].plot(p[:, i_lev], u, '.')
            ax[i_lev, 1].plot(p[:, i_lev], v, '.')
            ax[i_lev, 2].plot(p[:, i_lev], speed, '.')
            ax[i_lev, 0].set_ylabel(str(f.variables['levels'][i_lev])+" m")
    else:
        ax[0].plot(p[:, 0], u, '.')
        ax[1].plot(p[:, 0], v, '.')
        ax[2].plot(p[:, 0], speed, '.')
        #ax[0].set_ylabel(str(f.variables['levels'][i_lev])+" m")

    fig.set_size_inches(6, p.shape[1]*2)

    if grad:
        varname = "press grad"
    else:
        varname = "pressure"

    if resid:
        varname += " resid"

    if diff:
        if p.shape[1] > 1:
            ax[0, 0].set_title('U diff')
            ax[0, 1].set_title('V diff')
            ax[0, 2].set_title('speed diff')
        else:
            ax[0].set_title('U diff')
            ax[1].set_title('V diff')
            ax[2].set_title('speed diff')
        fig.suptitle(run_name+" minus ctrl\n x: "+varname+" diff "+pt+", y: wind diff solano")
        fig.savefig(os.path.join(plot_dir, varname.replace(" ", "")+"_vs_wind_DIFF_"+pt+".png")) 
    else:
        if p.shape[1] > 1:
            ax[0, 0].set_title('U')
            ax[0, 1].set_title('V')
            ax[0, 2].set_title('speed')
        else:
            ax[0].set_title('U')
            ax[1].set_title('V')
            ax[2].set_title('speed')
        fig.suptitle(run_name+" x: "+varname+" "+pt+", y: wind solano")
        fig.savefig(os.path.join(plot_dir, varname.replace(" ", "")+"_vs_wind_"+pt+".png"))
        plt.close('all')

def plot_regression(p, u, v, speed, slope, intercept, f, run_name, pt, plot_dir, diff=False, grad=False):

    prange = np.array([min(p[np.isfinite(p)]), max(p[np.isfinite(p)])])

    fig, ax = plt.subplots(nrows=1, ncols=3)
    
    ax[0].plot(p, u, '.')
    ax[0].plot(prange, slope['u'][0]*prange+intercept['u'][0], 'b')
    ax[0].plot(prange, slope['u'][1]*prange+intercept['u'][1], 'g')
    ax[0].set_xlabel('slope0: '+str(round(slope['u'][0]))+' int0: '+str(round(intercept['u'][0],2))+
                     '\nslope1: '+str(round(slope['u'][1]))+' int1: '+str(round(intercept['u'][1],2)))

    ax[1].plot(p, v, '.')
    ax[1].plot(prange, slope['v'][0]*prange+intercept['v'][0], 'b')
    ax[1].plot(prange, slope['v'][1]*prange+intercept['v'][1], 'g')
    ax[1].set_xlabel('slope0: '+str(round(slope['v'][0]))+' int0: '+str(round(intercept['v'][0],2))+
                     '\nslope1: '+str(round(slope['v'][1]))+' int1: '+str(round(intercept['v'][1],2)))

    ax[2].plot(p, speed, '.')
    ax[2].plot(prange, slope['sp'][0]*prange+intercept['sp'][0], 'b')
    ax[2].plot(prange, slope['sp'][1]*prange+intercept['sp'][1], 'g')
    ax[2].set_xlabel('slope0: '+str(round(slope['sp'][0]))+' int0: '+str(round(intercept['sp'][0],2))+
                     '\nslope1: '+str(round(slope['sp'][1]))+' int1: '+str(round(intercept['sp'][1],2)))

    fig.set_size_inches(12, 4)

    if grad:
        varname = "press grad 250m"
    else:
        varname = "pressure 250m"

    if diff:
        ax[0].set_ylabel('U diff')
        ax[1].set_ylabel('V diff')
        ax[2].set_ylabel('speed diff')
        fig.suptitle(run_name+" minus ctrl\n x: "+varname+" diff "+pt+", y: wind diff solano")
        fig.savefig(os.path.join(plot_dir, "regress_"+varname.replace(" ", "")+"_vs_wind_DIFF_"+pt+".png")) 
    else:
        ax[0].set_title('U')
        ax[1].set_title('V')
        ax[2].set_title('speed')
        fig.suptitle(run_name+" x: "+varname+" "+pt+", y: wind solano")
        fig.savefig(os.path.join(plot_dir, "regress_"+varname.replace(" ", "")+"_vs_wind_"+pt+".png"))
    plt.close('all')   

def calc_regression(pressure, wind):
    nlevs = wind.shape[1]
    slope = np.zeros(nlevs)
    intercept = np.zeros(nlevs)
    r = np.zeros(nlevs)
    p = np.zeros(nlevs)
    residuals = np.zeros_like(wind)
    for i in xrange(nlevs):
        slope[i], intercept[i], r[i], p[i], dum = linregress(pressure, wind[:, i])
        wind_predicted = slope[i]*pressure + intercept[i]
        residuals[:, i] = wind[:, i]-wind_predicted
    return slope, intercept, r, p, residuals


if __name__ == "__main__":

    ctrl_name = "CA-0.24"
    test_names = ["CA-dryCR", "CA-dryCV", "CA-drySN"]
    domains = ["d01", "d02", "d03"]
    #test_names = ["CA-dryCR"]
    #domains = ["d01"]  #["d02", "d03"]
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

                # if this is a diff case, scatter diffs in u, v, speed vs diffs in pressure
                if run_name != ctrl_name:
                    pdiff = pdict[pt]-control_vals[domain]['p'][pt]
                    udiff = usolano-control_vals[domain]['u']
                    vdiff = vsolano-control_vals[domain]['v']
                    spdiff = speed-control_vals[domain]['speed']

                    # for each of the p locations, scatter diffs of u, v, speed vs each level of pressure
                    plot_p_vs_wind(pdiff, udiff, vdiff, spdiff, f, run_name, pt, plot_dir, diff=True)

                else:
                    control_vals[domain]['p'] = pdict
                    control_vals[domain]['u'] = usolano
                    control_vals[domain]['v'] = vsolano
                    control_vals[domain]['speed'] = speed

            # now that we know we have the values for sac, loop through
            # pts again to calculate gradients from pt to sac.
            for pt in grad_pts:

                # for each of the gradients, scatter u, v, speed vs each level of gradient
                if (pt != 'sac') and (pt in pdict):
                    
                    distance = calc_distance(pt, domain, f)
                    gradP = (pdict[pt] - pdict['sac']) / distance

                    plot_p_vs_wind(gradP, usolano, vsolano, speed, f, run_name, pt, plot_dir, grad=True)

                    # if this is a diff case, scatter u, v, speed vs diffs in gradient
                    if run_name != ctrl_name:
                        gradPc = (control_vals[domain]['p'][pt]-control_vals[domain]['p']['sac'])/distance
                        udiff = usolano-control_vals[domain]['u']
                        vdiff = vsolano-control_vals[domain]['v']
                        spdiff = speed-control_vals[domain]['speed']
                        gradPdiff = gradP-gradPc

                        # for each of the p locations, scatter diffs of u, v, speed vs each level of pressure
                        plot_p_vs_wind(gradPdiff, udiff, vdiff, spdiff, f, run_name, pt, plot_dir, grad=True, diff=True)

            # calculate regression lines of bay-sac press grad @ 250 m vs winds
            distance = calc_distance('bay', domain, f)
            mask250 = f.variables['levels'][:] == 250
            gradP_baysac = (pdict['bay'][:, mask250] - pdict['sac'][:, mask250]) / distance
            gradP_baysac = gradP_baysac.squeeze()

            slope = {}
            intercept = {}
            r = {}
            p = {}
            resid = {}
            slope['u'], intercept['u'], r['u'], p['u'], resid['u'] = calc_regression(gradP_baysac, usolano)
            slope['v'], intercept['v'], r['v'], p['v'], resid['v'] = calc_regression(gradP_baysac, vsolano)
            slope['sp'], intercept['sp'], r['sp'], p['sp'], resid['sp'] = calc_regression(gradP_baysac, speed)

            # plot bay-sac press grad vs winds
            plot_regression(gradP_baysac, usolano, vsolano, speed, slope, intercept, f, 
                            run_name, "bay_sac", plot_dir, grad=True)
            
            # plot other pressure grads vs residual winds
            for pt in grad_pts:
                if (pt not in ['sac', 'bay']) and (pt in pdict):
                    distance = calc_distance(pt, domain, f)
                    gradP = (pdict[pt][:, mask250] - pdict['sac'][:, mask250]) / distance
                    plot_p_vs_wind(gradP, resid['u'], resid['v'], resid['sp'], f, run_name, 
                                    pt, plot_dir, grad=True, resid=True)
