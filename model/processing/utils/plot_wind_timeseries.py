"""
PAL 9/8/2014
plot time series of wind at solano for each test
"""

from scipy.io import netcdf as nc
import matplotlib.pyplot as plt
import numpy as np
import os
import datetime


#coords_solano = {'d01': (49, 51), 'd02': (91, 96), 'd03': (127, 143)}

def get_box_index(f, lat, lon):
    dif_lat = f.variables['XLAT'][:, :]-lat
    dif_lon = f.variables['XLONG'][:, :]-lon
    distance = dif_lat**2 + dif_lon**2
    ixlat, ixlon = np.unravel_index(np.argmin(distance), distance.shape)
    return ixlat, ixlon

def get_winds(f, domain, coords_solano):
    usolano = f.variables['U'][:, :2, coords_solano[domain][1], coords_solano[domain][0]]
    vsolano = f.variables['V'][:, :2, coords_solano[domain][1], coords_solano[domain][0]]
    speed = (usolano**2 + vsolano**2)**0.5  # calculate speed
    return usolano, vsolano, speed

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
    return np.array(times)


if __name__ == "__main__":

    ctrl_name = "CA-0.1"
    #ctrl_name = "CA-0.2"
    #ctrl_name = "CA-0.25"

    test_names = ["CA-wetCR", "CA-wetCV", "CA-wetSN"]
    colors = {ctrl_name: 'b', test_names[0]: 'g', test_names[1]: 'r', test_names[2]: 'c'}
    test_id = "wet_regions"
    
    #test_names = ["CA-CV0.05", "CA-CV0.15", "CA-CV0.2", "CA-CV0.25", "CA-CV0.3", "CA-CV0.35"]
    #colors = {ctrl_name: 'b', test_names[0]: 'g', test_names[1]: 'r', test_names[2]: 'c',
    #          test_names[3]: 'm', test_names[4]: 'y', test_names[5]: 'k'}
    #test_id = "vary_CV"
    
    #test_names = ["CA-CV0.05w", "CA-CV0.1w", "CA-CV0.15w", "CA-CV0.25w", "CA-CV0.3w", "CA-CV0.35w"]
    #colors = {ctrl_name: 'b', test_names[0]: 'g', test_names[1]: 'r', test_names[2]: 'c',
    #          test_names[3]: 'm', test_names[4]: 'y', test_names[5]: 'k'}
    #test_id = "vary_CV_0.2"
    
    #test_names = ["CA-dryCR", "CA-dryCV", "CA-drySN"]
    #colors = {ctrl_name: 'b', test_names[0]: 'g', test_names[1]: 'r', test_names[2]: 'c'}
    #test_id = "dry_regions"
    
    domains = ["d01", "d02"]
    root_dir = "/scratch2/scratchdirs/plink/WRF/output"

    coords_solano = {}

    plot_dir = os.path.join(root_dir, "timeseries_wind")
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    for domain in domains:

        print domain
        fig1, ax1 = plt.subplots(nrows=2, ncols=1)
        fig2, ax2 = plt.subplots(nrows=2, ncols=1)

        for run_name in [ctrl_name]+test_names:

            print run_name
            plot_dir2 = os.path.join(root_dir, run_name, "plots", domain, "wind_tseries")
            if not os.path.exists(plot_dir2):
                os.makedirs(plot_dir2)
    
            fi = nc.netcdf_file(os.path.join(root_dir, run_name, "wrfout_interp_"+domain+".nc"))
            #fr = nc.netcdf_file(os.path.join(root_dir, run_name, "wrfout_"+domain+"_2009-07-01_00:00:00"))

            if run_name == ctrl_name:
                ixlat, ixlon = get_box_index(fi, 38.166, -121.817)
                coords_solano[domain] = (ixlon, ixlat)

            # get u, v at solano
            usolano, vsolano, speed = get_winds(fi, domain, coords_solano)
            if run_name == ctrl_name:
                speed_c = speed

            # get times
            times = convert_times(fi)

            ax1[0].plot(times, speed[:, 0], color=colors[run_name], label=run_name)
            if run_name != ctrl_name:
                ax1[1].plot(times, speed[:, 0]-speed_c[:, 0], color=colors[run_name], label=run_name)

            ax2[0].plot(times, speed[:, 1], color=colors[run_name], label=run_name)
            if run_name != ctrl_name:
                ax2[1].plot(times, speed[:, 1]-speed_c[:, 1], color=colors[run_name], label=run_name)

            # plot the lower two levels next to each other
            fig3, ax3 = plt.subplots(nrows=3, ncols=1)
            ax3[0].plot(times, usolano)
            ax3[1].plot(times, vsolano)
            ax3[2].plot(times, speed)
            ax3[0].set_ylabel('u')
            ax3[1].set_ylabel('v')
            ax3[2].set_ylabel('speed')
            ax3[0].set_title(run_name+" "+domain)
            fig3.savefig(os.path.join(plot_dir2, "wind_tseries_"+run_name+"_"+domain+".png"))
            plt.close(fig3)

        ax1[0].legend(prop={"size": 8})
        ax1[0].set_title("level 0")
        ax1[0].set_ylabel("wind speed")
        ax1[1].set_ylabel("speed minus ctrl")
        fig1.set_size_inches(15, 8)
        fig1.savefig(os.path.join(plot_dir, "solano_wind_"+test_id+"_"+domain+"_level0.png"))

        ax2[0].legend(prop={"size": 8})
        ax2[0].set_title("level 1")
        ax2[0].set_ylabel("wind speed")
        ax2[1].set_ylabel("speed minus ctrl")
        fig2.set_size_inches(15, 8)
        fig2.savefig(os.path.join(plot_dir, "solano_wind_"+test_id+"_"+domain+"_level1.png"))

        plt.close("all")
        


