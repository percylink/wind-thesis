"""
PAL 7/29/2014
driver script to run plotting of WRF output - vertical slices
"""

from utils.plot_uw_vertslice import VerticalPlotter


if __name__ == "__main__":

    run_name = 'DK-Aug11-mod'
    plot_interval = 3
    file_in = 'wrfout_d01_2013-08-11_00:00:00'
    reg_diff = 'reg'
    
    files = ['wrfout_d01_2013-08-11_00:00:00', 'wrfout_d02_2013-08-11_00:00:00', 'wrfout_d03_2013-08-11_00:00:00']
    domains = ['d01', 'd02', 'd03']

    for i in xrange(len(files)):

        print run_name, files[i], reg_diff, 'interval', plot_interval, domains[i]

        p = VerticalPlotter(file_in=files[i], plot_interval=plot_interval, run_name=run_name, \
                reg_diff=reg_diff, lat_locs=[], lon_locs=[11.7], plot_type='quiver', domain=domains[i])
        p.run()

    #root_dir='/Users/percy/Documents/Research/Projects/wind/wind-thesis/model/output'