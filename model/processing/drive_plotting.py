"""
PAL 7/29/2014
driver script to run plotting of WRF output
"""

import os
os.environ['PYTHONPATH'] = '/scratch2/scratchdirs/plink/WRF:'+os.environ['PYTHONPATH']

from utils.plot_WRF import Plotter


if __name__ == "__main__":

    run_name = 'DK-Aug11-mod'
    plot_interval = 3
    reg_diff = 'reg'
    variables = ['smois', 'winds', 'PH', 'T_q', 'sfc_flx']

    files = ['wrfout_d01_2013-08-11_00:00:00', 'wrfout_d02_2013-08-11_00:00:00', 'wrfout_d03_2013-08-11_00:00:00']
    domains = ['d01', 'd02', 'd03']

    for i in xrange(len(files)):

        #print run_name, file_in, reg_diff, 'interval', plot_interval, variables
        #file_in = 'wrfout_d01_2009-07-01_00:00:00'
    
        p = Plotter(file_in=files[i], plot_interval=plot_interval, run_name=run_name, reg_diff=reg_diff, domain=domains[i])
        p.run(variables)