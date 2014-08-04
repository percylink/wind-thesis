"""
PAL 7/29/2014
driver script to run plotting of WRF output
"""

import os
os.environ['PYTHONPATH'] = '/scratch2/scratchdirs/plink/WRF:'+os.environ['PYTHONPATH']

from utils.plot_WRF import Plotter


if __name__ == "__main__":

    run_name = 'CA-0.08'
    plot_interval = 3
    file_in = 'wrfout_d01_2009-07-01_00:00:00'
    reg_diff = 'reg'
    variables = ['smois']

    print run_name, file_in, reg_diff, 'interval', plot_interval, variables

    p = Plotter(file_in=file_in, plot_interval=plot_interval, run_name=run_name, reg_diff=reg_diff)
    p.run(variables)