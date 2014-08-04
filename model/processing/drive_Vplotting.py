"""
PAL 7/29/2014
driver script to run plotting of WRF output - vertical slices
"""

from utils.plot_uw_vertslice import VerticalPlotter


if __name__ == "__main__":

    run_name = 'CA-0.08'
    plot_interval = 3
    file_in = 'wrfout_diff_CTRL'
    reg_diff = 'diff'
    
    print run_name, file_in, reg_diff, 'interval', plot_interval

    p = VerticalPlotter(file_in=file_in, plot_interval=plot_interval, run_name=run_name, reg_diff=reg_diff)
    p.run()