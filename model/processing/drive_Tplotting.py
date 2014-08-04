"""
PAL 7/31/2014
driver script to run plotting of WRF time series output
"""

from utils.plot_timeseries import TimePlotter


if __name__ == "__main__":

    files = ['KS-CTRL/wrfout_d01_2010-05-27_12:00:00',
             'KS-left_wet/wrfout_d01_2010-05-27_12:00:00',
             'KS-right_wet/wrfout_d01_2010-05-27_12:00:00',
             'KS-all_wet/wrfout_d01_2010-05-27_12:00:00',
             'KS-all_dry/wrfout_d01_2010-05-27_12:00:00']

    p = TimePlotter(files_in=files, lat=38.5, lon=-98., nlevs=3)
    p.run(['wind', 'bkgd'])  #, 'PH', 'T_q'])
