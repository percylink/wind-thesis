'''
PAL 7/30/2014
make plots of time series at the grid box closest to user-defined lat, lon
'''

import cPickle
import datetime
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.dates import date2num
from mpl_toolkits.basemap import Basemap
import numpy as np
import os
from scipy.interpolate import interp2d
from scipy.io import netcdf as nc
import shutil
import sys


class TimeHeightPlotter:

    def __init__(self, file_in=None, lat=37.75, lon=-121.65, run_name=None, reg_diff='reg', \
                 file2=None, domain='d01', nlevs=5):

        # instance parameters
        self.run_name = run_name
        if (reg_diff != 'reg') and (reg_diff != 'diff'):
            raise Exception('reg_diff must be "reg" or "diff"')
        self.reg_diff = reg_diff  # whether netcdf file is regular model output, or diff of 2 runs
        if root_dir is None:
            self.root_dir = os.path.join(os.environ['SCRATCH'], 'WRF', 'output', self.run_name)
        else:
            self.root_dir = os.path.join(root_dir, self.run_name)
        self.plot_dir = os.path.join(self.root_dir, 'plots', domain)
        self.domain = domain
        if not os.path.exists(self.plot_dir):
            os.mkdir(self.plot_dir)
        self.lat = lat
        self.lon = lon
        self.nlevs = nlevs

        # open netcdf file
        if file_in is None:
            raise Exception('file_in is a required parameter')
        self.file_in = file_in
        self.f = nc.netcdf_file(os.path.join(self.root_dir, file_in), 'r')
        if self.reg_diff == 'diff':
            self.f_abs = nc.netcdf_file(os.path.join(self.root_dir, file2), 'r')

        # convert string times to datetime
        self.convert_times()

        # mask for grid box of interest
        self.get_box_index()

        # get heights of staggered grid at grid box of interest
        self.get_heights()

    def convert_times(self):

        times = []
        for ii in xrange(self.f.variables['Times'].shape[0]):
            time_arr = self.f.variables['Times'][ii, :]
            year = int(''.join(time_arr[:4]))
            month = int(''.join(time_arr[5:7]))
            day = int(''.join(time_arr[8:10]))
            hour = int(''.join(time_arr[11:13]))
            minute = int(''.join(time_arr[14:16]))
            second = int(''.join(time_arr[17]))
            times.append(datetime.datetime(year, month, day, hour, minute, second))
        self.times = np.array(times)

    def get_box_index(self):

        dif_lat = self.f.variables['XLAT'][0, :, :]-self.lat
        dif_lon = self.f.variables['XLONG'][0, :, :]-self.lon
        distance = dif_lat**2 + dif_lon**2
        ixlat, ixlon = np.unravel_index(np.argmin(distance), distance.shape)
        self.ix_lat = ixlat
        self.ix_lon = ixlon

    def get_heights(self):
        if self.reg_diff == 'reg':
            f = self.f
        elif self.reg_diff == 'diff':
            f = self.f_abs

        PHB = f.variables['PHB'][0, :self.nlevs+1, self.ix_lat, self.ix_lon]
        self.heights = PHB/9.81  # convert from geopotential to height by dividing by grav. accel.

    def invert_wind_dir(self, wind_dir=None):
        # change wind dir from "direction to which vector points" to "direction from which vector points"
        wind_dir -= 180
        wind_dir[wind_dir < 0] = wind_dir[wind_dir < 0] + 360.
        return wind_dir

    def calc_wind_speed_dir(self, ix=None):

        f = self.f
        ix_lat = self.ix_lat
        ix_lon = self.ix_lon

        if ix is None:
            u_east = f.variables['U'][:, :self.nlevs, ix_lat, ix_lon]
            u_west = f.variables['U'][:, :self.nlevs, ix_lat, ix_lon+1]
            v_south = f.variables['V'][:, :self.nlevs, ix_lat, ix_lon]
            v_north = f.variables['V'][:, :self.nlevs, ix_lat+1, ix_lon]
        else:
            u_east = f.variables['U'][:, ix, ix_lat, ix_lon]
            u_west = f.variables['U'][:, ix, ix_lat, ix_lon+1]
            v_south = f.variables['V'][:, ix, ix_lat, ix_lon]
            v_north = f.variables['V'][:, ix, ix_lat+1, ix_lon]
     
        u_av = np.mean([u_east, u_west], axis=0)
        v_av = np.mean([v_south, v_north], axis=0)

        wind_speed = (u_av**2 + v_av**2)**0.5
        wind_dir = np.arctan2(u_av, v_av) * 180./np.pi
        wind_dir[wind_dir < 0] = wind_dir[wind_dir < 0] + 360.
        wind_dir = self.invert_wind_dir(wind_dir=wind_dir)

        return wind_speed, wind_dir

    def plot_wind(self):

        # get time x vertical wind speed and direction
        wind_speed, wind_dir = self.calc_wind_speed_dir()

        # plot with imshow
        fig, ax = plt.subplots(nrows=2, ncols=1)
        h_0 = ax[0].imshow(wind_speed, interpolation='none', aspect='auto', 
            extent=[self.times[0], self.times[-1], self.heights[0], self.heights[-1]],
            origin='lower')
        h_1 = ax[1].imshow(wind_dir, interpolation='none', aspect='auto', 
            extent=[self.times[0], self.times[-1], self.heights[0], self.heights[-1]],
            origin='lower')
        plt.colorbar(h_0, ax=ax[0])
        plt.colorbar(h_1, ax=ax[1])

        # format tick labels
        ax[0].set_yticks(range(0, self.nlevs+1))
        ax[1].set_yticks(range(0, self.nlevs+1))
        ax[0].set_yticklabels(self.heights)
        ax[1].set_yticklabels(self.heights)
        ax[0].xaxis_date()
        ax[1].xaxis_date()
        fig.autofmt_xdate()

        # save plot
        ax[0].set_title('speed, m/s')
        ax[0].set_ylabel('height, m')
        ax[1].set_title('direction, deg')
        ax[1].set_ylabel('height, m')
        fig.suptitle('wind speed, domain '+self.domain+' lat '+str(self.lat)+' lon '+str(self.lon))
        fig.savefig(os.path.join(self.plot_dir, 'timeheight', 'wind_speed.png'))
        plt.close()

    def run(self, plot_list=[]):

        if not isinstance(plot_list, list):
            plot_list = [plot_list]

        if not os.path.exists(os.path.join(self.plot_dir, 'timeheight')):
            os.makedirs(os.path.join(self.plot_dir, 'timeheight'))

        for vname in plot_list:

            if vname == 'wind':
                self.plot_wind()
            elif vname == 'pot_temp':
                self.plot_pot_temp()
            else:
                raise Exception('unrecognized plot variable: '+vname)


if __name__ == "__main__":

    files = ['CA-CTRL/wrfout_d01_2009-07-01_00:00:00', 
             'CA-0.08/wrfout_d01_2009-07-01_00:00:00', 
             'CA-0.3/wrfout_d01_2009-07-01_00:00:00']

    p = TimePlotter(files_in=files, lat=37.75, lon=-121.65, nlevs=3)
    p.run(['wind'])  #, 'PH', 'T_q'])
