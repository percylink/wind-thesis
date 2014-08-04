'''
PAL 7/30/2014
make plots of time series at the grid box closest to user-defined lat, lon
'''

import cPickle
import datetime
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.basemap import Basemap
import numpy as np
import os
from scipy.interpolate import interp2d
from scipy.io import netcdf as nc
import shutil
import sys


class TimePlotter:

    def __init__(self, files_in=None, lat=37.75, lon=-121.65, nlevs=3):

        # instance parameters
        self.root_dir = os.path.join(os.environ['SCRATCH'], 'WRF', 'output')
        self.plot_dir = os.path.join(self.root_dir, 'timeseries_plots')
        if not os.path.exists(self.plot_dir):
            os.mkdir(self.plot_dir)
        self.lat = lat
        self.lon = lon
        self.nlevs = nlevs

        # open netcdf file
        if files_in is None:
            raise Exception('files_in is a required parameter')
        self.files_in = files_in
        self.f = {}
        for filename in files_in:
            self.f[filename] = nc.netcdf_file(os.path.join(self.root_dir, filename), 'r')

        # convert string times to datetime
        self.convert_times()

        # mask for grid box of interest
        self.get_box_index()

    def convert_times(self):

        self.times = {}
        for filename in self.files_in:
            times = []
            for ii in xrange(self.f[filename].variables['Times'].shape[0]):
                time_arr = self.f[filename].variables['Times'][ii, :]
                year = int(''.join(time_arr[:4]))
                month = int(''.join(time_arr[5:7]))
                day = int(''.join(time_arr[8:10]))
                hour = int(''.join(time_arr[11:13]))
                minute = int(''.join(time_arr[14:16]))
                second = int(''.join(time_arr[17]))
                times.append(datetime.datetime(year, month, day, hour, minute, second))
            self.times[filename] = np.array(times)

    def get_box_index(self):

        self.ix_lats = {}
        self.ix_lons = {}

        for filename in self.files_in:
            dif_lat = self.f[filename].variables['XLAT'][0, :, :]-self.lat
            dif_lon = self.f[filename].variables['XLONG'][0, :, :]-self.lon
            distance = dif_lat**2 + dif_lon**2
            ixlat, ixlon = np.unravel_index(np.argmin(distance), distance.shape)
            self.ix_lats[filename] = ixlat
            self.ix_lons[filename] = ixlon

    def make_domain_map(self):

        f = self.f[self.files_in[0]]
        ix_lat = self.ix_lats[self.files_in[0]]
        ix_lon = self.ix_lons[self.files_in[0]]

        fig, ax = plt.subplots()
        m = Basemap(width=f.DX*1.2*getattr(f,'WEST-EAST_GRID_DIMENSION'),\
            height=f.DY*1.2*getattr(f,'SOUTH-NORTH_GRID_DIMENSION'),resolution='l',\
            projection='lcc',lat_1=f.TRUELAT1,lat_2=f.TRUELAT2,lat_0=f.CEN_LAT,lon_0=f.CEN_LON)

        x, y = m(f.variables['XLONG'][0,:,:],f.variables['XLAT'][0,:,:])
        m.plot(x, y, '+k')
        m.plot(x[ix_lat, ix_lon], y[ix_lat, ix_lon], 'or')

        minlon = np.floor(np.min(f.variables['XLONG'][0,:,:]/10.))*10.
        maxlon = np.ceil(np.max(f.variables['XLONG'][0,:,:]/10.))*10.
        minlat = np.floor(np.min(f.variables['XLAT'][0,:,:]/10.))*10.
        maxlat = np.ceil(np.max(f.variables['XLAT'][0,:,:]/10.))*10.

        m.drawparallels(np.arange(minlat, maxlat, 5))
        m.drawmeridians(np.arange(minlon, maxlon, 5))
        m.drawcoastlines()
        m.drawstates()

        run_name = self.files_in[0].split('/')[0]
        fig.savefig(os.path.join(self.plot_dir, 'domain_'+run_name+'.png'))
        plt.close()

    def invert_wind_dir(self, wind_dir=None):
        # change wind dir from "direction to which vector points" to "direction from which vector points"
        wind_dir -= 180
        wind_dir[wind_dir < 0] = wind_dir[wind_dir < 0] + 360.
        return wind_dir

    def calc_wind_speed_dir(self, filename=None, ix=None):

        f = self.f[filename]
        ix_lat = self.ix_lats[filename]
        ix_lon = self.ix_lons[filename]

        if ix is None:
            u_left = f.variables['U'][:, :self.nlevs, ix_lat, ix_lon]
            u_right = f.variables['U'][:, :self.nlevs, ix_lat, ix_lon+1]
            v_bot = f.variables['V'][:, :self.nlevs, ix_lat, ix_lon]
            v_top = f.variables['V'][:, :self.nlevs, ix_lat+1, ix_lon]
        else:
            u_left = f.variables['U'][:, ix, ix_lat, ix_lon]
            u_right = f.variables['U'][:, ix, ix_lat, ix_lon+1]
            v_bot = f.variables['V'][:, ix, ix_lat, ix_lon]
            v_top = f.variables['V'][:, ix, ix_lat+1, ix_lon]
     
        u_av = np.mean([u_left, u_right], axis=0)
        v_av = np.mean([v_bot, v_top], axis=0)

        wind_speed = (u_av**2 + v_av**2)**0.5
        wind_dir = np.arctan2(u_av, v_av) * 180./np.pi
        wind_dir[wind_dir < 0] = wind_dir[wind_dir < 0] + 360.
        wind_dir = self.invert_wind_dir(wind_dir=wind_dir)

        return wind_speed, wind_dir

    def plot_wind(self):

        for ilev in xrange(self.nlevs):  # new plot for each level

            print "model level", ilev
            fig, ax = plt.subplots(nrows=3, ncols=1)

            for fname in self.files_in:

                wind_speed, wind_dir = self.calc_wind_speed_dir(filename=fname)
                run_name = fname.split('/')[0]

                # wind speed
                ax[0].plot(self.times[fname], wind_speed[:, ilev], label=run_name)

                # wind direction
                ax[1].plot(self.times[fname], wind_dir[:, ilev], label=run_name)

                # ustar
                ustar = self.f[fname].variables['UST'][:, self.ix_lats[fname], self.ix_lons[fname]]
                ax[2].plot(self.times[fname], ustar, label=run_name)

            ax[0].set_ylabel('speed, m/s')
            ax[1].set_ylabel('direction, deg')
            ax[2].set_ylabel('u*, m/s')

            ax[2].legend()

            fig.suptitle('winds, model level '+str(ilev))
            fig.savefig(os.path.join(self.plot_dir, 'winds_lev'+str(ilev)+'.png'))
            plt.close()

    def plot_bkgd(self):
        
        f = self.f[self.files_in[0]]
        eta = f.variables['ZNU'][0, :]
        ix_steering = np.argmin(abs(eta-0.5))  # steering level, ~500 hPa
        wind_speed, wind_dir = self.calc_wind_speed_dir(filename=self.files_in[0], ix=ix_steering)

        fig, ax = plt.subplots(nrows=2, ncols=1)
        ax[0].plot(self.times[self.files_in[0]], wind_speed)
        ax[1].plot(self.times[self.files_in[0]], wind_dir)
        ax[0].set_ylabel('speed, m/s')
        ax[1].set_ylabel('direction, deg')

        fig.suptitle('winds, steering level (eta~0.5)')
        fig.savefig(os.path.join(self.plot_dir, 'winds_steering.png'), transparent=True)
        plt.close()

    def plot_slices(self):
        pass

    def run(self, plot_list=[]):

        if not isinstance(plot_list, list):
            plot_list = [plot_list]

        self.make_domain_map()

        for vname in plot_list:

            if vname == 'wind':
                self.plot_wind()
            elif vname == 'bkgd':
                self.plot_bkgd()
            elif vname == 'slices':
                self.plot_slices()
            else:
                raise Exception('unrecognized plot variable: '+vname)


if __name__ == "__main__":

    files = ['CA-CTRL/wrfout_d01_2009-07-01_00:00:00', 
             'CA-0.08/wrfout_d01_2009-07-01_00:00:00', 
             'CA-0.3/wrfout_d01_2009-07-01_00:00:00']

    p = TimePlotter(files_in=files, lat=37.75, lon=-121.65, nlevs=3)
    p.run(['wind'])  #, 'PH', 'T_q'])
