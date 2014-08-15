'''
PAL 5/21/2013 mod 8/2/2013 mod 7/29/2014
plot vertical--E-W slice of u-wind and w-wind
'''

import datetime
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import os
from scipy.interpolate import griddata
from scipy.io import netcdf as nc
import shutil


class VerticalPlotter:

    def __init__(self, file_in=None, plot_interval=1, run_name=None, reg_diff='reg', file2=None, zoom=False, 
            lat_locs=None, lon_locs=[], plot_type='quiver', root_dir=None, domain='d01'):
 
        # instance parameters
        self.plot_interval = plot_interval
        if run_name is None:
            raise Exception('run_name is a required parameter')
        self.run_name = run_name
        if (reg_diff != 'reg') and (reg_diff != 'diff'):
            raise Exception('reg_diff must be "reg" or "diff"')
        self.reg_diff = reg_diff  # whether netcdf file is regular model output, or diff of 2 runs
        if root_dir is None:
            self.root_dir = os.path.join(os.environ['SCRATCH'], 'WRF', 'output', self.run_name)
        else:
            self.root_dir = os.path.join(root_dir, self.run_name)
        self.plot_dir = os.path.join(self.root_dir, 'plots', domain, 'vertical_uw')
        self.zoom = zoom
        self.plot_type = plot_type
        self.w_multiplier = 20.
        if run_name[:2] == 'DK':
            self.utc_minus_hr = -1
        elif run_name[:2] == 'CA':
            self.utc_minus_hr = 8

        # open netcdf file
        if file_in is None:
            raise Exception('file_in is a required parameter')
        self.f = nc.netcdf_file(os.path.join(self.root_dir, file_in), 'r')
        if self.reg_diff == 'diff':
            self.f_abs = nc.netcdf_file(os.path.join(self.root_dir, file2), 'r')  # file with absolute magnitudes for diff case

        # list of string times
        self.times = [''.join(self.f.variables['Times'][ii,:]) \
                      for ii in xrange(self.f.variables['Times'].shape[0])]

        # list of approx lats and lons (not exact because not cartesian projection)
        self.approxlat = self.f.variables['XLAT'][0,:,self.f.dimensions['west_east']/2]
        self.approxlon = self.f.variables['XLONG'][0,self.f.dimensions['south_north']/2,:]
        if lat_locs is None:
            self.latproflocs = [36.5,38.,39.,39.7]
        else:
            self.latproflocs = lat_locs
        self.lonproflocs = lon_locs

        if self.reg_diff == 'reg':
            self.vmin = {'u': -7, 'w': -0.2}
            self.vmax = {'u': 7, 'w': 0.2}
        elif self.reg_diff == 'diff':
            self.vmin = {'u': -3, 'w': -0.1}
            self.vmax = {'u': 3, 'w': 0.1}

        self.pbottom = 1050
        if self.zoom:
            self.ptop = 700
        else:
            self.ptop = 300

    def get_hour(self, time=None):

        hour = int(time[11:13])-self.utc_minus_hr
        if hour<0: 
            hour = hour+24
        return hour

    def get_data(self, wind='u', ii=None):

        if self.slice_direction == 'xz':
            xzslice = self.f.variables[wind.upper()][ii, :, self.ixlat, :]
            if wind == 'u':
                sliceav = np.mean([xzslice[:, 1:], xzslice[:, :-1]], axis=0)  # interpolate staggered x
            elif wind == 'w':
                sliceav = np.mean([xzslice[1:, :], xzslice[:-1, :]], axis=0)  # interpolate staggered z
        elif self.slice_direction == 'yz':
            yzslice = self.f.variables[wind.upper()][ii, :, :, self.ixlon]
            if wind == 'v':
                sliceav = np.mean([yzslice[:, 1:], yzslice[:, :-1]], axis=0)  # interpolate staggered y
            elif wind == 'w':
                sliceav = np.mean([yzslice[1:, :], yzslice[:-1, :]], axis=0)  # interpolate staggered z

        return sliceav

    def plot_sub(self, wind='u', fig=None, ax=None, ii=None, pslice=None, psurf=None):

        hour = self.get_hour(self.times[ii])
        sliceav = self.get_data(wind=wind, ii=ii)

        if self.slice_direction == 'xz':
            h = ax.scatter(self.longrid, pslice, s=20, c=sliceav,\
                vmin=self.vmin[wind], vmax=self.vmax[wind], cmap=plt.cm.jet, edgecolors='none')
            ax.set_title(wind.upper()+' profile along '+str(self.lat)+'N, local time '+str(hour)+':00 hr')
            ax.set_xlabel('longitude')
            ax.fill_between(self.lonprof, psurf, self.pbottom, color='k', alpha=0.4)
            ax.set_xlim(np.min(self.lonprof), np.max(self.lonprof))
        elif self.slice_direction == 'yz':
            h = ax.scatter(self.latgrid, pslice, s=20, c=sliceav,\
                vmin=self.vmin[wind], vmax=self.vmax[wind], cmap=plt.cm.jet, edgecolors='none')
            ax.set_title(wind.upper()+' profile along '+str(self.lon)+'E, local time '+str(hour)+':00 hr')
            ax.set_xlabel('latitude')
            ax.fill_between(self.latprof, psurf, self.pbottom, color='k', alpha=0.4)
            ax.set_xlim(np.min(self.latprof), np.max(self.latprof))
        
        ax.set_ylabel('pressure, hPa')
        ax.set_ylim(self.pbottom, self.ptop)
        fig.colorbar(h, ax=ax, orientation='horizontal')
        
    def plot_quiver(self, fig=None, ax=None, ii=None, pslice=None, psurf=None):

        hour = self.get_hour(self.times[ii])
        wav = self.get_data(wind='w', ii=ii)

        if self.slice_direction == 'xz':
            uav = self.get_data(wind='u', ii=ii)
            speed = np.sqrt(uav**2 + wav**2)
            Q = ax.quiver(self.longrid, pslice, uav, wav*self.w_multiplier, speed)
            ax.set_title('wind profile along '+str(self.lat)+'N, local time '+str(hour)+':00 hr\n'+ \
                         'w magnified by '+str(self.w_multiplier))
            ax.set_xlabel('longitude')
            ax.fill_between(self.lonprof, psurf, self.pbottom, color='k', alpha=0.4)
            ax.set_xlim(np.min(self.lonprof), np.max(self.lonprof))
            ax.plot(self.lonprof[self.landmask.astype(bool)], 
                    np.zeros(sum(self.landmask))+self.pbottom, 'yo')
            ax.plot(self.lonprof[~self.landmask.astype(bool)], 
                    np.zeros(len(self.landmask)-sum(self.landmask))+self.pbottom, 'bo')
        elif self.slice_direction == 'yz':
            vav = self.get_data(wind='v', ii=ii)
            speed = np.sqrt(vav**2 + wav**2)
            Q = ax.quiver(self.latgrid, pslice, vav, wav*self.w_multiplier, speed)
            ax.set_title('wind profile along '+str(self.lon)+'E, local time '+str(hour)+':00 hr\n'+ \
                         'w magnified by '+str(self.w_multiplier))
            ax.set_xlabel('latitude')
            ax.fill_between(self.latprof, psurf, self.pbottom, color='k', alpha=0.4)
            ax.set_xlim(np.min(self.latprof), np.max(self.latprof))
            ax.plot(self.latprof[self.landmask.astype(bool)], 
                    np.zeros(sum(self.landmask))+self.pbottom, 'yo')
            ax.plot(self.latprof[~self.landmask.astype(bool)], 
                    np.zeros(len(self.landmask)-sum(self.landmask))+self.pbottom, 'bo')
        
        qk = plt.quiverkey(Q, 0.9, 0.95, 2, r'$2 \frac{m}{s}$',
                        labelpos='E',
                        coordinates='figure',
                        fontproperties={'weight': 'bold'})

        ax.set_ylabel('pressure, hPa')
        ax.set_ylim(self.ptop, self.pbottom)
        ax.invert_yaxis()

        #plt.show()
        #1/0
        
    def plot_master(self, ii=None):

        print 'plot type:', self.plot_type

        # set up figure
        if self.plot_type == 'scatter':
            fig, ax = plt.subplots(nrows=2, ncols=1)
        elif self.plot_type == 'quiver':
            fig, ax = plt.subplots(nrows=1, ncols=1)
    
        if self.reg_diff == 'reg':
            fpres = self.f
        elif self.reg_diff == 'diff':
            fpres = self.f_abs

        if self.slice_direction == 'xz':
            pslice = (fpres.variables['PB'][ii, :, self.ixlat, :] + \
                        fpres.variables['P'][ii, :, self.ixlat, :])/100. # convert to hPa
            psurf = fpres.variables['PSFC'][ii, self.ixlat, :]/100.
        elif self.slice_direction == 'yz':
            pslice = (fpres.variables['PB'][ii, :, :, self.ixlon] + \
                        fpres.variables['P'][ii, :, :, self.ixlon])/100. # convert to hPa
            psurf = fpres.variables['PSFC'][ii, :, self.ixlon]/100.

        if self.plot_type == 'scatter':
            self.plot_sub(wind='u', fig=fig, ax=ax[0], ii=ii, pslice=pslice, psurf=psurf)
            self.plot_sub(wind='w', fig=fig, ax=ax[1], ii=ii, pslice=pslice, psurf=psurf)
        elif self.plot_type == 'quiver':
            self.plot_quiver(fig=fig, ax=ax, ii=ii, pslice=pslice, psurf=psurf)

        if self.plot_type == 'scatter':
            fig.set_size_inches(8, 11)
        else:
            fig.set_size_inches(11, 8)
        if self.slice_direction == 'xz':
            fig.suptitle(self.run_name+' '+str(self.lat)+'N '+self.times[ii])
            filename = 'UW_vertslice_'+str(self.lat)+'N_'+self.times[ii]+'.png'
        elif self.slice_direction == 'yz':
            fig.suptitle(self.run_name+' '+str(self.lon)+'E '+self.times[ii])
            filename = 'VW_vertslice_'+str(self.lon)+'E_'+self.times[ii]+'.png'
        if self.zoom:
            filename = 'zoom_'+filename
        if self.reg_diff == 'diff':
            filename = 'diffCTRL_'+filename
        fig.savefig(os.path.join(self.plot_dir, filename))#, transparent=True)
        plt.close(fig)
        #1/0

    def run(self):

        if not os.path.exists(self.plot_dir):
            os.makedirs(self.plot_dir)

        if self.reg_diff == 'reg':
            f_land = self.f
        elif self.reg_diff == 'diff':
            f_land = self.f_abs

        for lat in self.latproflocs:

            # get index of lat profile
            self.slice_direction = 'xz'
            self.ixlat = np.argmin(abs(self.approxlat-lat))
            self.ixlon = None

            # get topography profile
            self.lat = lat
            self.lon = None
            self.lonprof = self.f.variables['XLONG'][0, self.ixlat, :]
            self.longrid = np.tile(self.lonprof, (self.f.dimensions['bottom_top'], 1))
            self.landmask = f_land.variables['LANDMASK'][0, self.ixlat, :]

            for ii in range(len(self.times))[0::self.plot_interval]:
                
                print lat, self.times[ii]
                self.plot_master(ii=ii)

        for lon in self.lonproflocs:

            # get index of lat profile
            self.slice_direction = 'yz'
            self.ixlon = np.argmin(abs(self.approxlon-lon))
            self.ixlat = None

            # get topography profile
            self.lon = lon
            self.lat = None
            self.latprof = self.f.variables['XLAT'][0, :, self.ixlon]
            self.latgrid = np.tile(self.latprof, (self.f.dimensions['bottom_top'], 1))
            self.landmask = f_land.variables['LANDMASK'][0, :, self.ixlon]

            for ii in range(len(self.times))[0::self.plot_interval]:
                
                print lon, self.times[ii]
                self.plot_master(ii=ii)


if __name__ == "__main__":

    p = VerticalPlotter(file_in='wrfout_d01_2009-07-01_00', plot_interval=3, run_name='CA-0.24', reg_diff='reg')
    p.run()
           
