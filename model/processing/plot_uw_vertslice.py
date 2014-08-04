'''
PAL 5/21/2013 mod 8/2/2013 mod 7/29/2014
plot vertical--E-W slice of u-wind and w-wind
'''

import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import os
from scipy.interpolate import griddata
from scipy.io import netcdf as nc
import shutil


class VerticalPlotter:

    def __init__(self, file_in=None, plot_interval=1, run_name=None, reg_diff='reg', file2=None, zoom=False, lat_locs=None):
 
        # instance parameters
        self.plot_interval = plot_interval
        if run_name is None:
            raise Exception('run_name is a required parameter')
        self.run_name = run_name
        if (reg_diff != 'reg') and (reg_diff != 'diff'):
            raise Exception('reg_diff must be "reg" or "diff"')
        self.reg_diff = reg_diff  # whether netcdf file is regular model output, or diff of 2 runs
        self.root_dir = os.path.join(os.environ['SCRATCH'], 'WRF', 'output', self.run_name)
        self.plot_dir = os.path.join(self.root_dir, 'plots', 'vertical_uw')
        self.zoom = zoom

        # open netcdf file
        if file_in is None:
            raise Exception('file_in is a required parameter')
        self.f = nc.netcdf_file(os.path.join(self.root_dir, file_in), 'r')
        if self.reg_diff == 'diff':
            self.f_abs = nc.netcdf_file(os.path.join(self.root_dir, file2), 'r')  # file with absolute magnitudes for diff case

        # list of string times
        self.times = [''.join(self.f.variables['Times'][ii,:]) \
                      for ii in xrange(self.f.variables['Times'].shape[0])]

        # list of approx latitudes (not exact because not cartesian projection)
        self.approxlat = self.f.variables['XLAT'][0,:,self.f.dimensions['west_east']/2]
        if lat_locs is None:
            self.latproflocs = [36.5,38.,39.,39.7]
        else:
            self.latproflocs = lat_locs

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

        hour = int(time[11:13])-8
        if hour<0: 
            hour = hour+24
        return hour

    def get_data(self, wind='u', ii=None):

        xzslice = self.f.variables[wind.upper()][ii, :, self.ixlat, :]
        if wind == 'u':
            sliceav = np.mean([xzslice[:, 1:], xzslice[:, :-1]], axis=0)  # interpolate staggered x
        elif wind == 'w':
            sliceav = np.mean([xzslice[1:, :], xzslice[:-1, :]], axis=0)  # interpolate staggered z

        return sliceav

    def plot_sub(self, wind='u', fig=None, ax=None, ii=None, pslice=None, psurf=None):

        hour = self.get_hour(self.times[ii])
        sliceav = self.get_data(wind=wind, ii=ii)

        h = ax.scatter(self.longrid, pslice, s=20, c=sliceav,\
            vmin=self.vmin[wind], vmax=self.vmax[wind], cmap=plt.cm.jet, edgecolors='none')
        
        ax.set_title(wind.upper()+' profile along '+str(self.lat)+'N, local time '+str(hour)+':00 hr')
        ax.set_xlabel('longitude')
        ax.set_ylabel('pressure, hPa')
        #if self.reg_diff == 'reg':
        ax.fill_between(self.lonprof, psurf, self.pbottom, color='k', alpha=0.4)
        ax.set_xlim(np.min(self.lonprof), np.max(self.lonprof))
        ax.set_ylim(self.pbottom, self.ptop)
        fig.colorbar(h, ax=ax, orientation='horizontal')
        
    def plot_master(self, ii=None):

        # set up figure
        fig, ax = plt.subplots(nrows=2, ncols=1)
    
        if self.reg_diff == 'reg':
            pslice = (self.f.variables['PB'][ii, :, self.ixlat, :] + \
                        self.f.variables['P'][ii, :, self.ixlat, :])/100. # convert to hPa
            psurf = self.f.variables['PSFC'][ii, self.ixlat, :]/100.
        elif self.reg_diff == 'diff':
            pslice = (self.f_abs.variables['PB'][ii, :, self.ixlat, :] + \
                        self.f_abs.variables['P'][ii, :, self.ixlat, :])/100. # convert to hPa
            psurf = self.f_abs.variables['PSFC'][ii, self.ixlat, :]/100.

        self.plot_sub(wind='u', fig=fig, ax=ax[0], ii=ii, pslice=pslice, psurf=psurf)
        self.plot_sub(wind='w', fig=fig, ax=ax[1], ii=ii, pslice=pslice, psurf=psurf)

        fig.suptitle(self.run_name+' '+str(self.lat)+'N '+self.times[ii])
        fig.set_size_inches(8, 11)
        filename = 'UW_vertslice_'+str(self.lat)+'N_'+self.times[ii]+'.png'
        if self.zoom:
            filename = 'zoom_'+filename
        if self.reg_diff == 'diff':
            filename = 'diffCTRL_'+filename
        fig.savefig(os.path.join(self.plot_dir, filename), transparent=True)
        plt.close(fig)

    def run(self):

        if os.path.exists(self.plot_dir):
            shutil.rmtree(self.plot_dir)
        os.makedirs(self.plot_dir)

        for lat in self.latproflocs:

             # get index of lat profile
            self.ixlat = np.argmin(abs(self.approxlat-lat))

            # get topography profile
            self.lat = lat
            self.topoprof = self.f.variables['HGT'][0, self.ixlat, :]
            self.lonprof = self.f.variables['XLONG'][0, self.ixlat, :]
            self.longrid = np.tile(self.lonprof, (self.f.dimensions['bottom_top'], 1))

            for ii in range(len(self.times))[0::self.plot_interval]:
                
                print lat, self.times[ii]
                self.plot_master(ii=ii)


if __name__ == "__main__":

    p = VerticalPlotter(file_in='wrfout_d01_2009-07-01_00', plot_interval=3, run_name='CA-0.24', reg_diff='reg')
    p.run()
           

# # open netcdf file
# outfolder = 'plots/vertical_uw'
# if not os.path.exists(outfolder):
#     os.mkdir(outfolder)
# fname2 = 'wrfout_d01_2009-07-01_00:00:00'
# f2 = nc.netcdf_file(fname2,'r')

# # list of string times
# times = [''.join(f2.variables['Times'][ii,:]) for ii in xrange(f2.variables['Times'].shape[0])]
# print times

# # list of approx latitudes (not exact because not cartesian projection)
# approxlat = f2.variables['XLAT'][0,:,f2.dimensions['west_east']/2]
# latproflocs = [36.5,38.,39.,39.7]

# eta = f2.variables['ZNU'][0,:]
# vmin = -7
# vmax = 7
# pbottom = 1050
# ptop = 300

# # loop through times
# for ii,tt in enumerate(times):
    
#     print tt
#     hour = int(tt[11:13])-8
#     if hour<0: 
#         hour = hour+24

#     for lat in latproflocs:

#         # get index of lat profile
#         ixlat = np.argmin(abs(approxlat-lat))

#         # get u x-z slice
#         Uxzslice = f2.variables['U'][ii,:,ixlat,:].copy()
#         Wxzslice = f2.variables['W'][ii,:,ixlat,:].copy()
#         Vxzslice = f2.variables['V'][ii,:,ixlat,:].copy()
#         Usliceav = np.mean([Uxzslice[:,1:],Uxzslice[:,:-1]],axis=0)
#         Wsliceav = np.mean([Wxzslice[1:,:],Wxzslice[:-1,:]],axis=0)
#         pslice = (f2.variables['PB'][ii,:,ixlat,:] + f2.variables['P'][ii,:,ixlat,:])/100. # convert to hPa
#         psurf = f2.variables['PSFC'][ii,ixlat,:]/100.

#         # get topography profile
#         topoprof = f2.variables['HGT'][0,ixlat,:]
#         lonprof = f2.variables['XLONG'][0,ixlat,:]
#         longrid = np.tile(lonprof,(Uxzslice.shape[0],1))

#         # uslice vs pressure slice, scatter
#         fig,ax = plt.subplots(nrows=2,ncols=1)
        
#         # plot
#         h = ax[0].scatter(longrid,pslice,s=20,c=Usliceav,\
#             vmin=vmin,vmax=vmax,cmap=plt.cm.jet,edgecolors='none')
#         #ax[0].scatter(longrid,pslice,s=20,facecolors='none')
#         ax[0].set_xlabel('longitude')
#         ax[0].set_ylabel('pressure, hPa')
#         ax[0].fill_between(lonprof,psurf,1050,color='k',alpha=0.4)
#         #ax[0].invert_yaxis()
#         ax[0].set_xlim(np.min(lonprof),np.max(lonprof))
#         ax[0].set_ylim(pbottom, ptop)
#         ax[0].set_title('U profile along '+str(lat)+'N, local time '+str(hour)+':00 hr')
#         fig.colorbar(h,ax=ax[0],orientation='horizontal')

#         h = ax[1].scatter(longrid,pslice,s=20,c=Wsliceav,\
#             vmin=-0.2,vmax=0.2,cmap=plt.cm.jet,edgecolors='none')
#         #ax[0].scatter(longrid,pslice,s=20,facecolors='none')
#         ax[1].set_xlabel('longitude')
#         ax[1].set_ylabel('pressure, hPa')
#         ax[1].fill_between(lonprof,psurf,1050,color='k',alpha=0.4)
#         #ax[0].invert_yaxis()
#         ax[1].set_xlim(np.min(lonprof),np.max(lonprof))
#         ax[1].set_ylim(pbottom, ptop)
#         ax[1].set_title('W profile along '+str(lat)+'N, local time '+str(hour)+':00 hr')
#         fig.colorbar(h,ax=ax[1],orientation='horizontal')
        
#         #plt.show()
#         #1/0    
        
#         # save
#         fig.suptitle(outfolder+'/'+fname2+' '+tt)
#         fig.set_size_inches(8,11)
#         fig.savefig(os.path.join(outfolder, 'UW_vertslice_'+str(lat)+'N_'+fname2+'_'+str(ii)+'.png'))
#         plt.close(fig)
# #   if ii==4:
# #       print ii
# #       1/0
# f2.close()
