'''
PAL 5/8/2013
make WRF plots for CS267 final project
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


class Plotter:

    def __init__(self, file_in=None, file2=None, plot_interval=1, run_name=None, reg_diff='reg', domain='d01'):

        # instance parameters
        self.plot_interval = plot_interval
        self.domain = domain
        if run_name is None:
            raise Exception('run_name is a required parameter')
        self.run_name = run_name
        if (reg_diff != 'reg') and (reg_diff != 'diff'):
            raise Exception('reg_diff must be "reg" or "diff"')
        self.reg_diff = reg_diff  # whether netcdf file is regular model output, or diff of 2 runs
        self.root_dir = os.path.join(os.environ['SCRATCH'], 'WRF', 'output', self.run_name)
        if run_name[:2] == 'CA':
            self.hr_diff = 8
        elif run_name[:2] == 'KS':
            self.hr_diff = 7

        print "file:", self.file_in, ", plot_interval:", self.plot_interval, ", run_name:", self.run_name, ", reg_diff:", reg_diff, ", domain:", self.domain

        # physical parameters
        self.latent_heat_evap = 2.27e6  # J/kg

        # open netcdf file
        if file_in is None:
            raise Exception('file_in is a required parameter')
        self.f = nc.netcdf_file(os.path.join(self.root_dir, file_in), 'r')
        if self.reg_diff == 'diff':
            self.f_abs = nc.netcdf_file(os.path.join(self.root_dir, file2), 'r')  # file with absolute magnitudes for diff case

        # get min/max lat and lon
        self.minlon = np.floor(np.min(self.f.variables['XLONG'][0,:,:]/10.))*10.
        self.maxlon = np.ceil(np.max(self.f.variables['XLONG'][0,:,:]/10.))*10.
        self.minlat = np.floor(np.min(self.f.variables['XLAT'][0,:,:]/10.))*10.
        self.maxlat = np.ceil(np.max(self.f.variables['XLAT'][0,:,:]/10.))*10.

        # list of string times
        self.times = [''.join(self.f.variables['Times'][ii,:]) \
                      for ii in xrange(self.f.variables['Times'].shape[0])]

        # subsample masks for winds
        self.subx = np.arange(0,getattr(self.f,'WEST-EAST_GRID_DIMENSION')-1,3)
        self.suby = np.arange(0,getattr(self.f,'SOUTH-NORTH_GRID_DIMENSION')-1,3)

        self.get_pressure_levels()

    def get_pressure_levels(self):

        if self.reg_diff == 'reg':
            self.pressures = self.f.variables['PB'][0,:,0,0] # background pressure profile
        elif self.reg_diff == 'diff':
            self.pressures = self.f_abs.variables['PB'][0,:,0,0] # background pressure profile
        self.maskP = {}
        self.P = {}
        
        # get index for ~950 hPa, ~850 hPa, ~500 hPa on unstaggered grid
        dif = abs(self.pressures-95000.)
        self.maskP['950'] = dif.argmin() # index of minimum difference
        self.P['950'] = self.pressures[self.maskP['950']]/100. # actual pressure near 950, in hPa
        dif = abs(self.pressures-85000.)
        self.maskP['850'] = dif.argmin() # index of minimum difference
        self.P['850'] = self.pressures[self.maskP['850']]/100. # actual pressure near 850, in hPa
        dif = abs(self.pressures-50000.)
        self.maskP['500'] = dif.argmin() # index of minimum difference
        self.P['500'] = self.pressures[self.maskP['500']]/100. # actual pressure near 500, in hPa

        # get index for ~950 hPa, ~850 hPa, ~500 hPa on vertically staggered grid
        if self.reg_diff == 'reg':
            self.eta = self.f.variables['ZNW'][0,:]
        elif self.reg_diff == 'diff':
            self.eta = self.f_abs.variables['ZNW'][0,:]
        dif = abs(self.eta-0.95)
        self.maskP['950s'] = dif.argmin() # index of minimum difference
        self.P['950s'] = self.eta[self.maskP['950s']]*1000. # actual pressure near 950, in hPa
        dif = abs(self.eta-0.85)
        self.maskP['850s'] = dif.argmin() # index of minimum difference
        self.P['850s'] = self.eta[self.maskP['850s']]*1000. # actual pressure near 850, in hPa
        dif = abs(self.eta-0.5)
        self.maskP['500s'] = dif.argmin() # index of minimum difference
        self.P['500s'] = self.eta[self.maskP['500s']]*1000. # actual pressure near 500, in hPa

    def setup_map(self, subplot_index=None, fig=None):

        f = self.f
        a = fig.add_subplot(subplot_index)
        m = Basemap(width=f.DX*1.2*getattr(f,'WEST-EAST_GRID_DIMENSION'),\
            height=f.DY*1.2*getattr(f,'SOUTH-NORTH_GRID_DIMENSION'),resolution='l',\
            projection='lcc',lat_1=f.TRUELAT1,lat_2=f.TRUELAT2,lat_0=f.CEN_LAT,lon_0=f.CEN_LON)
        return a, m

    def decorate_map(self, m=None, a=None, h=None, fig=None):

        m.drawparallels(np.arange(self.minlat, self.maxlat, 5))
        m.drawmeridians(np.arange(self.minlon, self.maxlon, 5))
        m.drawcoastlines()
        m.drawstates()
        fig.colorbar(h, ax=a)

    def make_filename(self, ii=None, variable=None):

        if variable == 'smois':
            fname = variable+'.png'
        else:
            fname = variable+'_'+self.times[ii]+'.png'

        if self.reg_diff == 'reg':
            file_out = os.path.join(self.plot_dir, variable, fname)
        elif self.reg_diff == 'diff':
            file_out = os.path.join(self.plot_dir, variable, 'diffCTRL_'+fname)

        return file_out

    def plot_PH(self, fig=None, subplot_index=None, pressure_level=None, contour_int=None, ii=None):

        a, m = self.setup_map(subplot_index=subplot_index, fig=fig)
        f = self.f

        if type(pressure_level) is not str:
            pressure_level = str(int(pressure_level))
        
        x,y = m(f.variables['XLONG'][0, :, :],f.variables['XLAT'][0, :, :])
        tmp = f.variables['PH'][ii, self.maskP[pressure_level], :, :] + f.variables['PHB'][ii, self.maskP[pressure_level], :, :]
        if contour_int is not None:
            h = m.contourf(x, y, tmp, contour_int)
        else:
            h = m.contourf(x, y, tmp)
        a.set_ylabel(str(round(self.P[pressure_level]))+' hPa')
        a.set_title('geopotential')
        
        self.decorate_map(m=m, a=a, h=h, fig=fig)

    def plot_wind(self, fig=None, subplot_index=None, pressure_level=None, contour_int=None, ii=None):

        a, m = self.setup_map(subplot_index=subplot_index, fig=fig)
        f = self.f

        if type(pressure_level) is not str:
            pressure_level = str(int(pressure_level))
        
        x,y = m(f.variables['XLONG'][0,:,:],f.variables['XLAT'][0,:,:])
        xu,yu = m(f.variables['XLONG_U'][0,:,:],f.variables['XLAT_U'][0,:,:])
        xv,yv = m(f.variables['XLONG_V'][0,:,:],f.variables['XLAT_V'][0,:,:])
        ui = np.mean([f.variables['U'][ii, self.maskP[pressure_level], :, :-1],
            f.variables['U'][ii, self.maskP[pressure_level], :, 1:]], axis=0)
        vi = np.mean([f.variables['V'][ii, self.maskP[pressure_level], :-1, :],
            f.variables['V'][ii, self.maskP[pressure_level], 1:, :]], axis=0)
        #wspeed = (ui**2+vi**2)**0.5
        #h = m.contourf(x, y, wspeed, contour_int)
        limits = [np.min(ui)]+range(-8, 9)+[np.max(ui)]
        h = m.contourf(x, y, ui, limits)
        ur,vr = m.rotate_vector(ui, vi, f.variables['XLONG'][0,:,:], f.variables['XLAT'][0,:,:])
        m.quiver(x[self.suby, :][:, self.subx],y[self.suby, :][:, self.subx], 
            ur[self.suby, :][:, self.subx], vr[self.suby, :][:, self.subx])
        a.set_title('wind, m/s '+str(round(self.P[pressure_level]))+' hPa')

        self.decorate_map(m=m, a=a, h=h, fig=fig)

    def plot_T(self, fig=None, subplot_index=None, pressure_level=None, contour_int=None, ii=None):

        a, m = self.setup_map(subplot_index=subplot_index, fig=fig)
        f = self.f

        if type(pressure_level) is not str:
            pressure_level = str(int(pressure_level))
        
        x,y = m(f.variables['XLONG'][0,:,:],f.variables['XLAT'][0,:,:])
        tmp = f.variables['T'][ii, self.maskP[pressure_level], :, :] + f.variables['T00'][ii]
        h = m.contourf(x,y,tmp)
        a.set_title('pot. temp. at '+str(round(self.P[pressure_level]))+' hPa')

        self.decorate_map(m=m, a=a, h=h, fig=fig)

    def plot_q(self, fig=None, subplot_index=None, pressure_level=None, contour_int=None, ii=None):
    
        a, m = self.setup_map(subplot_index=subplot_index, fig=fig)
        f = self.f

        if type(pressure_level) is not str:
            pressure_level = str(int(pressure_level))
        
        x,y = m(f.variables['XLONG'][0,:,:],f.variables['XLAT'][0,:,:])
        tmp = f.variables['QVAPOR'][ii, self.maskP[pressure_level], :, :]
        h = m.contourf(x,y,tmp)
        a.set_title('q at '+str(round(self.P[pressure_level]))+' hPa')

        self.decorate_map(m=m, a=a, h=h, fig=fig)

    def plot_latent(self, fig=None, subplot_index=None, contour_int=None, ii=None):
    
        a, m = self.setup_map(subplot_index=subplot_index, fig=fig)
        f = self.f

        x,y = m(f.variables['XLONG'][0,:,:],f.variables['XLAT'][0,:,:])
        tmp = f.variables['QFX'][ii, :, :] * self.latent_heat_evap
        if contour_int is not None:
            h = m.contourf(x,y,tmp,contour_int)
        else:
            h = m.contourf(x,y,tmp)
        a.set_title('LH, W/m2')

        self.decorate_map(m=m, a=a, h=h, fig=fig)

    def plot_sensible(self, fig=None, subplot_index=None, contour_int=None, ii=None):
    
        a, m = self.setup_map(subplot_index=subplot_index, fig=fig)
        f = self.f

        x,y = m(f.variables['XLONG'][0,:,:],f.variables['XLAT'][0,:,:])
        tmp = f.variables['HFX'][ii, :, :]
        if contour_int is not None:
            h = m.contourf(x,y,tmp,contour_int)
        else:
            h = m.contourf(x,y,tmp)
        a.set_title('SH, W/m2')

        self.decorate_map(m=m, a=a, h=h, fig=fig)

    def plot_smois_diff(self, fig=None, subplot_index=None, integrated=True, contour_int=None):
    
        a, m = self.setup_map(subplot_index=subplot_index, fig=fig)
        f = self.f

        x,y = m(f.variables['XLONG'][0,:,:],f.variables['XLAT'][0,:,:])
        smois_diff = f.variables['SMOIS'][-1,:,:,:]-f.variables['SMOIS'][0,:,:,:]

        if integrated:
            tmp = smois_diff.sum(axis=0)
        else:
            tmp = smois_diff[0, :, :]

        if contour_int is not None:
            h = m.contourf(x,y,tmp,contour_int)
        else:
            h = m.contourf(x,y,tmp)

        if integrated:
            a.set_title('d(SMOIS), column')
        else:
            a.set_title('d(SMOIS), surface')

        self.decorate_map(m=m, a=a, h=h, fig=fig)

    def plot_PH_winds(self, ii=None):

        fig = plt.figure()

        self.plot_PH(fig=fig, subplot_index=321, pressure_level='950', contour_int=20, ii=ii)
        self.plot_PH(fig=fig, subplot_index=323, pressure_level='850', contour_int=20, ii=ii)
        self.plot_PH(fig=fig, subplot_index=325, pressure_level='500', contour_int=20, ii=ii)

        self.plot_wind(fig=fig, subplot_index=322, pressure_level='950', ii=ii)
        self.plot_wind(fig=fig, subplot_index=324, pressure_level='850', ii=ii)
        self.plot_wind(fig=fig, subplot_index=326, pressure_level='500', ii=ii)

        fig.suptitle(self.times[ii])
        file_out = self.make_filename(ii=ii, variable='PH_winds')
        print file_out
        fig.savefig(file_out)
        plt.close(fig)

    def plot_wind_only(self, ii=None):

        fig = plt.figure()

        self.plot_wind(fig=fig, subplot_index=131, pressure_level='950', ii=ii)
        self.plot_wind(fig=fig, subplot_index=132, pressure_level='850', ii=ii)
        self.plot_wind(fig=fig, subplot_index=133, pressure_level='500', ii=ii)

        local_hour = int(self.times[ii][11:13])-self.hr_diff
        if local_hour < 0:
            local_hour += 24

        fig.set_size_inches(15,5)
        fig.suptitle('color = u-wind '+self.times[ii]+', local hour '+str(local_hour))
        file_out = self.make_filename(ii=ii, variable='winds')
        print file_out
        fig.savefig(file_out)
        plt.close(fig)      

    def plot_PH_only(self, ii=None):

        fig = plt.figure()

        self.plot_PH(fig=fig, subplot_index=131, pressure_level='950', ii=ii)
        self.plot_PH(fig=fig, subplot_index=132, pressure_level='850', ii=ii)
        self.plot_PH(fig=fig, subplot_index=133, pressure_level='500', ii=ii)

        fig.suptitle(self.times[ii])
        file_out = self.make_filename(ii=ii, variable='PH')
        print file_out
        fig.savefig(file_out)
        plt.close(fig)      

    def plot_T_q(self, ii=None):

        fig = plt.figure()

        self.plot_T(fig=fig, subplot_index=221, pressure_level='950', ii=ii)
        self.plot_T(fig=fig, subplot_index=223, pressure_level='850', ii=ii)
        
        self.plot_q(fig=fig, subplot_index=222, pressure_level='950', ii=ii)
        self.plot_q(fig=fig, subplot_index=224, pressure_level='850', ii=ii)
        
        fig.suptitle(self.times[ii])
        file_out = self.make_filename(ii=ii, variable='T_q')
        print file_out
        fig.savefig(file_out)
        plt.close(fig)

    def plot_sfc_flx(self, ii=None):

        fig = plt.figure()

        self.plot_latent(fig=fig, subplot_index=121, ii=ii)
        self.plot_sensible(fig=fig, subplot_index=122, ii=ii)
        
        fig.set_size_inches(15,7)
        fig.suptitle(self.times[ii])
        file_out = self.make_filename(ii=ii, variable='sfc_flx')
        print file_out
        fig.savefig(file_out)
        plt.close(fig)

    def plot_smois_diff_driver(self):

        fig = plt.figure()

        self.plot_smois_diff(fig=fig, subplot_index=121, integrated=False)
        self.plot_smois_diff(fig=fig, subplot_index=122, integrated=True)
        
        fig.set_size_inches(15,7)
        file_out = self.make_filename(variable='smois')
        print file_out
        fig.savefig(file_out)
        plt.close(fig)

    def plot_smois_init(self):

        fig = plt.figure()
        a, m = self.setup_map(subplot_index=111, fig=fig)
        f = self.f

        x,y = m(f.variables['XLONG'][0,:,:],f.variables['XLAT'][0,:,:])
        smois_init = f.variables['SMOIS'][0,0,:,:]  # surface

        contours = list(np.arange(0.08, 0.4, 0.04))
        if smois_init.min() < contours[0]:
            contours = [smois_init.min()]+contours
        if smois_init[smois_init < 1].max() > contours[-1]:  # don't use 1 as max
            contours = contours+[smois_init.max()]
        
        h = m.contourf(x, y, smois_init, contours)

        a.set_title('initial SMOIS, surface')
        self.decorate_map(m=m, a=a, h=h, fig=fig)

        #1/0
        if not os.path.exists(os.path.join(self.plot_dir, 'smois')):
            os.mkdir(os.path.join(self.plot_dir, 'smois'))
        fig.savefig(os.path.join(self.plot_dir, 'smois', 'initial_smois.png'))

    def run(self, plot_list=[]):

        if not isinstance(plot_list, list):
            plot_list = [plot_list]

        # make output directories and call plotting functions
        plot_dir = os.path.join(self.root_dir, 'plots', self.domain)
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        self.plot_dir = plot_dir

        if self.reg_diff == 'reg':
           self.plot_smois_init()

        for plot in plot_list:
            #if os.path.exists(os.path.join(plot_dir, plot)):
            #    shutil.rmtree(os.path.join(plot_dir, plot))
            #os.mkdir(os.path.join(plot_dir, plot))
            if not os.path.exists(os.path.join(plot_dir, plot)):
                os.mkdir(os.path.join(plot_dir, plot))

            if plot == 'smois':
                self.plot_smois_diff_driver()

            for ii in range(len(self.times))[0::self.plot_interval]:
                print self.times[ii]
                if plot == 'PH_winds':
                    self.plot_PH_winds(ii=ii)
                elif plot == 'T_q':
                    self.plot_T_q(ii=ii)
                elif plot == 'winds':
                    self.plot_wind_only(ii=ii)
                elif plot == 'PH':
                    self.plot_PH_only(ii=ii)
                elif plot == 'sfc_flx':
                    self.plot_sfc_flx(ii=ii)


if __name__ == "__main__":

    if len(sys.argv) >= 2:
        file_in = sys.argv[1]
    else:
        raise Exception("Usage: python plot_WRF.py <wrfout file name>")

    if len(sys.argv) >2:
        plot_interval = int(sys.argv[2])
    else:
        plot_interval = 1

    p = Plotter(file_in=file_in, plot_interval=plot_interval)
    p.run(['winds', 'PH', 'T_q'])



        
    #   """
    #   ######
    #   ## cloud and precip plots
    #   ######
    #   fig = plt.figure()
        
    #   # plot total precipitable water (integrate QCLOUD, QRAIN, QVAPOR) - all kg/kg
    #   a1 = fig.add_subplot(221)
    #   m = Basemap(width=f.DX*1.2*getattr(f,'WEST-EAST_GRID_DIMENSION'),\
    #       height=f.DY*1.2*getattr(f,'SOUTH-NORTH_GRID_DIMENSION'),resolution='l',\
    #       projection='lcc',lat_1=f.TRUELAT1,lat_2=f.TRUELAT2,lat_0=f.CEN_LAT,lon_0=f.CEN_LON)
    #   x,y = m(f.variables['XLONG'][0,:,:],f.variables['XLAT'][0,:,:])
    #   tmp = f.variables['QCLOUD'][ii,:,:,:]+f.variables['QRAIN'][ii,:,:,:]+f.variables['QVAPOR'][ii,:,:,:] # sum all forms of water
    #   tmpa = np.mean([tmp[1:,:,:],tmp[:-1,:,:]],axis=0) # average for layers between grid points
    #   ps = f.variables['PSFC'][ii,:,:]
    #   ptop = f.variables['P_TOP'][ii]
    #   eta = f.variables['ZNU'][ii,:]
    #   tpw = np.zeros_like(ps)
    #   for jj in xrange(np.shape(tmpa)[0]): # calculate sum(dp*qavg/g)
    #       dp = (ps-ptop)*(eta[jj]-eta[jj+1])
    #       dq = tmpa[jj,:,:]*dp/9.8
    #       tpw += dq
    #   h = m.contourf(x,y,tpw,20)
    #   m.drawparallels(np.arange(minlat,maxlat,5))
    #   m.drawmeridians(np.arange(minlon,maxlon,5))
    #   m.drawcoastlines()
    #   m.drawstates()
    #   fig.colorbar(h,ax=a1)
    #   a1.set_title('total precipitable water, kg/m2')
        
    #   # plot OLR
    #   a2 = fig.add_subplot(222)
    #   m = Basemap(width=f.DX*1.2*getattr(f,'WEST-EAST_GRID_DIMENSION'),\
    #       height=f.DY*1.2*getattr(f,'SOUTH-NORTH_GRID_DIMENSION'),resolution='l',\
    #       projection='lcc',lat_1=f.TRUELAT1,lat_2=f.TRUELAT2,lat_0=f.CEN_LAT,lon_0=f.CEN_LON)
    #   x,y = m(f.variables['XLONG'][0,:,:],f.variables['XLAT'][0,:,:])
    #   tmp = f.variables['OLR'][ii,:,:]
    #   h = m.contourf(x,y,tmp,20)
    #   m.drawparallels(np.arange(minlat,maxlat,5))
    #   m.drawmeridians(np.arange(minlon,maxlon,5))
    #   m.drawcoastlines()
    #   m.drawstates()
    #   fig.colorbar(h,ax=a2)
    #   a2.set_title('OLR, W/m2')
        
    #   # plot integrated QCLOUD below 500 hPa
    #   a3 = fig.add_subplot(223)
    #   m = Basemap(width=f.DX*1.2*getattr(f,'WEST-EAST_GRID_DIMENSION'),\
    #       height=f.DY*1.2*getattr(f,'SOUTH-NORTH_GRID_DIMENSION'),resolution='l',\
    #       projection='lcc',lat_1=f.TRUELAT1,lat_2=f.TRUELAT2,lat_0=f.CEN_LAT,lon_0=f.CEN_LON)
    #   x,y = m(f.variables['XLONG'][0,:,:],f.variables['XLAT'][0,:,:])
    #   tmp = f.variables['QCLOUD'][ii,:,:,:]
    #   tmpa = np.mean([tmp[1:,:,:],tmp[:-1,:,:]],axis=0) # average for layers between grid points
    #   ps = f.variables['PSFC'][ii,:,:]
    #   ptop = f.variables['P_TOP'][ii]
    #   eta = f.variables['ZNU'][ii,:]
    #   icld = np.zeros_like(ps)
    #   for jj in xrange(np.shape(tmpa)[0]): # calculate sum(dp*qavg/g)
    #       if 1000.*eta[jj+1]>=500.:
    #           dp = (ps-ptop)*(eta[jj]-eta[jj+1])
    #           dq = tmpa[jj,:,:]*dp/9.8
    #           icld += dq
    #   h = m.contourf(x,y,icld,20)
    #   m.drawparallels(np.arange(minlat,maxlat,5))
    #   m.drawmeridians(np.arange(minlon,maxlon,5))
    #   m.drawcoastlines()
    #   m.drawstates()
    #   fig.colorbar(h,ax=a3)
    #   a3.set_title('cloud water below 500 hPa, kg/m2')
        
    #   # plot integrated QCLOUD above 500 hPa
    #   a4 = fig.add_subplot(224)
    #   m = Basemap(width=f.DX*1.2*getattr(f,'WEST-EAST_GRID_DIMENSION'),\
    #       height=f.DY*1.2*getattr(f,'SOUTH-NORTH_GRID_DIMENSION'),resolution='l',\
    #       projection='lcc',lat_1=f.TRUELAT1,lat_2=f.TRUELAT2,lat_0=f.CEN_LAT,lon_0=f.CEN_LON)
    #   x,y = m(f.variables['XLONG'][0,:,:],f.variables['XLAT'][0,:,:])
    #   tmp = f.variables['QCLOUD'][ii,:,:,:]
    #   tmpa = np.mean([tmp[1:,:,:],tmp[:-1,:,:]],axis=0) # average for layers between grid points
    #   ps = f.variables['PSFC'][ii,:,:]
    #   ptop = f.variables['P_TOP'][ii]
    #   eta = f.variables['ZNU'][ii,:]
    #   icld = np.zeros_like(ps)
    #   for jj in xrange(np.shape(tmpa)[0]): # calculate sum(dp*qavg/g)
    #       if 1000.*eta[jj+1]<500.:
    #           dp = (ps-ptop)*(eta[jj]-eta[jj+1])
    #           dq = tmpa[jj,:,:]*dp/9.8
    #           icld += dq
    #   h = m.contourf(x,y,icld,20)
    #   m.drawparallels(np.arange(minlat,maxlat,5))
    #   m.drawmeridians(np.arange(minlon,maxlon,5))
    #   m.drawcoastlines()
    #   m.drawstates()
    #   fig.colorbar(h,ax=a4)
    #   a4.set_title('cloud water above 500 hPa, kg/m2')
        
    #   fig.suptitle(tt)
    #   pcloud.savefig(fig)
    #   plt.close(fig)
        
    #   ######
    #   ## rain, OLR, PBL, w plots
    #   ######
    #   fig = plt.figure()
        
    #   # plot rain since last timestep (RAINC+RAINNC) --> mm
    #   a1 = fig.add_subplot(221)
    #   m = Basemap(width=f.DX*1.2*getattr(f,'WEST-EAST_GRID_DIMENSION'),\
    #       height=f.DY*1.2*getattr(f,'SOUTH-NORTH_GRID_DIMENSION'),resolution='l',\
    #       projection='lcc',lat_1=f.TRUELAT1,lat_2=f.TRUELAT2,lat_0=f.CEN_LAT,lon_0=f.CEN_LON)
    #   x,y = m(f.variables['XLONG'][0,:,:],f.variables['XLAT'][0,:,:])
    #   if ii==0:
    #       tmp = f.variables['RAINNC'][ii,:,:]+f.variables['RAINC'][ii,:,:]
    #   else:
    #       tmp = (f.variables['RAINNC'][ii,:,:]-f.variables['RAINNC'][ii-1,:,:])+\
    #           (f.variables['RAINC'][ii,:,:]-f.variables['RAINC'][ii-1,:,:])
    #   h = m.contourf(x,y,tmp,20)
    #   m.drawparallels(np.arange(minlat,maxlat,5))
    #   m.drawmeridians(np.arange(minlon,maxlon,5))
    #   m.drawcoastlines()
    #   m.drawstates()
    #   fig.colorbar(h,ax=a1)
    #   a1.set_title('total precip, mm since last output')
        
    #   # PBL height
    #   a2 = fig.add_subplot(222)
    #   m = Basemap(width=f.DX*1.2*getattr(f,'WEST-EAST_GRID_DIMENSION'),\
    #       height=f.DY*1.2*getattr(f,'SOUTH-NORTH_GRID_DIMENSION'),resolution='l',\
    #       projection='lcc',lat_1=f.TRUELAT1,lat_2=f.TRUELAT2,lat_0=f.CEN_LAT,lon_0=f.CEN_LON)
    #   x,y = m(f.variables['XLONG'][0,:,:],f.variables['XLAT'][0,:,:])
    #   tmp = f.variables['PBLH'][ii,:,:]/1000.
    #   h = m.contourf(x,y,tmp,20)
    #   m.drawparallels(np.arange(minlat,maxlat,5))
    #   m.drawmeridians(np.arange(minlon,maxlon,5))
    #   m.drawcoastlines()
    #   m.drawstates()
    #   fig.colorbar(h,ax=a2)
    #   a2.set_title('PBL height, km')
        
    #   # w velocity near 850 hPa
    #   a3 = fig.add_subplot(223)
    #   m = Basemap(width=f.DX*1.2*getattr(f,'WEST-EAST_GRID_DIMENSION'),\
    #       height=f.DY*1.2*getattr(f,'SOUTH-NORTH_GRID_DIMENSION'),resolution='l',\
    #       projection='lcc',lat_1=f.TRUELAT1,lat_2=f.TRUELAT2,lat_0=f.CEN_LAT,lon_0=f.CEN_LON)
    #   x,y = m(f.variables['XLONG'][0,:,:],f.variables['XLAT'][0,:,:])
    #   tmp = f.variables['W'][ii,mask850s,:,:]
    #   h = m.contourf(x,y,tmp,20)
    #   m.drawparallels(np.arange(minlat,maxlat,5))
    #   m.drawmeridians(np.arange(minlon,maxlon,5))
    #   m.drawcoastlines()
    #   m.drawstates()
    #   fig.colorbar(h,ax=a3)
    #   a3.set_title('w velocity at '+str(round(p850s))+' hPa')
        
    #   # w velocity near 500 hPa
    #   a4 = fig.add_subplot(224)
    #   m = Basemap(width=f.DX*1.2*getattr(f,'WEST-EAST_GRID_DIMENSION'),\
    #       height=f.DY*1.2*getattr(f,'SOUTH-NORTH_GRID_DIMENSION'),resolution='l',\
    #       projection='lcc',lat_1=f.TRUELAT1,lat_2=f.TRUELAT2,lat_0=f.CEN_LAT,lon_0=f.CEN_LON)
    #   x,y = m(f.variables['XLONG'][0,:,:],f.variables['XLAT'][0,:,:])
    #   tmp = f.variables['W'][ii,mask500s,:,:]
    #   h = m.contourf(x,y,tmp,20)
    #   m.drawparallels(np.arange(minlat,maxlat,5))
    #   m.drawmeridians(np.arange(minlon,maxlon,5))
    #   m.drawcoastlines()
    #   m.drawstates()
    #   fig.colorbar(h,ax=a4)
    #   a4.set_title('w velocity at '+str(round(p500s))+' hPa')
        
    #   fig.suptitle(tt)
    #   pmisc.savefig(fig)
    #   plt.close(fig)
    #   """

