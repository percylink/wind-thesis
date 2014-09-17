'''
PAL 9/10/14
plot the domains
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


def get_box_index(f, lat, lon):
    dif_lat = f.variables['XLAT'][0, :, :]-lat
    dif_lon = f.variables['XLONG'][0, :, :]-lon
    distance = dif_lat**2 + dif_lon**2
    ixlat, ixlon = np.unravel_index(np.argmin(distance), distance.shape)
    return ixlat, ixlon

def get_smois_regions(fCR, fCV, fSN):
    smois_dry = 0.1
    smois_init_cr = fCR.variables['SMOIS'][0, 0, :, :]
    smois_init_cv = fCV.variables['SMOIS'][0, 0, :, :]
    smois_init_sn = fSN.variables['SMOIS'][0, 0, :, :]
    im = np.zeros_like(smois_init_cr) + np.nan
    im[smois_init_cr == 0.1] = 1
    im[smois_init_cv == 0.1] = 2
    im[smois_init_sn == 0.1] = 3
    return im

def decorate_map(f, m):
    minlon = np.floor(np.min(f.variables['XLONG'][0,:,:]/10.))*10.
    maxlon = np.ceil(np.max(f.variables['XLONG'][0,:,:]/10.))*10.
    minlat = np.floor(np.min(f.variables['XLAT'][0,:,:]/10.))*10.
    maxlat = np.ceil(np.max(f.variables['XLAT'][0,:,:]/10.))*10.
    m.drawparallels(np.arange(minlat, maxlat, 5), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(minlon, maxlon, 5), labels=[0, 0, 0, 1])
    m.drawcoastlines()
    m.drawstates()
    m.drawcoastlines()
    m.drawstates()

def plot_domain_box(f, m, domain_name):
    lons = f.variables['XLONG'][0,:,:]
    lats = f.variables['XLAT'][0,:,:]
    llcnr = m(lons[0, 0], lats[0, 0])
    lrcnr = m(lons[0, -1], lats[0, -1])
    urcnr = m(lons[-1, -1], lats[-1, -1])
    ulcnr = m(lons[-1, 0], lats[-1, 0])
    m.plot(*zip(llcnr, lrcnr), color='gray', lw=1.5)
    m.plot(*zip(lrcnr, urcnr), color='gray', lw=1.5)
    m.plot(*zip(urcnr, ulcnr), color='gray', lw=1.5)
    m.plot(*zip(ulcnr, llcnr), color='gray', lw=1.5)
    m.ax.annotate(domain_name, llcnr, xytext=(5, 5), textcoords='offset points', color='gray')

def plot_box(f, m, ix_llcnr, ix_urcnr, color='white'):
    # ix order: lon, lat
    lons = f.variables['XLONG'][0,:,:]
    lats = f.variables['XLAT'][0,:,:]
    llcnr = m(lons[ix_llcnr[1], ix_llcnr[0]], lats[ix_llcnr[1], ix_llcnr[0]])
    lrcnr = m(lons[ix_llcnr[1], ix_urcnr[0]], lats[ix_llcnr[1], ix_urcnr[0]])
    urcnr = m(lons[ix_urcnr[1], ix_urcnr[0]], lats[ix_urcnr[1], ix_urcnr[0]])
    ulcnr = m(lons[ix_urcnr[1], ix_llcnr[0]], lats[ix_urcnr[1], ix_llcnr[0]])
    m.plot(*zip(llcnr, lrcnr), color=color, lw=1.5)
    m.plot(*zip(lrcnr, urcnr), color=color, lw=1.5)
    m.plot(*zip(urcnr, ulcnr), color=color, lw=1.5)
    m.plot(*zip(ulcnr, llcnr), color=color, lw=1.5)


if __name__=="__main__":

    root_dir = "/scratch2/scratchdirs/plink/WRF/output"

    fCR1 = nc.netcdf_file(os.path.join(root_dir, "CA-dryCR", "wrfout_d01_2009-07-01_00:00:00"))
    fCR2 = nc.netcdf_file(os.path.join(root_dir, "CA-dryCR", "wrfout_d02_2009-07-01_00:00:00"))
    fCV1 = nc.netcdf_file(os.path.join(root_dir, "CA-dryCV", "wrfout_d01_2009-07-01_00:00:00"))
    fSN1 = nc.netcdf_file(os.path.join(root_dir, "CA-drySN", "wrfout_d01_2009-07-01_00:00:00"))

    lat_solano = 38.166
    lon_solano = -121.817

    # get index for wind farm
    ixlat, ixlon = get_box_index(fCR1, lat_solano, lon_solano)

    # get smois change regions

    fig, ax = plt.subplots(nrows=1, ncols=2)
    m1 = Basemap(width=fCR1.DX*1.2*getattr(fCR1,'WEST-EAST_GRID_DIMENSION'),\
        height=fCR1.DY*1.2*getattr(fCR1,'SOUTH-NORTH_GRID_DIMENSION'),resolution='l',\
        projection='lcc',lat_1=fCR1.TRUELAT1,lat_2=fCR1.TRUELAT2,lat_0=fCR1.CEN_LAT,lon_0=fCR1.CEN_LON,\
        ax=ax[0])
    m2 = Basemap(width=fCR1.DX*1.2*getattr(fCR1,'WEST-EAST_GRID_DIMENSION'),\
        height=fCR1.DY*1.2*getattr(fCR1,'SOUTH-NORTH_GRID_DIMENSION'),resolution='l',\
        projection='lcc',lat_1=fCR1.TRUELAT1,lat_2=fCR1.TRUELAT2,lat_0=fCR1.CEN_LAT,lon_0=fCR1.CEN_LON,\
        ax=ax[1])

    x, y = m1(fCR1.variables['XLONG'][0,:,:], fCR1.variables['XLAT'][0,:,:])

    # plot topography on left panel
    im = m1.pcolormesh(x, y, fCR1.variables['HGT'][0,:,:])
    fig.colorbar(im, ax=ax[0])

    # plot smois change regions on right panel
    im_smois = get_smois_regions(fCR1, fCV1, fSN1)
    i2 = m2.pcolormesh(x, y, np.ma.masked_array(im_smois, np.isnan(im_smois)))
    cax = fig.colorbar(i2, ax=ax[1])
    #cax.set_visible(False)
    m2.ax.annotate("CR", m2(-123.5, 39), color='gray')
    m2.ax.annotate("CV", m2(-121.8, 37.8), color='gray')
    m2.ax.annotate("SN", m2(-120.5, 38.5), color='gray')

    # draw coastlines and grid
    decorate_map(fCR1, m1)
    decorate_map(fCR1, m2)
    # ANNOTATE

    # plot domain boxes on both panels and annotate
    plot_domain_box(fCR1, m1, "d01")
    plot_domain_box(fCR2, m1, "d02")
    plot_domain_box(fCR1, m2, "d01")
    plot_domain_box(fCR2, m2, "d02")

    # plot bay box
    plot_box(fCR1, m1, [42, 48], [44, 52], color='white')
    plot_box(fCR2, m1, [70, 88], [76, 100], color='red')

    fig.set_size_inches(12, 6)
    #fig.savefig(os.path.join(root_dir, "domain_map.pdf"))
    plt.show()


    # plot solano point on both panels

    # m.plot(x, y, '+k')
    # m.plot(x[ix_lat, ix_lon], y[ix_lat, ix_lon], 'or')

    # minlon = np.floor(np.min(f.variables['XLONG'][0,:,:]/10.))*10.
    # maxlon = np.ceil(np.max(f.variables['XLONG'][0,:,:]/10.))*10.
    # minlat = np.floor(np.min(f.variables['XLAT'][0,:,:]/10.))*10.
    # maxlat = np.ceil(np.max(f.variables['XLAT'][0,:,:]/10.))*10.

    # m.drawparallels(np.arange(minlat, maxlat, 5))
    # m.drawmeridians(np.arange(minlon, maxlon, 5))
    # m.drawcoastlines()
    # m.drawstates()

    # run_name, tmp = self.files_in[i].split('/')
    # domain = tmp.split('_')[1]
    # print run_name, domain
    # fig.savefig(os.path.join(self.root_dir, 'domain_'+run_name+'_'+domain+'.png'))
    # plt.close()

