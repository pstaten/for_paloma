import numpy as np
import ps_ccmp_mpsi
from netCDF4 import Dataset
from subprocess import call
from os.path import isfile
from matplotlib import pyplot as plt

# calculate the divergent wind using NCL
ncl_script = "calc_div_aqua.ncl"
inname = 'small_slice.nc'
divname = 'small_div.nc'
if not isfile(divname):
    callstr = [
        "ncl",
        'infile="%s"' % inname,
        'divfile="%s"' % divname,
        ncl_script,
    ]
    call(callstr)

# read in some data from the input file...
# there's some abnoxious 'transpose' stuff because of how the fortran code needs the data arranged
# the squeeze stuff just gets rid of degenerate dimensions.
# (I really should get more comfortable with xarray instead.)
indat = Dataset(inname)
lat = np.squeeze(indat.variables["lat"][:])
lon = np.squeeze(indat.variables["lon"][:])
pp = np.squeeze(indat.variables['pfull'][::-1])*100.
ps = np.transpose(np.squeeze(indat.variables['ps'][:]))
indat.close()

# ...and from the divergent wind file.
vdat = Dataset(divname)
vd = np.transpose(np.squeeze(vdat.variables['vd'])[::-1,:,:])
vdat.close()

# calculate the local hadley circulation using fortran
missing_val = -5e35
psi_d_full = ps_ccmp_mpsi.ps_ccmp_mpsi(vd, lat, pp, ps, missing_val)

# save a filled contour plot as a pdf
plt.set_cmap('RdBu')
plt.contourf(lon, lat, np.transpose(psi_d_full[:,:,10]))
plt.xlabel('$^\circ$lon')
plt.ylabel('$^\circ$lat')
plt.title('A sample regionally varying "Hadley cell"')
plt.savefig('test_plot.pdf')
plt.close()
