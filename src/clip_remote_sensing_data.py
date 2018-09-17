"""

regrid data sets
- this routine clips and regrids the data layers to be fed into the RFR analysis to estimate potential biomass
- it predominately makes use of gdal
- the 30 arcsec worldclim2 bio data are used as the template onto which other data are regridded
- regridding modes are as follows:
  - AGB (Avitabile) - nearest neighbour
  - ESA CCI landcover - mode
  - SOILGRIDS - nearest neighbour
- the extent is specified based on the bounding box (lat long limits)
- if required a mask is generated to clip the output to the specified national boundaries

- the regridding routines utilise the following packages:
  - GDAL
  - NCO

14/09/2018
David T. Milodowski

"""

# This script just calls GDAL and NCO to clip existing datasets to a specified bounding box
import os
import glob

# Bounding Box
W = 92.
E = 154.
N = 8.
S = -12.

# create some directories to host the processed data_io
prefix = 'INDO'
outdir = '/disk/scratch/local.2/dmilodow/PotentialBiomass/processed/'
os.system('mkdir %s' % outdir)
os.system('mkdir %s%s' % (outdir,prefix))
os.system('mkdir %s%s/wc2' % (outdir,prefix))
os.system('mkdir %s%s/agb' % (outdir,prefix))
os.system('mkdir %s%s/esacci' % (outdir,prefix))
os.system('mkdir %s%s/soilgrids' % (outdir,prefix))


# Start with the worldclim2 data. This should be straightforward using gdal as are just clipping
# to the target extent
wc2dir = '/disk/scratch/local.2/worldclim2/'
wc2files = glob.glob('%swc2*tif' % wc2dir);wc2files.sort()
wc2subset = []
wc2vars = []
for ff,fname in enumerate(wc2files):
    variable = fname.split('_')[-1].split('.')[0]
    outfname = fname.split('/')[-1][:-4]
    #if int(variable) in ['01','04','05','06','12','13','14','15']:
    if int(variable) in range(1,20):
        wc2vars.append(variable)
        wc2subset.append(fname)
        os.system("gdalwarp -overwrite -te %f %f %f %f %s %s%s/wc2/%s_%s.tif" % (W,S,E,N,wc2files[ff],outdir,prefix,outfname,prefix))


"""
# ERA Interim
era_vars = ['prcp','t2m','u10w','v10w','d2m']
era_path = '/disk/scratch/local.2/dmilodow/ERAinterim/source_files/0.25deg_Mexico/'
era_savepath = '/exports/csce/datastore/geos/users/dmilodow/FOREST2020/hazardsINLAndscapes/fireINLAndscapes/data/external/jalisco/era_interim/'

for yy in range(start_year,end_year+1):
    for mm in range(1,13):
        for vv in range(0,len(era_vars)):
            src_file = '%s%s_%04i%02i.nc' % (era_path,era_vars[vv],yy,mm)
            save_file = '%s%s_%04i%02i.nc' % (era_savepath,era_vars[vv],yy,mm)
            os.system("ncks -O -d latitude,%f,%f -d longitude,%f,%f %s %s" % (S,N,W_,E_,src_file,save_file))

# MODIS Burned Area
modis_path = '/disk/scratch/local.2/MCD64A1.006/processed/'
modis_savepath = '/exports/csce/datastore/geos/users/dmilodow/FOREST2020/hazardsINLAndscapes/fireINLAndscapes/data/external/jalisco/modis_mcd64a1/'
for yy in range(start_year,end_year+1):
    for mm in range(1,13):
        src_file = '%s%04i.%02i.01_mex.nc' % (modis_path,yy,mm)
        save_file= '%s%04i.%02i.01.nc' % (modis_savepath,yy,mm)
        os.system("ncks -O -d lat,%f,%f -d lon,%f,%f %s %s" % (S,N,W,E,src_file,save_file))

# Land cover
lc_path = '/home/dmilodow/DataStore_GCEL/ESA_CCI_landcover/'
lc_savepath = '/exports/csce/datastore/geos/users/dmilodow/FOREST2020/hazardsINLAndscapes/fireINLAndscapes/data/external/jalisco/esa_cci/'
for yy in range(start_year,end_year+1):
    src_file = '%sESACCI-LC-L4-LCCS-Map-300m-P1Y-%04i-v2.0.7.nc' % (lc_path,yy)
    save_file= '%sESACCI-LC-L4-LCCS-Map-300m-P1Y-%04i-v2.0.7.nc' % (lc_savepath,yy)
    os.system("ncks -O -d lat,%f,%f -d lon,%f,%f %s %s" % (S,N,W,E,src_file,save_file))

# worldpop
years = [2000,2005,2010,2015]
pop_path = '/home/dmilodow/DataStore_GCEL/WorldPop/Population/'
pop_savepath = '/exports/csce/datastore/geos/users/dmilodow/FOREST2020/hazardsINLAndscapes/fireINLAndscapes/data/external/jalisco/worldpop/'
#'LAC_PPP_2015_adj_v2.tif'
for yy in range(0,4):
    os.system("gdalwarp -overwrite -te %f %f %f %f -r mode %sLAC_PPP_%04i_adj_v2.tif %sLAC_PPP_%04i_adj_v2.tif" % (W,S,E,N,pop_path,years[yy],pop_savepath,years[yy]))

# topography
srtm_path = '/disk/scratch/local.2/dmilodow/SRTM90m/Mexico/'
srtm_savepath = '/exports/csce/datastore/geos/users/dmilodow/FOREST2020/hazardsINLAndscapes/fireINLAndscapes/data/external/jalisco/srtm/'
os.system("gdalwarp -overwrite -te %f %f %f %f -r mode %smexico_srtm90.tif %smexico_srtm90.tif" % (W,S,E,N,srtm_path,srtm_savepath))

"""
