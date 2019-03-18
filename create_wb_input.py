# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 21:16:13 2019

@author: N.B.Tran
    
    Create netCDF file for simple water balance model
"""
#%% Import libs & Define parameters and paths
import numpy as np
import os
import tempfile
import pandas as pd
from WA_Hyperloop import becgis
import netCDF4
from matplotlib import pyplot as plt
from WA_Hyperloop.sheet4_functions.sheet4_functions import fractions,update_irrigation_fractions
from WA_Hyperloop import get_dictionaries as gd

start_date= '2010-01-01'
years=7
lu_ras=r"E:\Litani_WA\Data\LUCWA\LUCWA_Litan_v3i.tif"
p_path=r'E:\Litani_WA\Data\WaPOR\PCP_M\0_L1_monthly_PCP_{yyyy}{mm}.tif'
et_path=r'E:\Litani_WA\Data\WaPOR\AETI_M\AETI_{yyyy}{mm}.tif'
lai_path=r'E:\Litani_WA\Data\LAI\monthly_8daily\LAI_{yyyy}{mm}.tif'
swi_path=r'E:\\Litani_WA\\Data\\Soil_Moisture\\SWI\\040\\Mean\\SWI_{yyyy}{mm}.tif' 
swio_path=r'E:\\Litani_WA\\Data\\Soil_Moisture\\SWI\\040\\First\\SWI_{yyyy}{mm}.tif'
swix_path=r'E:\Litani_WA\Data\Soil_Moisture\SWI\040\Last\SWI_{yyyy}{mm}.tif'
qratio_y_path=r'E:\Litani_WA\Data\Runoff_GLDAS\runoff_ratio_yearly\runoff_ratio_{yyyy}.tif'
qratio_m_path=r'E:\Litani_and_Jordan\Data\Runoff_GLDAS\runoff_ratio_resampled\runoff_ratio_{yyyy}{mm}.tif'
etg_path=r'E:\Litani_WA\Hyperloop\Old_ETgbsplit\data\etg\ETgreen_{yyyy}{mm}.tif'
etb_path=r'E:\Litani_WA\Hyperloop\Old_ETgbsplit\data\etb\ETblue_{yyyy}{mm}.tif'
thetasat_ras=r"E:\Litani_WA\Data\HiHydroSoils\thetasat_topsoil.tif"
rootdepth_ras=r"E:\Litani_WA\Data\GlobCover\rd_mm.tif"
i_path= r'E:\Litani_WA\Data\WaPOR\Interception_M\I_{yyyy}{mm}.tif'
input_nc=r'E:\Litani_WA\Waterbal\input_nc\test_input_7years.nc'
logfile=True
epsg=4326

#Calculate  sw_supply_fraction_tif
lucs = gd.get_sheet4_6_classes() 
equiped_sw_irrigation_tif= r'E:\Litani_WA\Data\global_data\gmia_v5_aeisw_pct_aei.tif'
sw_supply_fractions = gd.get_sheet4_6_fractions()
sw_supply_fraction_tif = fractions(lu_ras, sw_supply_fractions, lucs, r'E:\Litani_WA\Waterbal', filename = 'sw_supply_fraction.tif')
sw_supply_fraction_tif = update_irrigation_fractions(lu_ras, sw_supply_fraction_tif, lucs, equiped_sw_irrigation_tif)
#%% Define netCDF file dimension parameters

srs, ts, te, ndv=becgis.GetGdalWarpInfo(lu_ras)
driver, NDV, xsize, ysize, GeoT, Projection = becgis.GetGeoInfo(lu_ras)
xmin, ymin, xmax, ymax = [float(t) for t in te.split(' ')]
latlim=[ymin,ymax]
lonlim=[xmin,xmax]

time_range = pd.date_range(start_date, periods=12*years, freq='MS') #DatetimeIndex
time_ls = [d.strftime('%Y%m') for d in time_range] #list of str: 'yyyymm'
time_dt = [pd.to_datetime(i, format='%Y%m') for i in time_ls] #list of timestamps

time_n = len(time_ls) #size of Time dimension

years_ls=set()
years_ls = [i.year for i in time_dt
            if i.year not in years_ls and not years_ls.add(i.year)]#[:-1]

time_indeces = {}

for j, item in enumerate(years_ls):
    temp_ls = [int(i.strftime('%Y%m')) for i in
               pd.date_range(str(item) + '0101',
                             str(item) + '1231', freq='MS')]
    time_indeces[item] = [time_ls.index(str(i)) for i in temp_ls]
    
lat_ls = pd.np.arange(latlim[0] + 0.5*GeoT[1], latlim[1],GeoT[1])
lat_ls = lat_ls[::-1]  # ArcGIS numpy
lon_ls = pd.np.arange(lonlim[0] - 0.5*GeoT[5], lonlim[1] - 0.5*GeoT[5],-GeoT[5])
lat_n= ysize
lon_n= xsize
#%% Listing variables and info    
monthly_var_dict={'p':  [p_path,'Precipitation_M', 'Precipitation','mm/month'],
                   'et': [et_path,'Evapotranspiration_M','Evapotranspiration','mm/month'],
                   'lai': [lai_path,'LeafAreaIndex_M', 'Leaf Area Index','m2/m2'],
                   'swi': [swi_path,'SWI_M','Soil Water Index - Monthly mean','%'],
                   'swio': [swio_path,'SWIo_M','Soil water index - First day of the month','%'],
                   'swix': [swix_path,'SWIx_M','Soil water index - Last day of the month','%'],
                   'etg': [etg_path,'ETgreen_M','ET Green','mm/month'],
                   'etb': [etb_path,'ETblue_M','ET Blue', 'mm/month'],
                   'i': [i_path, 'Interception_M', 'Interception', 'mm/month'],
                   'qratiom':[qratio_m_path,'RunoffRatio_M', 'Monthly Runoff Ratio','-']
                   }
yearly_var_dict={'qratio': [qratio_y_path,'RunoffRatio_Y', 'Yearly Runoff Ratio','-'],
                  }
stat_var_dict={'thetasat': [thetasat_ras,'SaturatedWaterContent','Saturated water content (top soil)','cm3/cm3'],
                'rootdepth': [rootdepth_ras,'RootDepth','Root depth','mm'],
                'lu': [lu_ras,'LandUse','LandUse Classes WA','cm3/cm3'],
                'sw_supply_frac':[sw_supply_fraction_tif,'SWSupplyFrac','Surfacewater Supply Fraction','-']
                }
######################################
###Create netCDF4
#%% Open netCDF file
nc_file = netCDF4.Dataset(input_nc, 'w', format="NETCDF4")
nc_file.set_fill_on()
#%% Create dimension variables
#Ccreate dimensons
lat_dim = nc_file.createDimension('latitude', lat_n)
lon_dim = nc_file.createDimension('longitude', lon_n)
month_dim = nc_file.createDimension('time_yyyymm', time_n)
year_dim = nc_file.createDimension('time_yyyy', len(years_ls))

#Define dimension Variables
crs_var = nc_file.createVariable('crs', 'i', (), fill_value=-9999)
crs_var.standard_name = 'crs'
crs_var.grid_mapping_name = 'latitude_longitude'
crs_var.crs_wkt = unicode(str(Projection),'utf-8')

lat_var = nc_file.createVariable('latitude', 'f8', ('latitude'),fill_value=-9999)
lat_var.units = 'degrees_north'
lat_var.standard_name = 'latitude'
lat_var[:] = lat_ls

lon_var = nc_file.createVariable('longitude', 'f8', ('longitude'),fill_value=-9999)
lon_var.units = 'degrees_east'
lon_var.standard_name = 'longitude'
lon_var[:] = lon_ls

month_var = nc_file.createVariable('time_yyyymm', 'l', ('time_yyyymm'),fill_value=-9999)
month_var.standard_name = 'time'
month_var.format = 'YYYYMM'
month_var[:] = time_ls

year_var = nc_file.createVariable('time_yyyy', 'l', ('time_yyyy'),fill_value=-9999)
year_var.standard_name = 'time'
year_var.format = 'YYYY'
year_var[:] = years_ls

#temporary folder to load data
temp_dir = tempfile.mkdtemp()
#%% Create variables and load data
#static variables
for var in stat_var_dict:    
    path,shortname,longname,unit=stat_var_dict[var]
    print var,path,shortname,longname,unit
    #Create var
    exec('{0}_var = nc_file.createVariable(\'{1}\',\'f8\',(\'latitude\',\'longitude\'),fill_value=-9999)'.format(var,shortname))
    exec('{0}_var.long_name=\'{1}\''.format(var,longname))
    exec('{0}_var.units=\'{1}\''.format(var,unit))
    exec('{0}_var.grid_mapping=\'crs\''.format(var))
    #Load data to var    
    exec('{0}_temp=becgis.MatchProjResNDV(lu_ras,[path],temp_dir)'.format(var))
    exec('{0}_array=becgis.OpenAsArray({0}_temp[0],nan_values=True)'.format(var))
    exec('{0}_var[:,:]={0}_array[:,:]'.format(var))
#monthly variables
for var in monthly_var_dict:    
    path,shortname,longname,unit=monthly_var_dict[var]
    print var,path,shortname,longname,unit
    #Create var
    exec('{0}_var = nc_file.createVariable(\'{1}\',\'f8\',(\'time_yyyymm\',\'latitude\',\'longitude\'),fill_value=-9999)'.format(var,shortname))
    exec('{0}_var.long_name=\'{1}\''.format(var,longname))
    exec('{0}_var.units=\'{1}\''.format(var,unit))
    exec('{0}_var.grid_mapping=\'crs\''.format(var))
    #Load data to var
    for yyyymm in time_ls:
        yyyy = yyyymm[:4]
        mm = yyyymm[-2:]
        ras = path.format(yyyy=yyyy, mm=mm)
        print "{0}\t{1}".format(var, ras)
        exec('{0}_temp=becgis.MatchProjResNDV(lu_ras,[ras],temp_dir)'.format(var))
        exec('{0}_array=becgis.OpenAsArray({0}_temp[0],nan_values=True)'.format(var))
        t_index = time_ls.index(yyyymm)
        exec('{0}_var[t_index,:,:]={0}_array[:,:]'.format(var))        
#yearly variables
for var in yearly_var_dict:    
    path,shortname,longname,unit=yearly_var_dict[var]
    print var,path,shortname,longname,unit
    #Create var
    exec('{0}_var = nc_file.createVariable(\'{1}\',\'f8\',(\'time_yyyy\',\'latitude\',\'longitude\'),fill_value=-9999)'.format(var,shortname))
    exec('{0}_var.long_name=\'{1}\''.format(var,longname))
    exec('{0}_var.units=\'{1}\''.format(var,unit))
    exec('{0}_var.grid_mapping=\'crs\''.format(var))
    #load data to var
    for yyyy in years_ls:
        ras = path.format(yyyy=yyyy)
        print "{0}\t{1}".format(var, ras)
        exec('{0}_temp=becgis.MatchProjResNDV(lu_ras,[ras],temp_dir)'.format(var))
        exec('{0}_array=becgis.OpenAsArray({0}_temp[0],nan_values=True)'.format(var))
        y_index = years_ls.index(yyyy)
        exec('{0}_var[y_index,:,:]={0}_array[:,:]'.format(var))         

#%% Close file
nc_file.close()

