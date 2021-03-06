# -*- coding: utf-8 -*-
"""
Authors: Elga Salvadore
         IHE Delft 2019
Contact: e.salvadore@un-ihe.org
Repository: 
Module: WAWB
The code is a simplified version of WaterPix model (Authors: Gonzalo Espinoza and Claire Michailovsky)
"""
import os
from __future__ import division
import datetime as dt
import numpy as np
import pandas as pd
import netCDF4
from wa_wb.functions import (lai_and_soil_calculations, SCS_surface_runoff, 
                             baseflow_calculation, SCS_surface_runoff_gr, total_supply)

def run(input_nc, output_nc, rootdepth_par = 1.1,
        wateryear = ['0101','1231'], filter_par=0.5, min_qratio=0.40, log=True):
    if log:
        fn=output_nc.replace('.nc','.txt')
        f=open(fn,'w')
        f.write('input_nc: {0} \n'.format(input_nc))
        f.write('output_nc: {0} \n'.format(output_nc))
        f.write('rootdepth_par: {0} \n'.format(rootdepth_par))
        f.write('wateryear: {0} \n'.format(wateryear))
        f.write('filter_par: {0} \n'.format(filter_par))
        f.write('min_qratio: {0} \n'.format(min_qratio))        
        f.close()
    '''
    Executes the main module of WAWB
    '''
    # Read inputs
    started = dt.datetime.now()
    print 'Reading input netcdf ...'
    inp_nc = netCDF4.Dataset(input_nc, 'r')
    ncv = inp_nc.variables
    inp_crs = ncv['crs']
    inp_lat = ncv['latitude']
    inp_lon = ncv['longitude']
    inp_time = ncv['time_yyyymm']
    lat_ls = list(inp_lat[:])
    lon_ls = list(inp_lon[:])
    lat_n = len(lat_ls)
    lon_n = len(lon_ls)
    time_ls = list(inp_time[:])
    time_dt = [pd.to_datetime(i, format='%Y%m')
               for i in time_ls]
    time_n = len(time_ls)
    years_ls = set()
    if wateryear == ['0101','1231']:
        years_ls = [i.year for i in time_dt
                    if i.year not in years_ls and not years_ls.add(i.year)]        
    else:
        years_ls = [i.year for i in time_dt
                    if i.year not in years_ls and not years_ls.add(i.year)][:-1]
    years_n = len(years_ls)
    time_indeces = {}
    for j in range(years_n):
        temp_ls = [int(i.strftime('%Y%m')) for i in
                   pd.date_range(start=str(years_ls[j]) + wateryear[0],
                                 periods=12, freq='MS')]
        time_indeces[years_ls[j]] = [time_ls.index(i) for i in temp_ls]

    for key in time_indeces.keys():
        if time_indeces[key] != range(time_indeces[key][0],
                                      time_indeces[key][-1] + 1):
            raise Exception('The year {0} in the netcdf file is incomplete'
                            ' or the dates are non-consecutive')
    # Create ouput NetCDF
    print 'Creating output netcdf ...'
    out_nc = netCDF4.Dataset(output_nc, 'w', format="NETCDF4")
    std_fv = -9999
    
    out_nc.createDimension(inp_lat.standard_name, lat_n)
    out_nc.createDimension(inp_lon.standard_name, lon_n)
    out_nc.createDimension('time_yyyymm', time_n)
    out_nc.createDimension('time_yyyy', years_n)
    
    # Reference system
    crs_var = out_nc.createVariable(inp_crs.standard_name, 'i', (),
                                    fill_value=std_fv)
    crs_var.standard_name = inp_crs.standard_name
    crs_var.grid_mapping_name = inp_crs.grid_mapping_name
    crs_var.crs_wkt = inp_crs.crs_wkt
    # Latitude
    lat_var = out_nc.createVariable(inp_lat.standard_name, 'f8',
                                    (inp_lat.standard_name),
                                    fill_value=inp_lat._FillValue)
    lat_var.units = inp_lat.units
    lat_var.standard_name = inp_lat.standard_name
    # Longitude
    lon_var = out_nc.createVariable(inp_lon.standard_name, 'f8',
                                    (inp_lon.standard_name),
                                    fill_value=inp_lon._FillValue)
    lon_var.units = inp_lon.units
    lon_var.standard_name = inp_lon.standard_name
    # Time (month)
    time_var = out_nc.createVariable('time_yyyymm', 'l', ('time_yyyymm'),
                                     fill_value=inp_time._FillValue)
    time_var.standard_name = inp_time.standard_name
    time_var.format = inp_time.format
    # Time (year)
    year_var = out_nc.createVariable('time_yyyy', 'l', ('time_yyyy'),
                                     fill_value=std_fv)
    year_var.standard_name = 'time_yyyy'
    year_var.format = 'yyyy'
    # FillValues
    lu_fv = ncv['LandUse']._FillValue
    p_fv = ncv['Precipitation_M']._FillValue
    et_fv = ncv['Evapotranspiration_M']._FillValue
    lai_fv = ncv['LeafAreaIndex_M']._FillValue
    swi_fv = ncv['SWI_M']._FillValue
    swio_fv = ncv['SWIo_M']._FillValue
    swix_fv = ncv['SWIx_M']._FillValue
    qratio_y_fv = ncv['RunoffRatio_Y']._FillValue
    thetasat_fv = ncv['SaturatedWaterContent']._FillValue
    rootdepth_fv = ncv['RootDepth']._FillValue
    etg_fv = ncv['ETgreen_M']._FillValue
    etb_fv = ncv['ETblue_M']._FillValue
    interception_fv = ncv['Interception_M']._FillValue
    # Copy data
    lat_var[:] = lat_ls
    lon_var[:] = lon_ls
    time_var[:] = time_ls
    year_var[:] = years_ls
    # Create output NetCDF variables:
    # Surface runoff (monthly)
    ss_var = out_nc.createVariable('SurfaceRunoff_M', 'f8',
                                   ('time_yyyymm', 'latitude', 'longitude'),
                                   fill_value=std_fv)
    ss_var.long_name = 'Surface runoff (fast)'
    ss_var.units = 'mm/month'
    ss_var.grid_mapping = 'crs'
    # Surface runoff (yearly)
    ssy_var = out_nc.createVariable('SurfaceRunoff_Y', 'f8',
                                    ('time_yyyy', 'latitude', 'longitude'),
                                    fill_value=std_fv)
    ssy_var.long_name = 'Surface runoff (fast)'
    ssy_var.units = 'mm/year'
    ssy_var.grid_mapping = 'crs'
    # Baseflow (monthly)
    bf_var = out_nc.createVariable('Baseflow_M', 'f8',
                                   ('time_yyyymm', 'latitude', 'longitude'),
                                   fill_value=std_fv)
    bf_var.long_name = 'Baseflow (slow)'
    bf_var.units = 'mm/month'
    bf_var.grid_mapping = 'crs'
    # Baseflow (yearly)
    bfy_var = out_nc.createVariable('Baseflow_Y', 'f8',
                                    ('time_yyyy', 'latitude', 'longitude'),
                                    fill_value=std_fv)
    bfy_var.long_name = 'Baseflow (slow)'
    bfy_var.units = 'mm/year'
    bfy_var.grid_mapping = 'crs'
    # Total runoff (monthly)
    sr_var = out_nc.createVariable('TotalRunoff_M', 'f8',
                                   ('time_yyyymm', 'latitude', 'longitude'),
                                   fill_value=std_fv)
    sr_var.long_name = 'Total runoff'
    sr_var.units = 'mm/month'
    sr_var.grid_mapping = 'crs'
    # Total runoff (yearly)
    sry_var = out_nc.createVariable('TotalRunoff_Y', 'f8',
                                    ('time_yyyy', 'latitude', 'longitude'),
                                    fill_value=std_fv)
    sry_var.long_name = 'Total runoff'
    sry_var.units = 'mm/year'
    sry_var.grid_mapping = 'crs'
    # Storage change - soil moisture (monthly)
    dsm_var = out_nc.createVariable('StorageChange_M', 'f8',
                                    ('time_yyyymm', 'latitude', 'longitude'),
                                    fill_value=std_fv)
    dsm_var.long_name = 'Change in soil moisture storage'
    dsm_var.units = 'mm/month'
    dsm_var.grid_mapping = 'crs'
    # Storage change - soil moisture (yearly)
    dsmy_var = out_nc.createVariable('StorageChange_Y', 'f8',
                                     ('time_yyyy', 'latitude', 'longitude'),
                                     fill_value=std_fv)
    dsmy_var.long_name = 'Change in soil moisture storage'
    dsmy_var.units = 'mm/year'
    dsmy_var.grid_mapping = 'crs'
    # Percolation (monthly)
    per_var = out_nc.createVariable('Percolation_M', 'f8',
                                    ('time_yyyymm', 'latitude', 'longitude'),
                                    fill_value=std_fv)
    per_var.long_name = 'Percolation'
    per_var.units = 'mm/month'
    per_var.grid_mapping = 'crs'
    # Percolation (yearly)
    pery_var = out_nc.createVariable('Percolation_Y', 'f8',
                                     ('time_yyyy', 'latitude', 'longitude'),
                                     fill_value=std_fv)
    pery_var.long_name = 'Percolation'
    pery_var.units = 'mm/year'
    pery_var.grid_mapping = 'crs'
    # Supply (monthly)
    sup_var = out_nc.createVariable('Supply_M', 'f8',
                                    ('time_yyyymm', 'latitude', 'longitude'),
                                    fill_value=std_fv)
    sup_var.long_name = 'Supply'
    sup_var.units = 'mm/month'
    sup_var.grid_mapping = 'crs'
    # Supply (yearly)
    supy_var = out_nc.createVariable('Supply_Y', 'f8',
                                     ('time_yyyy', 'latitude', 'longitude'),
                                     fill_value=std_fv)
    supy_var.long_name = 'Supply'
    supy_var.units = 'mm/year'
    supy_var.grid_mapping = 'crs'
    
    # Root depth soil moisture (monthly)
    rdsm_var = out_nc.createVariable('RootDepthSoilMoisture_M', 'f8',
                                     ('time_yyyymm', 'latitude', 'longitude'),
                                     fill_value=rootdepth_fv)
    rdsm_var.long_name = 'Root depth soil moisture'
    rdsm_var.units = 'cm3/cm3'
    rdsm_var.grid_mapping = 'crs'
    
    # Incremental surface runoff (monthly)
    incss_var = out_nc.createVariable('IncrementalRunoff_M', 'f8',
                                      ('time_yyyymm', 'latitude', 'longitude'),
                                      fill_value=std_fv)
    incss_var.long_name = 'Incremental runoff'
    incss_var.units = 'mm/month'
    incss_var.grid_mapping = 'crs'
    # Incremental surface runoff (yearly)
    incssy_var = out_nc.createVariable('IncrementalRunoff_Y', 'f8',
                                       ('time_yyyy', 'latitude', 'longitude'),
                                       fill_value=std_fv)
    incssy_var.long_name = 'Incremental runoff'
    incssy_var.units = 'mm/year'
    incssy_var.grid_mapping = 'crs'
    # Incremental percolation (monthly)
    incper_var = out_nc.createVariable('IncrementalPercolation_M', 'f8',
                                       ('time_yyyymm',
                                        'latitude', 'longitude'),
                                       fill_value=std_fv)
    incper_var.long_name = 'Incremental Percolation'
    incper_var.units = 'mm/month'
    incper_var.grid_mapping = 'crs'
    # Incremental percolation (yearly)
    incpery_var = out_nc.createVariable('IncrementalPercolation_Y', 'f8',
                                        ('time_yyyy', 'latitude', 'longitude'),
                                        fill_value=std_fv)
    incpery_var.long_name = 'Incremental Percolation'
    incpery_var.units = 'mm/year'
    incpery_var.grid_mapping = 'crs'
    
    # corrected rootdepth
    rootdepth_var = out_nc.createVariable('RootDepth_cor', 'f8',
                                           ('latitude', 'longitude'),
                                           fill_value=-9999)
    rootdepth_var.long_name = 'Root depth corrected'
    rootdepth_var.units = 'mm'
    rootdepth_var.grid_mapping = 'crs'
    
    # corrected qratio_y
    qratio_y_corr_var = out_nc.createVariable('Qratio_cor', 'f8',
                                           ('latitude', 'longitude'),
                                           fill_value=-9999)
    rootdepth_var.long_name = 'Q ratio corrected'
    rootdepth_var.units = 'mm'
    rootdepth_var.grid_mapping = 'crs'



    for yyyy in years_ls:
        print '\tyear: {0}'.format(yyyy)
        yyyyi = years_ls.index(yyyy)
        ti1 = time_indeces[yyyy][0]
        ti2 = time_indeces[yyyy][-1] + 1
 
        # Read data
        lu = np.array(ncv['LandUse'][:, :])
        p = np.array(ncv['Precipitation_M'][ti1:ti2, :, :])
        et = np.array(ncv['Evapotranspiration_M'][ti1:ti2, :, :])
        lai = np.array(ncv['LeafAreaIndex_M'][ti1:ti2, :, :])
        swi = np.array(ncv['SWI_M'][ti1:ti2, :, :])
        swio = np.array(ncv['SWIo_M'][ti1:ti2, :, :])
        swix = np.array(ncv['SWIx_M'][ti1:ti2, :, :])
        qratio_y = np.array(ncv['RunoffRatio_Y'][yyyyi, :, :])
        etg = np.array(ncv['ETgreen_M'][ti1:ti2, :, :])
        etb = np.array(ncv['ETblue_M'][ti1:ti2, :, :])
        thetasat = np.array(ncv['SaturatedWaterContent'][:, :])
        rootdepth = np.array(ncv['RootDepth'][:, :])
        # calibration parameter rootdepth_par
        rootdepth = rootdepth*rootdepth_par
        
        interception = np.array(ncv['Interception_M'][ti1:ti2, :, :])
        # Check for NoData values
        lu[np.isclose(lu, lu_fv)] = np.nan
        p[np.isclose(p, p_fv)] = np.nan
        et[np.isclose(et, et_fv)] = np.nan
        lai[np.isclose(lai, lai_fv)] = np.nan
        etg[np.isclose(etg, etg_fv)] = np.nan
        etb[np.isclose(etb, etb_fv)] = np.nan
        # filling gaps in LAI
        lai = np.nan_to_num(lai)

        swi[np.isclose(swi, swi_fv)] = np.nan
        swio[np.isclose(swio, swio_fv)] = np.nan
        swix[np.isclose(swix, swix_fv)] = np.nan
        
        # qratio limiting Qgw:
        qratio_y[np.isclose(qratio_y, qratio_y_fv)] = np.nan
        qratio_y[qratio_y<min_qratio]=min_qratio
        
        thetasat[np.isclose(thetasat, thetasat_fv)] = np.nan
        rootdepth[np.isclose(rootdepth, rootdepth_fv)] = np.nan
        interception[np.isclose(interception, interception_fv)] = np.nan
        
        #rootdepth soil moisture
        thetarz, thetarzo, thetarzx, dsm = lai_and_soil_calculations(thetasat, lai, swi, swio, swix, rootdepth)
        
        # Check P-ET-dsm
        p_et_dsm = np.sum(p, axis=0) - np.sum(et, axis=0) - np.sum(dsm, axis=0)
        
        # simple supply comuptation based on blue ET and LU
        supply = total_supply(etb, lu)
        
        Qsw_gr = np.zeros(np.shape(et))
        Qsw = np.zeros(np.shape(et))
        Qsw_gr= SCS_surface_runoff_gr(p, interception, rootdepth, thetasat, thetarz)
        Effective_supply = supply-etb
        Qsw= SCS_surface_runoff(p, Effective_supply, interception, rootdepth, thetasat, thetarz)
#        Qsw= SCS_surface_runoff(p, supply, interception, rootdepth, thetasat, thetarz)    
        # if rain = 0 and supply = to ETb --> Qsw = 0
        Qsw = np.where(np.logical_and(p==0, supply==etb), 0, Qsw)
        
        # check for green pixels
        Qsw_gr = np.where(p==0, 0, Qsw_gr)
        # otherwise I mess up delta Qsw
        Qsw_gr = np.where(etb==0, Qsw_gr, 0)
                
        incr_Qsw = Qsw-Qsw_gr
        incr_Qsw[incr_Qsw<0]=0
        
        perc_gr = p-etg-dsm-Qsw_gr
        # otherwise I mess up delta perc
        perc_gr = np.where(etb==0, perc_gr, 0)
        perc = p+supply-et-dsm-Qsw
        incr_perc = perc-perc_gr
        incr_perc[incr_perc<0]=0 # I'm not sure this should only be positive...
        
        Qgw_gr = baseflow_calculation(Qsw_gr, filter_par, qratio_y)
        Qgw = baseflow_calculation(Qsw, filter_par, qratio_y)
        incr_Qgw = Qgw - Qgw_gr
        incr_Qgw[incr_Qgw<0]=0
        incr_Q = incr_Qsw 
        Qtot = Qsw+Qgw
        
        # Store values in output NetCDF
        
        ss_var[ti1:ti2, :, :] = Qsw
        bf_var[ti1:ti2, :, :] = Qgw
        sr_var[ti1:ti2, :, :] = Qtot
        dsm_var[ti1:ti2, :, :] = dsm
        per_var[ti1:ti2, :, :] = perc
        rdsm_var[ti1:ti2, :, :] = thetarz
        sup_var[ti1:ti2, :, :] = supply
        incss_var[ti1:ti2, :, :] = incr_Q
        incper_var[ti1:ti2,:, :] = incr_perc
        rootdepth_var[:, :] = rootdepth
        qratio_y_corr_var[:,:] = qratio_y
        
    

    # Calculate yearly variables
    print 'Calculating values per year...'
    for yyyy in years_ls:
        # Time indeces
        yyyyi = years_ls.index(yyyy)
        ti1 = time_indeces[yyyy][0]
        ti2 = time_indeces[yyyy][-1] + 1
        # Sums used in efficiency calculation
        supply_yearly_val = np.sum(sup_var[ti1:ti2, :, :], axis=0)
        inc_ss_yearly_val = np.sum(incss_var[ti1:ti2, :, :], axis=0)
        inc_per_yearly_val = np.sum(incper_var[ti1:ti2, :, :], axis=0)
        # Store values
        ssy_var[yyyyi, :, :] = np.sum(ss_var[ti1:ti2, :, :], axis=0)
        incssy_var[yyyyi, :, :] = inc_ss_yearly_val
        bfy_var[yyyyi, :, :] = np.sum(bf_var[ti1:ti2, :, :], axis=0)
        sry_var[yyyyi, :, :] = np.sum(sr_var[ti1:ti2, :, :], axis=0)
        dsmy_var[yyyyi, :, :] = np.sum(dsm_var[ti1:ti2, :, :], axis=0)
        pery_var[yyyyi, :, :] = np.sum(per_var[ti1:ti2, :, :], axis=0)
        incpery_var[yyyyi, :, :] = inc_per_yearly_val
        supy_var[yyyyi, :, :] = supply_yearly_val
        
    # Finishing
    print 'Closing netcdf...'
    out_nc.close()
    ended = dt.datetime.now()
    print 'Time elapsed: {0}'.format(ended - started)
