# -*- coding: utf-8 -*-
"""
Waterbal

"""
from __future__ import division
import os
import datetime as dt
import numpy as np
import pandas as pd
import netCDF4
import get_dictionaries as gd
import numpy.ma as ma


def run(input_nc, output_nc, rootdepth_par = 1,
        wateryear = ['0101','1231'], dS_GRACE=None, log=True):
    if log:
        fn=output_nc.replace('.nc','.txt')
        f=open(fn,'w')
        f.write('input_nc: {0} \n'.format(input_nc))
        f.write('output_nc: {0} \n'.format(output_nc))
        f.write('rootdepth_par: {0} \n'.format(rootdepth_par))
        f.write('wateryear: {0} \n'.format(wateryear))
        f.write('dS_GRACE: {0} \n'.format(dS_GRACE))
       
        f.close()
    '''
    Executes the main module of WAWB
    '''
### Read inputs dimension
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
### Create dimension variables of ouput NetCDF 
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
    lat_var = out_nc.createVariable(inp_lat.standard_name, 'f4',
                                    (inp_lat.standard_name),
                                    fill_value=inp_lat._FillValue)
    lat_var.units = inp_lat.units
    lat_var.standard_name = inp_lat.standard_name
    # Longitude
    lon_var = out_nc.createVariable(inp_lon.standard_name, 'f4',
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

    # Copy data
    lat_var[:] = lat_ls
    lon_var[:] = lon_ls
    time_var[:] = time_ls
    year_var[:] = years_ls

### Read input variables
    # FillValues
    lu_fv = ncv['LandUse']._FillValue
    p_fv = ncv['Precipitation_M']._FillValue
    et_fv = ncv['Evapotranspiration_M']._FillValue
    qratio_y_fv = ncv['RunoffRatio_Y']._FillValue
    thetasat_fv = ncv['SaturatedWaterContent']._FillValue
    rootdepth_fv = ncv['RootDepth']._FillValue
    interception_fv = ncv['Interception_M']._FillValue
    SWSupplyFrac_fv= ncv['SWSupplyFrac']._FillValue

    #Data
    lu = np.array(ncv['LandUse'][:, :])
    p = np.array(ncv['Precipitation_M'][:, :, :])
    et = np.array(ncv['Evapotranspiration_M'][:, :, :])
    qratio_y = np.array(ncv['RunoffRatio_Y'][:, :, :])
    thetasat = np.array(ncv['SaturatedWaterContent'][:, :])
    rootdepth = np.array(ncv['RootDepth'][:, :])
    interception = np.array(ncv['Interception_M'][:, :, :])
    SWSupplyFrac= np.array(ncv['SWSupplyFrac'][:,:])

    # Check for NoData values
    lu[np.isclose(lu, lu_fv)] = np.nan
    p[np.isclose(p, p_fv)] = np.nan
    et[np.isclose(et, et_fv)] = np.nan
    qratio_y[np.isclose(qratio_y, qratio_y_fv)] = np.nan
    thetasat[np.isclose(thetasat, thetasat_fv)] = np.nan
    rootdepth[np.isclose(rootdepth, rootdepth_fv)] = np.nan
    interception[np.isclose(interception, interception_fv)] = np.nan
    SWSupplyFrac[np.isclose(SWSupplyFrac, SWSupplyFrac_fv)] = np.nan

### Create output NetCDF variables:

    #Blue ET (monthly)
    etb_var = out_nc.createVariable('ETBlue_M', 'f4',
                                   ('time_yyyymm', 'latitude', 'longitude'),
                                   fill_value=std_fv)
    etb_var.long_name = 'Monthly Blue ET'
    etb_var.units = 'mm/month'
    etb_var.grid_mapping = 'crs'
    
    #Green ET (monthly)
    etg_var = out_nc.createVariable('ETGreen_M', 'f4',
                                   ('time_yyyymm', 'latitude', 'longitude'),
                                   fill_value=std_fv)
    etg_var.long_name = 'Monthly Green ET'
    etg_var.units = 'mm/month'
    etg_var.grid_mapping = 'crs'   
    
   
    # Surface runoff (monthly)
    sr_var = out_nc.createVariable('SurfaceRunoff_M', 'f4',
                                   ('time_yyyymm', 'latitude', 'longitude'),
                                   fill_value=std_fv)
    sr_var.long_name = 'Surface runoff (fast)'
    sr_var.units = 'mm/month'
    sr_var.grid_mapping = 'crs'
    

    # Baseflow (monthly)
    bf_var = out_nc.createVariable('Baseflow_M', 'f4',
                                   ('time_yyyymm', 'latitude', 'longitude'),
                                   fill_value=std_fv)
    bf_var.long_name = 'Baseflow (slow)'
    bf_var.units = 'mm/month'
    bf_var.grid_mapping = 'crs'

    # Total runoff (monthly)
    tr_var = out_nc.createVariable('TotalRunoff_M', 'f4',
                                   ('time_yyyymm', 'latitude', 'longitude'),
                                   fill_value=std_fv)
    tr_var.long_name = 'Total runoff'
    tr_var.units = 'mm/month'
    tr_var.grid_mapping = 'crs'

    # Storage change - soil moisture (monthly)
    dsm_var = out_nc.createVariable('StorageChange_M', 'f4',
                                    ('time_yyyymm', 'latitude', 'longitude'),
                                    fill_value=std_fv)
    dsm_var.long_name = 'Change in soil moisture storage'
    dsm_var.units = 'mm/month'
    dsm_var.grid_mapping = 'crs'

    # Percolation (monthly)
    per_var = out_nc.createVariable('Percolation_M', 'f4',
                                    ('time_yyyymm', 'latitude', 'longitude'),
                                    fill_value=std_fv)
    per_var.long_name = 'Percolation'
    per_var.units = 'mm/month'
    per_var.grid_mapping = 'crs'

    # Supply (monthly)
    sup_var = out_nc.createVariable('Supply_M', 'f4',
                                    ('time_yyyymm', 'latitude', 'longitude'),
                                    fill_value=std_fv)
    sup_var.long_name = 'Supply'
    sup_var.units = 'mm/month'
    sup_var.grid_mapping = 'crs'
 
    
        # Supply (monthly)
    supsw_var = out_nc.createVariable('SupplySW_M', 'f4',
                                    ('time_yyyymm', 'latitude', 'longitude'),
                                    fill_value=std_fv)
    supsw_var.long_name = 'SupplySW'
    supsw_var.units = 'mm/month'
    supsw_var.grid_mapping = 'crs'

    
        # Supply (monthly)
    supgw_var = out_nc.createVariable('SupplyGW_M', 'f4',
                                    ('time_yyyymm', 'latitude', 'longitude'),
                                    fill_value=std_fv)
    supgw_var.long_name = 'SupplyGW'
    supgw_var.units = 'mm/month'
    supgw_var.grid_mapping = 'crs'

    
    # Root depth soil moisture (monthly)
    rdsm_var = out_nc.createVariable('RootDepthSoilMoisture_M', 'f4',
                                     ('time_yyyymm', 'latitude', 'longitude'),
                                     fill_value=rootdepth_fv)
    rdsm_var.long_name = 'Root depth soil moisture'
    rdsm_var.units = 'cm3/cm3'
    rdsm_var.grid_mapping = 'crs'
    
    # Incremental surface runoff (monthly)
    incsr_var = out_nc.createVariable('IncrementalRunoff_M', 'f4',
                                      ('time_yyyymm', 'latitude', 'longitude'),
                                      fill_value=std_fv)
    incsr_var.long_name = 'Incremental runoff'
    incsr_var.units = 'mm/month'
    incsr_var.grid_mapping = 'crs'
 
    # Incremental percolation (monthly)
    incper_var = out_nc.createVariable('IncrementalPercolation_M', 'f4',
                                       ('time_yyyymm',
                                        'latitude', 'longitude'),
                                       fill_value=std_fv)
    incper_var.long_name = 'Incremental Percolation'
    incper_var.units = 'mm/month'
    incper_var.grid_mapping = 'crs'
  
    
    # corrected rootdepth
    rootdepth_var = out_nc.createVariable('RootDepth_cor', 'f4',
                                           ('latitude', 'longitude'),
                                           fill_value=-9999)
    rootdepth_var.long_name = 'Root depth corrected'
    rootdepth_var.units = 'mm'
    rootdepth_var.grid_mapping = 'crs'
    
### Read total storage change    
    if dS_GRACE is not None:
        dfS=pd.read_csv(dS_GRACE,sep=';')
    else:
        d=dict()
        d['dS[mm]']=list(np.zeros(len(years_ls)))
        d['year']=years_ls
        dfS=pd.DataFrame(d)              
    
### Calculation
    # calibrating parameter rootdepth_par
    Rd = rootdepth*rootdepth_par
    rootdepth_var[:,:]=Rd
    #maximum storage in unsaturated zone
    SMmax=thetasat*Rd    
    sm_g=et*0
    sm_b=et*0
    f_consumed = Consumed_fraction(lu)
    for yyyy in years_ls:
        print '\tyear: {0}'.format(yyyy)
        yyyyi = years_ls.index(yyyy)
        ti1 = time_indeces[yyyy][0]
        ti2 = time_indeces[yyyy][-1] + 1
        for t in time_indeces[yyyy]:
            print '\tMonth: {0}'.format(time_ls[t])
### Step 0: Get data of the month and previous month
            P = p[t,:,:]
            ETa=et[t,:,:]
            I=interception[t,:,:]
            if t==0:
                SMgt_1=ETa*0
                Qsupply=ETa*0
                SMincrt_1=ETa*0
            else:
                SMgt_1=sm_g[t-1,:,:]  
                SMincrt_1=sm_b[t-1,:,:]            
### Step 1: Soil moisture
            SMg,SMincr,SM,dsm,Qsupply,ETincr,ETg=SM_bucket(P,ETa,I,SMmax,SMgt_1,SMincrt_1,f_consumed)    
### Step 2: SRO
            SRO,SROincr=SCS_calc_SRO(P,I,SMmax,SM,Qsupply)
### Step 3: Percolation
            Qperc=np.where((P+Qsupply)>(ETa+dsm+SRO),(P+Qsupply)-(ETa+dsm+SRO),0)
            Qperc_incr=np.where(Qsupply>(SROincr+ETincr),Qsupply-(SROincr+ETincr),0)
### Step 4: Split Qsupply 
            QsupplySW = Qsupply*SWSupplyFrac
            QsupplyGW = Qsupply - QsupplySW 
### Step 5: Store monthly data of the year
            sr_var[t,:,:]=SRO
            incsr_var[t,:,:]=SROincr
            dsm_var[t,:,:]=dsm
            per_var[t,:,:]=Qperc
            incper_var[t,:,:]=Qperc_incr            
            sup_var[t,:,:]=Qsupply
            supsw_var[t,:,:]=QsupplySW
            supgw_var[t,:,:]=QsupplyGW
            rootdepth_var[:,:]=Rd
            etb_var[t,:,:]=ETincr
            etg_var[t,:,:]=ETg
### Step 6: Estimate qratio and BF from yearly SRO
        dS=float(dfS.loc[dfS['year']==yyyy]['dS[mm]'])
        P=p[ti1:ti2,:,:]
        ETa=et[ti1:ti2,:,:]
        SRO=sr_var[ti1:ti2,:,:]
        QsupplySW=supsw_var[ti1:ti2,:,:]        
        BFratio=BFratio_y(P,ETa,QsupplySW,SRO,dS) 
        print('BF/SRO: {0}'.format(BFratio))
        BF=SRO*BFratio
        TR=SRO+BF
        bf_var[ti1:ti2,:,:]=BF
        tr_var[ti1:ti2,:,:]=TR        
        
### Finish
    print 'Closing netcdf...'
    out_nc.close()
    ended=dt.datetime.now()
    print 'Time elapsed: {0}'.format(ended-started)
#%% Functions        
def BFratio_y(P,ETa,Qsupply_sw,SRO,dS):
    '''
    P, ETa, SRO, Qsupply_sw (12xNxM) (mm)
    dS (1) [mm/year]
    '''            
    Pavg=np.nanmean(12*np.nanmean(P,axis=0))
    ETavg=np.nanmean(12*np.nanmean(ETa,axis=0))
    SROavg=np.nanmean(12*np.nanmean(SRO,axis=0))
    Supply_SWavg=np.nanmean(12*np.nanmean(Qsupply_sw,axis=0))
    BFratio=((Pavg-ETavg-dS)+Supply_SWavg-SROavg)/SROavg
    print 'check shape P: {0}, ET: {1}, Qsro: {2}, Qsupply_sw: {3}'.format(P.shape,ETa.shape,SRO.shape,Qsupply_sw.shape)
    print 'P: {0}, ET: {1}, dS: {2}, Qsupply_sw: {3}, Qsro: {4}'.format(Pavg,ETavg,dS,Supply_SWavg,SROavg)
    return BFratio            
            
def SM_bucket(P,ETa,I,SMmax,SMgt_1,SMincrt_1,f_consumed):
    SMt_1=SMgt_1+SMincrt_1
    ETr=ETa-I
    SMtemp=SMgt_1+np.maximum(P-I,P*0)
    SMg=np.where(SMtemp-ETr>0,SMtemp-ETr,0)
    ETincr=np.where(SMtemp-ETr>0,0,ETr-SMtemp)
    Qsupply=np.where(SMtemp-ETr>0,0,(ETincr-SMincrt_1)/f_consumed)            
    SMincr=np.where(ETincr>SMincrt_1,SMincrt_1+Qsupply-ETincr,SMincrt_1-ETincr)
    SMsum=SMg+SMincr
    SM=np.where(SMsum>SMmax,SMmax,SMsum)
    SMg=np.where(SMsum>SMmax,SMg-(SMmax-SMsum)*(SMg/SMsum),SMg)
    SMincr=np.where(SMsum>SMmax,SMincr-(SMmax-SMsum)*(SMincr/SMsum),SMincr)
    dsm=SM-SMt_1
    ETg=ETa-ETincr
    return SMg,SMincr,SM,dsm,Qsupply,ETincr,ETg
        
def SCS_calc_SRO(P,I,SMmax,SM,Qsupply):    
    SM=np.where(SM>SMmax,SMmax,SM)        
    SRO= ((P-I+Qsupply)**2)/(P-I+Qsupply+(SMmax-SM))
    SROg= ((P-I)**2)/(P-I+(SMmax-SM))
    SROincr=SRO-SROg
    return SRO,SROincr

def Consumed_fraction(lu):
    f_consumed=np.copy(lu)
    consumed_fractions, sw_supply_fractions, sw_return_fractions=gd.get_sheet4_6_fractions(version = '1.0')
    lu_categories = gd.get_sheet4_6_classes(version = '1.0')
    for key in consumed_fractions.keys():
        classes = lu_categories[key]
        mask = np.logical_or.reduce([lu == value for value in classes])
        consumed_fraction = consumed_fractions[key]
        f_consumed[mask] = consumed_fraction
    return f_consumed            

   
