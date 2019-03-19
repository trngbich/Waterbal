# -*- coding: utf-8 -*-
"""
Authors: Elga Salvadore
         IHE Delft 2019
Contact: e.salvadore@un-ihe.org
Repository: 
Module: wa_wb

Simplified version of WaterPix
"""

from __future__ import division
import numpy as np
from scipy.optimize import fsolve
import netCDF4
import get_dictionaries as gd
import becgis
import os
from scipy import interpolate
from wa_wb import davgis

def interpolate_nan(array2d,method='cubic',positive=True):   
    '''
    fill nan values in 2-d array by interpolating
    '''
    x=np.arange(0,array2d.shape[1])
    y=np.arange(0,array2d.shape[0])
    array2d = np.ma.masked_invalid(array2d)
    xx, yy = np.meshgrid(x, y)
    #get only the valid values
    x1 = xx[~array2d.mask]
    y1 = yy[~array2d.mask]
    newarr = array2d[~array2d.mask]
    
    filled = interpolate.griddata((x1, y1), newarr.ravel(),(xx, yy),method=method)
    if positive:
        filled=np.where(filled<0,0,filled)
    return filled

def lai_and_soil_calculations(thetasat, lai, swi, swio, swix, rootdepth):
    '''
    Calculate thetasat, lai, and soil moisture parameters
    '''
    # Soil moisture calculations
    theta0 = thetasat * swi/100.0
    thetao= thetasat * swio/100.0
    thetax = thetasat * swix/100.0
    # Argument for exponential function
    exp_arg = exp_arg_func(theta0,thetasat,lai)
    exp_argo = exp_arg_func(thetao,thetasat,lai)
    exp_argx = exp_arg_func(thetax,thetasat,lai)
    # Soil moisture values - root zone
    thetarz = (0.1*lai +(1-0.1*lai) * (exp_arg))*thetasat
    thetarzo = (0.1*lai +(1-0.1*lai) * (exp_argo))*thetasat
    thetarzx = (0.1*lai +(1-0.1*lai) * (exp_argx))*thetasat
    # Change in storage
    dsm = rootdepth*(thetarzx - thetarzo)
    return thetarz, thetarzo, thetarzx, dsm


def exp_arg_func(theta0, thetasat, lai):
    '''
    Compute the exponential argument in the calculation of the root depth
    soil moisture.
    '''
    value = 1 - np.exp((theta0/thetasat) * (-0.5*lai - 1))
    return value

def SCS_surface_runoff(p, supply, interception, rootdepth, thetasat, thetarz):
    P_I_supply = p-interception+supply
    Qsw = P_I_supply**2 / (P_I_supply + rootdepth * (thetasat - thetarz))
    return Qsw


def SCS_surface_runoff_gr(p, interception, rootdepth, thetasat, thetarz):
    P_I = p-interception
    Qsw_gr = P_I**2 / (P_I + rootdepth * (thetasat - thetarz))
    return Qsw_gr


def baseflow_calculation(Qsw, filter_par, qratio_y):
    '''
    Calculate the baseflow using the runoff ratio and the surface runoff
    '''
    Qgw_tot = np.nansum(Qsw/qratio_y-Qsw)
    q0 = fsolve(baseflow_function, 0.0, [Qsw, filter_par,
                                         qratio_y, Qgw_tot, True])
    
    Qgw = baseflow_function(q0, [Qsw, filter_par,
                                        qratio_y, Qgw_tot, False])
    # fixing negative values when Qsw is 0
    Qgw[Qgw<0]=0
    
    return Qgw



def baseflow_function(q0, args):
    '''
    Balance baseflow values in a yearly basis
    '''
    Qsw, filter_par, qratio, Qgw_tot, return_error = args
    q_temp = np.zeros(np.shape(Qsw)) * np.nan
    q_temp[0] = filter_par*q0 + 0.5*(1 + filter_par)*(
        Qsw[0] - Qsw[-1])
    for i in range(1, 12):
        q_temp[i] = filter_par*q_temp[i-1] + 0.5*(1 + filter_par)*(
            Qsw[i] - Qsw[i-1])
    # Baseflow
    Qgw = (1-qratio)/qratio*(Qsw - q_temp)
    # Error calculation
    error = Qgw_tot - np.nansum(Qgw)
    # Return error or vector
    if return_error:
        return error
    else:
        return Qgw


def total_supply(etb, lu):
    """
    Apply a consumed fraction to groups of landuse classes (categories) to aqcuire 
    maps of blue water supplies.
    """
    consumed_fractions, sw_supply_fractions, sw_return_fractions=gd.get_sheet4_6_fractions(version = '1.0')
    lu_categories = gd.get_sheet4_6_classes(version = '1.0')
    supply = np.copy(etb)
   
    for key in consumed_fractions.keys():
        classes = lu_categories[key]
        mask = np.logical_or.reduce([lu == value for value in classes])
        consumed_fraction = consumed_fractions[key]
        for i in range(np.shape(supply)[0]):
            supply[i,mask] /= consumed_fraction
    return supply



def output_nc_to_tiffs(output_nc, output_path):
    """
    Create raster files from the variables in the output netcdf file
    """
    # Output folders
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    path_y = os.path.join(output_path, 'yearly')
    path_m = os.path.join(output_path, 'monthly')
    path_a = os.path.join(output_path, 'additional')
    if not os.path.isdir(path_y):
        os.mkdir(path_y)
    if not os.path.isdir(path_m):
        os.mkdir(path_m)
    if not os.path.isdir(path_a):
        os.mkdir(path_a)

    # Read netcdf file
    nc_file = netCDF4.Dataset(output_nc, 'r')
    variables_ls = nc_file.variables.keys()
    time_y = nc_file.variables['time_yyyy'][:]
    time_m = nc_file.variables['time_yyyymm'][:]
    nc_file.close()

    # Remove variables
    for variable in ['latitude', 'longitude', 'time_yyyy', 'time_yyyymm', 'crs']:
        variables_ls.remove(variable)

    # Add sub-folders
    for variable in variables_ls:
        if '_Y' in variable:
            if not os.path.exists(os.path.join(path_y, variable)):
                os.mkdir(os.path.join(path_y, variable))
        elif '_M' in variable:
            if not os.path.exists(os.path.join(path_m, variable)):
                os.mkdir(os.path.join(path_m, variable))
        else:
            if not os.path.exists(os.path.join(path_a, variable)):
                os.mkdir(os.path.join(path_a, variable))

    # Main Loop
    for variable in variables_ls:
        # Yearly rasters
        if '_Y' in variable:
            for time in time_y:
                print '{0}\t{1}'.format(variable, time)
                file_name = variable[:-1] + '{0}.tif'.format(time)
                output_tiff = os.path.join(path_y, variable, file_name)
                davgis.NetCDF_to_Raster(output_nc, output_tiff, variable,
                                 x_variable='longitude',
                                 y_variable='latitude',
                                 time={'variable': 'time_yyyy',
                                       'value': time})
        # Monthly rasters
        elif '_M' in variable:
            for time in time_m:
                print '{0}\t{1}'.format(variable, time)
                file_name = variable[:-1] + '{0}.tif'.format(time)
                output_tiff = os.path.join(path_m, variable, file_name)
                print output_tiff
                davgis.NetCDF_to_Raster(output_nc, output_tiff, variable,
                                 x_variable='longitude',
                                 y_variable='latitude',
                                 time={'variable': 'time_yyyymm',
                                       'value': time})
        # Additional rasters
        else:
            print '{0}'.format(variable)
            file_name = variable[:-1] + '.tif'
            output_tiff = os.path.join(path_a, variable, file_name)
            davgis.NetCDF_to_Raster(output_nc, output_tiff, variable,
                             x_variable='longitude',
                             y_variable='latitude')
            
        




