# -*- coding: utf-8 -*-
"""
Created on Mon Oct 03 13:36:07 2016

@author: Bert Coerver, b.coerver@unesco-ihe.org
"""

def get_lulcs(lulc_version = '4.0'):
    
    lulc = dict()
        
    lulc['4.0'] = {
    'legend': ['Code', 'Landuse', 'Description', 'Beneficial T [%]', 'Beneficial E [%]', 'Beneficial I [%]', 'Agriculture [%]', 'Environment [%]', 'Economic [%]', 'Energy [%]', 'Leisure [%]'],
    0: ['X','X','X',0.,0.,0.,0.,0.,0.,0.,0.],
    1: ['PLU1', 'Protected', 'Protected forests', 100.0, 100.0, 0.0, 0.0, 85.0, 0.0, 0.0, 15.0], 
    2: ['PLU2', 'Protected', 'Protected shrubland', 100.0, 100.0, 0.0, 0.0, 85.0, 0.0, 0.0, 15.0], 
    3: ['PLU3', 'Protected', 'Protected natural grasslands', 100.0, 100.0, 0.0, 0.0, 85.0, 0.0, 0.0, 15.0], 
    4: ['PLU4', 'Protected', 'Protected natural waterbodies', 100.0, 100.0, 0.0, 0.0, 85.0, 0.0, 0.0, 15.0], 
    5: ['PLU5', 'Protected', 'Protected wetlands', 100.0, 100.0, 0.0, 0.0, 85.0, 0.0, 0.0, 15.0], 
    6: ['PLU6', 'Protected', 'Glaciers', 0.0, 100.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0], 
    7: ['PLU7', 'Protected', 'Protected other', 100.0, 100.0, 0.0, 0.0, 85.0, 0.0, 0.0, 15.0], 
    8: ['ULU1', 'Utilized', 'Closed deciduous forest', 100.0, 0.0, 0.0, 5.0, 90.0, 0.0, 0.0, 5.0], 
    9: ['ULU2', 'Utilized', 'Open deciduous forest', 100.0, 0.0, 0.0, 5.0, 90.0, 0.0, 0.0, 5.0], 
    10: ['ULU3', 'Utilized', 'Closed evergreen forest', 100.0, 0.0, 0.0, 5.0, 90.0, 0.0, 0.0, 5.0], 
    11: ['ULU4', 'Utilized', 'Open evergreen forest', 100.0, 0.0, 0.0, 5.0, 90.0, 0.0, 0.0, 5.0], 
    12: ['ULU5', 'Utilized', 'Closed savanna', 100.0, 0.0, 0.0, 5.0, 80.0, 0.0, 10.0, 5.0], 
    13: ['ULU6', 'Utilized', 'Open savanna', 100.0, 0.0, 0.0, 10.0, 80.0, 0.0, 5.0, 5.0], 
    14: ['ULU7', 'Utilized', 'Shrub land & mesquite', 100.0, 0.0, 0.0, 5.0, 85.0, 0.0, 10.0, 0.0], 
    15: ['ULU8', 'Utilized', ' Herbaceous cover', 100.0, 0.0, 0.0, 5.0, 95.0, 0.0, 0.0, 0.0], 
    16: ['ULU9', 'Utilized', 'Meadows & open grassland', 100.0, 0.0, 0.0, 60.0, 30.0, 0.0, 0.0, 10.0], 
    17: ['ULU10', 'Utilized', 'Riparian corridor', 100.0, 0.0, 0.0, 10.0, 60.0, 10.0, 0.0, 20.0], 
    18: ['ULU11', 'Utilized', 'Deserts', 100.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0], 
    19: ['ULU12', 'Utilized', 'Wadis', 100.0, 0.0, 0.0, 15.0, 80.0, 0.0, 0.0, 5.0], 
    20: ['ULU13', 'Utilized', 'Natural alpine pastures', 100.0, 0.0, 0.0, 70.0, 20.0, 0.0, 0.0, 10.0], 
    21: ['ULU14', 'Utilized', 'Rocks & gravel & stones & boulders', 100.0, 0.0, 0.0, 0.0, 95.0, 0.0, 0.0, 5.0], 
    22: ['ULU15', 'Utilized', 'Permafrosts', 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0], 
    23: ['ULU16', 'Utilized', 'Brooks & rivers & waterfalls', 0.0, 50.0, 0.0, 25.0, 55.0, 5.0, 0.0, 15.0], 
    24: ['ULU17', 'Utilized', 'Natural lakes\xa0', 0.0, 50.0, 0.0, 25.0, 40.0, 5.0, 0.0, 30.0], 
    25: ['ULU18', 'Utilized', 'Flood plains & mudflats', 100.0, 50.0, 0.0, 40.0, 60.0, 0.0, 0.0, 0.0], 
    26: ['ULU19', 'Utilized', 'Saline sinks & playas & salinized soil', 100.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0], 
    27: ['ULU20', 'Utilized', 'Bare soil', 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0], 
    28: ['ULU21', 'Utilized', 'Waste land', 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0], 
    29: ['ULU22', 'Utilized', 'Moorland', 100.0, 0.0, 0.0, 5.0, 80.0, 0.0, 0.0, 15.0], 
    30: ['ULU23', 'Utilized', 'Wetland', 100.0, 50.0, 0.0, 5.0, 80.0, 0.0, 5.0, 10.0], 
    31: ['ULU24', 'Utilized', 'Mangroves', 100.0, 50.0, 0.0, 5.0, 80.0, 0.0, 5.0, 10.0], 
    32: ['ULU25', 'Utilized', 'Alien invasive species', 0.0, 0.0, 0.0, 0.0, 60.0, 0.0, 10.0, 30.0], 
    33: ['MLU1', 'Modified', 'Forest plantations', 100.0, 0.0, 0.0, 90.0, 0.0, 10.0, 0.0, 0.0], 
    34: ['MLU2', 'Modified', 'Rainfed production pastures', 100.0, 0.0, 0.0, 90.0, 0.0, 10.0, 0.0, 0.0], 
    35: ['MLU3', 'Modified', 'Rainfed crops - cereals', 100.0, 0.0, 0.0, 90.0, 0.0, 10.0, 0.0, 0.0], 
    36: ['MLU4', 'Modified', 'Rainfed crops - root/tuber', 100.0, 0.0, 0.0, 90.0, 0.0, 10.0, 0.0, 0.0], 
    37: ['MLU5', 'Modified', 'Rainfed crops - legumious', 100.0, 0.0, 0.0, 90.0, 0.0, 10.0, 0.0, 0.0], 
    38: ['MLU6', 'Modified', 'Rainfed crops - sugar', 100.0, 0.0, 0.0, 90.0, 0.0, 10.0, 0.0, 0.0], 
    39: ['MLU7', 'Modified', 'Rainfed crops - fruit and nuts', 100.0, 0.0, 0.0, 90.0, 0.0, 10.0, 0.0, 0.0], 
    40: ['MLU8', 'Modified', 'Rainfed crops - vegetables and melons', 100.0, 0.0, 0.0, 90.0, 0.0, 10.0, 0.0, 0.0], 
    41: ['MLU9', 'Modified', 'Rainfed crops - oilseed', 100.0, 0.0, 0.0, 45.0, 0.0, 15.0, 40.0, 0.0], 
    42: ['MLU10', 'Modified', 'Rainfed crops - beverage and spice', 100.0, 0.0, 0.0, 90.0, 0.0, 10.0, 0.0, 0.0], 
    43: ['MLU11', 'Modified', 'Rainfed crops - other ', 100.0, 0.0, 0.0, 80.0, 0.0, 20.0, 0.0, 0.0], 
    44: ['MLU12', 'Modified', 'Mixed species agro-forestry', 100.0, 0.0, 0.0, 90.0, 0.0, 10.0, 0.0, 0.0], 
    45: ['MLU13', 'Modified', 'Fallow & idle land', 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0], 
    46: ['MLU14', 'Modified', 'Dump sites & deposits', 0.0, 0.0, 0.0, 0.0, 60.0, 40.0, 0.0, 0.0], 
    47: ['MLU15', 'Modified', 'Rainfed homesteads and gardens (urban cities) - outdoor', 100.0, 0.0, 0.0, 0.0, 0.0, 35.0, 0.0, 65.0], 
    48: ['MLU16', 'Modified', 'Rainfed homesteads and gardens (rural villages) - outdoor', 100.0, 0.0, 0.0, 0.0, 0.0, 35.0, 0.0, 65.0], 
    49: ['MLU17', 'Modified', 'Rainfed industry parks - outdoor', 100.0, 0.0, 0.0, 0.0, 0.0, 50.0, 0.0, 50.0], 
    50: ['MLU18', 'Modified', 'Rainfed parks (leisure & sports)', 100.0, 0.0, 0.0, 0.0, 15.0, 0.0, 0.0, 85.0], 
    51: ['MLU19', 'Modified', 'Rural paved surfaces (lots, roads, lanes)', 100.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0], 
    52: ['MWU1', 'Managed', 'Irrigated forest plantations', 100.0, 0.0, 0.0, 90.0, 0.0, 10.0, 0.0, 0.0], 
    53: ['MWU2', 'Managed', 'Irrigated production pastures', 100.0, 0.0, 0.0, 90.0, 0.0, 10.0, 0.0, 0.0], 
    54: ['MWU3', 'Managed', 'Irrigated crops - cereals', 100.0, 0.0, 0.0, 90.0, 0.0, 10.0, 0.0, 0.0], 
    55: ['MWU4', 'Managed', 'Irrigated crops - root/tubers', 100.0, 0.0, 0.0, 90.0, 0.0, 10.0, 0.0, 0.0], 
    56: ['MWU5', 'Managed', 'Irrigated crops - legumious', 100.0, 0.0, 0.0, 90.0, 0.0, 10.0, 0.0, 0.0], 
    57: ['MWU6', 'Managed', 'Irrigated crops - sugar', 100.0, 0.0, 0.0, 90.0, 0.0, 10.0, 0.0, 0.0], 
    58: ['MWU7', 'Managed', 'Irrigated crops - fruit and nuts', 100.0, 0.0, 0.0, 90.0, 0.0, 10.0, 0.0, 0.0], 
    59: ['MWU8', 'Managed', 'Irrigated crops - vegetables and melons', 100.0, 0.0, 0.0, 90.0, 0.0, 10.0, 0.0, 0.0], 
    60: ['MWU9', 'Managed', 'Irrigated crops - Oilseed', 100.0, 0.0, 0.0, 65.0, 0.0, 10.0, 25.0, 0.0], 
    61: ['MWU10', 'Managed', 'Irrigated crops - beverage and spice', 100.0, 0.0, 0.0, 90.0, 0.0, 10.0, 0.0, 0.0], 
    62: ['MWU11', 'Managed', 'Irrigated crops - other', 100.0, 0.0, 0.0, 80.0, 0.0, 20.0, 0.0, 0.0], 
    63: ['MWU12', 'Managed', 'Managed water bodies (reservoirs, canals, harbors, tanks)', 0.0, 100.0, 0.0, 35.0, 5.0, 30.0, 20.0, 10.0], 
    64: ['MWU13', 'Managed', 'Greenhouses - indoor', 100.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0], 
    65: ['MWU14', 'Managed', 'Aquaculture', 0.0, 100.0, 0.0, 90.0, 0.0, 10.0, 0.0, 0.0], 
    66: ['MWU15', 'Managed', 'Domestic households - indoor (sanitation)', 0.0, 100.0, 0.0, 0.0, 0.0, 35.0, 0.0, 65.0], 
    67: ['MWU16', 'Managed', 'Manufacturing & commercial industry - indoor', 0.0, 100.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0], 
    68: ['MWU17', 'Managed', 'Irrigated homesteads and gardens (urban cities) - outdoor', 100.0, 0.0, 0.0, 30.0, 5.0, 15.0, 0.0, 50.0], 
    69: ['MWU18', 'Managed', 'Irrigated homesteads and gardens (rural villages) - outdoor', 100.0, 0.0, 0.0, 30.0, 5.0, 15.0, 0.0, 50.0], 
    70: ['MWU19', 'Managed', 'Irrigated industry parks - outdoor', 100.0, 0.0, 0.0, 0.0, 15.0, 35.0, 0.0, 50.0], 
    71: ['MWU20', 'Managed', 'Irrigated parks (leisure, sports)', 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0], 
    72: ['MWU21', 'Managed', 'Urban paved Surface (lots, roads, lanes)', 100.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0], 
    73: ['MWU22', 'Managed', 'Livestock and domestic husbandry', 100.0, 0.0, 0.0, 90.0, 0.0, 10.0, 0.0, 0.0], 
    74: ['MWU23', 'Managed', 'Managed wetlands & swamps', 100.0, 50.0, 0.0, 0.0, 65.0, 10.0, 0.0, 25.0], 
    75: ['MWU24', 'Managed', 'Managed other inundation areas', 100.0, 50.0, 0.0, 0.0, 55.0, 20.0, 0.0, 25.0], 
    76: ['MWU25', 'Managed', 'Mining/ quarry & shale exploiration', 100.0, 50.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0], 
    77: ['MWU26', 'Managed', 'Evaporation ponds', 0.0, 100.0, 0.0, 0.0, 75.0, 25.0, 0.0, 0.0], 
    78: ['MWU27', 'Managed', 'Waste water treatment plants', 0.0, 100.0, 0.0, 0.0, 55.0, 45.0, 0.0, 0.0], 
    79: ['MWU28', 'Managed', 'Hydropower plants', 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 97.5, 2.5], 
    80: ['MWU29', 'Managed', 'Thermal power plants', 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0]
    }
    
    return lulc[lulc_version]


def get_bluegreen_classes(version = '1.0'):
    
    gb_cats = dict()
    mvg_avg_len = dict()
    
    gb_cats['1.0'] = {
    'crops':                [53,54,55,56,57,58,59,60,61,62, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 33, 44],
    'perennial crops':      [52],
    'savanna':              [12,13],
    'water':                [63,74,75,77,4, 19, 23, 24],
    'forests':              [1, 8, 9, 10, 11, 17],
    'grass':                [3, 16, 20, 2, 14, 15],
    'other':                [68,69,70,71,72,76,78,73,67,65,66,64,79,80,6, 7, 18, 21, 22, 26, 27, 28, 29, 32, 45, 46, 47, 48, 49, 50, 51, 5, 25, 30, 31],
    }

    mvg_avg_len['1.0'] = {
    'crops':                2,
    'perennial crops':      3,
    'savanna':              4,
    'water':                1,
    'forests':              5,
    'grass':                1,
    'other':                1,
    }
    
    return gb_cats[version], mvg_avg_len[version]
    
def get_sheet4_6_classes(version = '1.0'):
    lucs = dict()
    
    lucs['1.0'] = {
    'Forests':              [1, 8, 9, 10, 11, 17],
    'Shrubland':            [2, 12, 14, 15],
    'Rainfed Crops':        [34, 35, 36, 37, 38, 39, 40, 41, 42, 43],
    'Forest Plantations':   [33, 44],
    'Natural Water Bodies': [4, 19, 23, 24],
    'Wetlands':             [5, 25, 30, 31],
    'Natural Grasslands':   [3, 13, 16, 20],
    'Other (Non-Manmade)':  [6, 7, 18, 21, 22, 26, 27, 28, 29, 32, 45, 46, 48, 49, 50, 51],
    'Irrigated crops':      [52,53,54,55,56,57,58,59,60,61,62],
    'Managed water bodies': [63,74,75,77],
    'Aquaculture':          [65],
    'Residential':          [47, 66, 68, 72],
    'Greenhouses':          [64],
    'Other':                [68,69,70,71,76,78]}
    
    return lucs[version]

def get_sheet4_6_fractions(version = '1.0'):
    consumed_fractions = dict()
    sw_supply_fractions = dict()
    sw_return_fractions = dict()
    
    consumed_fractions['1.0'] = {
    'Forests':              1.00,
    'Shrubland':            1.00,
    'Rainfed Crops':        1.00,
    'Forest Plantations':   1.00,
    'Natural Water Bodies': 0.15,
    'Wetlands':             0.15,
    'Natural Grasslands':   0.70,
    'Other (Non-Manmade)':  1.0,#0.4
    'Irrigated crops':      0.80,#0.8
    'Managed water bodies': 0.4,#0.4
    'Other':                0.40,#0.4
    'Residential':          1.00,
    'Greenhouses':          0.95,
    'Aquaculture':          0.20}
    
    sw_supply_fractions['1.0'] = {
    'Forests':              0.005,
    'Shrubland':            0.10,
    'Rainfed Crops':        0.05,
    'Forest Plantations':   0.005,
    'Natural Water Bodies': 0.95,
    'Wetlands':             0.95,
    'Natural Grasslands':   0.30,
    'Other (Non-Manmade)':  0.50,
    'Irrigated crops':      9999,
    'Managed water bodies': 0.95,
    'Other':                0.50,
    'Residential':          0.90,
    'Greenhouses':          0.50,
    'Aquaculture':          0.95}
    
    sw_return_fractions['1.0'] = {
    'Forests':              9999,
    'Shrubland':            9999,
    'Rainfed Crops':        9999,
    'Forest Plantations':   9999,
    'Natural Water Bodies': 0.95,
    'Wetlands':             0.95,
    'Natural Grasslands':   0.10,
    'Other (Non-Manmade)':  0.50,
    'Irrigated crops':      0.50,
    'Managed water bodies': 0.95,
    'Other':                0.50,
    'Residential':          0.70,
    'Greenhouses':          0.50,
    'Aquaculture':          0.95}
    
    return consumed_fractions[version], sw_supply_fractions[version], sw_return_fractions[version]
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    