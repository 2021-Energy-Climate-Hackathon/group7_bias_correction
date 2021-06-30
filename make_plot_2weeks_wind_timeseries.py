""" Make a plot of 1 month of solar, wind and T2m data in ERA5 and CPDN data """

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from read_gridded_data_functions import load_country_weather_data,read_spatial_dimensions
from matplotlib import gridspec

month = 3
country_names = ['United Kingdom','Ireland']
variable_names = ['t2m','u10','ssrd']
#title = [r'2m Temperature January 2018 $(^{\circ}C)$',r'10m wind speed January 2018 $(ms^{-1})$',r'Incoming solar radiation January 2018 $(Wm^{-2})$']
title = [r'10m wind speed $(^{\circ}C)$']
ylim = [ [-2,11], [2,12],[0,200] ]
var_name = 'u10'

# set up figure for plotting
plt.figure(figsize=(16,11))
gs = gridspec.GridSpec(2,2)

def convert_to_windpower(wind_speed_data,power_curve_file):

    """
    This function takes the ERA5 reanalysis data, loads it and applied a
    country mask (ready for conversion to energy) it then returns
    the array (of original size) with all irrelvelant gridpoints
    set to zeros.

    You will need the shpreader.natural_earth data downloaded
    to find the shapefiles.

    Args:

        gridded_wind_power (array): wind power capacity factor data, dimensions
            [time,lat,lon]. Capacity factors range between 0 and 1.

        power_curve_file (str): The filename of a .csv file
            containing the wind speeds (column 0) and capacity factors
            (column 2) of the chosen wind turbine.

    Returns:

        wind_power_cf (array): Gridded wind Power capacity factor
            data, dimensions [time,lat,lon]. Values vary between 0 and 1.

    """

    # first load in the power curve data
    pc_w = []
    pc_p = []

    with open(power_curve_file) as f:
        for line in f:
            columns = line.split()
            #print columns[0]
            pc_p.append(np.float(columns[1]))
            pc_w.append(np.float(columns[0]))  # get power curve output (CF)

    # convert to an array
    power_curve_w = np.array(pc_w)
    power_curve_p = np.array(pc_p)

    #interpolate to fine resolution.
    pc_winds = np.linspace(0,50,501) # make it finer resolution
    pc_power = np.interp(pc_winds,power_curve_w,power_curve_p)

    reshaped_speed = wind_speed_data.flatten()
    test = np.digitize(reshaped_speed,pc_winds,right=False) # indexing starts
    #from 1 so needs -1: 0 in the next bit to start from the lowest bin.
    test[test ==len(pc_winds)] = 500 # make sure the bins don't go off the
    #end (power is zero by then anyway)
    wind_power_flattened = 0.5*(pc_power[test-1]+pc_power[test])

    wind_power_cf = np.reshape(wind_power_flattened,(np.shape(wind_speed_data)))

    return(wind_power_cf)


def make_plot(t,variable,ylim,title='',xlabel='',ylabel='',fontsize=18,color='r'):
    print(variable.shape)
    plt.title(title,fontsize=fontsize)
    plt.plot(t,variable,color=color,linewidth=2)
    plt.xlabel(xlabel,fontsize=fontsize)
    plt.ylabel('Wind capacity factor',fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    #plt.xlim(xlim[0],xlim[1])
    #print(ylim)
    plt.ylim(0,1)
    plt.xlim(0,15)


for i, country_name in enumerate(country_names): # np.arange(8):    

    power_curve_file = '/gws/pw/j05/cop26_hackathons/oxford/Group_folders/group_7/energy_conversion/group_7/inputs/power_onshore.csv'

    # read in ERA5 data
    ERA5_dir = '/gws/pw/j05/cop26_hackathons/oxford/Data/ERA5_data_EU_domain/field_set_1/'
    file_name = 'ERA5_1hr_field_set_1_2018_'+str(month).zfill(2)+'.nc'
    var_ERA5, country_mask_ERA5 = load_country_weather_data(country_name,ERA5_dir,file_name,var_name)

    # if looking at wind speed, read in v wind too
    if var_name == 'u10':
        V10_ERA5, _ = load_country_weather_data(country_name,ERA5_dir,file_name,'v10')
        var_ERA5 = np.sqrt(var_ERA5**2 + V10_ERA5**2)

    # calculate hourly and monthly country means
    ERA5_var_mean = np.nanmean(np.nanmean(var_ERA5,axis=2),axis=1)
    ERA5_var_month = np.mean(ERA5_var_mean)
    ERA5_wp = convert_to_windpower(ERA5_var_mean,power_curve_file)
    
    # Read in CPDN
    CPDN_dir = '/gws/pw/j05/cop26_hackathons/oxford/Group_folders/group_7/cpdn_example_data/remap_ERA5/'
    if var_name == 't2m':
        file_name = 'tas_1hrly_mean_a000_2018-01_2018-12.nc'
        var_name_CPDN = 'item3236_1hrly_mean'
    elif var_name == 'u10':
        file_name = '10mwd_1hrly_mean_a000_2018-01_2018-12.nc'
        var_name_CPDN = 'item3249_1hrly_mean'
    elif var_name == 'ssrd':
        file_name = 'rsds_1hrly_mean_a000_2018-01_2018-12.nc'
        var_name_CPDN = 'field203'

    if var_name == 'ssrd': time_name='time0'
    else: time_name = 'time1'
    var_CPDN, country_mask_CPDN = load_country_weather_data(country_name,CPDN_dir,file_name,var_name_CPDN,chosen_months=(month),time_name=time_name)
    if var_name == 't2m': var_CPDN = var_CPDN - 273.15 # convert to deg Celsius

    # calculate hourly and monthly country means
    CPDN_var_mean = np.nanmean(np.nanmean(var_CPDN[:,0,:,:],axis=2),axis=1)
    CPDN_var_month = np.mean(CPDN_var_mean)
    CPDN_wp = convert_to_windpower(CPDN_var_mean,power_curve_file)
    
    title = [country_name+' ERA5 wind capacity factor',country_name+' CPDN wind capacity factor']

    ax = plt.subplot(gs[i,0])
    #if i==0: plt.text(15,12,'ERA5',fontsize=25) #,transform=ax.transAxes)
    t = np.linspace(0,31,744)
    if country_name == 'Ireland': xlabel = 'Day'
    else: xlabel = ''
    make_plot(t,ERA5_wp,ylim[i],title=title[0],color='royalblue',xlabel=xlabel)

    plt.subplot(gs[i,1])
    #if i==0: plt.text(15,12,'CPDN',fontsize=25)
    if country_name == 'Ireland': 
        xlabel = 'Day'
    else: xlabel = ''
    t = np.linspace(0,30,720)
    make_plot(t,CPDN_wp,ylim[i],title=title[1],color='red',xlabel=xlabel)


figure_name = 'daily_variability_Ire_GB_10m_wind_speed.png'
print('saving to %s' % (figure_name))
plt.savefig(figure_name,bbox_inches='tight')
                                               
