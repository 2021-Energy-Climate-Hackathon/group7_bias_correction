""" Make a plot of 1 month of solar, wind and T2m data in ERA5 and CPDN data """

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from read_gridded_data_functions import load_country_weather_data,read_spatial_dimensions
from matplotlib import gridspec

month = 1
country_name = 'United Kingdom'
variable_names = ['t2m','u10','ssrd']
title = [r'2m Temperature January 2018 $(^{\circ}C)$',r'10m wind speed January 2018 $(ms^{-1})$',r'Incoming solar radiation January 2018 $(Wm^{-2})$']
ylim = [ [-2,11], [2,12],[0,200] ]

# set up figure for plotting
plt.figure(figsize=(18,15))
gs = gridspec.GridSpec(3,2)

def make_plot(t,variable,ylim,title='',xlabel='',ylabel='',fontsize=20,color='r'):
    print(variable.shape)
    plt.title(title,fontsize=fontsize)
    plt.plot(t,variable,color=color)
    plt.xlabel(xlabel,fontsize=fontsize)
    plt.ylabel(ylabel,fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    #plt.xlim(xlim[0],xlim[1])
    #print(ylim)
    plt.ylim(ylim[0],ylim[1])
    plt.xlim(0,30)


for i, var_name in enumerate(variable_names): # np.arange(8):    

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
    

    ax = plt.subplot(gs[i,0])
    if i==0: plt.text(13,14,'ERA5',fontsize=25) #,transform=ax.transAxes)
    t = np.linspace(0,31,744)
    if var_name == 'ssrd': xlabel = 'Day'
    else: xlabel = ''
    make_plot(t,ERA5_var_mean,ylim[i],title=title[i],color='royalblue',xlabel=xlabel)

    plt.subplot(gs[i,1])
    if i==0: plt.text(13,14,'CPDN',fontsize=25)
    if var_name == 'ssrd': 
        t = np.linspace(0,30,240)
        xlabel = 'Day'
    else: 
        t = np.linspace(0,30,720)
        xlabel = ''
    make_plot(t,CPDN_var_mean,ylim[i],title=title[i],color='gold',xlabel=xlabel)


figure_name = 'daily_variability_ssrd_t2m_10m_wind_speed.png'
print('saving to %s' % (figure_name))
plt.savefig(figure_name,bbox_inches='tight')
                                               
