""" A set of functions to read gridded climate data """

import numpy as np
from netCDF4 import Dataset, num2date, date2num
import cartopy.io.shapereader as shpreader
import shapely.geometry


def read_spatial_dimensions(file_name,lat_name=None,lon_name=None,lev_name=None):
    """ Returns the latitude, longitude and pressure (if available) dimensions """
    nc = Dataset(file_name,'r')
    lats=None
    lons=None
    levs=None
    common_latitude_names = ['lat','lats','latitude','latitude0','latitude1']
    common_longitude_names = ['lon','lons','longitude','longitude0','longitude1']
    common_pressure_names = ['level','levels','levs','p','plevs','plev','z0','z5','z6','z7','z8']
    if lat_name is not None: lats = nc.variables[lat_name][:]
    else:
        for name in common_latitude_names:
            try:
                lats = nc.variables[name][:]
                break
            except: pass
    if lon_name is not None: lons = nc.variables[lon_name][:]
    else:
        for name in common_longitude_names:
            try:
                lons = nc.variables[name][:]
                break
            except: pass
    if lev_name is not None: levs = nc.variables[lev_name][:]
    else:
        for name in common_pressure_names:
            try:
                levs = nc.variables[name][:]
                break
            except: pass
    return levs,lats,lons


def read_time_dimension(file_name,time_name=None):
    """ Reads in data on the time dimension and returns it along with units
    and calendar """
    nc = Dataset(file_name)
    common_time_names = ['t','time','times']
    times = None
    calendar = None
    units = None
    if time_name is not None: times = nc.variables[time_name][:]
    else:
        for time_name in common_time_names:
            try:
                times = nc.variables[time_name][:]
                break
            except: pass
    try: calendar = nc.variables[time_name].calendar
    except: print('Warning: calendar not defined')
    try: units = nc.variables[time_name].units
    except: print('Warning: time units not defined')
    return times,calendar,units


def get_variable(file_name,var_name,lat_name=None,lon_name=None,lev_name=None,chosen_lev=None):
    """ Read in a given variable on a chosen pressure level """
    nc = Dataset(file_name,'r')
    levs,lats,lons = get_dimensions(file_name,lat_name=lat_name,lon_name=lon_name,lev_name=lev_name)
    lev_idx = np.argmin(np.abs(levs - chosen_lev))
    if chosen_lev is not None: var = nc.variables[var_name][:,lev_idx,:,:]
    else: var = nc.variables[var_name][:]

    return var


def load_country_weather_data(COUNTRY,data_dir,filename,nc_key,chosen_months='all',lon_name='longitude',lat_name='latitude',time_name='times'):

    """
    This function takes the ERA5 reanalysis data, loads it and applied a 
    country mask (ready for conversion to energy) it then returns
    the array (of original size) with all irrelvelant gridpoints 
    set to zeros.
    You will need the shpreader.natural_earth data downloaded 
    to find the shapefiles.
    Args:
        COUNTRY (str): This must be a name of a country (or set of) e.g. 
            'United Kingdom','France','Czech Republic'
        data_dir (str): The parth for where the data is stored.
            e.g '/home/users/zd907959/'
        filename (str): The filename of a .netcdf file
            e.g. 'ERA5_1979_01.nc'
        nc_key (str): The string you need to load the .nc data 
            e.g. 't2m','rsds'
    Returns:
        country_masked_data (array): Country-masked weather data, dimensions 
            [time,lat,lon] where there are 0's in locations where the data is 
            not within the country border.
        MASK_MATRIX_RESHAPE (array): Dimensions [lat,lon] where there are 1's if 
           the data is within a country border and zeros if data is outside a 
           country border. 
    """


    # first loop through the countries and extract the appropraite shapefile
    countries_shp = shpreader.natural_earth(resolution='10m',category='cultural',
                                            name='admin_0_countries')
    country_shapely = []
    for country in shpreader.Reader(countries_shp).records():
        if country.attributes['NAME_LONG'] == COUNTRY:
            print('Found country')
            country_shapely.append(country.geometry)

    # load in the data you wish to mask
    file_str = data_dir + filename
    dataset = Dataset(file_str,mode='r')
    lons = dataset.variables[lon_name][:]
    lats = dataset.variables[lat_name][:]
    if chosen_months=='all': data = dataset.variables[nc_key][:] # data in shape [time,lat,lon]
    else: 
        # only read in chosen month, months is now an index for that month
        times = dataset.variables[time_name][:]
        calendar, units = dataset.variables[time_name].calendar, dataset.variables[time_name].units
        dates = num2date(times,calendar=calendar,units=units)
        months = np.array([])
        for d in dates:
            months = np.append(months,d.timetuple()[1])
        data = dataset.variables[nc_key][months==chosen_months]

    dataset.close()

    # get data in appropriate units for models
    if nc_key == 't2m':
        data = data-273.15 # convert to Kelvin from Celsius
    if nc_key == 'ssrd':
        data = data/3600. # convert Jh-1m-2 to Wm-2

    LONS, LATS = np.meshgrid(lons,lats) # make grids of the lat and lon data
    x, y = LONS.flatten(), LATS.flatten() # flatten these to make it easier to 
    #loop over.
    points = np.vstack((x,y)).T
    MASK_MATRIX = np.nan*np.zeros((len(x),1))
    # loop through all the lat/lon combinations to get the masked points
    for i in range(0,len(x)):
        my_point = shapely.geometry.Point(x[i],y[i]) 
        if country_shapely[0].contains(my_point) == True: 
            MASK_MATRIX[i,0] = 1.0 # creates 1s and 0s where the country is
    
    MASK_MATRIX_RESHAPE = np.reshape(MASK_MATRIX,(len(lats),len(lons)))

    # now apply the mask to the data that has been loaded in:

    country_masked_data = data*MASK_MATRIX_RESHAPE
                                     


    return(country_masked_data,MASK_MATRIX_RESHAPE)
