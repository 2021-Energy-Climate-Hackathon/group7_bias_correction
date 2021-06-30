""" Script containing function to perform bias correction of output from 
CPDN simulations. 

CPDN simulations are typically output on a rotated grid, hence prior to
bias correction the variables must be regridded.

The bias correction consists of a spatially varying, constant offset, 
which is calculated using a climatology of ERA5 data. The correction
is different for each month of the year. """

import numpy as np
from netCDF4 import Dataset

def check_dimensions(x,y):
    """ Check that the dimensions of two arrays (x,y) are the same """

    return x.shape == y.shape



def read_bias_correction_file(file_str,offset_name,chosen_month_index):
    """ Read in the bias correction offset for a particular month.

    Args:
        file_str (str) = file name and path of bias correction file.
        offset_name (str) = the name of the offset variable in the bias correction file.
        chosen_month_index = number 1...12 indicating the corresponding month (Jan...Dec)
    Returns:
        offset = the bias correction offset, but for a particular month only
  """
    
    # read in bias correction variable
    dataset = Dataset(file_str,'r')
    offset = dataset.variables[offset_name][:]
    month_index = dataset.variables[month_index][:] # 1...12

    return offset[chosen_month_index-1,:,:]



def apply_bias_correction(variable_biased,file_str_correction,offset_name='offset'):
    """ Function to apply the bias correction. This reads in a bias correction which is a constant
    offset for a given month with shape: (latitude, longitude). 

    Args:
        variable_biased = the variable which is to be bias corrected. 
        file_str_correction (str) = file name and path of bias correction file.
        offset_name (str) = the name of the offset variable in the bias correction file.

    Returns:
        variable_corrected = the bias corrected variable"""
    
    offset = read_bias_correction_file(file_str_correction,chosen_month_index)
        
    # check that the bias file and variable have the same spatial dimensions
    if check_dimensions(offset,variable_biased[0,:,:]) == False:
        print('Spatial dimensions of variable and bias correction offset are different')
        raise ValueError      
    
    # bias correction
    variable_corrected = variable_biased + offset

    return variable_corrected 
    



