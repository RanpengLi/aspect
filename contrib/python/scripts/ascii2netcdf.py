'''
This function creates netcdf files from an input ascii file 
formatted following the ASPECT convention for structured ascii data files.
The generated netcdf file can be used for visualization in Paraview or 
as an input file in ASPECT.

see https://unidata.github.io/netcdf4-python/ for more information

'''

import netCDF4 as nc
import numpy as np
import sys

def main():

    msg = "Usage            : ascii2netcdf.py input_ascii_filename  [dimensions, coordinate_system, output_format] \n \n" + \
          "Optional arguments: \n" + \
          "dimensions       : '2' or '3', the dimensions of the data in the ascii file \n" +  \
          "coordinate_system: 'spherical' or 'cartesian' \n" +  \
          "output_format    : 'aspect' or 'paraview'\n" 
          
    if (sys.argv[1] == '-h'):
        print(msg)
        sys.exit(0)

    # Load command line arguments
    args = sys.argv[0:]
    argc = len(sys.argv[0:])
    print(args, argc)

    if (argc < 2):
        raise Exception("Not enough arguments. Please provide the name of the input ascii file for conversion.")

    # The first argument is the name of the file, which is required.
    current_argument = 2
    file_name        = sys.argv[1]

    # default argument values
    dim               = 2
    coordinate_system = 'cartesian'
    output_format     = 'aspect'
    
    while (current_argument < argc):
        if (args[current_argument] == '3' or args[current_argument] == '2'):
            dim = sys.argv[2]
            current_argument += 1
        elif (args[current_argument] == 'spherical' or args[current_argument] == 'cartesian'):
            coordinate_system = sys.argv[3]
            current_argument += 1
        elif (args[current_argument] == 'paraview' or args[current_argument] == 'aspect'):
            output_format     = sys.argv[4]
            current_argument += 1
        else:
            raise Exception("Invalid argument list. Please check the usage of the script using -h flag.")
        
    ascii_data = np.loadtxt(file_name)

    # TODO: make this an optional argument.
    ofilename  = file_name.split('.')[0] + '.nc'

    if (coordinate_system == 'cartesian'):
        cartesian_coordinate_system(ascii_data, int(dim), ofilename)
    elif (coordinate_system == 'spherical'):
        spherical_coordinate_system (ascii_data, int(dim), output_format, ofilename)


def cartesian_coordinate_system(ascii_data, dim, ofile) :

    ds = nc.Dataset(ofile, 'w', format='NETCDF4')

    coord_x = np.unique(ascii_data[:, 0])
    coord_y = np.unique(ascii_data[:, 1])
    n_x     = len(coord_x)
    n_y     = len(coord_y)
    coords  = [ds.createDimension('X', len(coord_x)),
               ds.createDimension('Y', len(coord_y))]
    
    # variables or dimensions:
    data_x = ds.createVariable('X', np.float64, ('X'))
    data_x.units = 'm'
    data_y = ds.createVariable('Y', np.float64, ('Y'))
    data_y.units = 'm'

    if (dim == 3):
        coord_z = np.unique(ascii_data[:, 2])
        n_z     = len(coord_z)
        coords.append(ds.createDimension('Z', len(coord_z)))

        data_z = ds.createVariable('Z', np.float64, ('Z'))
        data_z.units = 'm'

    # load all the columns after the dimension columns
    field_data = ascii_data[:, dim:]

    for i in range(field_data.shape[1]):
        data = ds.createVariable(('field_'+ str(i)), np.float64, ('Z','Y','X') if dim==3 else ('Y','X'))
        # np.reshape takes the indices in reverse order
        data[:] = np.reshape(ascii_data[:, dim+i], (n_z, n_y, n_x)) if dim==3 else np.reshape(ascii_data[:, dim+i], (n_y, n_x))

   # coordinate data
    data_x[:] = coord_x
    data_y[:] = coord_y
    if (dim == 3):
        data_z[:] = coord_z

    ds.close()


def spherical_coordinate_system(ascii_data, dim, output_format, ofile):

    ds = nc.Dataset(ofile, 'w', format='NETCDF4')

    coords_r     = np.unique(ascii_data[:, 0])
    coords_phi   = np.unique(ascii_data[:, 1])
    n_r          = len(coords_r)
    n_phi        = len(coords_phi)
    
    coords = [ds.createDimension('radius', n_r),
              ds.createDimension('longitude', n_phi)]
    
    # It is important to follow the convention here in Paraview for
    # spherical coordinate system visualization.
    # For more details see http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/cf-conventions.html#latitude-coordinate
    radius           = ds.createVariable('radius', np.float64, ('radius'), fill_value=np.nan)
    radius.units     = 'meters'
    longitudes       = ds.createVariable('longitude', np.float64, ('longitude'), fill_value=np.nan)
    longitudes.units = 'degrees_east'

    if (dim == 3):
        coords_theta     = (np.unique(ascii_data[:, 2]))
        n_theta          = len(coords_theta)
        coords.append(ds.createDimension('latitude',  n_theta))
        latitudes        = ds.createVariable('latitude', np.float64, ('latitude'), fill_value=np.nan)
        latitudes.units  = 'degrees_north'
    
    # load all the columns after the dimension columns
    field_data = ascii_data[:, dim:]
    
    for i in range(field_data.shape[1]):
        data = ds.createVariable(('field_'+ str(i)), np.float64, ('latitude', 'longitude', 'radius') if dim==3 else \
                                 ('longitude', 'radius'))
        # np.reshape takes the indices in reverse order
        data[:] = np.reshape(ascii_data[:, dim+i], (n_theta, n_phi, n_r)) if dim==3 else np.reshape(ascii_data[:, dim+i], (n_phi, n_r))

        if (output_format == 'paraview' and dim == 3):
            data_temp    = np.reshape(ascii_data[:, dim+i], (n_theta, n_phi, n_r))
            data[:]      = np.flip(data_temp, axis=0)

    # coordinate data
    radius[:]     = coords_r
    longitudes[:] = coords_phi        
    if (dim == 3):
        latitudes[:] = coords_theta
    
    # convert phi and theta from spherical coordinate system to geographical
    if (output_format == 'paraview'):
        longitudes[:] = np.degrees(coords_phi)
        if (dim == 3):
            latitudes[:] = 90 - np.degrees(coords_theta)
    
    ds.close()


if __name__ == "__main__":
    main()
