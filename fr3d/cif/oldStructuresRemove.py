#In order for these structures to be annotated in the old way, change the name of this file to "oldStructures.py".


# Processing of complex symmetries were not available in the Matlab Version of FR3D or in early versions of fr3d-python.
# As a result, the production server has some structures from the PDB with poorly named symmetries. For the sake of our database, we need to be able 
# to emulate the old naming convention for unit ids, this file will be used to read those files. 
# This file is just an old version of the fr3d python cif reader that has the old functionality maintained. The list of 'oldStructures' below can be appended if more 
# legacy structures are found that seem to have issues. 

# Some old structures need to have processing for the sake of the server and naming conventions. Default to this for these structures.
# Updated 11/29/2022
# SQL Statement to get list: 
# select distinct pdb_id from unit_info where sym_op='P_1';
# some structures from the old structures are appended to the end of each row as well that were not returned from the sql statement. 
oldStructures = ['1A34', '1AQ3', '1AQ4', '1BMV', '1CWP','1DDL','1F8V','1KUO',  '1LAJ','1PGL','1RMV','1U1Y','1VTM','1ZDH','1ZDI','1ZDJ','1ZDK','1ZSE', '1FJF', '1VS9',
                 '2B2D','2B2E','2B2G','2BBV','2BNY','2BQ5', '2BS0','2BS1','2BU1','2C4Q','2C4Y','2C4Z','2C50','2C51','2FZ2','2IZ8','2IZ9','2IZM','2IZN','2Q23','2Q25','2Q26','2QQP','2TMV','2XPJ', '2Z2Q', '2I1C',
                 '3LOB', '3S4G', 
                 '4ANG', '4FSJ', '4FTB', '4FTE', '4FTS', '4OQ8', '4OQ9', '4WR6', '4WRO', '4WZD','4Z92',
                 '5A79', '5A7A', '5A9Z', '5AA0', '5AFI', '5APO', '5L7Q','5L8Q','5M74','5MJV','5MSF','5MUP','5MV5','5MV6','5UF6','5UHY','5VA1','5VA2','5VA3','5WSN','5Z9W','5FN1',
                 '6GV4', '6H5Q', '6I2N', '6MSF','6NUT','6THN',
                 '7EWQ','7MSF','7NUN','7NUQ','7OI3']