1) Copy Template.ipynb to your folder.
2) You will also need a .xslx or .csv of your ['eblock', 'name', 'Destination Plate', 'Destination Well', 'Vector']
    *** This file needs to have the same number of proteins as you have SEC traces***
    *** The plate and well position of each protein should match up to what you put in the AKTA/HPLC ***
    -Name: the name of your design. (can be any name you give it)
    -Destination Well: the destination well of the transfer. (i.e. A1, A2, etc)
    -Destination Plate Name: 96w expression plate identifier. (i.e. 1, 2, etc) This should be an int. if not the script will strip down to anything after a '-'
    -eblock: the eBlock sequence you ordered. NOT what it should be after cloning
    -Vector: the GGA entry vector ID, e.g. LM0627.
    -n_chains: if your proteins are oligomers, this is the oligomeric state. not necessary if all your proteins are monomers
3) You need a folder containing 
    -Your .zip files if you used the AKTA
    -Folders containing .CSV files of your runs if you used the HPLC
        -In NanoES/MolES the 'cowboy' protocol should have an export file that automatically exports in the correct format
        -use 00_export_csv.pmx
        -If you used the HPLC you will also need to export a .csv file containing information on the fractions. 
        -use 00_export_fractions.rdl
        -make sure 'fractions' is somewhere in the name
        -place this in the same folder as your data file
4) make sure you provide a path to your CSV in the 2nd block of the .ipynb
    -comment out the vector line if you included your vector in your .csv
5) in the 3rd data block make sure to provide a path to the correct calibration file for the column you used and a path to the folder with your sec traces
