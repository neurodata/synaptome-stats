# Anish's synapse detection locations in the chessboard dataset

#### Note that Anish's locations are in `(row, col, z) = (y,x,z)` indexed from 1 with y = 0 as the top.

- Each folder is named according to the token of the dataset. 
- resultVol\_[1-6] correspond to Excitatory Queries 1 to 6.   
- resultVol\_[7-9] correspond to Inhibitory Queries 1, 2, 3

Query details should be given in tables 3-4
[here](http://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1005493&type=printable)


## Reformatted to XYZ

The Rscript `reformatCSV2xyz_m1.R` takes Anish's locations and reformats
the columns to be in `(x,y,z)` order and subtracts 1 from every element so
as to match up with python/ndio indexing from 0.

The reformatted location files are in the folder [XYZ](./XYZ/)

The script [XYZ/aggregateLocationsXYZ.R](./XYZ/aggretateLocationsXYZ.R)
aggregates all ot the locations and joins them into one data table
keeping the token and query for each.

