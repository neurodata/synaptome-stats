

## F0 feature calculation in Docker

### Usage 

TOKEN: The ndstore token for your dataset of interest.  
INPUT: The name of the file of locations with columns x, y, z.  
 NAME: The string to use as output.  
   BF: An integer to use as buffer around each location  

Given a csv file of locations in (x,y,z) column order


`docker run -v data/:/data $TOKEN1 $INPUT $NAME $BF`


