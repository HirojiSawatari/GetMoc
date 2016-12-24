#!/bin/bash

echo "-------------------------------------------------------------------------------------------------------"
echo "|                                            GetMoc V1.0                                              |"
echo "| Get meridional overturning circulation (MOC) from MASK and calculate the annual mean series of MOC. |"
echo "| See https://github.com/HirojiSawatari/GetMoc for more details.                                      |"
echo "| Copyright (C) 2016 Tou Ka All Rights Reserved                                                       |"
echo "-------------------------------------------------------------------------------------------------------"

path="/gpfsES/geo/the/MocArchieve"

# Create loop path
if [ -d "$path" ]; then
        rm -r "$path"
        mkdir "$path"
else
        mkdir "$path"
fi
echo "Create loop path complete!"

# Start create AMOC of ALL
bash all-moc.sh

# Start create AMOC of TSI
bash tsi-moc.sh

# Start create AMOC of VOL
bash vol-moc.sh

# Start create AMOC of ANTH
bash anth-moc.sh

# Draw volcanic series
echo "Start volseries.ncl"
ncl volseries.ncl
echo "Draw volcanic series complete!"
