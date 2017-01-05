#!/bin/bash

echo "Start create SST of ALL"

path="/gpfsES/geo/the/MocArchieve/ALL-SST"

# Create ALL path
if [ -d "$path" ]; then
	rm -r "$path"
	mkdir "$path"
else
	mkdir "$path"
fi
echo "Create ALL-SST path complete!"

# Create .nc to save monthly data and annual mean data of ALL-SST
mkdir "$path/monthly"
mkdir "$path/annual"
echo "Start month2year_ALL-SST.ncl"
ncl month2year_ALL-SST.ncl
echo "Create .nc to save monthly data and annual mean data of ALL-SST complete!"

# Draw series of ALL-SST
echo "Start drawseries_ALL-SST.ncl"
ncl drawseries_ALL-SST.ncl
echo "Draw series of ALL-SST complete!"
