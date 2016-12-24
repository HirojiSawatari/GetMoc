#!/bin/bash

echo "Start create AMOC of TSI"

data="/gpfsES/geo/zywang/CESM/Archive/shapiroMAX/ocn/hist/months"
path="/gpfsES/geo/the/MocArchieve/TSI"

# Create TSI path
if [ -d "$path" ]; then
	rm -r "$path"
	mkdir "$path"
else
	mkdir "$path"
fi
echo "Create TSI path complete!"

# Get AMOC from archive data of TSI
mkdir "$path/month"
for ((y=1; y<=2000; y++))
do
	ystr=$y
	if [ $y -lt 1000 ]; then
		if [ $y -lt 10 ]; then
			ystr="000"$ystr
		elif [ $y -lt 100 ]; then
			ystr="00"$ystr
		else
			ystr="0"$ystr
		fi
	fi
	for ((m=1; m<=12; m++))
	do
		mstr=$m
		if [ $m -lt 10 ]; then
			mstr="0"$m	
		fi
		ncks -v MOC $data/shapiroMAX.pop.h.$ystr-$mstr.nc $path/month/Moc.TSI.$ystr-$mstr.nc
	done
done
echo "Get AMOC data of TSI complete!"

# Create .nc to save monthly data and annual mean data of TSI
mkdir "$path/monthly"
mkdir "$path/annual"
echo "Start month2year_TSI.ncl"
ncl month2year_TSI.ncl
echo "Create .nc to save monthly data and annual mean data of TSI complete!"

# Draw maximum series and do wavelet transform of TSI
echo "Start drawseries_TSI.ncl"
ncl drawseries_TSI.ncl
echo "Draw maximum series and do wavelet transform of TSI complete!"

# Draw volcanic series
echo "Start volseries.ncl"
ncl volseries.ncl
echo "Draw volcanic series complete!"
