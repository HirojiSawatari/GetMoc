#!/bin/bash

echo "Start create AMOC of VOL"

data="/gpfsES/geo/zywang/CESM/Archive/VOL/ocn/hist/months"
path="/gpfsES/geo/the/MocArchieve/VOL"

# Create VOL path
if [ -d "$path" ]; then
	rm -r "$path"
	mkdir "$path"
else
	mkdir "$path"
fi
echo "Create VOL path complete!"

# Get AMOC from archive data of VOL (VOL HAS ONLY 1500 YEARS)
mkdir "$path/month"
for ((y=501; y<=2000; y++))
do
	ystr=$y
	if [ $y -lt 1000 ]; then
		ystr="0"$ystr
	fi
	for ((m=1; m<=12; m++))
	do
		mstr=$m
		if [ $m -lt 10 ]; then
			mstr="0"$m	
		fi
		ncks -v MOC $data/VOL.pop.h.$ystr-$mstr.nc $path/month/Moc.VOL.$ystr-$mstr.nc
	done
done
echo "Get AMOC data of VOL complete!"

# Create .nc to save monthly data and annual mean data of VOL
mkdir "$path/monthly"
mkdir "$path/annual"
echo "Start month2year_VOL.ncl"
ncl month2year_VOL.ncl
echo "Create .nc to save monthly data and annual mean data of VOL complete!"

# Draw maximum series and do wavelet transform of VOL
echo "Start drawseries_VOL.ncl"
ncl drawseries_VOL.ncl
echo "Draw maximum series and do wavelet transform of VOL complete!"
