#!/bin/bash

echo "Start create AMOC of ANTH"

data="/gpfsES/geo/zywang/CESM/Archive/anthropic/ocn/hist/months"
path="/gpfsES/geo/the/MocArchieve/ANTH"

# Create ANTH path
if [ -d "$path" ]; then
	rm -r "$path"
	mkdir "$path"
else
	mkdir "$path"
fi
echo "Create ANTH path complete!"

# Get AMOC from archive data of ANTH
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
		ncks -v MOC $data/anthropic.pop.h.$ystr-$mstr.nc $path/month/Moc.ANTH.$ystr-$mstr.nc
	done
done
echo "Get AMOC data of ANTH complete!"

# Create .nc to save monthly data and annual mean data of ANTH
mkdir "$path/monthly"
mkdir "$path/annual"
echo "Start month2year_ANTH.ncl"
ncl month2year_ANTH.ncl
echo "Create .nc to save monthly data and annual mean data of ANTH complete!"

# Draw maximum series and do wavelet transform of ANTH
echo "Start drawseries_ANTH.ncl"
ncl drawseries_ANTH.ncl
echo "Draw maximum series and do wavelet transform of ANTH complete!"

# Draw volcanic series
echo "Start volseries.ncl"
ncl volseries.ncl
echo "Draw volcanic series complete!"
