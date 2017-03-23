#!/bin/bash

echo "Start create AMOC of Ctrl"

data="/data/CESM/CTRL/xiao845_2400/ocn/hist/months"
path="/gpfsES/geo/the/MocArchieve/Ctrl"

# Create Ctrl path
if [ -d "$path" ]; then
	rm -r "$path"
	mkdir "$path"
else
	mkdir "$path"
fi
echo "Create Ctrl path complete!"

# Get AMOC from archive data of Ctrl
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
		ncks -v MOC $data/xiao845.pop.h.$ystr-$mstr.nc $path/month/Moc.Ctrl.$ystr-$mstr.nc
	done
done
echo "Get AMOC data of Ctrl complete!"

# Create .nc to save monthly data and annual mean data of Ctrl
mkdir "$path/monthly"
mkdir "$path/annual"
echo "Start month2year_Ctrl.ncl"
ncl month2year_Ctrl.ncl
echo "Create .nc to save monthly data and annual mean data of Ctrl complete!"

# Draw maximum series and do wavelet transform of Ctrl
echo "Start drawseries_Ctrl.ncl"
ncl drawseries_Ctrl.ncl
echo "Draw maximum series and do wavelet transform of Ctrl complete!"
