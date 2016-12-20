#!/bin/bash

data="/gpfsES/geo/zywang/CESM/Archive/MASK/ocn/hist/months"
path="/gpfsES/geo/the/MocArchieve"

# Create path
if [ -d "$path" ]; then
	rm -r "$path"
	mkdir "$path"
else
	mkdir "$path"
fi
echo "Create path complete!"

# Get moc from archive data of CESM
mkdir "$path/month"
for ((y=1; y<=50; y++))
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
		ncks -v MOC $data/MASK.pop.h.$ystr-$mstr.nc $path/month/Moc.$ystr-$mstr.nc
	done
done
echo "Get moc data complete!"

# Create .nc to save monthly data and annual mean data
mkdir "$path/monthly"
mkdir "$path/annual"
ncl month2year.ncl
echo "Create .nc to save monthly data and annual mean data complete!"


