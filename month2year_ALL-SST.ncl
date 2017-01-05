load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

	; Open monthly files
	f = addfile("/gpfsES/geo/zywang/data_CESM/MASK_CESM/data/sst_new.nc", "r")

	; TEMP(time=24000, z_t=1, lat=48, lon=96)
	TEMP = f->TEMP
	lat = f->lat
	lon = f->lon

	; Calculate annual mean of NH Atlantic Ocean
	do i = 0, 23999, 12
		; Create new varible to save annual mean data (0, lat, lon)
		annual = new((/1, dimsizes(lat), dimsizes(lon)/), float, "No_FillValue")
		annual@_FillValue = -999

		; Get annual mean
        	do j = 24, 40
        		do k = 75, 95
				annual(0, j, k) = avg(TEMP(i:(i + 11), 0, j, k))
			end do
        	end do

		annual!0 = "time"
                annual!1 = "lat"
                annual!2 = "lon"

		; Create new variables for each dimensions of annual
                annual&time = new((/1/), float, "No_FillValue")
                annual&lat  = new((/dimsizes(lat)/), double, "No_FillValue")
                annual&lon  = new((/dimsizes(lon)/), double, "No_FillValue")

		; Assign dimensions
                annual&time = ispan(1, 1, 1)
                annual&lat = f->lat
                annual&lon = f->lon

		; Name the variable and unit
                annual@long_name = "Potential Temperature"
                annual@units = "degC"

		; Create a new .nc to save monthly file
		y = (i + 12) / 12
		ystr = sprinti("%0.4i", y)
                system("rm -f /gpfsES/geo/the/MocArchieve/ALL-SST/annual/SST.ALL.annual." + ystr + ".nc")
                out = addfile("/gpfsES/geo/the/MocArchieve/ALL-SST/annual/SST.ALL.annual." + ystr + ".nc", "c")

                ; Assign the value to out
                out->TEMP = annual
                out->lat  = lat
                out->lon  = lon
    	end do
end
