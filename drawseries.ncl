;load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
                
begin 

	; Create new varible to save annual maximum data
	annual_max = new((/2000/), float, "No_FillValue")

	; Traversal 2000 years
        do y = 1, 2000
		; Change int to string
                ystr = sprinti("%0.4i", y)

		; Open annual mean file
		f = addfile("/gpfsES/geo/the/MocArchieve/annual/Moc.annual." + ystr + ".nc", "r")
		
		MOC = f->MOC
		
		; Get the maximum from MOC (Under 500m)
		annual_max((y - 1)) = max(MOC(0, 33:, :))
		
	end do	

	annual_max!0 = "time"
	annual_max&time = new((/2000/), float, "No_FillValue")
	annual_max&time = ispan(1, 2000, 1)
	annual_max@long_name = "Maximum of annual mean of Meridional Overturning Circulation"
        annual_max@units = "Sverdrups"

	; Create a new .nc to save maximum file
        system("rm -f /gpfsES/geo/the/MocArchieve/Moc.maximum.nc")
        out = addfile("/gpfsES/geo/the/MocArchieve/Moc.maximum.nc", "c")
	
	; Assign the value to out
	out->MOCMax = annual_max
	
	; Smooth
	annual_sm = runave(annual_max, 31, 0)	

	; Draw series
	time = new((/2000/), float, "No_FillValue")
	time = ispan(1, 2000, 1)
	
	wks = gsn_open_wks("png","/gpfsES/geo/the/MocArchieve/Max_Series")

	res = True
	res@tiMainString = "Maximum of MOC"
	res@tiYAxisString = "Meridional Overturning Circulation (Sverdrups)"
	res@tiXAxisString = "Year"
	res@xyLineColors = "red"
	res@vpHeightF = 0.43
	res@vpWidthF = 0.70
	
	plot = gsn_xy(wks, time, annual_sm, res)
	
end
