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
	
	; Set time
        time = new((/2000/), float, "No_FillValue")
        time = ispan(1, 2000, 1)

	; Draw original series
        wks = gsn_open_wks("png","/gpfsES/geo/the/MocArchieve/Original_Series")

        res = True
        res@tiMainString = "Maximum of MOC (Original)"
        res@tiYAxisString = "Meridional Overturning Circulation (Sverdrups)"
        res@tiXAxisString = "Year"
        res@xyLineColors = "red"
        res@vpHeightF = 0.43
        res@vpWidthF = 0.70

        plot = gsn_xy(wks, time, annual_max, res)

	; Smooth
	annual_sm = runave(annual_max, 31, 0)	

	; Get volcanic series
	vol = addfile("/gpfsES/geo/zywang/Rec_Obs_For_DATA/forcings/volc/IVI2.5_Gao_Sigl_1-2000.nc", "r")

	; Draw smooth series
	wks2 = gsn_open_wks("png","/gpfsES/geo/the/MocArchieve/Smooth_Series")

	res2 = True
	res2@tiMainString = "Maximum of MOC (Smooth)"
	res2@tiYAxisString = "Meridional Overturning Circulation (Sverdrups)"
	res2@tiXAxisString = "Year"
	res2@xyLineColors = "red"
	res2@vpHeightF = 0.43
	res2@vpWidthF = 0.70
	
	plot2 = gsn_xy(wks2, time, annual_sm, res2)

	; Calculate wavelet
	w = wavelet_default(annual_max, 0)

	; Draw wavelet series
        wks3 = gsn_open_wks("png","/gpfsES/geo/the/MocArchieve/Wavelet_Transform")

        res3 = True
        res3@tiMainString = "Wavelet Transform"
	res3@cnFillOn = True
	res3@trYReverse = True
	res3@gsnSpreadColors = True
	res3@gsnSpreadColorStart = 24
	res3@gsnSpreadColorEnd = -26
        res3@vpHeightF = 0.43
        res3@vpWidthF = 0.70

	plot3 = gsn_csm_contour(wks3,w(0,:,:),res3) 
        
end
