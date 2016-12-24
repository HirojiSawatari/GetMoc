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
		f = addfile("/gpfsES/geo/the/MocArchieve/ANTH/annual/Moc.ANTH.annual." + ystr + ".nc", "r")
		
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
        system("rm -f /gpfsES/geo/the/MocArchieve/ANTH/Moc.maximum.nc")
        out = addfile("/gpfsES/geo/the/MocArchieve/ANTH/Moc.maximum.nc", "c")
	
	; Assign the value to out
	out->MOCMax = annual_max
	
	; Set time
        time = new((/2000/), float, "No_FillValue")
        time = ispan(1, 2000, 1)

        ; Open GHGs file
        g = addfile("/gpfsES/geo/zywang/Rec_Obs_For_DATA/forcings/ghgs/GHG.nc", "r")

	; Create new varible to save GHGs data
        ghg = new((/2003/), double, "No_FillValue")

        ; ANTH(time=2003)
        ch4 = g->CH4
	co2 = g->CO2
	n2o = g->N2O
        time2 = g->time

	do i = 0, 2002
		ghg(i) = (ch4(i) * 1000) + co2(i) + (n2o(i) * 1000)
	end do

	; Set workspace
	output = "png"
	output@wkWidth = 1500
        output@wkHeight = 1080
	
	; Draw original series
        wks = gsn_open_wks(output, "/gpfsES/geo/the/MocArchieve/ANTH/Original_Series")

        resL = True
        resL@tiMainString = "Maximum of MOC (Original)"
        resL@tiYAxisString = "Meridional Overturning Circulation (Sverdrups)"
        resL@tiXAxisString = "Year"
        resL@xyLineColors = "blue"
        resL@vpHeightF = 0.43
        resL@vpWidthF = 0.65
	resL@trXMinF = 0
        resL@trXMaxF = 2000
	resL@vpXF = 0.15

	resR = True
	resR@tiYAxisString = "Greenhouse Gases (1e-6 mol/mol)"
        resR@trXMinF = 0
        resR@trXMaxF = 2000
        resR@xyLineColors = "red"

        plot = gsn_csm_xy2(wks, time, annual_max, ghg(:1999), resL, resR)

	; Smooth
	annual_sm = runave(annual_max, 31, 0)	

	; Draw smooth series
	wks2 = gsn_open_wks(output, "/gpfsES/geo/the/MocArchieve/ANTH/Smooth_Series")

	resL2 = True
	resL2@tiMainString = "Maximum of MOC (Smooth)"
	resL2@tiYAxisString = "Meridional Overturning Circulation (Sverdrups)"
	resL2@tiXAxisString = "Year"
	resL2@xyLineColors = "blue"
	resL2@vpHeightF = 0.43
	resL2@vpWidthF = 0.65
	resL2@trXMinF = 0
        resL2@trXMaxF = 2000
	resL2@vpXF = 0.15

	resR2 = True
	resR2@tiYAxisString = "Greenhouse Gases (1e-6 mol/mol)"
        resR2@trXMinF = 0
        resR2@trXMaxF = 2000
        resR2@xyLineColors = "red"
	
	plot2 = gsn_csm_xy2(wks2, time, annual_sm, ghg(:1999), resL2, resR2)

	; Calculate wavelet
	w = wavelet_default(annual_max, 0)

	; Draw wavelet transform
	map = new(2,graphic)
	
        wks3 = gsn_open_wks(output, "/gpfsES/geo/the/MocArchieve/ANTH/Wavelet_Transform")

        res3 = True
        res3@tiMainString = "Wavelet Transform"
	res3@cnFillOn = True
	res3@trYReverse = True
	res3@gsnSpreadColors = True
	res3@gsnSpreadColorStart = 24
	res3@gsnSpreadColorEnd = -26
        res3@vpHeightF = 0.43
        res3@vpWidthF = 0.65
	res3@vpXF = 0.15

	map(0) = gsn_csm_contour(wks3,w(0,:,:),res3) 

	res4 = True
        res4@xyLineColors = "red"
	res4@vpHeightF = 0.1
        res4@vpWidthF = 0.65
        res4@trXMinF = 0
        res4@trXMaxF = 2000
	res4@vpXF = 0.15

	map(1) = gsn_xy(wks3, time(:1999), ghg(:1999), res4)
        
	yres = True

	gsn_panel(wks3, map, (/2, 1/), yres)

end
