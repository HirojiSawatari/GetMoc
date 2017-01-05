load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
                
begin 
	; 1, do bias correction to the wavelet powers
	bias_correct = 1
	
	; Create new varible to save annual maximum data
	annual_ave = new((/2000/), double, "No_FillValue")

	; Create new varible to save annual maximum data of AMOC
        annual_max = new((/2000/), float, "No_FillValue")

	; Traversal 2000 years
        do y = 1, 2000
		; Change int to string
                ystr = sprinti("%0.4i", y)

		; Open annual mean file
		f = addfile("/gpfsES/geo/the/MocArchieve/ALL-SST/annual/SST.ALL.annual." + ystr + ".nc", "r")
		f2 = addfile("/gpfsES/geo/the/MocArchieve/ALL/annual/Moc.ALL.annual." + ystr + ".nc", "r")		

		TEMP = f->TEMP
		MOC = f2->MOC		

		; Get the average from SST of NH Atlantic Ocean
		annual_ave((y - 1)) = sum(TEMP(0, 24:40, 75:95)) / 250
		
		; Get the maximum from MOC (Under 500m)
                annual_max((y - 1)) = max(MOC(0, 33:, :))

	end do	

	annual_ave!0 = "time"
	annual_ave&time = new((/2000/), float, "No_FillValue")
	annual_ave&time = ispan(1, 2000, 1)
	annual_ave@long_name = "AMO"
        annual_ave@units = "degC"

	; Create a new .nc to save AMO file
        system("rm -f /gpfsES/geo/the/MocArchieve/ALL-SST/SST.average.nc")
        out = addfile("/gpfsES/geo/the/MocArchieve/ALL-SST/SST.average.nc", "c")
	
	; Assign the value to out
	out->AMO = annual_ave

        annual_max!0 = "time"
        annual_max&time = new((/2000/), float, "No_FillValue")
        annual_max&time = ispan(1, 2000, 1)
        annual_max@long_name = "Maximum of annual mean of Meridional Overturning Circulation"
        annual_max@units = "Sverdrups"
	
	; Set time
        time = new((/2000/), float, "No_FillValue")
        time = ispan(1, 2000, 1)

	; Set workspace
	output = "png"
	output@wkWidth = 1500
        output@wkHeight = 1080
	
	; Draw original series
        wks = gsn_open_wks(output, "/gpfsES/geo/the/MocArchieve/ALL-SST/Original_Series")

        resL = True
        resL@tiMainString = "AMO (Original)"
        resL@tiYAxisString = "Potential Temperature (degC)"
        resL@tiXAxisString = "Year"
        resL@xyLineColors = "blue"
        resL@vpHeightF = 0.43
        resL@vpWidthF = 0.65
	resL@trXMinF = 0
        resL@trXMaxF = 2000
        ;resL@trYMinF = 14
        ;resL@trYMaxF = 23
	resL@vpXF = 0.15

	resR = True
        resR@tiYAxisString = "Meridional Overturning Circulation (Sverdrups)"
        resR@tiXAxisString = "Year"
        resR@xyLineColors = "red"
        resR@trYMinF = 14
        resR@trYMaxF = 23

        plot = gsn_csm_xy2(wks, time, annual_ave, annual_max, resL, resR)

	; Smooth
	annual_sm = runave(annual_ave, 31, 0)	
	annual_sm2 = runave(annual_max, 31, 0)

	; Draw smooth series
	wks2 = gsn_open_wks(output, "/gpfsES/geo/the/MocArchieve/ALL-SST/Smooth_Series")

	resL2 = True
	resL2@tiMainString = "AMO (Smooth)"
	resL2@tiYAxisString = "Potential Temperature (degC)"
	resL2@tiXAxisString = "Year"
	resL2@xyLineColors = "blue"
	resL2@vpHeightF = 0.43
	resL2@vpWidthF = 0.65
        resL2@trXMinF = 0
        resL2@trXMaxF = 2000
        ;resL2@trYMinF = 15
        ;resL2@trYMaxF = 21
	resL2@vpXF = 0.15

	resR2 = True
        resR2@tiYAxisString = "Meridional Overturning Circulation (Sverdrups)"
        resR2@tiXAxisString = "Year"
        resR2@xyLineColors = "red"
        resR2@trYMinF = 15
        resR2@trYMaxF = 21
	
	plot2 = gsn_csm_xy2(wks2, time, annual_sm, annual_sm2, resL2, resR2)

end
