load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
                
begin 
	; 1, do bias correction to the wavelet powers
        bias_correct = 1

	; Create new varible to save annual maximum data
	annual_max = new((/2000/), float, "No_FillValue")

	; Traversal 1500 years (VOL HAS ONLY 1500 YEARS)
        do y = 501, 2000
		; Change int to string
                ystr = sprinti("%0.4i", y)

		; Open annual mean file
		f = addfile("/gpfsES/geo/the/MocArchieve/VOL/annual/Moc.VOL.annual." + ystr + ".nc", "r")
		
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
        system("rm -f /gpfsES/geo/the/MocArchieve/VOL/Moc.maximum.nc")
        out = addfile("/gpfsES/geo/the/MocArchieve/VOL/Moc.maximum.nc", "c")
	
	; Assign the value to out
	out->MOCMax = annual_max
	
	; Set time
        time = new((/2000/), float, "No_FillValue")
        time = ispan(1, 2000, 1)

	; Create new varible to save volcanic data
        volc_global = new((/24000/), double, "No_FillValue")
        volc_north = new((/24000/), double, "No_FillValue")
        volc_south = new((/24000/), double, "No_FillValue")

        ; Open volcanic file
        v = addfile("/gpfsES/geo/zywang/Rec_Obs_For_DATA/forcings/volc/IVI2.5_Gao_Sigl_1-2000.nc", "r")

        ; volc(time=24002, lev=18, lat=64)
        volc = v->MMRVOLC
        time2 = v->time
        lat = v->lat

        ; Create temp
        temp_latdata = new((/64/), double, "No_FillValue")

        ; Get the series
        do i = 0, 23999
                do j = 0, 63
                        temp_latdata(j) = avg(volc(i, :, j))
                end do
                volc_global(i) = avg(temp_latdata(:))
                volc_south(i) = avg(temp_latdata(0:31))
                volc_north(i) = avg(temp_latdata(32:63))
        end do

	; Set workspace
	output = "png"
	output@wkWidth = 1500
        output@wkHeight = 1080
	
	; Draw original series
        wks = gsn_open_wks(output, "/gpfsES/geo/the/MocArchieve/VOL/Original_Series")

        resL = True
        resL@tiMainString = "Maximum of MOC (Original)"
        resL@tiYAxisString = "Meridional Overturning Circulation (Sverdrups)"
        resL@tiXAxisString = "Year"
        resL@xyLineColors = "blue"
        resL@vpHeightF = 0.43
        resL@vpWidthF = 0.65
        resL@trXMinF = 0
        resL@trXMaxF = 2000
        resL@trYMinF = 14
        resL@trYMaxF = 23
	resL@vpXF = 0.15

	resR = True
	resR@tiYAxisString = "Volcanic Aerosol Mass Mixing Ratio"
        resR@trYMinF = -0.00000002
        resR@trYMaxF = 0.00000008
        resR@trXMinF = 0
        resR@trXMaxF = 2000
        resR@xyLineColors = "red"

        plot = gsn_csm_x2y2(wks, time(500:), time2, annual_max(500:), volc_global, resL, resR)

	; Smooth
	annual_sm = runave(annual_max(500:), 31, 0)	

	; Draw smooth series
	wks2 = gsn_open_wks(output, "/gpfsES/geo/the/MocArchieve/VOL/Smooth_Series")

	resL2 = True
	resL2@tiMainString = "Maximum of MOC (Smooth)"
	resL2@tiYAxisString = "Meridional Overturning Circulation (Sverdrups)"
	resL2@tiXAxisString = "Year"
	resL2@xyLineColors = "blue"
	resL2@vpHeightF = 0.43
	resL2@vpWidthF = 0.65
        resL2@trXMinF = 0
        resL2@trXMaxF = 2000
        resL2@trYMinF = 15
        resL2@trYMaxF = 21
	resL2@vpXF = 0.15

	resR2 = True
	resR2@tiYAxisString = "Volcanic Aerosol Mass Mixing Ratio"
        resR2@trYMinF = -0.00000002
        resR2@trYMaxF = 0.00000008
	resR2@trXMinF = 0
        resR2@trXMaxF = 2000
        resR2@xyLineColors = "red"
	
	plot2 = gsn_csm_x2y2(wks2, time(500:), time2, annual_sm, volc_global, resL2, resR2)

	; Calculate wavelet
	w = wavelet_default(annual_max(500:), 0)

	; Create coordinate arrays for plot
        dt = 1
        s0 = 2 * dt
        dj = 0.25
        N = dimsizes(annual_max(500:))
        jtot = 1 + floattointeger(((log(N * dt / s0)) / dj) / log(2.))
        power = onedtond(w@power, (/jtot, N/))
        power!0 = "period"                        ; Y axis
        power&period = w@period
        power!1 = "time"                          ; X axis
        power&time = time(500:)
        power@long_name = "Power Spectrum"
        power@units = "Sv^2"

        ; Correct bias
        if (bias_correct .eq. 1) then
                power = power / conform(power, w@scale, 0)
        end if

        ; Compute significance ( >= 1 is significant)
        SIG = power
        SIG = power / conform(power, w@signif, 0)

	; Draw wavelet transform
	wks3 = gsn_open_wks(output, "/gpfsES/geo/the/MocArchieve/VOL/Wavelet_Transform")

        res3 = True
        res3@gsnDraw = False
        res3@gsnFrame = False
        res3@tiMainString = "Wavelet Transform"
        res3@cnFillOn = True
        res3@cnLinesOn = False
        ; res3@trYReverse = True
        res3@gsnSpreadColors = True
        res3@gsnSpreadColorStart = 24
        res3@gsnSpreadColorEnd = -26
        res3@cnMinLevelValF = 0
        res3@cnMaxLevelValF = 0.6
        res3@cnLevelSpacingF = 0.05
        res3@vpHeightF = 0.43
        res3@vpWidthF = 0.65
        res3@vpXF = 0.15
        res3@tiYAxisString = "Period (year)"
        ; res3@trYMinF = 1
        res3@trYMaxF = 500
        ; res3@trGridType = "TriangularMesh"    ; This option has to be set to make the grid appear regular
        res3@cnLevelSelectionMode = "ManualLevels"

        plot3 = gsn_csm_contour(wks3, power, res3)

        res33 = True
        plot3 = ShadeCOI(wks3, plot3, w, time(500:), res33)
        draw(plot3)
        frame(wks3)

end
