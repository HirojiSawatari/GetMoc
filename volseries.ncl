;load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

	; Create new varible to save volcanic data
	volc_global = new((/24000/), double, "No_FillValue")
	volc_north = new((/24000/), double, "No_FillValue")
	volc_south = new((/24000/), double, "No_FillValue")

	; Open volcanic file
        v = addfile("/gpfsES/geo/zywang/Rec_Obs_For_DATA/forcings/volc/IVI2.5_Gao_Sigl_1-2000.nc", "r")
	
	; volc(time=24002, lev=18, lat=64)
	volc = v->MMRVOLC
	time = v->time
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

	; Draw global series
	wks = gsn_open_wks("png","/gpfsES/geo/the/MocArchieve/Volc_Global_Series")
            
        res = True
	res@tiMainString = "Global"
        res@tiYAxisString = "Volcanic Aerosol Mass Mixing Ratio"
        res@tiXAxisString = "Year"
        res@trYMinF = 0
        res@trYMaxF = 0.0000001
        res@xyLineColors = "red"
        res@vpHeightF = 0.43
        res@vpWidthF = 0.70

	plot = gsn_xy(wks, time(0:23999), volc_global, res)
	
	; Draw north series
        wks = gsn_open_wks("png","/gpfsES/geo/the/MocArchieve/Volc_North_Series")

        res = True
	res@tiMainString = "Northern Hemisphere"
        res@tiYAxisString = "Volcanic Aerosol Mass Mixing Ratio"
        res@tiXAxisString = "Year"
        res@trYMinF = 0
        res@trYMaxF = 0.0000001
        res@xyLineColors = "black"
        res@vpHeightF = 0.43
        res@vpWidthF = 0.70

        plot = gsn_xy(wks, time(0:23999), volc_north, res)

	; Draw south series
        wks = gsn_open_wks("png","/gpfsES/geo/the/MocArchieve/Volc_South_Series")

        res = True
	res@tiMainString = "Southern Hemisphere"
        res@tiYAxisString = "Volcanic Aerosol Mass Mixing Ratio"
        res@tiXAxisString = "Year"
        res@trYMinF = 0
        res@trYMaxF = 0.0000001
        res@xyLineColors = "blue"
        res@vpHeightF = 0.43
        res@vpWidthF = 0.70

        plot = gsn_xy(wks, time(0:23999), volc_south, res)

end
