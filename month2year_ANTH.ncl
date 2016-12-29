load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

	; Traversal 2000 years
	do y = 1, 2000
		; Change int to string
		ystr = sprinti("%0.4i", y)

		; Open monthly files
		f1 = addfile("/gpfsES/geo/the/MocArchieve/ANTH/month/Moc.ANTH." + ystr + "-01.nc", "r")
		f2 = addfile("/gpfsES/geo/the/MocArchieve/ANTH/month/Moc.ANTH." + ystr + "-02.nc", "r")
	        f3 = addfile("/gpfsES/geo/the/MocArchieve/ANTH/month/Moc.ANTH." + ystr + "-03.nc", "r")
	        f4 = addfile("/gpfsES/geo/the/MocArchieve/ANTH/month/Moc.ANTH." + ystr + "-04.nc", "r")
	        f5 = addfile("/gpfsES/geo/the/MocArchieve/ANTH/month/Moc.ANTH." + ystr + "-05.nc", "r")
	        f6 = addfile("/gpfsES/geo/the/MocArchieve/ANTH/month/Moc.ANTH." + ystr + "-06.nc", "r")
	        f7 = addfile("/gpfsES/geo/the/MocArchieve/ANTH/month/Moc.ANTH." + ystr + "-07.nc", "r")
	        f8 = addfile("/gpfsES/geo/the/MocArchieve/ANTH/month/Moc.ANTH." + ystr + "-08.nc", "r")
	        f9 = addfile("/gpfsES/geo/the/MocArchieve/ANTH/month/Moc.ANTH." + ystr + "-09.nc", "r")
	        f10 = addfile("/gpfsES/geo/the/MocArchieve/ANTH/month/Moc.ANTH." + ystr + "-10.nc", "r")
       		f11 = addfile("/gpfsES/geo/the/MocArchieve/ANTH/month/Moc.ANTH." + ystr + "-11.nc", "r")
       		f12 = addfile("/gpfsES/geo/the/MocArchieve/ANTH/month/Moc.ANTH." + ystr + "-12.nc", "r")

		MOC = f1->MOC
		moc_z = f1->moc_z
		lat_aux_grid = f1->lat_aux_grid
	
	        MOC2 = f2->MOC
	        moc_z2 = f1->moc_z
	        lat_aux_grid2 = f1->lat_aux_grid

	        MOC3 = f3->MOC
	        moc_z3 = f1->moc_z
	        lat_aux_grid3 = f1->lat_aux_grid

	        MOC4 = f4->MOC
	        moc_z4 = f1->moc_z
	        lat_aux_grid4 = f1->lat_aux_grid

	        MOC5 = f5->MOC
	        moc_z5 = f1->moc_z
	        lat_aux_grid5 = f1->lat_aux_grid

	        MOC6 = f6->MOC
	        moc_z6 = f1->moc_z
	        lat_aux_grid6 = f1->lat_aux_grid

	        MOC7 = f7->MOC
	        moc_z7 = f1->moc_z
	        lat_aux_grid7 = f1->lat_aux_grid

	        MOC8 = f8->MOC
	        moc_z8 = f1->moc_z
	        lat_aux_grid8 = f1->lat_aux_grid

        	MOC9 = f9->MOC
        	moc_z9 = f1->moc_z
        	lat_aux_grid9 = f1->lat_aux_grid

        	MOC10 = f10->MOC
        	moc_z10 = f1->moc_z
        	lat_aux_grid10 = f1->lat_aux_grid

        	MOC11 = f11->MOC
        	moc_z11 = f1->moc_z
	        lat_aux_grid11 = f1->lat_aux_grid

        	MOC12 = f12->MOC
        	moc_z12 = f1->moc_z
        	lat_aux_grid12 = f1->lat_aux_grid

		; Create new varible to save monthly data (month, moc_z, lat_aux_grid)
		new_monthly = new((/12, dimsizes(moc_z), dimsizes(lat_aux_grid)/), float, "No_FillValue")
	
		new_monthly(0, :, :) = 	MOC(0, 1, 0, : , :)
        	new_monthly(1, :, :) =  MOC2(0, 1, 0, : , :)
        	new_monthly(2, :, :) =  MOC3(0, 1, 0, : , :)
        	new_monthly(3, :, :) =  MOC4(0, 1, 0, : , :)
        	new_monthly(4, :, :) =  MOC5(0, 1, 0, : , :)
        	new_monthly(5, :, :) =  MOC6(0, 1, 0, : , :)
        	new_monthly(6, :, :) =  MOC7(0, 1, 0, : , :)
        	new_monthly(7, :, :) =  MOC8(0, 1, 0, : , :)
        	new_monthly(8, :, :) =  MOC9(0, 1, 0, : , :)
        	new_monthly(9, :, :) =  MOC10(0, 1, 0, : , :)
        	new_monthly(10, :, :) =  MOC11(0, 1, 0, : , :)
        	new_monthly(11, :, :) =  MOC12(0, 1, 0, : , :)

		new_monthly!0 = "time"
		new_monthly!1 = "moc_z"
		new_monthly!2 = "lat_aux_grid"
	
		; Create new variables for each dimensions of new_monthly
		new_monthly&time = new((/12/), float, "No_FillValue")
		new_monthly&moc_z  = new((/dimsizes(moc_z)/), float, "No_FillValue")
		new_monthly&lat_aux_grid  = new((/dimsizes(lat_aux_grid)/), float, "No_FillValue")
	
		; Assign dimensions
		new_monthly&time = ispan(1, 12, 1)
		new_monthly&moc_z = f1->moc_z
		new_monthly&lat_aux_grid = f1->lat_aux_grid
		
		; Name the variable and unit
		new_monthly@long_name = "Meridional Overturning Circulation"
		new_monthly@units = "Sverdrups"
	
		; Create a new .nc to save monthly file
		system("rm -f /gpfsES/geo/the/MocArchieve/ANTH/monthly/Moc.ANTH.monthly." + ystr + ".nc")
		out = addfile("/gpfsES/geo/the/MocArchieve/ANTH/monthly/Moc.ANTH.monthly." + ystr + ".nc", "c")
	
		; Assign the value to out
		out->MOC = new_monthly
		out->moc_z  = moc_z
		out->lat_aux_grid  = lat_aux_grid
		
		; Create new varible to save annual mean data (0, moc_z, lat_aux_grid)
		annual = new((/1, dimsizes(moc_z), dimsizes(lat_aux_grid)/), float, "No_FillValue")
		
		; Calculate annual mean
		do i = 0, 60
			do j = 0, 104
				annual(0, i, j) = avg(new_monthly(:, i, j))
			end do
		end do

		annual!0 = "time"
                annual!1 = "moc_z"
                annual!2 = "lat_aux_grid"

                ; Create new variables for each dimensions of annual
                annual&time = new((/1/), float, "No_FillValue")
                annual&moc_z  = new((/dimsizes(moc_z)/), float, "No_FillValue")
                annual&lat_aux_grid  = new((/dimsizes(lat_aux_grid)/), float, "No_FillValue")

                ; Assign dimensions
                annual&time = ispan(1, 1, 1)
                annual&moc_z = f1->moc_z
                annual&lat_aux_grid = f1->lat_aux_grid

                ; Name the variable and unit
                annual@long_name = "Annual Mean of Meridional Overturning Circulation"
                annual@units = "Sverdrups"

		; Create a new .nc to save monthly file
                system("rm -f /gpfsES/geo/the/MocArchieve/ANTH/annual/Moc.ANTH.annual." + ystr + ".nc")
                out2 = addfile("/gpfsES/geo/the/MocArchieve/ANTH/annual/Moc.ANTH.annual." + ystr + ".nc", "c")

		; Assign the value to out
                out2->MOC = annual
                out2->moc_z  = moc_z
                out2->lat_aux_grid  = lat_aux_grid

	end do
end
