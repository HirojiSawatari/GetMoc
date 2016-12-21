;**************************************************************
; thcal.ncl
;**************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin
;**************************************************************
; read in data and average over 50 yrs
;**************************************************************
yrs=(/"00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49"/)
maindir=(/"/raid2/MSS/GUIDO/csm/modtvm3-ocn/Modtvm3_O_AVG00"/)
pc3=(/".nc"/)
name=maindir+yrs(0)+pc3
in = addfile(name,"r")
v  = in->V(0,:,:,{0:359})
do i=1,49
  name=maindir+yrs(i)+pc3
  in = addfile(name,"r")
  vtemp  = in->V(0,:,:,{0:359})     ; remove cyclic points
  v=v+vtemp
end do
v=v/50.
  msk = in->rmaska
  dz = in->dz
  dx = in->dxu(0)
  lat= in->lat_t
  lon= in->lon_t
  kmt= in->kmt
;**************************************************************
; some parameters
;**************************************************************
  nzt = getfilevardimsizes(in,"z_t")   ; z_t grid (25)
  nzw = nzt+1                         ; z_w grid (26)
  nlat = dimsizes(lat)
  nlon = dimsizes(lon)
  
  d2rad = 0.017453         ; degrees to radians
;**************************************************************
; masking for basins goes here
msk2=new((/25,116,100/),float)
do i=0,115
 do j=0,99
   if(msk(i,j).eq.6.or.msk(i,j).eq.8.or.msk(i,j).eq.9.or.msk(i,j).eq.10) then
    do depth=0,24
     msk2(depth,i,j)=1.0
    end do
   else
    do depth=0,24
     msk2(depth,i,j)=0.0
    end do
   end if
 end do
end do
;printVarSummary(v)
;printVarSummary(msk2)
v=v*msk2
;**************************************************************  
; calculate first intergral 
; int[lon1:lon2]v*cos(lat)*dx*dz
; this calculation is done on the z_t grid
;**************************************************************
  zone_int = new((/nlat,nzt/),typeof(v))     ; allocate space
  do k = 0, nzt-1
    do j = 0, nlat-1
      zone_int(j,k) = dim_sum(v(k,j,:)*cos(lat(j)*d2rad)*dx*dz(k))
    end do
  end do
;************************************************************** 
; calculate second integral (partial summation) over levels on z_w grid
; psi(k,y)=int[k:0]zone_int 
;**************************************************************
  mocm      = new((/nzw,nlat/),typeof(v))      ; allocate space
  zone_int = zone_int(:,::-1)                 ; rearrange so bottom to top
  mocm(0,:) = 0.                               ; bottom is zero 
  do k=1,nzw-2
    mocm(k+1,:) = -1.0 * dim_sum(zone_int(:,0:k))
  end do
  mocm(1,:) = mocm(2,:) 
  mocm = mocm(::-1,:)                           ; put back in original order
;************************************************************** 
; scale to sverdrups
;**************************************************************
  mocm = mocm/1e12                              
;************************************************************** 
; assign meta data
;**************************************************************
  mocm!0         = "depth"
  mocm!1         = "lat"
  mocm&depth     = in->z_w
  mocm&lat       = lat
  mocm@long_name = "Atl. basin MOTSF - MOD"
  mocm@units     = "SV"
;*************************************************************
;Block out the topography
;printVarSummary(mocm)
;print(mocm@_FillValue)
do lati=0,nlat-1
 ;print("lati = "+lati)
 mxd=0
 numd=0
 do loni=0,nlon-1
   if(msk(lati,loni).eq.6.or.msk(lati,loni).eq.8.or.msk(lati,loni).eq.9.or.msk(lati,loni).eq.10) then
     if(kmt(lati,loni).gt.mxd) then
      mxd=kmt(lati,loni)
     end if
     numd=numd+1
   end if
 end do
 if(numd.eq.0)
  mxd=0
 end if
 ;print("lati = "+lati+" av depth grid pts" +mxd+" num. depths "+numd)
 if(mxd.lt.(nzw-1)) then
  do zz=mxd,nzw-1
   mocm(zz,lati)=mocm@_FillValue
;  mocm(zz,lati)=-999.
  end do
 end if
end do
  
  
;************************************************************** 
; assign meta data and take out the southern ocean
;**************************************************************
lat2=lat({-30:90})
mocm2=mocm(:,{-30:90})
  mocm2!0         = "depth"
  mocm2!1         = "lat"
  mocm2&depth     = in->z_w
  mocm2&lat       = lat2
  mocm2@long_name = "Eulerian MOTSF - MOD"
  mocm2@units     = "SV"
 
;**************************************************************
; read in data and average over 50 yrs
;**************************************************************
yrs=(/"00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49"/)
maindir=(/"/raid2/MSS/GUIDO/csm/gitvm3-ocn/gitvm3_O_AVG00"/)
pc3=(/".nc"/)
name=maindir+yrs(0)+pc3
in = addfile(name,"r")
v  = in->V(0,:,:,{0:359})
do i=1,49
  name=maindir+yrs(i)+pc3
  in = addfile(name,"r")
  vtemp  = in->V(0,:,:,{0:359})     ; remove cyclic points
  v=v+vtemp
end do
v=v/50.
  msk = in->rmaska
  dz = in->dz
  dx = in->dxu(0)
  lat= in->lat_t
  lon= in->lon_t
  kmt= in->kmt
;**************************************************************
; some parameters
;**************************************************************
  nzt = getfilevardimsizes(in,"z_t")   ; z_t grid (25)
  nzw = nzt+1                         ; z_w grid (26)
  nlat = dimsizes(lat)
  nlon = dimsizes(lon)
  
  d2rad = 0.017453         ; degrees to radians
;**************************************************************
; masking for basins goes here
msk2=new((/25,116,100/),float)
do i=0,115
 do j=0,99
   if(msk(i,j).eq.6.or.msk(i,j).eq.8.or.msk(i,j).eq.9.or.msk(i,j).eq.10) then
    do depth=0,24
     msk2(depth,i,j)=1.0
    end do
   else
    do depth=0,24
     msk2(depth,i,j)=0.0
    end do
   end if
 end do
end do
;printVarSummary(v)
;printVarSummary(msk2)
v=v*msk2
;**************************************************************  
; calculate first intergral 
; int[lon1:lon2]v*cos(lat)*dx*dz
; this calculation is done on the z_t grid
;**************************************************************
  zone_int = new((/nlat,nzt/),typeof(v))     ; allocate space
  do k = 0, nzt-1
    do j = 0, nlat-1
      zone_int(j,k) = dim_sum(v(k,j,:)*cos(lat(j)*d2rad)*dx*dz(k))
    end do
  end do
;************************************************************** 
; calculate second integral (partial summation) over levels on z_w grid
; psi(k,y)=int[k:0]zone_int 
;**************************************************************
  mocl      = new((/nzw,nlat/),typeof(v))      ; allocate space
  zone_int = zone_int(:,::-1)                 ; rearrange so bottom to top
  mocl(0,:) = 0.                               ; bottom is zero 
  do k=1,nzw-2
    mocl(k+1,:) = -1.0 * dim_sum(zone_int(:,0:k))
  end do
  mocl(1,:) = mocl(2,:) 
  mocl = mocl(::-1,:)                           ; put back in original order
;************************************************************** 
; scale to sverdrups
;**************************************************************
  mocl = mocl/1e12                              
;************************************************************** 
; assign meta data
;**************************************************************
  mocl!0         = "depth"
  mocl!1         = "lat"
  mocl&depth     = in->z_w
  mocl&lat       = lat
  mocl@long_name = "Atl. basin MOTSF - GI"
  mocl@units     = "SV"
;*************************************************************
;Block out the topography
;printVarSummary(mocl)
;print(mocl@_FillValue)
do lati=0,nlat-1
 ;print("lati = "+lati)
 mxd=0
; numd=0
 do loni=0,nlon-1
   if(msk(lati,loni).eq.6.or.msk(lati,loni).eq.8.or.msk(lati,loni).eq.9.or.msk(lati,loni).eq.10) then
    if(kmt(lati,loni).gt.mxd) then
     mxd=kmt(lati,loni)
    end if
     numd=numd+1
   end if
 end do
 if(numd.eq.0)
  mxd=0
 end if
 ;print("lati = "+lati+" av depth grid pts" +mxd+" num. depths "+numd)
 if(mxd.lt.(nzw-1)) then
  do zz=mxd,nzw-1
   mocl(zz,lati)=mocl@_FillValue
;  mocl(zz,lati)=-999.
  end do
 end if
end do
  
  
;************************************************************** 
; assign meta data and take out the southern ocean
;**************************************************************
lat2=lat({-30:90})
mocl2=mocl(:,{-30:90})
  mocl2!0         = "depth"
  mocl2!1         = "lat"
  mocl2&depth     = in->z_w
  mocl2&lat       = lat2
  mocl2@long_name = "Eulerian MOTSF - GI"
  mocl2@units     = "SV"
;************************************************************** 
; Create difference, assign meta data and take out the southern ocean
;**************************************************************

;mocdiff=mocm2-mocl2
mocdiff=mocl2-mocm2

  mocdiff!0         = "depth"
  mocdiff!1         = "lat"
  mocdiff&depth     = in->z_w
  mocdiff&lat       = lat2
  mocdiff@long_name = "Diff. in MOTSF (GI - MOD)"
  mocdiff@units     = "SV"
 
;*********************************
; create plot
;*********************************
  wks = gsn_open_wks("ps","mocatl")                 ; open a ps file
;  gsn_define_colormap(wks,"ViBlGrWhYeOrRe")      ; choose colormap
  gsn_define_colormap(wks,"BlGrYeOrReVi200")      ; choose colormap
  plot = new(3, graphic)
  res                         = True             ; plot mods desired
  res@gsnDraw = False 			         ; don't draw
  res@gsnFrame = False 				 ; don't dvance frame
  res@lbLabelBarOn = False			 ; turn off individual labels bars 

  res@cnFillOn                = True             ; turn on color fill
  res@cnLineLabelsOn          = False            ; turns off contour line labels
  res@cnLinesOn               = False            ; turn off contour lines
  res@cnInfoLabelOn           = False            ; turns off contour info label
  res@gsnSpreadColors         = True             ; use full colormap
  
  res@sfXArray                = mocm2&lat          ; uses lat_t as plot x-axis
  res@sfYArray                = mocm2&depth/100000 ; convert cm to km

  res@cnMissingValPerimOn     = True             ; turn on perimeter
  res@cnMissingValFillPattern = 3                ; choose a fill pattern
  res@cnMissingValFillColor   = "black"          ; choose a pattern color 
  res@gsnYAxisIrregular2Linear = True   ; converts irreg depth to linear
  res@gsnXAxisIrregular2Linear = True   ; converts irreg depth to linear

  res@cnLevelSelectionMode    = "ManualLevels"   ; manually set contour levels
  res@cnMinLevelValF          =  -4.             ; min level
  res@cnMaxLevelValF          =   28.             ; max level
  res@cnLevelSpacingF   = 2.                    ; interval spacing

  res@trYReverse              = True             ; reverses y-axis
  res@vpXF                   = 0.14        ; default is 0.2 (aspect ratio)
  res@vpYF                   = 0.3         ; default is 0.8
  res@vpHeightF              = 0.2         ; default is 0.6
  res@vpWidthF               = 0.75        ; default is 0.6


 

  plot(0)  = gsn_csm_contour(wks,mocm2,res)           ; create plot
  plot(1)  = gsn_csm_contour(wks,mocl2,res)           ; create plot

  res@cnMinLevelValF          =  -2.             ; min level
  res@cnMaxLevelValF          =   4.             ; max level
  res@cnLevelSpacingF         = 0.5                    ; interval spacing

  plot(2)  = gsn_csm_contour(wks,mocdiff,res)           ; create plot


;==============================================
;  create panel
;=============================================
 pres                 = True
 pres@gsnPanelRowSpec = True                   ; tell panel what order to plt
 pres@txString = "Ann. Av. MOTSF - Atlantic and Arctic basins"
 pres@gsnPanelLabelBar = True                   ; add common colorbar
 pres@lbOrientation           = "Vertical"       ; vertical label bar
 pres@lbLabelAutoStride       = True             ; nice label bar labels
 ;pres@lbLabelStride            = 2          ; every 9th label
 pres@gsnPanelXWhiteSpacePercent = 0.1

  gsn_panel(wks,plot,(/1,1,1/),pres)

end

