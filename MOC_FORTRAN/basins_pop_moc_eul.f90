program poptrans
  ! Calculate POP transports, ie, MOC and heat transports, from 
  ! PCM1/PCTM/CCSM2 POP history tapes
  !
  implicit none
  include 'netcdf.inc'
  !
  ! f90 poptrans.f90 -L/contrib/lib -I/contrib/include -lnetcdf
  ! modified by tc 3/1999 - psi computed from the w method
  ! modified by gs 6/1999 - f90 and output is now netCDF
  ! modified by gs 1/2000 - do each level of U,V,T one at a time to fit
  !                         into memory of T3E
  ! modified by gs 4/2001 - get grid info from pcm1_grid.nc file
  ! modified by ecb 5/2003 - add arctic, hudsons bay, med to amask
  !                        - added namelist input
  !                        - make gridfile = infile 
  !                        - use counts to put in spval in arrays
  !
  ! added Atlantic, Pacific, Indian MOC after S. Yaeger (CW, 5/2008)
  ! also using UET, VNT, UES, VNS if existing in input NetCDF file
  !
  !gx1v3:
  !integer,parameter::nphi=360
  !integer,parameter::nphi=91
  integer :: nphi=181
  !integer,parameter::nd =  3,nv = 5,ndv=nd+nv
  integer,parameter::nd =  3,nv = 20,ndv=nd+nv
  !real   ,parameter::phi_start=-89.75,phi_end=89.75
  real ::phi_start=-89.75,phi_end=89.75
  real   ,parameter::spval=1.e30
  real   ,parameter::cal2joule = 4.1868
  !
  integer::varxtype,varndims,varnatts,status,itime,imt,jmt,km
  integer::ocnncid,tranncid,ndims,nvars,ngatts,unlimdimid,timelen
  integer::latd,levd,timd,vdim,varid,gridid
  integer::i,j,k,l,m,n,jc,ip1,im1,jp1,jm1,iv
  integer,dimension(ndv)::varlev,vrdims
  integer,dimension(4)::istart,icount,vardimids
  integer,dimension(3)::ostart,ocount,vdims
  integer,dimension(4)::nc_dimlen
  !
  real::pi,r2d,dphi
  real::wtkb,wtk,fuew,fuww,fvnw,fvsw,work
  real::fuetw,fuwtw,fvntw,fvstw,fuesw,fuwsw,fvnsw,fvssw
  integer,dimension(:,:),allocatable::amaskt,imaskt,pmaskt,kmu,kmt,region_mask
  integer,dimension(:,:),allocatable::cnta,cntp,cnti,cntg
  real,dimension(:,:),allocatable::ulat,htn,hte,dxu,dyu,ulon,tlat,tarea
  real,dimension(:,:),allocatable::u,v,t,s,uet,vnt,ues,vns
  real,dimension(:,:),allocatable::gpsiw,apsiw,ipsiw,ppsiw
  real,dimension(:,:),allocatable::gmhtw,amhtw,imhtw,pmhtw
  real,dimension(:,:),allocatable::gmhsw,amhsw,imhsw,pmhsw
  real,dimension(:,:,:),allocatable::w,htw,fue,fvn,fuet,fvnt,fues,fvns,hsw
  real,dimension(:),allocatable::psia,mhta,psii,mhti,psip,mhtp,mhsa,mhsi,mhsp
  real,dimension(:),allocatable::dz,depth
  real,dimension(:),allocatable::latr,htrans
  real,dimension(120)::time
  !
  character(len=132)::varname,infile,title,ncfile,nc_dimname,gridfile
  character(len=132),dimension(ndv)::varnam,vrunit,vrlong,vrtype,vrgrid
  character(len=132),dimension(ndv)::varstr
  !
  logical::does_exist,uetexist,vntexist,uesexist,vnsexist
  !
  namelist /io_info/ infile
  namelist /grid_info/ nphi, phi_start, phi_end

  varstr = (/&
       !VARNAM     Lv VRLONG                             T  D  Gr VRUNITS
       'lat         1 Latitude                           F  1  XX degrees_north        ',&
       'lev         1 Level                              F  1  XX m                    ',&
       'time        1 Time                               F  1  XX days since 0000-01-01',&
       'gmsf        1 Global meridional streamfunction   F  3  TS Sv (10^6 m^3/s)      ',&
       'amsf        1 Atlantic meridional streamfunction F  3  TS Sv (10^6 m^3/s)      ',&
       'imsf        1 Indian meridional streamfunction   F  3  TS Sv (10^6 m^3/s)      ',&
       'pmsf        1 Pacific meridional streamfunction  F  3  TS Sv (10^6 m^3/s)      ',&
       'gmht        1 Global meridional heat transport   F  2  TS Pw                   ',&
       'amht        1 Atlantic meridional heat transport F  2  TS Pw                   ',&
       'imht        1 Indian meridional heat transport   F  2  TS Pw                   ',&
       'pmht        1 Pacific meridional heat transport  F  2  TS Pw                   ',&
       'gmhs        1 Global meridional salt transport   F  2  TS ppt*Sv               ',&
       'amhs        1 Atlantic meridional salt transport F  2  TS ppt*Sv               ',&
       'imhs        1 Indian meridional salt transport   F  2  TS ppt*Sv               ',&
       'pmhs        1 Pacific meridional salt transport  F  2  TS ppt*Sv               ',&
       'gmhtw       1 Global heat transport              F  3  TS Pw                   ',&
       'amhtw       1 Atlantic heat transpor             F  3  TS Pw                   ',&
       'imhtw       1 Indian heat transport              F  3  TS Pw                   ',&
       'pmhtw       1 Pacific heat transport             F  3  TS Pw                   ',&
       'gmhsw       1 Global salt transport              F  3  TS ppt*Sv               ',&
       'amhsw       1 Atlantic salt transport            F  3  TS ppt*Sv               ',&
       'imhsw       1 Indian salt transport              F  3  TS ppt*Sv               ',&
       'pmhsw       1 Pacific salt transport             F  3  TS ppt*Sv               '/)
  !
  !=======================================================================


  uesexist = .false.
  uetexist = .false.
  vnsexist = .false.
  vntexist = .false.

  read(*,io_info)
  write(*,io_info)
  read(*,grid_info)
  write(*,grid_info)
! write(*,'(''Input file: '',$)')
! read (*,'(a)') infile      
  write(*,*)
  inquire(file=trim(infile),exist=does_exist)
  if (.not.(does_exist)) then
     write(*,*) 'Ocean tape ',trim(infile),' missing. Stop.'
     stop
  endif
  !
  status = nf_open(infile,nf_nowrite,ocnncid)
  if (status /= nf_noerr) write(*,*) nf_strerror(status)
  status = nf_inq(ocnncid,ndims,nvars,ngatts,unlimdimid)
  if (status /= nf_noerr) write(*,*) nf_strerror(status)

  ! Find out I and J and K dimensions
  do i = 1,ndims
     status = nf_inq_dim(ocnncid,i,nc_dimname,nc_dimlen(i))
     if (status /= nf_noerr) write(*,*) nf_strerror(status)
     if (trim(nc_dimname) == 'lon')  imt = nc_dimlen(i)
     if (trim(nc_dimname) == 'lat')  jmt = nc_dimlen(i)
     if (trim(nc_dimname) == 'lev')   km = nc_dimlen(i)
     if (trim(nc_dimname) == 'nlon') imt = nc_dimlen(i)
     if (trim(nc_dimname) == 'nlat') jmt = nc_dimlen(i)
     if (trim(nc_dimname) == 'z_t')   km = nc_dimlen(i)
  enddo
  !
  allocate(amaskt(imt,jmt),imaskt(imt,jmt),pmaskt(imt,jmt))
  allocate(kmu(imt,jmt),kmt(imt,jmt),region_mask(imt,jmt))
  allocate(ulat(imt,jmt),htn(imt,jmt),hte(imt,jmt))
  allocate(dxu(imt,jmt),dyu(imt,jmt),ulon(imt,jmt),tlat(imt,jmt),tarea(imt,jmt))
  allocate(u(imt,jmt),v(imt,jmt),t(imt,jmt),s(imt,jmt))
  allocate(uet(imt,jmt),vnt(imt,jmt),ues(imt,jmt),vns(imt,jmt))
  allocate(gpsiw(km,nphi),apsiw(km,nphi),gmhtw(km,nphi),amhtw(km,nphi))
  allocate(ipsiw(km,nphi),ppsiw(km,nphi),imhtw(km,nphi),pmhtw(km,nphi))
  allocate(cntg(km,nphi),cnta(km,nphi),cntp(km,nphi),cnti(km,nphi))
  allocate(gmhsw(km,nphi),amhsw(km,nphi),imhsw(km,nphi),pmhsw(km,nphi))
  allocate(w(imt,jmt,km),htw(imt,jmt,km),fue(imt,jmt,km),fvn(imt,jmt,km),hsw(imt,jmt,km))
  allocate(fuet(imt,jmt,km),fvnt(imt,jmt,km),fues(imt,jmt,km),fvns(imt,jmt,km))
  allocate(psia(km),mhta(km),psii(km),mhti(km),psip(km),mhtp(km))
  allocate(mhsa(km),mhsi(km),mhsp(km))
  allocate(dz(km),depth(km))
  allocate(latr(nphi),htrans(nphi))
  gpsiw = 0.;apsiw = 0.;ipsiw = 0.;ppsiw = 0.
  cntg = 0;cnta= 0;cntp= 0;cnti= 0
  gmhtw = 0.;amhtw = 0.;imhtw = 0.;pmhtw = 0.
  gmhsw = 0.;amhsw = 0.;imhsw = 0.;pmhsw = 0.

  ! Find out how many time samples there are
  do i = 1,nvars
     status = nf_inq_var(ocnncid,i,varname,varxtype,varndims,vardimids,varnatts)
     if (status /= nf_noerr) write(*,*) nf_strerror(status)
     if ((trim(varname) == 'time').or.(trim(varname) ==  'TIME')) then
        status = nf_inq_dim(ocnncid,unlimdimid,varname,timelen)
        if (status /= nf_noerr) write(*,*) nf_strerror(status)
        status = nf_get_var_real(ocnncid,i,time)
        if (status /= nf_noerr) write(*,*) nf_strerror(status)
     endif
  enddo
  !
  istart =(/   1,  1,  1,  1/)
  icount =(/ imt,jmt,  1,  1/)
  ostart =(/   1,  1,  1/)
  ocount =(/nphi, km,  1/)
  !
  pi=atan(1.)*4.
  r2d=180./pi
  ! Get grid info
! if (jmt == 288) gridfile = 'pcm1_grid.nc'
! if (jmt == 320) gridfile = 'pcm2_grid.nc'
! if (jmt == 384) gridfile = 'ccsm2_ocn.nc'
  gridfile = infile 
  inquire(file=trim(gridfile),exist=does_exist)
  if (.not.(does_exist)) then
     write(*,*) 'Need to link in ',trim(gridfile),' as follows: '
     write(*,*) 'ln -s /fs/cgd/data0/strandwg/ferret/data/',trim(gridfile),' .'
     write(*,*) 'Stopping.'
     stop
  endif

  status = nf_open(trim(gridfile),nf_nowrite,gridid)
  if (status /= nf_noerr) write(*,*) nf_strerror(status)
  status = nf_inq(gridid,ndims,nvars,ngatts,unlimdimid)
  if (status /= nf_noerr) write(*,*) nf_strerror(status)
  do i = 1,nvars
     status = nf_inq_var(gridid,i,varname,varxtype,varndims,vardimids,varnatts)
     if (status /= nf_noerr) write(*,*) nf_strerror(status)
     select case (trim(varname))
     case ('kmt','KMT')
        status = nf_get_var_int(gridid,i,kmt)
        if (status /= nf_noerr) write(*,*) nf_strerror(status)
     case ('amaskt')
        status = nf_get_var_int(gridid,i,amaskt)
        if (status /= nf_noerr) write(*,*) nf_strerror(status)
     case ('imaskt')
        status = nf_get_var_int(gridid,i,imaskt)
        if (status /= nf_noerr) write(*,*) nf_strerror(status)
     case ('pmaskt')
        status = nf_get_var_int(gridid,i,pmaskt)
        if (status /= nf_noerr) write(*,*) nf_strerror(status)
     case ('REGION_MASK')
        status = nf_get_var_int(gridid,i,region_mask)
        if (status /= nf_noerr) write(*,*) nf_strerror(status)
        amaskt = 0 ; imaskt = 0 ; pmaskt = 0
        where (region_mask == 2) ; pmaskt = 1 ; endwhere
        where (region_mask == 3) ; imaskt = 1 ; endwhere
        where (region_mask == 6) ; amaskt = 1 ; endwhere
!        where (region_mask == 8) ; amaskt = 1 ; endwhere
!        where (region_mask == 9) ; amaskt = 1 ; endwhere
!        where (region_mask == 10) ; amaskt = 1 ; endwhere
!        where (region_mask == 11) ; amaskt = 1 ; endwhere
     case ('lat','ULAT')
        status = nf_get_var_real(gridid,i,ulat)
        if (status /= nf_noerr) write(*,*) nf_strerror(status)
        ulat = ulat / r2d
     case ('lon','ULON')
        status = nf_get_var_real(gridid,i,ulon)
        if (status /= nf_noerr) write(*,*) nf_strerror(status)
        ulon = ulon / r2d
     case ('htn','HTN')
        status = nf_get_var_real(gridid,i,htn)
        if (status /= nf_noerr) write(*,*) nf_strerror(status)
     case ('hte','HTE')
        status = nf_get_var_real(gridid,i,hte)
        if (status /= nf_noerr) write(*,*) nf_strerror(status)
     case ('dz')
        status = nf_get_var_real(gridid,i,dz)
        if (status /= nf_noerr) write(*,*) nf_strerror(status)
     case default
     end select
  enddo
  status = nf_close(gridid)
  if (status /= nf_noerr) write(*,*) nf_strerror(status)

  do j=1,jmt
     do i=1,imt
        im1=mod(i-2+imt,imt)+1
        ip1=mod(i,imt)+1
        jm1=mod(j-2+jmt,jmt)+1
        jp1=mod(j,jmt)+1
        kmu(i,j)=min(kmt(i,j),kmt(ip1,j),kmt(i,jp1),kmt(ip1,jp1))
        tlat(i,j)=0.25*(ulat(i,j)+ulat(im1,j)+ulat(i,jm1)+ulat(im1,jm1))*r2d
        tarea(i,j)=0.25*(htn(i,j)+htn(i,jm1))*(hte(i,j)+hte(im1,j))
        dxu(i,j)=0.5*(htn(i,j)+htn(ip1,j))
        dyu(i,j)=0.5*(hte(i,j)+hte(i,jp1))
     enddo
  enddo

  depth(1)=dz(1)
  do k=2,km
     depth(k)=(depth(k-1)+dz(k))
  enddo
  depth = depth/100 ! cm to m

  dphi = (phi_end - phi_start)/float(nphi - 1)
  do n = 1,nphi
     latr(n) = (phi_start+float(n-1)*dphi)
  enddo

  ! Create new netCDF file for transports
  ncfile = 'TVTS_'//trim(infile)
  do i = 1,ndv
     read(varstr(i),'(a10,i3,1x,a34,1x,a1,2x,i1,2x,a2,1x,a)') &
          varnam(i),varlev(i),vrlong(i),vrtype(i),vrdims(i),vrgrid(i),vrunit(i)
  enddo
  status = nf_create(trim(ncfile),nf_clobber,tranncid)
  !
  status = nf_def_dim(tranncid,varnam(1),nphi,latd)
  status = nf_def_dim(tranncid,varnam(2),km,levd)
  status = nf_def_dim(tranncid,varnam(3),NF_UNLIMITED,timd)
  !
  do i = 1,nd
     if (i == 1) vdim = latd
     if (i == 2) vdim = levd
     if (i == 3) vdim = timd
     status = nf_def_var(tranncid,varnam(i),ncfloat,vrdims(i),vdim,varid)
     status = nf_put_att_text(tranncid,varid,'units',len_trim(vrunit(i)),vrunit(i))
     status = nf_put_att_text(tranncid,varid,'long_name',len_trim(vrlong(i)),vrlong(i))
     if (i == 2) status = nf_put_att_text(tranncid,varid,'positive',4,'down')
  enddo
  !
  do i = 1,nv
     j = (i+nd)
     vdims(1) = latd
     if (vrdims(j) == 3) then
        vdims(2) = levd
        vdims(3) = timd
     else
        vdims(2) = timd
     endif
     status = nf_def_var(tranncid,varnam(j),ncfloat,vrdims(j),vdims,varid)
     status = nf_put_att_text(tranncid,varid,'units',len_trim(vrunit(j)),vrunit(j))
     status = nf_put_att_text(tranncid,varid,'long_name',len_trim(vrlong(j)),vrlong(j))
     status = nf_put_att_real(tranncid,varid,'missing_value',ncfloat,1,spval)
     status = nf_put_att_real(tranncid,varid,'_FillValue',ncfloat,1,spval)
  enddo
  !
  title  = 'Transports from '//trim(infile)
  status = nf_put_att_text(tranncid,nf_global,'title',len_trim(title),title)
  status = nf_enddef(tranncid)
  !
  status = nf_inq_varid(tranncid,trim(varnam(1)),varid)
  status = nf_put_var_real(tranncid,varid,latr)
  if (status /= nf_noerr) write(*,*) nf_strerror(status)
  status = nf_inq_varid(tranncid,trim(varnam(2)),varid)
  status = nf_put_var_real(tranncid,varid,depth)
  if (status /= nf_noerr) write(*,*) nf_strerror(status)
  ! Now do the calculation 
  status = nf_inq(ocnncid,ndims,nvars,ngatts,unlimdimid)
  if (status /= nf_noerr) write(*,*) nf_strerror(status)
  do itime = 1,timelen
     do k=1,km
        ! get U, V, and T
        u = spval ; v = spval ; t = spval; s = spval
        istart(3) = k ; istart(4) = itime
        !
        do iv = 1,nvars
           status = nf_inq_var(ocnncid,iv,varname,varxtype,varndims,vardimids,varnatts)
           if (status /= nf_noerr) write(*,*) nf_strerror(status)
           select case(trim(varname))
              case ('temp','TEMP')
                 status = nf_get_vara_real(ocnncid,iv,istart,icount,t)
                 if (status /= nf_noerr) write(*,*) nf_strerror(status)
              case ('salt','SALT')
                 status = nf_get_vara_real(ocnncid,iv,istart,icount,s)
                 if (status /= nf_noerr) write(*,*) nf_strerror(status)
              case ('u','U','UVEL')
                 status = nf_get_vara_real(ocnncid,iv,istart,icount,u)
                 if (status /= nf_noerr) write(*,*) nf_strerror(status)
              case ('v','V','VVEL')
                 status = nf_get_vara_real(ocnncid,iv,istart,icount,v)
                 if (status /= nf_noerr) write(*,*) nf_strerror(status)
              case ('vnt','VNT')
                 vntexist = .true.
                 status = nf_get_vara_real(ocnncid,iv,istart,icount,vnt)
                 if (status /= nf_noerr) write(*,*) nf_strerror(status)
                 vnt(:,:) = vnt(:,:)*tarea*dz(k)
              case ('uet','UET')
                 uetexist = .true.
                 status = nf_get_vara_real(ocnncid,iv,istart,icount,uet)
                 if (status /= nf_noerr) write(*,*) nf_strerror(status)
                 uet(:,:) = uet(:,:)*tarea*dz(k)
              case ('vns','VNS')
                 vnsexist = .true.
                 status = nf_get_vara_real(ocnncid,iv,istart,icount,vns)
                 if (status /= nf_noerr) write(*,*) nf_strerror(status)
                 vns(:,:) = vns(:,:)*tarea*dz(k)
              case ('ues','UES')
                 uesexist = .true.
                 status = nf_get_vara_real(ocnncid,iv,istart,icount,ues)
                 if (status /= nf_noerr) write(*,*) nf_strerror(status)
                 ues(:,:) = ues(:,:)*tarea*dz(k)
              case default
           end select
        enddo

        ! Clean them up a wee bit
        where (abs(u)>1.e10) ; u = 0. ; end where
        where (abs(v)>1.e10) ; v = 0. ; end where
        where (abs(t)>1.e10) ; t = 0. ; end where
        where (abs(s)>1.e10) ; s = 0. ; end where
        where (kmu<k) ; u(:,:) = 0. ; v(:,:) = 0. ; end where
        where (kmt<k) ; t(:,:) = 0. ; s(:,:) = 0. ; end where

        ! precompute flux information
        ! compute dxu,dyu,fue,fvn like wcalc in pop
        ! fue, fvn is the flux through the east and north edge, u*dy or v*dx
        do j=1,jmt
           do i=1,imt
              im1=mod(i-2+imt,imt)+1
              ip1=mod(i,imt)+1
              jm1=mod(j-2+jmt,jmt)+1
              jp1=mod(j,jmt)+1
              !
              fue(i,j,k)=0.5*(u(i,j)*dyu(i,j)+u(i,jm1)*dyu(i,jm1))
              fvn(i,j,k)=0.5*(v(i,j)*dxu(i,j)+v(im1,j)*dxu(im1,j))
              if (uetexist.and.vntexist) then
                fuet(i,j,k)=uet(i,j)
                fvnt(i,j,k)=vnt(i,j)
              else
                fuet(i,j,k)=fue(i,j,k)*0.5*(t(i,j)+t(ip1,j))*dz(k)
                fvnt(i,j,k)=fvn(i,j,k)*0.5*(t(i,j)+t(i,jp1))*dz(k)
              endif
              if (uesexist.and.vnsexist) then
                fues(i,j,k)=ues(i,j)
                fvns(i,j,k)=vns(i,j)
              else
                fues(i,j,k)=fue(i,j,k)*0.5*(s(i,j)+s(ip1,j))*dz(k)
                fvns(i,j,k)=fvn(i,j,k)*0.5*(s(i,j)+s(i,jp1))*dz(k)
              endif
           enddo
        enddo
     enddo ! do k = 1,km
     ! ------------------------------
     ! new w method
     ! ------------------------------
     ! calc psia,mhta for amaskt,imaskt,pmaskt
     mhta=0.;psia=0.;mhti=0.;psii=0.;mhtp=0.;psip=0.
     mhsa=0.;mhsi=0.;mhsp=0.
     do j=1,jmt
        do i=1,imt
           if (amaskt(i,j)>0) goto 720
        enddo
     enddo
720  continue
     jc=j-1
     do i=1,imt
        im1=mod(i-2+imt,imt)+1
        if (amaskt(i,jc+1)>0) then
           do k=1,km
              fvnw=fvn(i,jc,k)
              fvntw=fvnt(i,jc,k)
	      fvnsw=fvns(i,jc,k)
              psia(k)=psia(k)+fvnw
              mhta(k)=mhta(k)+fvntw
	      mhsa(k)=mhsa(k)+fvnsw
           enddo
        endif
        if (imaskt(i,jc+1)>0) then
           do k=1,km
              fvnw=fvn(i,jc,k)
              fvntw=fvnt(i,jc,k)
              fvnsw=fvns(i,jc,k)
              psii(k)=psii(k)+fvnw
              mhti(k)=mhti(k)+fvntw
              mhsi(k)=mhsi(k)+fvnsw
           enddo
        endif
        if (pmaskt(i,jc+1)>0) then
           do k=1,km
              fvnw=fvn(i,jc,k)
              fvntw=fvnt(i,jc,k)
              fvnsw=fvns(i,jc,k)
              psip(k)=psip(k)+fvnw
              mhtp(k)=mhtp(k)+fvntw
              mhsp(k)=mhsp(k)+fvnsw
           enddo
        endif
     enddo
     psia(km)=-dz(km)*psia(km)
     psii(km)=-dz(km)*psii(km)
     psip(km)=-dz(km)*psip(km)
     do k=km-1,1,-1
        psia(k)=psia(k+1)-dz(k)*psia(k)
        psii(k)=psii(k+1)-dz(k)*psii(k)
        psip(k)=psip(k+1)-dz(k)*psip(k)
     enddo
     !
     do j=1,jmt
        do i=1,imt
           wtk=0.
           wtkb=0.
           im1=mod(i-2+imt,imt)+1
           ip1=mod(i,imt)+1
           jm1=mod(j-2+jmt,jmt)+1
           jp1=mod(j,jmt)+1
           do k=km,1,-1
              fuew=fue(i,j,k)
              fuww=fue(im1,j,k)
              fvnw=fvn(i,j,k)
              fvsw=fvn(i,jm1,k)
              fuetw=fuet(i,j,k)
              fuwtw=fuet(im1,j,k)
              fvntw=fvnt(i,j,k)
              fvstw=fvnt(i,jm1,k)
              fuesw=fues(i,j,k)
              fuwsw=fues(im1,j,k)
              fvnsw=fvns(i,j,k)
              fvssw=fvns(i,jm1,k)
              work=(fvnw-fvsw+fuew-fuww)/tarea(i,j)
              htw(i,j,k)=(fvntw-fvstw+fuetw-fuwtw)
	      hsw(i,j,k)=(fvnsw-fvssw+fuesw-fuwsw)
              if (k<=kmt(i,j)) then
                 wtk = wtkb - dz(k)*work
              else
                 wtk = 0.
              endif
              w(i,j,k)=wtk
              wtkb = wtk
           enddo
        enddo
     enddo
     !
     gpsiw(:,1)=0.;gmhtw(:,1)=0.;gmhsw(:,1)=0.
     apsiw(:,1)=0.;amhtw(:,1)=0.;amhsw(:,1)=0.
     ipsiw(:,1)=0.;imhtw(:,1)=0.;imhsw(:,1)=0.
     ppsiw(:,1)=0.;pmhtw(:,1)=0.;pmhsw(:,1)=0.            
     !
     do n = 2,nphi
        do k=1,km
           gpsiw(k,n)=gpsiw(k,n-1);gmhtw(k,n)=gmhtw(k,n-1);gmhsw(k,n)=gmhsw(k,n-1)
           apsiw(k,n)=apsiw(k,n-1);amhtw(k,n)=amhtw(k,n-1);amhsw(k,n)=gmhsw(k,n-1)
           ipsiw(k,n)=ipsiw(k,n-1);imhtw(k,n)=imhtw(k,n-1);imhsw(k,n)=gmhsw(k,n-1)
           ppsiw(k,n)=ppsiw(k,n-1);pmhtw(k,n)=pmhtw(k,n-1);pmhsw(k,n)=gmhsw(k,n-1)
        enddo
        !
        do j=1,jmt
           do i=1,imt
              if (tlat(i,j)<latr(n).and.tlat(i,j)>=latr(n-1)) then
                 do k=1,km
                    if (k<=kmt(i,j)) then
                       gpsiw(k,n)=gpsiw(k,n)+w(i,j,k)*tarea(i,j)
                       gmhtw(k,n)=gmhtw(k,n)+htw(i,j,k)
		       gmhsw(k,n)=gmhsw(k,n)+hsw(i,j,k)
                       cntg(k,n) =cntg(k,n) + 1
                       if (amaskt(i,j) == 1) then
                          apsiw(k,n)=apsiw(k,n)+w(i,j,k)*tarea(i,j)
                          amhtw(k,n)=amhtw(k,n)+htw(i,j,k)
                          amhsw(k,n)=amhsw(k,n)+hsw(i,j,k)
                          cnta(k,n) =cnta(k,n) + 1
                       endif
                       if (imaskt(i,j) == 1) then
                          ipsiw(k,n)=ipsiw(k,n)+w(i,j,k)*tarea(i,j)
                          imhtw(k,n)=imhtw(k,n)+htw(i,j,k)
                          imhsw(k,n)=imhsw(k,n)+hsw(i,j,k)
                          cnti(k,n) =cnti(k,n) + 1
!                          if (i == 222) then
!                             if ((j >= 87).and.(j <= 118)) then
!                                write(*,'(3i5,2e18.8)') j,k,n,ipsiw(k,n),imhtw(k,n)
!                             endif
!                          endif
                       endif
                       if (pmaskt(i,j) == 1) then
                          ppsiw(k,n)=ppsiw(k,n)+w(i,j,k)*tarea(i,j)
                          pmhtw(k,n)=pmhtw(k,n)+htw(i,j,k)
                          pmhsw(k,n)=pmhsw(k,n)+hsw(i,j,k)
                          cntp(k,n) =cntp(k,n) + 1
                       endif
                    endif
                 enddo ! do k=1,km
              endif
           enddo ! do i=1,imt
        enddo ! do j=1,jmt
     enddo ! do n=1,nphi
     ! Fix 'em up
     do n = 1,nphi
        do k = 1,km
           if ((apsiw(k,n)/=0.).or.(apsiw(1,n)/=0.).or.(apsiw(2,n)/=0.)) then
              apsiw(k,n)=apsiw(k,n)+psia(k)
              amhtw(k,n)=amhtw(k,n)+mhta(k)
              amhsw(k,n)=amhsw(k,n)+mhsa(k)
           endif
           if ((ipsiw(k,n)/=0.).or.(ipsiw(1,n)/=0.).or.(ipsiw(2,n)/=0.)) then
              ipsiw(k,n)=ipsiw(k,n)+psii(k)
              imhtw(k,n)=imhtw(k,n)+mhti(k)
              imhsw(k,n)=imhsw(k,n)+mhsi(k)
           endif
           if ((ppsiw(k,n)/=0.).or.(ppsiw(1,n)/=0.).or.(ppsiw(2,n)/=0.)) then
              ppsiw(k,n)=ppsiw(k,n)+psip(k)
              pmhtw(k,n)=pmhtw(k,n)+mhtp(k)
              pmhsw(k,n)=pmhsw(k,n)+mhsp(k)
           endif
        enddo
     enddo

     ! Transports cannot exist outside north & south bounds of basin masks
     !write(*,*)'max val of amaskt*tlat= ',maxval(tlat*amaskt)
     do n = 1,nphi
        if (latr(n)>maxval(tlat*amaskt)) then
           apsiw(:,n) = 0;amhtw(:,n) = 0;amhsw(:,n) = 0
        endif
        if (latr(n)<minval(tlat*amaskt)) then
           apsiw(:,n) = 0;amhtw(:,n) = 0;amhsw(:,n) = 0
        endif
        if (latr(n)>maxval(tlat*imaskt)) then 
           ipsiw(:,n) = 0;imhtw(:,n) = 0;imhsw(:,n) = 0
        endif
        if (latr(n)<minval(tlat*imaskt)) then 
           ipsiw(:,n) = 0;imhtw(:,n) = 0;imhsw(:,n) = 0
        endif
        if (latr(n)>maxval(tlat*pmaskt)) then 
           ppsiw(:,n) = 0;pmhtw(:,n) = 0;pmhsw(:,n) = 0
        endif
        if (latr(n)<minval(tlat*pmaskt)) then 
           ppsiw(:,n) = 0;pmhtw(:,n) = 0;pmhsw(:,n) = 0
        endif
     enddo
     ! Block tiny values
!    where (abs(gpsiw) <= 1.e8) ; gpsiw = 0 ; end where
!    where (abs(gmhtw) <= 1.e7) ; gmhtw = 0 ; end where
!    where (abs(gmhsw) <= 1.e7) ; gmhsw = 0 ; end where
!    where (abs(apsiw) <= 1.e9) ; apsiw = 0 ; end where
!    where (abs(amhtw) <= 1.e9) ; amhtw = 0 ; end where
!    where (abs(amhsw) <= 1.e9) ; amhsw = 0 ; end where
!    where (abs(ipsiw) <= 1.e8) ; ipsiw = 0 ; end where
!    where (abs(imhtw) <= 1.e8) ; imhtw = 0 ; end where
!    where (abs(imhsw) <= 1.e8) ; imhsw = 0 ; end where
!    where (abs(ppsiw) <= 1.e8) ; ppsiw = 0 ; end where
!    where (abs(pmhtw) <= 1.e8) ; pmhtw = 0 ; end where
!    where (abs(pmhsw) <= 1.e8) ; pmhsw = 0 ; end where

!    where (cntg == 0) ; gpsiw = 0 ; end where 
!    where (cntg == 0) ; gmhtw = 0 ; end where
!    where (cntg == 0) ; gmhsw = 0 ; end where

!    where (cnta == 0) ; apsiw = 0 ; end where
!    where (cnta == 0) ; amhtw = 0 ; end where
!    where (cnta == 0) ; amhsw = 0 ; end where

!    where (cnti == 0) ; ipsiw = 0 ; end where
!    where (cnti == 0) ; imhtw = 0 ; end where
!    where (cnti == 0) ; imhsw = 0 ; end where

!    where (cntp == 0) ; ppsiw = 0 ; end where
!    where (cntp == 0) ; pmhtw = 0 ; end where
!    where (cntp == 0) ; pmhsw = 0 ; end where

     ! Mask and convert to MKS
     ! cgs to Sverdrups
     where (cntg == 0) ; gpsiw = spval ; elsewhere ; gpsiw = gpsiw * 1e-12 ; end where
     where (cnta == 0) ; apsiw = spval ; elsewhere ; apsiw = apsiw * 1e-12 ; end where
     where (cnti == 0) ; ipsiw = spval ; elsewhere ; ipsiw = ipsiw * 1e-12 ; end where
     where (cntp == 0) ; ppsiw = spval ; elsewhere ; ppsiw = ppsiw * 1e-12 ; end where
     ! cgs to Petawatts
     where (cntg == 0) ; gmhtw = spval ; elsewhere ; gmhtw = gmhtw * cal2joule * 1e-15 ; end where
     where (cnta == 0) ; amhtw = spval ; elsewhere ; amhtw = amhtw * cal2joule * 1e-15 ; end where
     where (cnti == 0) ; imhtw = spval ; elsewhere ; imhtw = imhtw * cal2joule * 1e-15 ; end where
     where (cntp == 0) ; pmhtw = spval ; elsewhere ; pmhtw = pmhtw * cal2joule * 1e-15 ; end where
     ! cgs to ppt*Sv
     where (cntg == 0) ; gmhsw = spval ; elsewhere ; gmhsw = gmhsw * 1e-12 ; end where
     where (cnta == 0) ; amhsw = spval ; elsewhere ; amhsw = amhsw * 1e-12 ; end where
     where (cnti == 0) ; imhsw = spval ; elsewhere ; imhsw = imhsw * 1e-12 ; end where
     where (cntp == 0) ; pmhsw = spval ; elsewhere ; pmhsw = pmhsw * 1e-12 ; end where
     !
     ! Add variable data to netCDF file
     status = nf_inq_varid(tranncid,trim(varnam(3)),varid)
     if (status /= nf_noerr) write(*,*) nf_strerror(status)
     status = nf_put_vara_real(tranncid,varid,itime,1,time(itime))
     if (status /= nf_noerr) write(*,*) nf_strerror(status)
     do l = 1,nv
        m = l+nd
        ostart(1) = 1 ; ocount(1) = nphi
        status = nf_inq_varid(tranncid,varnam(m),varid)
        if (vrdims(m) == 3) then
           ostart(2) = 1     ; ocount(2) =   km
           ostart(3) = itime ; ocount(3) =    1
        else
           ostart(2) = itime ; ocount(2) =    1
        endif
        select case (trim(varnam(m)))
        case ('gmsf')
           status = nf_put_vara_real(tranncid,varid,ostart,ocount,transpose(gpsiw))
           if (status /= nf_noerr) write(*,*) nf_strerror(status)
        case ('amsf')
           status = nf_put_vara_real(tranncid,varid,ostart,ocount,transpose(apsiw))
           if (status /= nf_noerr) write(*,*) nf_strerror(status)
        case ('imsf')
           status = nf_put_vara_real(tranncid,varid,ostart,ocount,transpose(ipsiw))
           if (status /= nf_noerr) write(*,*) nf_strerror(status)
        case ('pmsf')
           status = nf_put_vara_real(tranncid,varid,ostart,ocount,transpose(ppsiw))
           if (status /= nf_noerr) write(*,*) nf_strerror(status)
        case ('gmht')
           htrans = sum(gmhtw,dim=1,mask=(gmhtw /= spval))
           where (htrans == 0) ; htrans = spval ; end where
           status = nf_put_vara_real(tranncid,varid,ostart,ocount,htrans)
           if (status /= nf_noerr) write(*,*) nf_strerror(status)
         case ('amht')
           htrans = sum(amhtw,dim=1,mask=(amhtw /= spval))
           where (htrans == 0) ; htrans = spval ; end where
           status = nf_put_vara_real(tranncid,varid,ostart,ocount,htrans)
           if (status /= nf_noerr) write(*,*) nf_strerror(status)
         case ('imht')
           htrans = sum(imhtw,dim=1,mask=(imhtw /= spval))
           where (htrans == 0) ; htrans = spval ; end where
           status = nf_put_vara_real(tranncid,varid,ostart,ocount,htrans)
           if (status /= nf_noerr) write(*,*) nf_strerror(status)
         case ('pmht')
           htrans = sum(pmhtw,dim=1,mask=(pmhtw /= spval))
           where (htrans == 0) ; htrans = spval ; end where
           status = nf_put_vara_real(tranncid,varid,ostart,ocount,htrans)
           if (status /= nf_noerr) write(*,*) nf_strerror(status)
        case ('gmhs')
           htrans = sum(gmhsw,dim=1,mask=(gmhsw /= spval))
           where (htrans == 0) ; htrans = spval ; end where
           status = nf_put_vara_real(tranncid,varid,ostart,ocount,htrans)
           if (status /= nf_noerr) write(*,*) nf_strerror(status)
        case ('amhs')
           htrans = sum(amhsw,dim=1,mask=(amhsw /= spval))
           where (htrans == 0) ; htrans = spval ; end where
           status = nf_put_vara_real(tranncid,varid,ostart,ocount,htrans)
           if (status /= nf_noerr) write(*,*) nf_strerror(status)
        case ('imhs')
           htrans = sum(imhsw,dim=1,mask=(imhsw /= spval))
           where (htrans == 0) ; htrans = spval ; end where
           status = nf_put_vara_real(tranncid,varid,ostart,ocount,htrans)
           if (status /= nf_noerr) write(*,*) nf_strerror(status)
        case ('pmhs')
           htrans = sum(pmhsw,dim=1,mask=(pmhsw /= spval))
           where (htrans == 0) ; htrans = spval ; end where
           status = nf_put_vara_real(tranncid,varid,ostart,ocount,htrans)
           if (status /= nf_noerr) write(*,*) nf_strerror(status)
        case ('gmhtw')
           status = nf_put_vara_real(tranncid,varid,ostart,ocount,transpose(gmhtw))
           if (status /= nf_noerr) write(*,*) nf_strerror(status)
        case ('amhtw')
           status = nf_put_vara_real(tranncid,varid,ostart,ocount,transpose(amhtw))
           if (status /= nf_noerr) write(*,*) nf_strerror(status)
        case ('imhtw')
           status = nf_put_vara_real(tranncid,varid,ostart,ocount,transpose(imhtw))
           if (status /= nf_noerr) write(*,*) nf_strerror(status)
        case ('pmhtw')
           status = nf_put_vara_real(tranncid,varid,ostart,ocount,transpose(pmhtw))
           if (status /= nf_noerr) write(*,*) nf_strerror(status)
        case ('gmhsw')
           status = nf_put_vara_real(tranncid,varid,ostart,ocount,transpose(gmhsw))
           if (status /= nf_noerr) write(*,*) nf_strerror(status)
        case ('amhsw')
           status = nf_put_vara_real(tranncid,varid,ostart,ocount,transpose(amhsw))
           if (status /= nf_noerr) write(*,*) nf_strerror(status)
        case ('imhsw')
           status = nf_put_vara_real(tranncid,varid,ostart,ocount,transpose(imhsw))
           if (status /= nf_noerr) write(*,*) nf_strerror(status)
        case ('pmhsw')
           status = nf_put_vara_real(tranncid,varid,ostart,ocount,transpose(pmhsw))
           if (status /= nf_noerr) write(*,*) nf_strerror(status)
        case default
           write(*,*) 'Cannot happen.'
        end select
     enddo
  enddo
  status = nf_close(tranncid)
end program poptrans

