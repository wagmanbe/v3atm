module surface_emis_mod

! This is a collective module for codes related to the reading surface
! emissivity data and updating of surface temperature to account for
! surface spectral emissivity

  use time_manager  , only: get_curr_date
! use radiation     , only: umich_surf_emis_file
  use radconstants  , only: nlwbands 
  use ppgrid        , only: pcols, begchunk, endchunk
  use shr_kind_mod  , only: r8 => shr_kind_r8
  use shr_const_mod,  only: shr_const_stebol
  use spmd_utils,     only: masterproc
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun
  use perf_mod,       only: t_startf, t_stopf
!
! !PUBLIC TYPES:
  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

! public :: surface_emis_init
  public :: surface_emis_intr

!--------------------------------------------------------------------------
! Private interfaces
!--------------------------------------------------------------------------

  integer, save :: month_last_read = -999
  real(r8), save, allocatable:: emis0(:,:)
  real(r8) :: v1_rrtmg_lw(nlwbands + 1)   ! RRTMG_LW band edges
  real(r8) :: desert_emis(nlwbands), water_emis(nlwbands),ice_emis(nlwbands), &
              grass_emis(nlwbands), snow_emis(nlwbands)

!
!================================================================================
CONTAINS
!================================================================================

  subroutine surface_emis_intr (cam_in, cam_out)

    use camsrfexch , only: cam_out_t, cam_in_t
    use phys_grid  , only: get_ncols_p, get_rlat_all_p, get_rlon_all_p
    use shr_const_mod,  only: shr_const_stebol
    use radconstants,   only: nlwbands, wavenumber1_longwave
    type(cam_in_t) , intent(INOUT) :: cam_in(begchunk:endchunk)
    type(cam_out_t), intent(IN) :: cam_out(begchunk:endchunk)

    ! This change is to define and initialize new variables
    integer :: i,j,c,ncols,sizebuf
    real(r8) :: ts_lw              ! surface temperature derived from longwave upward flux
    real(r8) :: lats(pcols)           ! array of chunk latitudes
    real(r8) :: lons(pcols)           ! array of chunk longitude
    real(r8) :: pi
    real(r8) :: lats2(pcols*(endchunk-begchunk+1)),lons2(pcols*(endchunk-begchunk+1))

    integer :: yr              ! CAM current year
    integer :: mon             ! CAM current month
    integer :: day             ! CAM current day
    integer :: tod             ! CAM current time of day (sec)



    v1_rrtmg_lw(1:nlwbands) = wavenumber1_longwave
    v1_rrtmg_lw(nlwbands+1) = 3250._r8

! This change is to compute surface skin temperature and pass it to cam_in
    call get_curr_date(yr, mon, day, tod)
    pi = 4._r8*atan(1._r8)
    sizebuf=0
    do c=begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, ncols, lats)
       call get_rlon_all_p(c, ncols, lons)
       do i=1,ncols
          sizebuf = sizebuf+1
          lats2(sizebuf) = lats(i)*180._r8/pi
          if (lons(i) .lt. 0._r8) then
             lons2(sizebuf) = lons(i)*180._r8/pi + 180._r8
          else
             lons2(sizebuf) = lons(i)*180._r8/pi
          endif
       enddo
    enddo

    !*** using realstic emissivity ***       
    if (cam_out(begchunk)%do_emis(1) .eq. 1) then
       if (month_last_read .le. 0) then
          ! never read
          allocate(emis0(pcols*(endchunk-begchunk+1),nlwbands))
       end if
       if (mon .ne. month_last_read) then
          call t_startf('SURF_EMIS_READ')
          call read_surface_emis(sizebuf,lats2,lons2,mon,emis0(1:sizebuf,:) )
                              
          call t_stopf('SURF_EMIS_READ')
           month_last_read = mon
           if(masterproc) then
             write(iulog,*)'Read surface emissivity data for month ',mon
           endif
       endif

       sizebuf=0
       do c=begchunk, endchunk
          ncols = get_ncols_p(c)

          do i=1,ncols
             sizebuf = sizebuf + 1

             if (cam_in(c)%landfrac(i).gt. 0.99_r8 .and. cam_in(c)%icefrac(i) .lt. 0.01_r8 .and. &
                     (cam_in(c)%snowhland(i) + cam_in(c)%snowhice(i)).lt. 0.001_r8) then
                ! if original emissivity over band 1080-1180 cm-1, ie.
                ! cam_in%srf_emis_spec(i,8) is
                ! smaller than 0.8, then the original surface type
                ! of this grid is desert, otherwise is non-desert
                 cam_in(c)%srf_emis_spec(i,:)= emis0(sizebuf,:)
  
                ! if the original surface type is non-desert but LAI is
                ! smaller than 0.001, change to desert type
                 if (cam_in(c)%tlai(i).lt. 0.001_r8 .and. emis0(sizebuf,8)>0.8_r8) then
                   cam_in(c)%srf_emis_spec(i,:)  = desert_emis
                 endif
                ! if the orignal surface type is desert but LAI larger than
                ! 2, change to grass type
                 if (cam_in(c)%tlai(i).gt. 2._r8 .and. emis0(sizebuf,8)<0.8_r8) then
                   cam_in(c)%srf_emis_spec(i,:)  = grass_emis
                 endif
             else
                 cam_in(c)%srf_emis_spec(i,:)  =  emis0(sizebuf,:) * cam_in(c)%landfrac(i) + &
                         ice_emis * cam_in(c)%icefrac(i) + water_emis * cam_in(c)%ocnfrac(i)
             endif  ! end if of cam_in(c)%landfrac, %icefrac, and %snowhland 
  
             call  get_Ts_from_LW_emis(v1_rrtmg_lw, cam_in(c)%srf_emis_spec(i,:), &
                     cam_in(c)%lwup(i),cam_out(c)%flwds_spec(i,:), nlwbands + 1, 3, ts_lw)
             cam_in(c)%ts_atm(i) = ts_lw
             !write(iulog,*) 'chen2',sum(cam_out(c)%flwds_spec(i,:)), real(cam_in(c)%srf_emis_spec(i,1)),&
             !               ts_lw-sqrt(sqrt((cam_in(c)%lwup(i)/shr_const_stebol)))
             cam_in(c)%ts(i) =  cam_in(c)%ts_atm(i)
  
           enddo ! end do of i, ncols
        enddo ! end do of c, chunk

    endif

  end subroutine surface_emis_intr

!-------------------------------------------
! U-MICH team added the following subroutines/functions -->  

  subroutine read_surface_emis(ncols,lats2,lons2,mn,surface_emis)
  !  This subroutine is added by U-MICH team on Dec.15, 2019
  !  This subroutine is to read surface emissivity from dataset

      use netcdf
  
      use time_manager, only: get_curr_date
      use ppgrid

      use error_messages, only : handle_ncerr
      use radconstants  , only : nlwbands  ! added by U-MICH team on Dec.15, 2019
      use cam_logfile   , only : iulog

      implicit none
      integer :: ncid, status, latid,lonid,bandid,timeid
      character(256) filename

      integer :: ntime, nlat, nlon, nband,i,mn
      real, allocatable :: band_emissivity(:)
      real :: water_emis_tmp(nlwbands), ice_emis_tmp(nlwbands), &
              desert_emis_tmp(nlwbands), grass_emis_tmp(nlwbands)
      character(len = nf90_max_name) :: RecordDimName
      integer :: lat_varID,emis_varID, lon_varID, waterID, iceID, desertID, grassID
      real, allocatable:: lat_tmp(:), lon_tmp(:)
      real(r8), allocatable:: lat(:), lon(:)
      !integer ::start(4),count(4)
      integer ::start(4),cnt(4)
      integer :: ncols,j
      integer ::ilats, ilons
      real(r8) :: lats2(ncols)           ! array of chunk latitudes
      real(r8) :: lons2(ncols)           ! array of chunk longitude
      real(r8) :: minvalue
      real(r8),intent(out) :: surface_emis(ncols, nlwbands)

      filename = "surface_emissivity_1x1_UMRad_53deg.nc"
      !filename = "/global/cscratch1/sd/xianwen/data/emis/surface_emissivity_1x1_RRTMGP_53deg.nc"
      status = nf90_open(trim(filename), nf90_nowrite, ncid)
      status = nf90_inq_dimid(ncid, "time", timeID)
      status = nf90_inq_dimid(ncid, "lat", latID)
      status = nf90_inq_dimid(ncid, "lon", lonID)
      status = nf90_inq_dimid(ncid, "band", bandID)
      
      if (month_last_read .le. 0) then
         status = nf90_inq_dimid(ncid, "water_emissivity", waterID)
         status = nf90_inq_dimid(ncid, "ice_emissivity", iceID)
         status = nf90_inq_dimid(ncid, "desert_emissivity", desertID)
         status = nf90_inq_dimid(ncid, "grass_emissivity", grassID)
      endif

      status =  nf90_inquire_dimension( ncid, timeID,len=ntime )
      status =  nf90_inquire_dimension( ncid, latID,len=nlat )
      status =  nf90_inquire_dimension( ncid, lonID,len=nlon )
      status =  nf90_inquire_dimension( ncid, bandID,len=nband )
      ! nband needs to be get only once. xianwen. 
      !status =  nf90_inquire_dimension( ncid, waterID,len=nband )
      !status =  nf90_inquire_dimension( ncid, iceID,len=nband )
      !status =  nf90_inquire_dimension( ncid, desertID,len=nband )
      !status =  nf90_inquire_dimension( ncid, grassID,len=nband )

      allocate(band_emissivity(nband))
      allocate(lat_tmp(nlat))
      allocate(lon_tmp(nlon))
      allocate(lat(nlat))
      allocate(lon(nlon))
      
       
      status = nf90_inq_varid (ncid, 'lat', lat_varID )
      !status = nf90_get_var (ncid, lat_varID, lat)
      status = nf90_get_var (ncid, lat_varID, lat_tmp)
      lat(:) = real(lat_tmp(:), r8)

      status = nf90_inq_varid (ncid, 'lon', lon_varID )
      !status = nf90_get_var (ncid, lon_varID, lon)
      status = nf90_get_var (ncid, lon_varID, lon_tmp)
      lon(:) = real(lon_tmp(:), r8)


      if (month_last_read .le. 0) then
         status = nf90_inq_varid (ncid, 'water_emissivity', waterID )
         status = nf90_get_var (ncid, waterID, water_emis_tmp)
         water_emis(:) = real(water_emis_tmp(:), r8)

         status = nf90_inq_varid (ncid, 'ice_emissivity', iceID )
         status = nf90_get_var (ncid, iceID, ice_emis_tmp)
         ice_emis(:) = real(ice_emis_tmp(:), r8)

         status = nf90_inq_varid (ncid, 'desert_emissivity', desertID )
         status = nf90_get_var (ncid, desertID, desert_emis_tmp)
         desert_emis(:) = real(desert_emis_tmp(:), r8)

         status = nf90_inq_varid (ncid, 'grass_emissivity', grassID )
         status = nf90_get_var (ncid, grassID, grass_emis_tmp)
         grass_emis(:) = real(desert_emis_tmp(:), r8)
      endif

      cnt =(/nlwbands,1,1,1/)
      do i = 1, ncols
         minvalue= 10000.0_r8
         do j=1, nlat
           if (abs(lat(j) - lats2(i)) .le. minvalue)  then
             ilats = j
             minvalue = abs(lat(j) - lats2(i))
           endif
         enddo
         minvalue= 10000.0_r8
         do j=1, nlon
           if (abs(lon(j) - lons2(i)) .le. minvalue)  then
             ilons = j
             minvalue = abs(lon(j) - lons2(i))
           endif
         enddo
         start =(/1,ilons,ilats,mn/)
         status = nf90_inq_varid (ncid, 'band_emissivity', emis_varID )
         status = nf90_get_var (ncid, emis_varID, band_emissivity,start = start,count = cnt)
         surface_emis(i,:) = real(band_emissivity(:),r8)
      enddo
       
      status = NF90_CLOSE( ncid )

  end subroutine read_surface_emis

  subroutine get_Ts_from_LW_emis(v1, emis, LW, LWdown, v1_num, nguass_point, Ts)
  !  This subroutine is made by U-MICH team on Dec.15, 2019
  !  This subroutine is to obtain the Ts that can give the right upward LW
  ! flux with given band-averaged surface emissivity
  ! Input variables:
  ! v1, band ranges.e.g. [v11, v12, v13] will indicates 2bands, with
  ! band1 from [v11, v12] and band2 from [v12, v13];
  ! emis, the band-average surface emissivity, its size should be len(v1)-1
  ! LW, the upward LW flux in Wm^{-2}

   use cam_logfile,     only: iulog

   !integer   i, j, Count, v1_num, nguass_point
   integer   i, j, cnt, v1_num, nguass_point
   real(r8)  v1(v1_num)
   real(r8)  emis(v1_num-1)
   real(r8)  LW, LW2, LW3
   real(r8)  LWdown(v1_num-1)
   real(r8)  T1, T2, T3, F1, F2, F3, A1, A2, A3
   real(r8)  x(nguass_point), w(nguass_point)
   real(r8)  rad1, rad2, rad3, pi
   real(r8)  Ts
   
   pi = 4._r8*atan(1._r8)

! T1 and T2 should be two values encompassing all possible Ts
   T1 = 150._r8;
   T2 = 400._r8;

   F1 = 0._r8;
   F2 = 0._r8;
   F3 = 0._r8;

   LW3 = LW
   do i = 1, v1_num - 1
        if (emis(i) .eq. 0._r8) then
            emis(i) = 1.0_r8
        endif
        LW3 = LW3 - (1._r8 - emis(i)) * LWdown(i)
   enddo

   LW2 = int(LW3*10._r8)/10.0_r8

   !Count = 0
   cnt = 0
   do while (abs(F3-LW2)>0.001_r8)

!       for count how many iteration needed

        !Count = Count + 1
        cnt = cnt + 1
        T3 = (T1 + T2)/2._r8

        do i = 1,v1_num - 1

              call gaulegf(v1(i), v1(i+1), nguass_point, x, w)

               A1 = 0._r8
               A2 = 0._r8
               A3 = 0._r8

               do j = 1,nguass_point

                  A1 = A1 + planck(x(j), T1) * w(j)
                  A2 = A2 + planck(x(j), T2) * w(j)
                  A3 = A3 + planck(x(j), T3) * w(j)

                enddo

!              don't mix up i and j subscripts

               F1 = F1 + emis(i) * A1;
               F2 = F2 + emis(i) * A2;
               F3 = F3 + emis(i) * A3;
        enddo

! covert to Wm-2
        F1 = F1 * pi * 1.e-3_r8
        F2 = F2 * pi * 1.e-3_r8
        F3 = F3 * pi * 1.e-3_r8


        !if (Count .eq.1 .and. (LW2 .lt. F1 .or. LW2 .gt. F2)) then
        if (cnt .eq.1 .and. (LW2 .lt. F1 .or. LW2 .gt. F2)) then
           !write(iulog,*) 'CCC',LW2,LW, LWdown(1),emis(1),LWdown(6),emis(6),Count
           !write(iulog, *)'the LW2 is too low or too high so it is beyond the reasonable range'
           !return
           write(iulog,*) 'CCC_1, LW=',LW,', LW2=',LW2,', LW3=',LW3
           write(iulog,*) 'CCC_2, LWdown: ', LWdown
           write(iulog,*) 'CCC_3, emis: ', emis
           write(iulog, *)'the LW2 is too low or too high so it is beyond the reasonable range. This is likely due to incorrect emissivity values'
           call endrun('STOP: subroutin get_Ts_from_LW_emis: LW2 is too low or too high')
        elseif (LW2 .lt. F1) then !the case due to numerical error
           T1 = T1 - 2._r8
        elseif (LW2 .gt. F2) then ! the case due to numberical error
           T2 = T2 + 2._r8
        elseif (LW2 .ge. F3) then
           T1 = T3
        elseif (LW2 .lt. F3) then
           T2 = T3
        endif
   enddo

!   Ts = T3
!   To correct the bias

    if (T3.lt.270._r8) then
         Ts = T3 + (T3-160._r8)*0.0008_r8 + 0.13_r8
    elseif (T3.lt.300._r8) then
         Ts = T3 + (T3-270._r8)*0.0006_r8 + 0.2171_r8
    elseif (T3.lt.320._r8) then
         Ts = T3 + (T3-300._r8)*0.0002_r8 + 0.2342_r8
    elseif (T3.lt.340._r8) then
         Ts = T3 - (T3-320._r8)*0.0004_r8 + 0.2386_r8
    else
         Ts = T3 - (T3-340._r8)*0.0013_r8+0.2316_r8
    endif

  end subroutine get_Ts_from_LW_emis


  subroutine gaulegf(x1, x2, n, x, w)
  !  This subroutine is added by U-MICH team on Dec.15, 2019
  !  This subroutine is for Gauss Legendre integration
 
  ! gauleg.f90     P145 Numerical Recipes in Fortran
  ! compute x(i) and w(i)  i=1,n  Legendre ordinates and weights
  ! on interval -1.0 to 1.0 (length is 2.0)
  ! use ordinates and weights for Gauss Legendre integration
  implicit none
  integer, intent(in) :: n
  real(r8), intent(in) :: x1, x2
  real(r8), dimension(n), intent(out) :: x, w
  integer :: i, j, m
  real(r8) :: p1, p2, p3, pp, xl, xm, z, z1
  real(r8), parameter :: eps=3.e-14_r8

  m = (n+1)/2
  xm = 0.5_r8*(x2+x1)
  xl = 0.5_r8*(x2-x1)

  do i=1,m
    !z = cos(3.141592654d0*(i-0.25d0)/(n+0.5d0))
    z = cos(3.141592654_r8*(real(i,r8)-0.25_r8)/(real(n,r8)+0.5_r8))
    z1 = 0.0_r8
    do while(abs(z-z1) .gt. eps)
      p1 = 1.0_r8
      p2 = 0.0_r8
      do j=1,n
        p3 = p2
        p2 = p1
        !p1 = ((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
        p1 = ((2.0_r8*real(j,r8)-1.0_r8)*z*p2-(real(j,r8)-1.0_r8)*p3)/real(j,r8)
      end do
      pp = real(n,r8)*(z*p1-p2)/(z*z-1.0_r8)
      z1 = z
      z = z1 - p1/pp
    end do
    x(i) = xm - xl*z
    x(n+1-i) = xm + xl*z
    w(i) = (2.0_r8*xl)/((1.0_r8-z*z)*pp*pp)
     w(n+1-i) = w(i)
  end do
 end subroutine gaulegf

 function  planck(freq,temp)
 !  This function is added by U-MICH team on Dec.15, 2019
 !  This function is to compute blackbody thermal emission based on a temperature
 ! freq in wavenumber, temp in Kelvin degree and
 ! radiance in 1e-3 W per square meter per sr per wavenumber

      !real(8):: freq
      real(r8):: freq
      real(r8) ::  temp, ca, cb, cof, arg, zeroind
      real(r8) ::  planck

!      ca    = 3.741832e-05 / 3.14159
      ca = 1.191043934e-05_r8

!      cb    = 1.438837896
      cb = 1.438769911_r8

      cof   = ca * (freq **3)

      arg   = cb * freq

      planck    = cof / ( exp( (arg/ temp)) - 1.0_r8 )

 end function planck
! <-- end of U-MICH add. 

end module surface_emis_mod
