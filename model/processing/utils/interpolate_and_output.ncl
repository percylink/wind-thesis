; downloaded 8/26/14 from http://www2.mmm.ucar.edu/wrf/OnLineTutorial/Graphics/NCL/Examples/LEVELS_INTERP/wrf_Height2.ncl


;   Example script to produce plots for a WRF real-data run,
;   with the ARW coordinate dynamics option.
;   Interpolating to specified height levels
;   In this example we obtain all the fields at once and interpolate to all height levels at once

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

  run_name = "CA-0.2"
  root_dir = "$SCRATCH/WRF/output"

  do i_d = 0,1  ; domain loop

    if (i_d .eq. 0) then
      domain = "d01"
    else if (i_d .eq. 1) then
      domain = "d02"
    else if (i_d .eq. 2) then
      domain = "d03"
    end if
    end if
    end if
    
    print(domain)

    run_dir = root_dir + "/" + run_name
    system("cd "+run_dir+"; ln -sf wrfout_"+domain+"_2009-07-01_00:00:00 wrfout_"+domain+"_2009-07-01_00:00:00.nc")
    a = addfile(run_dir+"/wrfout_"+domain+"_2009-07-01_00:00:00.nc","r")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ; What times and how many time steps are in the data set?
	times = wrf_user_getvar(a,"times",-1)  ; get all times in the file
	ntimes = dimsizes(times)         ; number of times in the file

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; First get the variables we will need        

    tc = wrf_user_getvar(a,"tc",-1)        ; T in C
    u  = wrf_user_getvar(a,"ua",-1)        ; u averaged to mass points
    v  = wrf_user_getvar(a,"va",-1)        ; v averaged to mass points
    p  = wrf_user_getvar(a, "pressure",-1) ; pressure is our vertical coordinate
    z  = wrf_user_getvar(a, "z",-1)        ; grid point height
    xlat  = wrf_user_getvar(a, "XLAT",0)
    xlon  = wrf_user_getvar(a, "XLONG",0)
    hgt  = wrf_user_getvar(a, "HGT",0)

    ; The specific height levels that we want the data interpolated to.
    ; And interpolate to these levels
    height_levels = (/ 110., 150., 250., 350., 550., 750., 1050./)    ; height levels to plot - in meter
    nlevels       = dimsizes(height_levels)     ; number of height levels

    p_plane  = wrf_user_intrp3d( p,z,"h",height_levels,0.,False)
    tc_plane = wrf_user_intrp3d(tc,z,"h",height_levels,0.,False)
    u_plane  = wrf_user_intrp3d( u,z,"h",height_levels,0.,False)
    v_plane  = wrf_user_intrp3d( v,z,"h",height_levels,0.,False)

    print(dimsizes(p_plane))
    print(p_plane!0)
    print(p_plane!1)
    print(p_plane!2)
    print(p_plane!3)
    
    system("/bin/rm -f "+run_dir+"/wrfout_interp_"+domain+".nc")   ; remove any pre-existing file
    ncdf     = addfile(run_dir+"/wrfout_interp_"+domain+".nc" ,"c")  ; open output netCDF file

    ; make time and UNLIMITED dimension      ; recommended  for most applications
	filedimdef(ncdf,"time",-1,True) 

    ; output variables directly
    ncdf->T  = tc_plane                          ; 4D               
	ncdf->P = p_plane
	ncdf->U = u_plane
	ncdf->V = v_plane
	
	ncdf->levels = height_levels
	ncdf->XLAT = xlat
	ncdf->XLONG = xlon
	ncdf->HGT = hgt
	ncdf->T2 = wrf_user_getvar(a, "T2",-1)
	ncdf->TSK = wrf_user_getvar(a, "TSK",-1)
	ncdf->Times = wrf_user_getvar(a,"Times",-1)
	ncdf->DX = a@DX
	ncdf->DY = a@DY

    delete(tc)
    delete(u)
    delete(v)
    delete(p)
    delete(z)
    delete(p_plane)
    delete(tc_plane)
    delete(u_plane)
    delete(v_plane)
    delete(xlat)
    delete(xlon)
    delete(hgt)
    
  end do   ; END DOMAIN LOOP
  
end