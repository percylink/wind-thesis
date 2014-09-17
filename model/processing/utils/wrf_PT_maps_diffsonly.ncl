; Plot maps of T, P at a few different levels, and diffs, for 

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

; ctrl run
  ctrl_name = "CA-CV0.05w"
  tstep_per_hr = 2

; We generate plots, but what kind do we prefer?
; type = "png"
; type = "x11"
  type = "pdf"
; type = "ps"
; type = "ncgm"


; Set some basic resources
  res = True
  res@MainTitle = "REAL-TIME WRF"
  ;res@Footer = False

  pltres = True
  mpres = True
  mpres@mpGeophysicalLineColor = "Black"
  mpres@mpNationalLineColor    = "Black"
  mpres@mpUSStateLineColor     = "Black"
  mpres@mpGridLineColor        = "Black"
  mpres@mpLimbLineColor        = "Black"
  mpres@mpPerimLineColor       = "Black"


do i_run = 1, 4  ; RUN LOOP
  
  if (i_run .eq. 1) then
    run_name = "CA-CV0.05a"
  end if
  if (i_run .eq. 2) then
    run_name = "CA-CV0.05d"
  end if
  if (i_run .eq. 3) then
    run_name = "CA-CV0.05s"
  end if
  if (i_run .eq. 4) then
    run_name = "CA-CV0.05n"
  end if
  
  print("run "+run_name)
  run_dir = "$SCRATCH/WRF/output/"+run_name

  do i_d = 0, 1  ; DOMAIN LOOP
  
    if (i_d .eq. 0) then
      domain = "d01"
	end if
	if (i_d .eq. 1) then
      domain = "d02"
	end if
	if (i_d .eq. 2) then
      domain = "d03"
	end if
	
	print("domain "+domain)
	plot_dir = "plots/"+domain+"/map_P_T"
    system("cd "+run_dir+"; mkdir plots; mkdir plots/"+domain+"; mkdir "+plot_dir+"; ln -sf wrfout_"+domain+"_2009-07-01_00:00:00 wrfout_"+domain+"_2009-07-01_00:00:00.nc")
    a = addfile(run_dir+"/wrfout_"+domain+"_2009-07-01_00:00:00.nc","r")
    a_c = addfile("$SCRATCH/WRF/output/"+ctrl_name+"/wrfout_"+domain+"_2009-07-01_00:00:00.nc","r")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; What times and how many time steps are in the data set?
    times = wrf_user_getvar(a,"times",-1)  ; get all times in the file
    ntimes = dimsizes(times)         ; number of times in the file

  ; The specific height levels that we want the data interpolated to.
    height_levels = (/ 250., 350., 550., 1050./)   ; height levels to plot - in meter
    nlevels       = dimsizes(height_levels)     ; number of height levels

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    do it = 42*tstep_per_hr, 64*tstep_per_hr, 2*tstep_per_hr  ; TIME LOOP - only 7-2 18:00 to 7-3 16:00

      print("Working on time: " + times(it) )
      ;wks = gsn_open_wks(type, run_dir+"/"+plot_dir+"/plt_HeightLevel_" + flt2string(round(height, 3)) + "_" + times(it))
      res@TimeLabel = times(it)   ; Set Valid time to use on plots
      
      ; get surface level p and T
      ps = wrf_user_getvar(a, "PSFC", it) / 100.  ; surface pressure, converted to hPa
      p  = wrf_user_getvar(a, "pressure",it)      ; pressure
      Ts = wrf_user_getvar(a, "TSK", it)-273.     ; surface skin temperature, converted to C
      T2 = wrf_user_getvar(a, "T2", it)-273.      ; 2m temperature, converted to C
      tc = wrf_user_getvar(a,"tc",it)             ; T in C
      z  = wrf_user_getvar(a, "z",it)             ; grid point height
      
	  ps_c = wrf_user_getvar(a_c, "PSFC", it) / 100.  ; surface pressure, converted to hPa
	  p_c  = wrf_user_getvar(a_c, "pressure",it)      ; pressure
	  Ts_c = wrf_user_getvar(a_c, "TSK", it)-273.     ; surface skin temperature, converted to C
	  T2_c = wrf_user_getvar(a_c, "T2", it)-273.      ; 2m temperature, converted to C
	  tc_c = wrf_user_getvar(a_c, "tc", it)           ; T in C
	  z_c  = wrf_user_getvar(a_c, "z",it)             ; grid point height

	  ; Plot surface pressure diff
	  wks = gsn_open_wks(type, run_dir+"/"+plot_dir+"/Ps_diff_" + times(it))
	  opts = res
	  opts@cnFillOn = True
	  opts@FieldTitle = "Surface pressure diff (test-control)"
	  opts@ContourParameters = (/ -1.2, 0.2, 0.1/)
	  contour = wrf_contour(a,wks,ps-ps_c,opts)
	  plot = wrf_map_overlays(a,wks,(/contour/),pltres,mpres)
	  delete(opts)
	  delete(wks)
	  delete(contour)

	  ; Plot surface temperature diff
	  wks = gsn_open_wks(type, run_dir+"/"+plot_dir+"/Tsk_diff_" + times(it))
	  opts = res
	  opts@cnFillOn = True
	  opts@FieldTitle = "Surface temperature diff (test-control)"
	  opts@ContourParameters = (/ -6., 14., 1./)
	  contour = wrf_contour(a,wks,Ts-Ts_c,opts)
	  plot = wrf_map_overlays(a,wks,(/contour/),pltres,mpres)
	  delete(opts)
	  delete(wks)
	  delete(contour)

	  ; Plot 2 m T diff
	  wks = gsn_open_wks(type, run_dir+"/"+plot_dir+"/T2m_diff_" + times(it))
	  opts = res
	  opts@cnFillOn = True
	  opts@FieldTitle = "2 m temperature diff (test-control)"
	  opts@ContourParameters = (/ -2., 6., 0.5/)
	  contour = wrf_contour(a,wks,T2-T2_c,opts)
	  plot = wrf_map_overlays(a,wks,(/contour/),pltres,mpres)
	  delete(opts)
	  delete(wks)
	  delete(contour)
		  
	  delete(ps_c)
	  delete(Ts_c)
	  delete(T2_c)
      delete(ps)
      delete(Ts)
      delete(T2)


 	  do level = 0,nlevels-1                 ; LOOP OVER LEVELS

        height = height_levels(level)

        p_plane  = wrf_user_intrp3d( p,z,"h",height,0.,False)
        tc_plane = wrf_user_intrp3d(tc,z,"h",height,0.,False)

        p_plane_c  = wrf_user_intrp3d( p_c,z_c,"h",height,0.,False)
        tc_plane_c = wrf_user_intrp3d(tc_c,z_c,"h",height,0.,False)
        
        if (height .eq. 250.) then
          Tcontours = (/ -3., 5.5, 0.5 /)
          Pcontours = (/ -1.5, 0., 0.1 /)
        end if
        if (height .eq. 350.) then
          Tcontours = (/ -2.5, 5.5, 0.5 /)
          Pcontours = (/ -1.4, 0., 0.1 /)
        end if
        if (height .eq. 550.) then
          Tcontours = (/ -2., 5., 0.5 /)
          Pcontours = (/ -1.2, 0., 0.1 /)
        end if
        if (height .eq. 1050.) then
          Tcontours = (/ -1., 5., 0.5 /)
          Pcontours = (/ -1.2, 0., 0.1 /)
        end if


		; Plot T diff
		wks = gsn_open_wks(type, run_dir+"/"+plot_dir+"/Tair_diff_" + flt2string(round(height, 3)) + "m_" + times(it))
		opts = res
		opts@cnFillOn = True
		opts@FieldTitle = "Air temperature diff (test-control)"
		opts@ContourParameters = Tcontours
	    contour = wrf_contour(a,wks,tc_plane-tc_plane_c,opts)
		plot = wrf_map_overlays(a,wks,(/contour/),pltres,mpres)
		delete(opts)
		delete(wks)
		delete(contour)
		  
		; Plot P diff
		wks = gsn_open_wks(type, run_dir+"/"+plot_dir+"/Pair_diff_" + flt2string(round(height, 3)) + "m_" + times(it))
		opts = res
		opts@cnFillOn = True
		opts@FieldTitle = "Air pressure diff (test-control)"
		opts@ContourParameters = Pcontours
	    contour = wrf_contour(a,wks,p_plane-p_plane_c,opts)
		plot = wrf_map_overlays(a,wks,(/contour/),pltres,mpres)
		delete(opts)
		delete(wks)
		delete(contour)
		  
	    delete(p_plane_c)
	    delete(tc_plane_c)
  	    delete(p_plane)
	    delete(tc_plane)
		
	  end do  ; end loop over levels
	  
	  delete(p)
	  delete(tc)
	  delete(z)
	  delete(p_c)
	  delete(tc_c)
	  delete(z_c)
	  
    end do      ; END OF TIME LOOP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  end do        ; END OF DOMAIN LOOP
  
  end do        ; END OF RUN LOOP

end
