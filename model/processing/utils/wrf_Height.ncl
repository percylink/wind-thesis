; Downloaded on 8/21/2014 from http://www2.mmm.ucar.edu/wrf/OnLineTutorial/Graphics/NCL/Examples/LEVELS_INTERP/wrf_Height.ncl


;   Example script to produce plots for a WRF real-data run,
;   with the ARW coordinate dynamics option.
;   Interpolating to specified height levels

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin
;
; The WRF ARW input file.  
; This needs to have a ".nc" appended, so just do it.
  run_name = "CA-0.25"
  root_dir = "$SCRATCH/WRF/output"
  tstep_per_hr = 2

; We generate plots, but what kind do we prefer?
; type = "png"
; type = "x11"
  type = "pdf"
; type = "ps"
; type = "ncgm"
; wks = gsn_open_wks(type,"plt_HeightLevel")


; Set some basic resources
  res = True
  res@MainTitle = "REAL-TIME WRF"
  res@Footer = False

  pltres = True
  mpres = True
  mpres@mpGeophysicalLineColor = "Black"
  mpres@mpNationalLineColor    = "Black"
  mpres@mpUSStateLineColor     = "Black"
  mpres@mpGridLineColor        = "Black"
  mpres@mpLimbLineColor        = "Black"
  mpres@mpPerimLineColor       = "Black"


  do i_d = 0,1   ; DOMAIN LOOP

    if (i_d .eq. 0) then      
      domain = "d01"
    else if (i_d .eq. 1) then
      domain = "d02"
    else if (i_d .eq. 2) then
      domain = "d03"
    end if
    end if
    end if

    run_dir = root_dir + "/" + run_name
    plot_dir = "plots/"+domain+"/map_P_T_uv_diff"
    system("cd "+run_dir+"; mkdir "+plot_dir+"; ln -sf wrfout_"+domain+"_2009-07-01_00:00:00 wrfout_"+domain+"_2009-07-01_00:00:00.nc")
    ;system("cd "+root_dir+"/"+ctrl_name+"; ln -sf wrfout_"+domain+"_2009-07-01_00:00:00 wrfout_"+domain+"_2009-07-01_00:00:00.nc")
    a = addfile(run_dir+"/wrfout_"+domain+"_2009-07-01_00:00:00.nc","r")
    ;a_ctrl = addfile(root_dir+"/"+ctrl_name+"/wrfout_"+domain+"_2009-07-01_00:00:00.nc","r")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; What times and how many time steps are in the data set?
  times = wrf_user_getvar(a,"times",-1)  ; get all times in the file
  ntimes = dimsizes(times)         ; number of times in the file

; The specific height levels that we want the data interpolated to.
  height_levels = (/ 100., 200., 300., 500., 750., 1000., 2000./)   ; height levels to plot - in meter
  nlevels       = dimsizes(height_levels)     ; number of height levels

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 do level = 0,nlevels-1                 ; LOOP OVER LEVELS

    height = height_levels(level)
    ;wks = gsn_open_wks(type, "plt_HeightLevel_" + flt2string(round(height, 3)))

    do it = 42*tstep_per_hr, 88*tstep_per_hr, 2*tstep_per_hr             ; TIME LOOP

      print("Working on level " + flt2string(round(height, 3)) + " time: " + times(it) )
      wks = gsn_open_wks(type, run_dir+"/"+plot_dir+"/plt_HeightLevel_" + flt2string(round(height, 3)) + "_" + times(it))
      res@TimeLabel = times(it)   ; Set Valid time to use on plots

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; First get the variables we will need        

      tc = wrf_user_getvar(a,"tc",it)        ; T in C
      u  = wrf_user_getvar(a,"ua",it)        ; u averaged to mass points
      v  = wrf_user_getvar(a,"va",it)        ; v averaged to mass points
      p  = wrf_user_getvar(a, "pressure",it) ; pressure is our vertical coordinate
      z  = wrf_user_getvar(a, "z",it)        ; grid point height
      rh = wrf_user_getvar(a,"rh",it)        ; relative humidity

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


      ;do level = 0,nlevels-1                 ; LOOP OVER LEVELS

        ;height = height_levels(level) 

      p_plane  = wrf_user_intrp3d( p,z,"h",height,0.,False)
      tc_plane = wrf_user_intrp3d(tc,z,"h",height,0.,False)
      rh_plane = wrf_user_intrp3d(rh,z,"h",height,0.,False)
      u_plane  = wrf_user_intrp3d( u,z,"h",height,0.,False)
      v_plane  = wrf_user_intrp3d( v,z,"h",height,0.,False)

      ;u_plane = u_plane*1.94386     ; kts
      ;v_plane = v_plane*1.94386     ; kts
      ;u_plane@units = "kts"
      ;v_plane@units = "kts"
      speed = sqrt(u_plane^2 + v_plane^2)


    ; Plotting options for T                
      opts = res                          
      opts@cnLineColor = "Red"
      opts@ContourParameters = (/ 1.0 /)
      opts@cnInfoLabelOrthogonalPosF = 0.07  ; offset second label information
      opts@gsnContourLineThicknessesScale = 2.0
      contour_tc = wrf_contour(a,wks,tc_plane,opts)
      delete(opts)


    ; Plotting options for Pressure
      opts = res                          
      opts@cnLineColor = "Blue"
      opts@ContourParameters = (/ 0.5 /)
      opts@gsnContourLineThicknessesScale = 3.0
      contour_p = wrf_contour(a,wks,p_plane,opts)
      delete(opts)


    ; Plotting options for RH                
    ;  opts = res                          
    ;  opts@cnFillOn = True  
    ;  opts@ContourParameters = (/ 10., 90., 10./)
    ;  opts@cnFillColors = (/"White","White","White", \
    ;                        "White","Chartreuse","Green",\
    ;                        "Green3","Green4", \
    ;                        "ForestGreen","PaleGreen4"/)
    ;  contour_rh = wrf_contour(a,wks,rh_plane,opts)
    ;  delete(opts)


    ; Plotting options for speed
      opts = res
      opts@cnFillOn = True
      opts@FieldTitle = "Wind speed"
      opts@ContourParameters = (/ -5., 5, 0.5 /)
      ;opts@cnFillColors = (/"White","White","White", \
      ;                      "White","Chartreuse","Green",\
      ;                      "Green3","Green4", \
      ;                      "ForestGreen","PaleGreen4"/)
      contour_s = wrf_contour(a,wks,speed,opts)
      delete(opts)


    ; Plotting options for Wind Vectors                 
      opts = res          
      opts@FieldTitle = "Wind"   ; overwrite Field Title
      opts@vcGlyphStyle = "CurlyVector"
      opts@NumVectors = 47       ; wind barb density
      vector = wrf_vector(a,wks,u_plane,v_plane,opts)
      delete(opts)


    ; MAKE PLOTS                                       
    ;  plot = wrf_map_overlays(a,wks,(/contour_rh,contour_tc,contour_p, \
    ;                            vector/),pltres,mpres)
      plot = wrf_map_overlays(a,wks,(/contour_s, contour_tc,contour_p, \
                                vector/),pltres,mpres)
      
      delete(plot)
      delete(wks)
      delete(contour_tc)
      delete(contour_p)
      delete(vector)

        delete(tc_plane)
        delete(u_plane)
        delete(v_plane)
        delete(p_plane)
        delete(rh_plane)
        delete(speed)

    end do      ; END OF LEVEL LOOP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  end do        ; END OF TIME LOOP

      delete(tc)
      delete(u)
      delete(v)
      delete(p)
      delete(z)
      delete(rh)

  end do        ; END OF DOMAIN LOOP

end
