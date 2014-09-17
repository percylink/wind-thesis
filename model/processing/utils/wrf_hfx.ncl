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
  ctrl_name = "CA-0.24"
  test_name = "CA-dryCV"
  root_dir = "$SCRATCH/WRF/output"

; We generate plots, but what kind do we prefer?
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


  do i_d = 0,2  ; domain loop

    if (i_d .eq. 0) then
      domain = "d01"
    else if (i_d .eq. 1) then
      domain = "d02"
    else if (i_d .eq. 2) then
      domain = "d03"
    end if
    end if
    end if
    run_dir_test = root_dir + "/" + test_name    
    run_dir_ctrl = root_dir + "/" + ctrl_name
    plot_dir = "plots/"+domain+"/HFX"
    
    system("cd "+run_dir_test+"; rm -r "+plot_dir+"; mkdir "+plot_dir+"; ln -sf wrfout_"+domain+"_2009-07-01_00:00:00 wrfout_"+domain+"_2009-07-01_00:00:00.nc")
    system("cd "+run_dir_ctrl+"; rm -r "+plot_dir+"; mkdir "+plot_dir+"; ln -sf wrfout_"+domain+"_2009-07-01_00:00:00 wrfout_"+domain+"_2009-07-01_00:00:00.nc")
    a_test = addfile(run_dir_test+"/wrfout_"+domain+"_2009-07-01_00:00:00.nc","r")
    a_ctrl = addfile(run_dir_ctrl+"/wrfout_"+domain+"_2009-07-01_00:00:00.nc","r")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; What times and how many time steps are in the data set?
  times = wrf_user_getvar(a_test,"times",-1)  ; get all times in the file
  ntimes = dimsizes(times)         ; number of times in the file

; The specific height levels that we want the data interpolated to.
  ;height_levels = (/ 250., 2000./)   ; height levels to plot - in meter
  ;nlevels       = dimsizes(height_levels)     ; number of height levels

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  do it = 4,ntimes-1,2             ; TIME LOOP

    print("Working on time: " + times(it) )
    res@TimeLabel = times(it)   ; Set Valid time to use on plots

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; First get the variables we will need        

	hfx_t = wrf_user_getvar(a_test,"HFX",it)
	hfx_c = wrf_user_getvar(a_ctrl,"HFX",it)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



    ; Plotting options for HFX ctrl
      wks = gsn_open_wks(type, run_dir_ctrl+"/"+plot_dir+"/HFXctrl_" + times(it))            
      opts = res                          
      opts@cnFillOn = True
      opts@ContourParameters = (/ 0., 500.0, 50. /)
      contour_ctrl = wrf_contour(a_ctrl,wks,hfx_c,opts)
      delete(opts)
      plot = wrf_map_overlays(a_ctrl,wks,(/contour_ctrl/),pltres,mpres)
	  delete(wks)
      

    ; Plotting options for HFX diff
      wks = gsn_open_wks(type, run_dir_test+"/"+plot_dir+"/HFXdiff_" + times(it))            
      opts = res                          
      opts@cnFillOn = True
      opts@ContourParameters = (/ -300., 300.0, 50. /)
      contour_test = wrf_contour(a_test,wks,hfx_t-hfx_c,opts)
      delete(opts)
      plot = wrf_map_overlays(a_test,wks,(/contour_test/),pltres,mpres)
      delete(wks)


      delete(hfx_t)
      delete(hfx_c)                                       


    end do      ; END OF TIME LOOP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  end do        ; END OF DOMAIN LOOP

end
