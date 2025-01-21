; NCL script
; convert_CLMhist_to_PLUMBER_v3.ncl ; Keith Oleson, May 2021
; This script converts a CLM history file to a PLUMBER history file with ALMA history fields
; Has "standard_name" and "y,x" variable dimensions
; Use this for Phase II
;**************************************

;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  print ("=========================================")
  print ("Start Time: "+systemfunc("date") )
  print ("=========================================")

; ==========================================
; SP
; ==========================================
; FSA                            -> SWnet        net shortwave radiation (W/m2)
; -FIRA                          -> LWnet        net longwave radiation (W/m2)
; FSDS                           -> SWdown       surface downwelling shortwave flux (W/m2)
; FLDS                           -> LWdown       surface downwelling longwave flux (W/m2)
; FSR                            -> SWup         surface upwelling shortwave flux (W/m2)
; FIRE                           -> LWup         surface upwelling longwave flux (W/m2)
; TBOT                           -> Tair         air temperature at forcing height (K)
; QBOT                           -> Qair         specific humidity at forcing height (kg/kg)
; PBOT                           -> PSurf        air pressure at forcing height (Pa)
; WIND                           -> Wind         wind speed at forcing height (m/s)
; EFLX_LH_TOT                    -> Qle          surface upward latent heat flux (W/m2)
; FSH                            -> Qh           surface upward sensible heat flux (W/m2)
; FGR+URBAN_HEAT-URBAN_AC      	 -> Qstor        net storage heat flux in all materials (W/m2)
; URBAN_HEAT+WASTEHEAT           -> Qanth        anthropogenic heat flux
; URBAN_HEAT+WASTEHEAT           -> Qanth_Qh     anthropogenic sensible heat flux
; 0                              -> Qanth_Qle    anthropogenic latent heat flux
; TSKIN                          -> AvgSurfT     surface temperature (K)
; FSR/FSDS                       -> Albedo       surface albedo (-)
; FSRSNOW/FSDS                   -> SAlbedo      snow albedo (-)
; FSNO                   	 -> SnowFrac     snow covered fraction (-)
; RAIN                           -> Rainf        rainfall rate (kg/m2/s)
; SNOW                           -> Snowf        snowfall rate (kg/m2/s)
; QVEGT+QVEGE+QSOIL              -> Evap         total evapotranspiration (kg/m2/s)
; QOVER+QH2OSFC          	 -> Qs           surface runoff (kg/m2/s)
; QDRAI+QDRAI_PERCH      	 -> Qsb          subsurface runoff (kg/m2/s)
; QSNOMELT               	 -> Qsm          snow melt (kg/m2/s)
; -TAUX                          -> Qtau         momentum flux (N/m2) - can use just TAUX when U=V
; H2OSNO                 	 -> SWE          snow water equivalent (kg/m2)
; Change in H2OSNO               -> DelSWE       change in snow water equivalent (kg/m2)
; TG(index 1 of 1d file)         -> RoofSurfT    roof surface temperature (K)
; TG(avg index 2&3 of 1d file    -> WallSurfT    wall surface temperature (K)
; TG(index 4 of 1d file)         -> RoadSurfT    road surface temperature (K)
; TSA                            -> TairCanyon   air temperature in street canyon (K)
; TBUILD                         -> TairBuilding air temperature in buildings (K)
; SNO_T(0)                       -> SnowT        snow surface temperature (K)
; TSOI(0)(index 5 of 1d file)    -> BaresoilT    temperature of bare soil (K)
; SOILLIQ+SOILICE(index 5 of 1d) -> SoilMoist    average layer soil moisture (kg/m2)
; TSOI(index 5 of 1d file)       -> SoilTemp     average layer soil temperature (K)
; Change in SOILLIQ+SOILICE      -> DelSoilMoist change in soil moisture (kg/m2)
; Change in roof/imprvrd water   -> DelIntercept change in interception storage (kg/m2)
;  and H2OCAN and H2OSFC if v4
; ==========================================

; Sites
  Psites    = (/"AU-Preston","AU-SurreyHills","CA-Sunset","FI-Kumpula","FI-Torni","FR-Capitole","GR-HECKOR","JP-Yoyogi","KR-Jungnang","KR-Ochang","MX-Escandon","NL-Amsterdam","PL-Lipowa","PL-Narutowicza","SG-TelokKurau","UK-KingsCollege","UK-Swindon","US-Baltimore","US-Minneapolis1","US-Minneapolis2","US-WestPhoenix"/)
  Msites    = (/"AU-Pre","AU-Sur","CA-Sun","FI-Kum","FI-Tor","FR-Cap","GR-HEC","JP-Yoy","KR-Jun","KR-Och","MX-Esc","NL-Ams","PL-Lip","PL-Nar","SG-Tel","UK-Kin","UK-Swi","US-Bal","US-Mi1","US-Mi2","US-Wes"/)

; Starting year
  startyear = (/"1993","1994","2002","2000","2000","1994","2009","2006","2007","2005","2001","2009","1998","1998","2005","2002","2001","1992","1996","1996","2001"/)

; ==========================================
; Specify version and output file

; version = "v1"  ; This corresponds to v1 of Phase 1 (AU-Preston only; V1_SUBMITTED_AUG_14_2021_NO_SWup_LWup_ALBEQ1)
; version = "v2"  ; This corresponds to v2 of Phase 1 (AU-Preston only: V2_SUBMITTED_FEB_16_2021)
  version = "v3"  ; This is the first version of Phase 2 submission (all sites; V3_SUBMITTED_JUN_17_2021)
                  ; Also used for V3 with Bowen Fang v1 and v2 emissivity (all sites; V3_EMISSIVITY_JUN_29_2023 and V3_EMISSIVITY_JUL_10_2023)
; version = "v4"  ; This is Phase 2 but with urb/veg lu mix (all sites; V4_SUBMITTED_JUL_07_2021)

; Write netcdf output file?
  do_netcdf = True
; Write csv output file?
  do_csv = True

; Compute landunit soil moisture instead of pervious road soil moisture
  lu_soilmoist = False  ; This appears to work for AU-Preston for v3, but doesn't work as well for v4.
; Add DelIntercept, DelSoilMoist, and DelSWE variables and make changes associated with water balance project
  add_WaterBalVars = True
; ==========================================

; Get experiment type
  expr_dir  = getenv("expr_dir")
  expr      = getenv("expr")

; Get site
  sel = stringtointeger(getenv("sel"))

  f_out_nc  = "CLMU5_"+Psites(sel)+"_"+expr+"_"+version+".nc" 
  f_out_csv = "CLMU5_"+Psites(sel)+"_"+expr+"_"+version+".csv"

; Specify input directory 
  f_in_dir  = "/glade/work/oleson/PLUMBER/PLUMBER/"+Psites(sel)+"/"+expr_dir+"/Figure5/"
  system("mkdir -p " + f_in_dir)

  f_out_dir = "/glade/work/oleson/PLUMBER/PLUMBER/"+Psites(sel)+"/"+expr_dir+"/Figure5/"

  f_in        = f_in_dir+"CLMU5_"+Msites(sel)+"_"+expr+".clm2.h2."+startyear(sel)+"-01-01-00000.nc"
  f_stupid    = f_in_dir+"CLMU5_"+Msites(sel)+"_"+expr+".stupid.clm2.h2."+startyear(sel)+"-01-01-00000.nc"
  f_in_1d     = f_in_dir+"CLMU5_"+Msites(sel)+"_"+expr+".clm2.h1."+startyear(sel)+"-01-01-00000.nc"
  f_stupid_1d = f_in_dir+"CLMU5_"+Msites(sel)+"_"+expr+".stupid.clm2.h1."+startyear(sel)+"-01-01-00000.nc"

  command = "ncks -O -d time,1, "
  command = command + f_in + " " + f_stupid 
  print((/command/))
  system(command)
  command = "ncks -O -d time,1, "
  command = command + f_in_1d + " " + f_stupid_1d 
  print((/command/))
  system(command)

  fls  = f_stupid

  fls_1d = f_stupid_1d

; Read in some time-invariant data
  fstfile  = addfile(fls,"r")
  levgrnd  = fstfile->levgrnd
  levdcmp  = fstfile->levdcmp
  hist_interval = 2
  string_length = 16
  nlevsoi       = 20
  nlevsno       = 12
  soil_layer    = ispan(1,20,1)
  soil_layer!0  = "soil_layer"
  soil_layer@long_name = "Soil layer number"
  print(soil_layer)
  y             = ispan(1,1,1)
  ny            = dimsizes(y)
  y!0           = "y"
  y&y           = ny
  y@long_name   = "y dimension"
  x             = ispan(1,1,1)
  nx            = dimsizes(x)
  x!0           = "x"
  x&x           = nx
  x@long_name   = "x dimension"

  data          = addfile(fls,"r")
; Read in some time-varying data
  time          = data->time
  ntim          = dimsizes(time)
; Convert to seconds since
; time          = ((/time/) * 24. * 60. * 60.) - 1800.
  time          = ((/time/) * 24. * 60. * 60.) - 1800.
  time@units    = "seconds since "+startyear(sel)+"-01-01 00:00:00"
  time@standard_name = "time"
  delete(time@bounds)
  time&time = (/time/)
  printVarSummary(time)
  seconds_since = time@units

; NOTE that all variables that have dimension lndgrid must be converted to x,y
  var2d    = new((/ny,nx/),float)
  var3d    = new((/ntim,ny,nx/),float)
  var4d    = new((/ntim,nlevsoi,ny,nx/),float)

  lat      = conform(var2d,fstfile->lat,(/0/))
  ; Coordinates now need to be more precise for modelevaluation analysis to work
  if (Psites(sel) .eq. "AU-Preston") then
    lat = -37.7306
  end if
  lat@standard_name = "latitude"
  lat@units = "degrees_north"
  lat!0 = "y"
  lat&y = ny
  lat!1 = "x"
  lat&x = nx
  printVarSummary(lat)
  lon      = conform(var2d,fstfile->lon,(/0/))
  if (Psites(sel) .eq. "AU-Preston") then
    lon = 145.0145
  end if
  lon@standard_name = "longitude"
  lon@units = "degrees_east"
  copy_VarCoords(lat,lon)
  printVarSummary(lon)
  area     = conform(var2d,fstfile->area,(/0/))
  copy_VarCoords(lat,area)
  landfrac = conform(var2d,fstfile->landfrac,(/0/))
  copy_VarCoords(lat,landfrac)
  landmask = conform(var2d,fstfile->landmask,(/0/))
  copy_VarCoords(lat,landmask)
  pftmask  = conform(var2d,fstfile->pftmask,(/0/))
  copy_VarCoords(lat,pftmask)
  nbedrock = conform(var2d,fstfile->nbedrock,(/0/))
  copy_VarCoords(lat,nbedrock)

  nlevgrnd = dimsizes(levgrnd)
  nlevdcmp = dimsizes(levdcmp)

  data_1d  = addfile(fls_1d,"r")
  data_init = addfile(f_in,"r")
  data_1d_init = addfile(f_in_1d,"r")

  mcdate = data_1d->mcdate
  mcsec  = data_1d->mcsec

; FSA                            -> SWnet        net shortwave surface radiation (W/m2)
  SWnet = conform(var3d,data->FSA,(/0,1/))
  SWnet!0 = "time"
  SWnet!1 = "y"
  SWnet!2 = "x"
  SWnet&time = time
  SWnet&y = y
  SWnet&x = x
  SWnet@long_name = "net shortwave radiation"
  SWnet@standard_name = "surface_net_downward_shortwave_flux"
  SWnet@units = "W/m2"
  printVarSummary(SWnet)

; FIRA                           -> LWnet        net longwave radiation (W/m2)
; Request is positive down, our net longwave is positive up
  LWnet   =  - conform(var3d,data->FIRA,(/0,1/))
  copy_VarCoords(SWnet,LWnet)
  LWnet@long_name = "net longwave radiation"
  LWnet@standard_name = "surface_net_downward_longwave_flux"
  LWnet@units = "W/m2"
  printVarSummary(LWnet)
  print((/"MAX(LWnet): "+max(LWnet)/))
  print((/"MIN(LWnet): "+min(LWnet)/))

; FSDS                           -> SWdown     surface downwelling shortwave flux (W/m2)
  SWdown= conform(var3d,data->FSDS,(/0,1/))
  SWdown@long_name = "downward shortwave radiation"
  SWdown@standard_name = "surface_downwelling_shortwave_flux_in_air"
  SWdown@units = "W/m2"
  copy_VarCoords(SWnet,SWdown)

; FLDS                           -> LWdown     surface downwelling shortwave flux (W/m2)
  LWdown= conform(var3d,data->FLDS,(/0,1/))
  LWdown@long_name = "downward longwave radiation"
  LWdown@standard_name = "surface_downwelling_longwave_flux_in_air"
  LWdown@units = "W/m2"
  copy_VarCoords(SWnet,LWdown)

; FSR                            -> SWup       surface upwelling shortwave flux (W/m2)
  SWup= conform(var3d,data->FSR,(/0,1/))
  SWup@long_name = "upward shortwave radiation"
  SWup@standard_name = "surface_upwelling_shortwave_flux_in_air"
  SWup@units = "W/m2"
  copy_VarCoords(SWnet,SWup)

; FIRE                           -> LWup       surface upwelling longwave flux (W/m2)
  LWup= conform(var3d,data->FIRE,(/0,1/))
  LWup@long_name = "upward longwave radiation"
  LWup@standard_name = "surface_upwelling_longwave_flux_in_air"
  LWup@units = "W/m2"
  copy_VarCoords(SWnet,LWup)

; TBOT                           -> Tair       air temperature at forcing height (K) 
  Tair  = conform(var3d,data->TBOT,(/0,1/))
  Tair@long_name = "air temperature"
  Tair@standard_name = "air_temperature"
  Tair@units = "K"
  copy_VarCoords(SWnet,Tair)

; QBOT                           -> Qair       specific humidity at forcing height (kg/kg) 
  Qair  = conform(var3d,data->QBOT,(/0,1/))
  Qair@long_name = "specific humidity"
  Qair@standard_name = "surface_specific_humidity"
  Qair@units = "kg/kg"
  copy_VarCoords(SWnet,Qair)

; PBOT                           -> PSurf      air pressure at forcing height (Pa) 
  PSurf = conform(var3d,data->PBOT,(/0,1/))
  PSurf@long_name = "air pressure"
  PSurf@standard_name = "surface_air_pressure"
  PSurf@units = "Pa"
  copy_VarCoords(SWnet,PSurf)

; WIND                           -> Wind       wind at forcing height (m/s) 
  Wind  = conform(var3d,data->WIND,(/0,1/))
  Wind@long_name = "wind speed"
  Wind@standard_name = "wind_speed"
  Wind@units = "m/s"
  copy_VarCoords(SWnet,Wind)

; EFLX_LH_TOT            	 -> Qle        surface upward latent heat flux (W/m2)
  Qle   = conform(var3d,data->EFLX_LH_TOT,(/0,1/))
  Qle@long_name = "latent heat flux"
  Qle@standard_name = "surface_upward_latent_heat_flux"
  Qle@units = "W/m2"
  copy_VarCoords(SWnet,Qle)

; FSH                    	 -> Qh         surface upward sensible heat flux (W/m2)
  Qh    = conform(var3d,data->FSH,(/0,1/))
  Qh@long_name = "sensible heat flux"
  Qh@standard_name = "surface_upward_sensible_heat_flux"
  Qh@units = "W/m2"
  copy_VarCoords(SWnet,Qh)

; See my notes of 07-31-20
; FGR+URBAN_HEAT-URBAN_AC      	 -> Qstor      net storage heat flux in all materials (W/m2)
  Qstor = conform(var3d,data->FGR + data->URBAN_HEAT - data->URBAN_AC,(/0,1/))
  Qstor@long_name = "net storage heat flux in all materials"
  Qstor@standard_name = "surface_thermal_storage_heat_flux"
  Qstor@units = "W/m2"
  Qstor@_FillValue = SWdown@_FillValue
  copy_VarCoords(SWdown,Qstor)
  printVarSummary(Qstor)

; See my notes of 07-31-20
; URBAN_HEAT+WASTEHEAT           -> Qanth      anthropogenic heat flux (W/m2)
  Qanth = conform(var3d,data->URBAN_HEAT + data->WASTEHEAT,(/0,1/))
  Qanth@long_name = "anthropogenic heat flux"
  Qanth@standard_name = "surface_upward_heat_flux_due_to_anthropogenic_energy_consumption"
  Qanth@units = "W/m2"
  Qanth@_FillValue = SWdown@_FillValue
  copy_VarCoords(SWdown,Qanth)
  printVarSummary(Qanth)

; URBAN_HEAT+WASTEHEAT           -> Qanth_Qh   anthropogenic sensible heat flux (W/m2)
  Qanth_Qh= conform(var3d,data->URBAN_HEAT + data->WASTEHEAT,(/0,1/))
  Qanth_Qh@long_name = "anthropogenic sensible heat flux"
  Qanth_Qh@standard_name = "surface_upward_sensible_heat_flux_due_to_anthropogenic_energy_consumption"
  Qanth_Qh@units = "W/m2"
  Qanth_Qh@_FillValue = SWdown@_FillValue
  copy_VarCoords(SWdown,Qanth_Qh)
  printVarSummary(Qanth_Qh)

; 0                              -> Qanth_Qle  anthropogenic latent heat flux (W/m2)
  Qanth_Qle= Qanth_Qh
  Qanth_Qle= 0. 
  Qanth_Qle@long_name = "anthropogenic latent heat flux"
  Qanth_Qle@standard_name = "surface_upward_latent_heat_flux_due_to_anthropogenic_energy_consumption"
  Qanth_Qle@units = "W/m2"
  Qanth_Qle@_FillValue = SWdown@_FillValue
  copy_VarCoords(SWdown,Qanth_Qle)
  printVarSummary(Qanth_Qle)

; TSKIN                          -> AvgSurfT   surface temperature (K)
  AvgSurfT = conform(var3d,data->TSKIN,(/0,1/))
  AvgSurfT@long_name = "average surface temperature"
  AvgSurfT@standard_name = "surface_temperature"
  AvgSurfT@units = "K"
  copy_VarCoords(SWnet,AvgSurfT)

; FSR/FSDS                       -> Albedo     surface albedo (-)
  fsr   = data->FSR
  fsds  = data->FSDS
  if (any(fsds .eq. 0.)) then
    fsds = mask(fsds,(fsds .eq. 0.),False)
  end if
  Albedo  = conform(var3d,fsr/fsds,(/0,1/))
  Albedo@long_name = "surface albedo"
  Albedo@standard_name = "surface_albedo"
  Albedo@units = "1"
  Albedo@_FillValue = SWdown@_FillValue
  copy_VarCoords(SWdown,Albedo)

; FSRSNOW/FSDS                   -> SAlbedo      snow albedo (-)
  fsrsnow = data->SNOFSRVD + data->SNOFSRND + data->SNOFSRVI + data->SNOFSRNI
  SAlbedo = conform(var3d,fsrsnow/fsds,(/0,1/))
  SAlbedo@long_name = "snow albedo"
  SAlbedo@standard_name = "snow_and_ice_albedo"
  SAlbedo@units = "1"
  SAlbedo@_FillValue = SWdown@_FillValue
  copy_VarCoords(SWdown,SAlbedo)

; FSNO                   	 -> SnowFrac   snow covered fraction (-)
  SnowFrac = conform(var3d,data->FSNO,(/0,1/))
  SnowFrac@long_name = "snow covered fraction"
  SnowFrac@standard_name = "surface_snow_area_fraction"
  SnowFrac@units = "1"
  copy_VarCoords(SWdown,SnowFrac)

; RAIN                           -> Rainf      rainfall rate (kg/m2/s)
  Rainf    = conform(var3d,data->RAIN,(/0,1/))
  Rainf@long_name = "rainfall rate"
  Rainf@standard_name = "rainfall_flux"
  Rainf@units = "kg/m2/s"
  Rainf@_FillValue = SWdown@_FillValue
  copy_VarCoords(SWdown,Rainf)

; SNOW                           -> Snowf      snowfall rate (kg/m2/s)
  Snowf    = conform(var3d,data->SNOW,(/0,1/))
  Snowf@long_name = "snowfall rate"
  Snowf@standard_name = "snowfall_flux"
  Snowf@units = "kg/m2/s"
  Snowf@_FillValue = SWdown@_FillValue
  copy_VarCoords(SWdown,Snowf)

; QVEGT+QVEGE+QSOIL      	 -> Evap       total evapotranspiration (kg/m2/s)
  Evap    = conform(var3d,data->QVEGT + data->QVEGE + data->QSOIL,(/0,1/))
  Evap@long_name = "total evapotranspiration"
  Evap@standard_name = "surface_evapotranspiration"
  Evap@units = "kg/m2/s"
  Evap@_FillValue = SWdown@_FillValue
  copy_VarCoords(SWdown,Evap)

; QOVER+QH2OSFC          	 -> Qs         surface runoff (kg/m2/s)
  Qs = conform(var3d,data->QOVER + data->QH2OSFC,(/0,1/))
  Qs@long_name = "surface runoff"
  Qs@standard_name = "surface_runoff_flux"
  Qs@units = "kg/m2/s"
  Qs@_FillValue = SWdown@_FillValue
  copy_VarCoords(SWdown,Qs)

; QDRAI+QDRAI_PERCH      	 -> Qsb        subsurface runoff (kg/m2/s)
  Qsb = conform(var3d,data->QDRAI + data->QDRAI_PERCH,(/0,1/))
  Qsb@long_name = "subsurface runoff"
  Qsb@standard_name = "subsurface_runoff_flux"
  Qsb@units = "kg/m2/s"
  Qsb@_FillValue = SWdown@_FillValue
  copy_VarCoords(SWdown,Qsb)

; QSNOMELT                       -> Qsm        snow melt (kg/m2/s)
  Qsm   = conform(var3d,data->QSNOMELT,(/0,1/))
  Qsm@long_name = "snowmelt"
  Qsm@standard_name = "surface_snow_and_ice_melt_flux"
  Qsm@units = "kg/m2/s"
  copy_VarCoords(SWdown,Qsm)

; TAUX                           -> Qtau       momentum flux (N/m2) - can use just TAUX when U=V
; Request is positive down, our TAUX is positive up
  Qtau  =  - conform(var3d,data->TAUX,(/0,1/))
  Qtau@long_name = "momentum flux"
  Qtau@standard_name = "magnitude_of_surface_downward_stress"
  Qtau@units = "N/m2"
  copy_VarCoords(SWdown,Qtau)
  
; H2OSNO                 	 -> SWE        snow water equivalent (kg/m2)
  SWE   = conform(var3d,data->H2OSNO,(/0,1/))
  SWE@long_name = "snow water equivalent"
  SWE@standard_name = "surface_snow_amount"
  SWE@units = "kg/m2"
  copy_VarCoords(SWdown,SWE)

; Change in H2OSNO               -> DelSWE     change in snow water equivalent (kg/m2)
  if (add_WaterBalVars) then
    DelSWE = new((/ntim,y,x/),"float")
;   ; Don't get initial snow water equivalent from original netcdf file (not stupid file) that has the initialization step.
;   ; The initialization step does some snow stuff. Snow is initialized at zero mm.
;   init_h2osno = 0.
    ; Get initial snow water equivalent from original netcdf file (not stupid file) that has the initialization step
    init_h2osno   = data_init->H2OSNO(0,0)
    do t = 0,ntim-1
      if (t .eq. 0) then
        DelSWE(t,0,0) = SWE(t,0,0) - init_h2osno
      else
        DelSWE(t,0,0) = SWE(t,0,0) - SWE(t-1,0,0)
      end if
    end do
    DelSWE@long_name = "Change in snow water equivalent (increase)"
    DelSWE@standard_name = "change_over_time_in_surface_snow_and_ice_amount"
    DelSWE@units = "kg/m2"
    DelSWE@_FillValue = SWdown@_FillValue
    DelSWE!0 = "time"
    DelSWE!1 = "y"
    DelSWE!2 = "x"
    printVarSummary(DelSWE)
    print((/"MAX(DelSWE): "+max(DelSWE)/))
    print((/"MIN(DelSWE): "+min(DelSWE)/))
  end if

; TG(index 1 of 1d file)         -> RoofSurfT  roof surface temperature (K)
  RoofSurfT = SWdown
  RoofSurfT(:,0,0) = data_1d->TG(:,1)
  RoofSurfT@long_name = "roof surface temperature"
  RoofSurfT@standard_name = "surface_roof_skin_temperature"
  RoofSurfT@units = "K"
  RoofSurfT@_FillValue = SWdown@_FillValue
  copy_VarCoords(SWdown,RoofSurfT)
  printVarSummary(RoofSurfT)

; TG(avg index 2&3 of 1d file    -> WallSurfT  wall surface temperature (K)
  WallSurfT = SWdown
  WallSurfT(:,0,0) = (data_1d->TG(:,2) + data_1d->TG(:,3))/2.
  WallSurfT@long_name = "wall surface temperature"
  WallSurfT@standard_name = "surface_wall_skin_temperature"
  WallSurfT@units = "K"
  WallSurfT@_FillValue = SWdown@_FillValue
  copy_VarCoords(SWdown,WallSurfT)

; TG(index 4 of 1d file)         -> RoadSurfT  road surface temperature (K)
  RoadSurfT = SWdown
  RoadSurfT(:,0,0) = data_1d->TG(:,4)
  RoadSurfT@long_name = "road surface temperature"
  RoadSurfT@standard_name = "surface_road_skin_temperature"
  RoadSurfT@units = "K"
  RoadSurfT@_FillValue = SWdown@_FillValue
  copy_VarCoords(SWdown,RoadSurfT)
  printVarSummary(RoadSurfT)

; TSA                            -> TairCanyon    air temperature in street canyon (K)
  if (version .eq. "v1" .or. version .eq. "v2" .or. version .eq. "v3") then
    TairCanyon = conform(var3d,data->TSA,(/0,1/))
    TairCanyon@long_name = "air temperature in canyon"
    TairCanyon@standard_name = "air_temperature_in_street_canyon"
    TairCanyon@units = "K"
    copy_VarCoords(SWdown,TairCanyon)
  else
    if (version .eq. "v4") then
      ; check for pfts1d_itype_lunit = 9 and use first instance (all values are the same at the pft level)
      pfts1d_itype_lunit = data_1d->pfts1d_itype_lunit
      indx = ind(pfts1d_itype_lunit .eq. 9)
      TairCanyon = SWdown
      TairCanyon = data_1d->TSA(:,indx(0))
      TairCanyon@long_name = "air temperature in canyon"
      TairCanyon@standard_name = "air_temperature_in_street_canyon"
      TairCanyon@units = "K"
      copy_VarCoords(SWdown,TairCanyon)
      delete(indx)
      delete(pfts1d_itype_lunit)
    end if
  end if

; TBUILD                         -> TairBuilding  air temperature in buildings (K)
  TairBuilding = conform(var3d,data->TBUILD,(/0,1/))
  TairBuilding@long_name = "air temperature in buildings"
  TairBuilding@standard_name = "air_temperature_in_buildings"
  TairBuilding@units = "K"
  copy_VarCoords(SWdown,TairBuilding)

; SNO_T(0)                       -> SnowT         snow surface temperature (K)
  tsno = data->SNO_T
  tsno = where(ismissing(tsno),0.,tsno)
  printVarSummary(tsno)
  sno_existence = data->SNO_EXISTENCE
  printVarSummary(sno_existence)
  SnowT = SWdown
  do i = 0,ntim-1
    indx1 = maxind(sno_existence(i,:,0))
    if (sno_existence(i,indx1,0) .eq. 1) then
      SnowT(i,0,0) = (/tsno(i,indx1,0)/)
    else
      SnowT(i,0,0) = SWdown@_FillValue
    end if
  end do
  SnowT@long_name = "snow surface temperature"
  SnowT@standard_name = "surface_snow_skin_temperature"
  SnowT@units = "K"
  SnowT@_FillValue = SWdown@_FillValue
  copy_VarCoords(SWdown,SnowT)
  printVarSummary(SnowT)

; TSOI(0)                        -> BaresoilT     temperature of bare soil (K)
  if (version .eq. "v1" .or. version .eq. "v2" .or. version .eq. "v3") then
    BaresoilT  = SWdown
    BaresoilT(:,0,0)  = data_1d->TSOI(:,0,5)
    BaresoilT@long_name = "temperature of bare soil"
    BaresoilT@standard_name = "surface_ground_skin_temperature"
    BaresoilT@units = "K"
    BaresoilT@_FillValue = SWdown@_FillValue
    copy_VarCoords(SWdown,BaresoilT)
    printVarSummary(BaresoilT)
  else
    if (version .eq. "v4") then
      ; check for cols1d_itype_lunit = 1 (should only be one instance)
      cols1d_itype_lunit = data_1d->cols1d_itype_lunit
      BaresoilT  = SWdown
      indx = ind(cols1d_itype_lunit .eq. 1)
      BaresoilT(:,0,0)  = data_1d->TSOI(:,0,indx(0))
      BaresoilT@long_name = "temperature of bare soil"
      BaresoilT@standard_name = "surface_ground_skin_temperature"
      BaresoilT@units = "K"
      BaresoilT@_FillValue = SWdown@_FillValue
      copy_VarCoords(SWdown,BaresoilT)
      printVarSummary(BaresoilT)
      delete(indx)
      delete(cols1d_itype_lunit)
    end if
  end if

; SOILLIQ+SOILICE                -> SoilMoist     average layer soil moisture (kg/m2)
  SoilMoist = new((/ntim,nlevsoi,y,x/),"float")
  if (version .eq. "v1" .or. version .eq. "v2" .or. version .eq. "v3") then
    cols1d_wtgcell        = dble2flt(data_1d->cols1d_wtgcell)
    if (lu_soilmoist) then
      ; Here we are adding the roof water and impervious road water, averaged appropriately, to the top layer of the pervious
      ; road soil moisture in an attempt to construct a water balance.
      ; NOTE that this only results in a reasonable water balance (average imbalance over the AU-Preston simulation
      ; is 1.24445e-11 mm/s) in places with no snow.
      soilice_roof          = (data_1d->SOILICE(:,0,1))*cols1d_wtgcell(1)                ; Roof (top layer only)
      soilice_imprvroad     = (1.-cols1d_wtgcell(1))*(data_1d->SOILICE(:,0,4)) * \
                              (cols1d_wtgcell(4)/(cols1d_wtgcell(4)+cols1d_wtgcell(5)))  ; Imprv Road (top layer only)
      soilice_pervroad      = (1.-cols1d_wtgcell(1))*(data_1d->SOILICE(:,:,5)) * \
                              (cols1d_wtgcell(5)/(cols1d_wtgcell(4)+cols1d_wtgcell(5)))  ; Perv Road (all layers)
      soilice               = soilice_pervroad
      soilice(:,0)          = soilice_pervroad(:,0) + soilice_roof + soilice_imprvroad
      soilliq_roof          = (data_1d->SOILLIQ(:,0,1))*cols1d_wtgcell(1)                ; Roof (top layer only)
      soilliq_imprvroad     = (1.-cols1d_wtgcell(1))*(data_1d->SOILLIQ(:,0,4)) * \
                              (cols1d_wtgcell(4)/(cols1d_wtgcell(4)+cols1d_wtgcell(5)))  ; Imprv Road (top layer only)
      soilliq_pervroad      = (1.-cols1d_wtgcell(1))*(data_1d->SOILLIQ(:,:,5)) * \
                              (cols1d_wtgcell(5)/(cols1d_wtgcell(4)+cols1d_wtgcell(5)))  ; Perv Road (all layers)
      soilliq               = soilliq_pervroad
      soilliq(:,0)          = soilliq_pervroad(:,0) + soilliq_roof + soilliq_imprvroad
      soilh2o               = soilice + soilliq
      SoilMoist(:,:,0,0)    = soilh2o
    else
      soilice = data_1d->SOILICE(:,:,5)
      soilliq = data_1d->SOILLIQ(:,:,5)
      ; This was reported initially for just unweighted pervious road
      ; For water balance project report this as weighted gridcell average
      if (add_WaterBalVars) then
        soilh2o = (1.-cols1d_wtgcell(1))*(cols1d_wtgcell(5)/(cols1d_wtgcell(4)+cols1d_wtgcell(5))) * \
                  (soilice + soilliq) 
      else
        soilh2o = soilice + soilliq
      end if
      SoilMoist(:,:,0,0) = soilh2o
    end if
    SoilMoist@long_name = "average layer soil moisture"
    SoilMoist@standard_name = "moisture_content_of_soil_layer"
    SoilMoist@units = "kg/m2"
    SoilMoist@_FillValue = SWdown@_FillValue
    SoilMoist!0 = "time"
    SoilMoist!1 = "soil_layer"
    SoilMoist!2 = "y"
    SoilMoist!3 = "x"
    printVarSummary(SoilMoist)
    print((/"MAX(SoilMoist): "+max(SoilMoist)/))
    print((/"MIN(SoilMoist): "+min(SoilMoist)/))
    print((/"MAX(sum(SoilMoist)): "+max(dim_sum_n(SoilMoist(:,:,0,0),1))/))
    print((/"MIN(sum(SoilMoist)): "+min(dim_sum_n(SoilMoist(:,:,0,0),1))/))
  else
    if (version .eq. "v4") then
      if (lu_soilmoist) then
        ; Here we are adding the roof water, impervious road water, averaged appropriately, to the top 
        ; layer of the pervious road soil moisture, and the canopy water to the total water, in an attempt to construct a water balance.
        ; NOTE that this results in a water imbalance of e-6 to e-9 at a given time step and an average imbalance of 6.6e-9 mm/s (AU-Preston). 
        ; This is not as good as the v3 simulations.  Not sure what if anything is missing.  Also doesn't work in the presence of snow.
        cols1d_wtgcell        = dble2flt(data_1d->cols1d_wtgcell)
        soilice_roof          = (1.-cols1d_wtgcell(0))*(data_1d->SOILICE(:,0,1))*cols1d_wtgcell(1)          ; Roof (top layer only)
        soilice_imprvroad     = (1.-cols1d_wtgcell(0))*(1.-cols1d_wtgcell(1))*(data_1d->SOILICE(:,0,4)) * \ ; Imprv Road (top layer only)
                                (cols1d_wtgcell(4)/(cols1d_wtgcell(4)+cols1d_wtgcell(5)))
        soilice_pervroad      = (1.-cols1d_wtgcell(0))*(1.-cols1d_wtgcell(1))*(data_1d->SOILICE(:,:,5)) * \ ; Prev Road (all layers)
                                (cols1d_wtgcell(5)/(cols1d_wtgcell(4)+cols1d_wtgcell(5)))
        soilice_pervroad(:,0) = soilice_pervroad(:,0) + soilice_roof + soilice_imprvroad
        soilice_veg           = (data_1d->SOILICE(:,:,0))*cols1d_wtgcell(0)                                 ; Veg (all layers)
        soilice               = soilice_pervroad + soilice_veg
        soilliq_roof          = (1.-cols1d_wtgcell(0))*(data_1d->SOILLIQ(:,0,1))*cols1d_wtgcell(1)          ; Roof (top layer only)
        soilliq_imprvroad     = (1.-cols1d_wtgcell(0))*(1.-cols1d_wtgcell(1))*(data_1d->SOILLIQ(:,0,4)) * \ ; Imprv Road (top layer only)
                                (cols1d_wtgcell(4)/(cols1d_wtgcell(4)+cols1d_wtgcell(5)))
        soilliq_pervroad      = (1.-cols1d_wtgcell(0))*(1.-cols1d_wtgcell(1))*(data_1d->SOILLIQ(:,:,5)) * \ ; Perv Road (all layers)
                                (cols1d_wtgcell(5)/(cols1d_wtgcell(4)+cols1d_wtgcell(5)))
        soilliq_pervroad(:,0) = soilliq_pervroad(:,0) + soilliq_roof + soilliq_imprvroad
        soilliq_veg           = (data_1d->SOILLIQ(:,:,0))*cols1d_wtgcell(0)                                 ; Veg (all layers)
        soilliq               = soilliq_pervroad + soilliq_veg
        soilh2o_veg           = soilice_veg + soilliq_veg
        soilh2o               = soilice + soilliq
        pfts1d_itype_col      = data_1d->pfts1d_itype_col
        indxpft               = ind(pfts1d_itype_col .eq. 1)
        pfts1d_wtcol          = data_1d->pfts1d_wtcol
        h2ocan                = cols1d_wtgcell(0)*dim_sum_n((data_1d->H2OCAN(:,indxpft))*conform_dims((/ntim,dimsizes(indxpft)/),dble2flt(pfts1d_wtcol(indxpft)),(/1/)),1)
        h2osfc                = cols1d_wtgcell(0)*(data_1d->H2OSFC(:,0))
        soilh2o(:,0)          = soilh2o(:,0) + h2ocan + h2osfc
        SoilMoist(:,:,0,0)    = soilh2o
        delete(indxpft)
        delete(pfts1d_wtcol)
      else
        ; This was reported initially for just unweighted veg landunit
        ; For water balance project report this as weighted gridcell average
        if (add_WaterBalVars) then
          cols1d_wtgcell     = dble2flt(data_1d->cols1d_wtgcell)
          soilice_veg        = (data_1d->SOILICE(:,:,0))*cols1d_wtgcell(0)                                 ; Veg (all layers)
          soilice_pervroad   = (1.-cols1d_wtgcell(0))*(1.-cols1d_wtgcell(1))*(data_1d->SOILICE(:,:,5)) * \ ; Prev Road (all layers)
                                  (cols1d_wtgcell(5)/(cols1d_wtgcell(4)+cols1d_wtgcell(5)))
          soilliq_veg        = (data_1d->SOILLIQ(:,:,0))*cols1d_wtgcell(0)                                 ; Veg (all layers)
          soilliq_pervroad   = (1.-cols1d_wtgcell(0))*(1.-cols1d_wtgcell(1))*(data_1d->SOILLIQ(:,:,5)) * \ ; Prev Road (all layers)
                                  (cols1d_wtgcell(5)/(cols1d_wtgcell(4)+cols1d_wtgcell(5)))
          soilice            = soilice_veg + soilice_pervroad
          soilliq            = soilliq_veg + soilliq_pervroad
          soilh2o            = soilice + soilliq
          SoilMoist(:,:,0,0) = soilh2o
        else
          ; check for cols1d_itype_lunit = 1 (should only be one instance)
          cols1d_itype_lunit = data_1d->cols1d_itype_lunit
          indx = ind(cols1d_itype_lunit .eq. 1)
          soilice = data_1d->SOILICE(:,:,indx(0))
          soilliq = data_1d->SOILLIQ(:,:,indx(0))
          soilh2o = soilice + soilliq
          SoilMoist(:,:,0,0) = soilh2o
          delete(indx)
          delete(cols1d_itype_lunit)
        end if
      end if
      SoilMoist@long_name = "average layer soil moisture"
      SoilMoist@standard_name = "moisture_content_of_soil_layer"
      SoilMoist@units = "kg/m2"
      SoilMoist@_FillValue = SWdown@_FillValue
      SoilMoist!0 = "time"
      SoilMoist!1 = "soil_layer"
      SoilMoist!2 = "y"
      SoilMoist!3 = "x"
      printVarSummary(SoilMoist)
      print((/"MAX(SoilMoist): "+max(SoilMoist)/))
      print((/"MIN(SoilMoist): "+min(SoilMoist)/))
    end if
  end if

; Change in SOILLIQ+SOILICE      -> DelSoilMoist change in soil moisture (kg/m2)
  if (add_WaterBalVars) then
    if (version .eq. "v1" .or. version .eq. "v2" .or. version .eq. "v3") then
      DelSoilMoist = new((/ntim,y,x/),"float")
      cols1d_wtgcell = dble2flt(data_1d->cols1d_wtgcell)
      ; Get initial pervious road soil moisture from original netcdf file (not stupid file) that has the initialization step
      init_soilliq = dim_sum(data_1d_init->SOILLIQ(0,:,5))
      init_soilice = dim_sum(data_1d_init->SOILICE(0,:,5))
      init_soilh2o = (1.-cols1d_wtgcell(1))*(cols1d_wtgcell(5)/(cols1d_wtgcell(4)+cols1d_wtgcell(5))) * \
                     (init_soilliq + init_soilice)
;     print((/"init_soilh2o: "+init_soilh2o/))
;     print((/"SoilMoist t=0: "+dim_sum(SoilMoist(0,:,0,0))/))
      do t = 0,ntim-1
        if (t .eq. 0) then
          DelSoilMoist(t,0,0) = dim_sum(SoilMoist(t,:,0,0)) - init_soilh2o
        else
          DelSoilMoist(t,0,0) = dim_sum(SoilMoist(t,:,0,0)) - dim_sum(SoilMoist(t-1,:,0,0))
        end if
;       print((/"DelSoilMoist: "+DelSoilMoist(t,0,0)/))
      end do
    else
      if (version .eq. "v4") then
        DelSoilMoist = new((/ntim,y,x/),"float")
        cols1d_wtgcell = dble2flt(data_1d->cols1d_wtgcell)
        ; Get initial veg and pervious road soil moisture from original netcdf file (not stupid file) that has the initialization step
        init_soilliq_veg      = dim_sum(data_1d_init->SOILLIQ(0,:,0))
        init_soilice_veg      = dim_sum(data_1d_init->SOILICE(0,:,0))
        init_soilh2o_veg      = cols1d_wtgcell(0)*(init_soilliq_veg + init_soilice_veg)
        init_soilliq_pervroad = dim_sum(data_1d_init->SOILLIQ(0,:,5))
        init_soilice_pervroad = dim_sum(data_1d_init->SOILICE(0,:,5))
        init_soilh2o_pervroad = (1.-cols1d_wtgcell(0))*(1.-cols1d_wtgcell(1)) * \
                                (cols1d_wtgcell(5)/(cols1d_wtgcell(4)+cols1d_wtgcell(5))) * \
                                (init_soilliq_pervroad + init_soilice_pervroad)
        init_soilh2o          = init_soilh2o_veg + init_soilh2o_pervroad
;       print((/"init_soilh2o: "+init_soilh2o/))
;       print((/"SoilMoist t=0: "+dim_sum(SoilMoist(0,:,0,0))/))
        do t = 0,ntim-1
          if (t .eq. 0) then
            DelSoilMoist(t,0,0) = dim_sum(SoilMoist(t,:,0,0)) - init_soilh2o
          else
            DelSoilMoist(t,0,0) = dim_sum(SoilMoist(t,:,0,0)) - dim_sum(SoilMoist(t-1,:,0,0))
          end if
;         print((/"DelSoilMoist: "+DelSoilMoist(t,0,0)/))
        end do
      end if
    end if
    DelSoilMoist@long_name = "change in soil moisture (increase)"
    DelSoilMoist@standard_name = "change_over_time_in_mass_content_of_water_in_soil"
    DelSoilMoist@units = "kg/m2"
    DelSoilMoist@_FillValue = SWdown@_FillValue
    DelSoilMoist!0 = "time"
    DelSoilMoist!1 = "y"
    DelSoilMoist!2 = "x"
    printVarSummary(DelSoilMoist)
    print((/"MAX(DelSoilMoist): "+max(DelSoilMoist)/))
    print((/"MIN(DelSoilMoist): "+min(DelSoilMoist)/))
  end if

; Change in roof/imprvrd water   -> DelIntercept change in interception storage (kg/m2)
  if (add_WaterBalVars) then
    if (version .eq. "v1" .or. version .eq. "v2" .or. version .eq. "v3") then
      DelIntercept = new((/ntim,y,x/),"float")
      cols1d_wtgcell = dble2flt(data_1d->cols1d_wtgcell)
      ; Get initial roof and impervious road ponding from original netcdf file (not stupid file) that has the initialization step
      init_soilice_roof      = (data_1d_init->SOILICE(0,0,1))*cols1d_wtgcell(1)                ; Roof (top layer only)
      init_soilice_imprvroad = (1.-cols1d_wtgcell(1))*(data_1d_init->SOILICE(0,0,4)) * \
                               (cols1d_wtgcell(4)/(cols1d_wtgcell(4)+cols1d_wtgcell(5)))  ; Imprv Road (top layer only)
      init_soilliq_roof      = (data_1d_init->SOILLIQ(0,0,1))*cols1d_wtgcell(1)                ; Roof (top layer only)
      init_soilliq_imprvroad = (1.-cols1d_wtgcell(1))*(data_1d_init->SOILLIQ(0,0,4)) * \
                               (cols1d_wtgcell(4)/(cols1d_wtgcell(4)+cols1d_wtgcell(5)))  ; Imprv Road (top layer only)
      init_soilh2o_roof_imprvroad = init_soilice_roof + init_soilice_imprvroad + init_soilliq_roof + init_soilliq_imprvroad 
      print((/"init_soilh2o_roof_imprvroad: "+init_soilh2o_roof_imprvroad/))
      soilice_roof      = (data_1d->SOILICE(:,0,1))*cols1d_wtgcell(1)                ; Roof (top layer only)
      soilice_imprvroad = (1.-cols1d_wtgcell(1))*(data_1d->SOILICE(:,0,4)) * \
                          (cols1d_wtgcell(4)/(cols1d_wtgcell(4)+cols1d_wtgcell(5)))  ; Imprv Road (top layer only)
      soilliq_roof      = (data_1d->SOILLIQ(:,0,1))*cols1d_wtgcell(1)                ; Roof (top layer only)
      soilliq_imprvroad = (1.-cols1d_wtgcell(1))*(data_1d->SOILLIQ(:,0,4)) * \
                          (cols1d_wtgcell(4)/(cols1d_wtgcell(4)+cols1d_wtgcell(5)))  ; Imprv Road (top layer only)
      soilh2o_roof_imprvroad = soilice_roof + soilice_imprvroad + soilliq_roof + soilliq_imprvroad
      do t = 0,ntim-1
        if (t .eq. 0) then
          DelIntercept(t,0,0) = soilh2o_roof_imprvroad(t) - init_soilh2o_roof_imprvroad
        else
          DelIntercept(t,0,0) = soilh2o_roof_imprvroad(t) - soilh2o_roof_imprvroad(t-1)
        end if
;       print((/"DelIntercept: "+DelIntercept(t,0,0)/))
      end do
    else
      if (version .eq. "v4") then
        DelIntercept = new((/ntim,y,x/),"float")
        cols1d_wtgcell = dble2flt(data_1d->cols1d_wtgcell)
        ; Get initial canopy, roof, and impervious road ponding from original netcdf file (not stupid file) that has the initialization step
        pfts1d_itype_col       = data_1d->pfts1d_itype_col
        indxpft                = ind(pfts1d_itype_col .eq. 1)
        pfts1d_wtcol           = data_1d->pfts1d_wtcol
        print(indxpft)
        print(pfts1d_wtcol)
        init_h2ocan            = cols1d_wtgcell(0)*dim_sum((data_1d_init->H2OCAN(0,indxpft))* \
                                 dble2flt(pfts1d_wtcol(indxpft)))
        init_h2osfc            = cols1d_wtgcell(0)*(data_1d_init->H2OCAN(0,0))
        init_soilice_roof      = (1.-cols1d_wtgcell(0))*(data_1d_init->SOILICE(0,0,1))*cols1d_wtgcell(1)            ; Roof (top layer only)
        init_soilice_imprvroad = (1.-cols1d_wtgcell(0))*(1.-cols1d_wtgcell(1))*(data_1d_init->SOILICE(0,0,4)) * \   ; Imprv Road (top layer only)
                                 (cols1d_wtgcell(4)/(cols1d_wtgcell(4)+cols1d_wtgcell(5)))
        init_soilliq_roof      = (1.-cols1d_wtgcell(0))*(data_1d_init->SOILLIQ(0,0,1))*cols1d_wtgcell(1)            ; Roof (top layer only)
        init_soilliq_imprvroad = (1.-cols1d_wtgcell(0))*(1.-cols1d_wtgcell(1))*(data_1d_init->SOILLIQ(0,0,4)) * \   ; Imprv Road (top layer only)
                                 (cols1d_wtgcell(4)/(cols1d_wtgcell(4)+cols1d_wtgcell(5)))   
        init_soilh2o_canopy_roof_imprvroad = init_h2ocan + init_h2osfc + init_soilice_roof + init_soilice_imprvroad + \
                                 init_soilliq_roof + init_soilliq_imprvroad 
        print((/"init_soilh2o_canopy_roof_imprvroad: "+init_soilh2o_canopy_roof_imprvroad/))
        h2ocan            = cols1d_wtgcell(0)*dim_sum_n((data_1d->H2OCAN(:,indxpft))* \
                            conform_dims((/ntim,dimsizes(indxpft)/),dble2flt(pfts1d_wtcol(indxpft)),(/1/)),1)
        h2osfc            = cols1d_wtgcell(0)*(data_1d->H2OSFC(:,0))
        soilice_roof      = (1.-cols1d_wtgcell(0))*(data_1d->SOILICE(:,0,1))*cols1d_wtgcell(1)                ; Roof (top layer only)
        soilice_imprvroad = (1.-cols1d_wtgcell(0))*(1.-cols1d_wtgcell(1))*(data_1d->SOILICE(:,0,4)) * \       ; Imprv Road (top layer only)
                            (cols1d_wtgcell(4)/(cols1d_wtgcell(4)+cols1d_wtgcell(5)))
        soilliq_roof      = (1.-cols1d_wtgcell(0))*(data_1d->SOILLIQ(:,0,1))*cols1d_wtgcell(1)                ; Roof (top layer only)
        soilliq_imprvroad = (1.-cols1d_wtgcell(0))*(1.-cols1d_wtgcell(1))*(data_1d->SOILLIQ(:,0,4)) * \       ; Imprv Road (top layer only)
                            (cols1d_wtgcell(4)/(cols1d_wtgcell(4)+cols1d_wtgcell(5)))
        soilh2o_canopy_roof_imprvroad = h2ocan + h2osfc + soilice_roof + soilice_imprvroad + soilliq_roof + soilliq_imprvroad
        do t = 0,ntim-1
          if (t .eq. 0) then
            DelIntercept(t,0,0) = soilh2o_canopy_roof_imprvroad(t) - init_soilh2o_canopy_roof_imprvroad
          else
            DelIntercept(t,0,0) = soilh2o_canopy_roof_imprvroad(t) - soilh2o_canopy_roof_imprvroad(t-1)
          end if
;         print((/"DelIntercept: "+DelIntercept(t,0,0)/))
        end do
      end if
    end if
    DelIntercept@long_name = "Change in interception storage (increase)"
    DelIntercept@standard_name = "change_over_time_in_canopy_water_amount"
    DelIntercept@units = "kg/m2"
    DelIntercept@_FillValue = SWdown@_FillValue
    DelIntercept!0 = "time"
    DelIntercept!1 = "y"
    DelIntercept!2 = "x"
    printVarSummary(DelIntercept)
    print((/"MAX(DelIntercept): "+max(DelIntercept)/))
    print((/"MIN(DelIntercept): "+min(DelIntercept)/))
  end if

; TSOI                           -> SoilTemp      average layer soil temperature (K)
  SoilTemp = new((/ntim,nlevsoi,y,x/),"float")
  if (version .eq. "v1" .or. version .eq. "v2" .or. version .eq. "v3") then
    SoilTemp(:,:,0,0) = data_1d->TSOI(:,0:19,5)
    SoilTemp@long_name = "average layer soil temperature"
    SoilTemp@standard_name = "soil_temperature"
    SoilTemp@units = "K"
    SoilTemp@_FillValue = SWdown@_FillValue
    copy_VarCoords(SoilMoist,SoilTemp)
    printVarSummary(SoilTemp)
    print((/"MAX(SoilTemp): "+max(SoilTemp)/))
    print((/"MIN(SoilTemp): "+min(SoilTemp)/))
  else
    if (version .eq. "v4") then
      ; check for cols1d_itype_lunit = 1 (should only be one instance)
      cols1d_itype_lunit = data_1d->cols1d_itype_lunit
      indx = ind(cols1d_itype_lunit .eq. 1)
      SoilTemp(:,:,0,0) = data_1d->TSOI(:,0:19,indx(0))
      SoilTemp@long_name = "average layer soil temperature"
      SoilTemp@standard_name = "soil_temperature"
      SoilTemp@units = "K"
      SoilTemp@_FillValue = SWdown@_FillValue
      copy_VarCoords(SoilMoist,SoilTemp)
      printVarSummary(SoilTemp)
      print((/"MAX(SoilTemp): "+max(SoilTemp)/))
      print((/"MIN(SoilTemp): "+min(SoilTemp)/))
      delete(indx)
      delete(cols1d_itype_lunit)
    end if
  end if

; Create netcdf output file
  if (do_netcdf) then

  cdf_filename = f_out_dir+f_out_nc
  system("/bin/rm -f " + cdf_filename)

; setfileoption("nc","Format","LargeFile")
  setfileoption("nc","Format","NetCDF4Classic")

  cdf_file_out = addfile(cdf_filename,"c")

  setfileoption(cdf_file_out,"DefineMode",True)

  fAtt = True
  if (version .eq. "v1" .or. version .eq. "v2" .or. version .eq. "v3") then
    fAtt@title = "CLMU5 model output for the Urban-PLUMBER project"
  else
    if (version .eq. "v4") then
      fAtt@title = "CLMU5VEG model output for the Urban-PLUMBER project"
    end if
  end if
  fAtt@site = Psites(sel)
  fAtt@experiment = expr
  fAtt@institution = "NCAR"
  fAtt@primary_contact = "Keith Oleson (oleson@ucar.edu)"
  fAtt@secondary_contact = "Keith Oleson (oleson@ucar.edu)"
  if (version .eq. "v1" .or. version .eq. "v2" .or. version .eq. "v3") then
    fAtt@model = "CLMU5"
  else
    if (version .eq. "v4") then
      fAtt@model = "CLMU5VEG"
    end if
  end if
  fAtt@source = "Community Land Model Urban version 5 - release-clm5.0.34"
  fAtt@references = "Oleson, K.W., and J. Feddema, 2019: Parameterization and surface data improvements and new capabilities for the Community Land Model Urban (CLMU), JAMES, 11, doi:10.1029/2018MS001586."
  fAtt@repository = "https://urldefense.com/v3/__https://github.com/ESCOMP/CTSM/releases/tag/release-clm5.0.34__;!!PDiH4ENfjr2_Jw!BZ9nqK13xFDelwkqWL8bgcyod381s7FoHRYHOLQKZs9762l9KahcY9LPIQF5HJgQ52CcDeRFmtRyj1AbPltxJDhMwYmIZ-Ry$ [github[.]com]"
  if (Msites(sel) .eq. "AU-Pre" .or. Msites(sel) .eq. "FI-Tor") then
    fAtt@site_experience = "Yes"
  else
    fAtt@site_experience = "No"
  end if
  fAtt@additional_data = "None"
  if (expr .eq. "baseline") then
    fAtt@comment = "Roughness length/Displacement height derived from Macdonald et al. (1998) using provided baseline input parameters"
  else
    if (expr .eq. "detailed") then
      fAtt@comment = "Roughness length/Displacement height from Macdonald et al. (1998) as provided in parameter input file (roughness_length_momentum_mac and displacement_height_mac)"
    end if
  end if
  fAtt@Conventions = fstfile@Conventions
  fAtt@history = "created on "+systemfunc("date")
  fAtt@hostname = fstfile@hostname
  fAtt@username = "oleson"
  fAtt@version = fstfile@version
  fAtt@Surface_dataset = fstfile@Surface_dataset
  fAtt@Initial_conditions_dataset = fstfile@Initial_conditions_dataset
  fAtt@PFT_physiological_constants_dataset = fstfile@PFT_physiological_constants_dataset
  fAtt@time_period_freq = fstfile@time_period_freq
  fileattdef(cdf_file_out,fAtt)

  dimNames = (/"string_length","hist_interval","levgrnd","soil_layer","levsno","levdcmp","y",  "x",  "time"/)
  dimSizes = (/string_length,  hist_interval,  nlevgrnd, nlevsoi,     nlevsno, nlevdcmp, y,    x,    1/)
  dimUnlim = (/False,          False,          False,    False,       False,   False,    False,False,True/)
  filedimdef(cdf_file_out,dimNames,dimSizes,dimUnlim)

  filevardef(cdf_file_out,"latitude",typeof(lat),getvardims(lat))
  filevardef(cdf_file_out,"longitude",typeof(lon),getvardims(lon))
  filevardef(cdf_file_out,"y",typeof(y),getvardims(y))
  filevardef(cdf_file_out,"x",typeof(x),getvardims(x))
  filevardef(cdf_file_out,"time",typeof(time),getvardims(time))
  filevardef(cdf_file_out,"area",typeof(area),getvardims(area))
  filevardef(cdf_file_out,"landfrac",typeof(landfrac),getvardims(landfrac))
  filevardef(cdf_file_out,"landmask",typeof(landmask),getvardims(landmask))
  filevardef(cdf_file_out,"pftmask",typeof(pftmask),getvardims(pftmask))
  filevardef(cdf_file_out,"levgrnd",typeof(levgrnd),getvardims(levgrnd))
  filevardef(cdf_file_out,"levdcmp",typeof(levdcmp),getvardims(levdcmp))
  filevardef(cdf_file_out,"nbedrock",typeof(nbedrock),getvardims(nbedrock))
  filevardef(cdf_file_out,"soil_layer",typeof(soil_layer),getvardims(soil_layer))
  filevardef(cdf_file_out,"SWnet",typeof(SWnet),getvardims(SWnet))
  filevardef(cdf_file_out,"LWnet",typeof(LWnet),getvardims(LWnet))
  filevardef(cdf_file_out,"SWdown",typeof(SWdown),getvardims(SWdown))
  filevardef(cdf_file_out,"LWdown",typeof(LWdown),getvardims(LWdown))
  filevardef(cdf_file_out,"SWup",typeof(SWup),getvardims(SWup))
  filevardef(cdf_file_out,"LWup",typeof(LWup),getvardims(LWup))
  filevardef(cdf_file_out,"Tair",typeof(Tair),getvardims(Tair))
  filevardef(cdf_file_out,"Qair",typeof(Qair),getvardims(Qair))
  filevardef(cdf_file_out,"PSurf",typeof(PSurf),getvardims(PSurf))
  filevardef(cdf_file_out,"Wind",typeof(Wind),getvardims(Wind))
  filevardef(cdf_file_out,"Qle",typeof(Qle),getvardims(Qle))
  filevardef(cdf_file_out,"Qh",typeof(Qh),getvardims(Qh))
  filevardef(cdf_file_out,"Qstor",typeof(Qstor),getvardims(Qstor))
  filevardef(cdf_file_out,"Qanth",typeof(Qanth),getvardims(Qanth))
  filevardef(cdf_file_out,"Qanth_Qh",typeof(Qanth_Qh),getvardims(Qanth_Qh))
  filevardef(cdf_file_out,"Qanth_Qle",typeof(Qanth_Qle),getvardims(Qanth_Qle))
  filevardef(cdf_file_out,"AvgSurfT",typeof(AvgSurfT),getvardims(AvgSurfT))
  filevardef(cdf_file_out,"Albedo",typeof(Albedo),getvardims(Albedo))
  filevardef(cdf_file_out,"SAlbedo",typeof(SAlbedo),getvardims(SAlbedo))
  filevardef(cdf_file_out,"SnowFrac",typeof(SnowFrac),getvardims(SnowFrac))
  filevardef(cdf_file_out,"Rainf",typeof(Rainf),getvardims(Rainf))
  filevardef(cdf_file_out,"Snowf",typeof(Snowf),getvardims(Snowf))
  filevardef(cdf_file_out,"Evap",typeof(Evap),getvardims(Evap))
  filevardef(cdf_file_out,"Qs",typeof(Qs),getvardims(Qs))
  filevardef(cdf_file_out,"Qsb",typeof(Qsb),getvardims(Qsb))
  filevardef(cdf_file_out,"Qsm",typeof(Qsm),getvardims(Qsm))
  filevardef(cdf_file_out,"Qtau",typeof(Qtau),getvardims(Qtau))
  filevardef(cdf_file_out,"SWE",typeof(SWE),getvardims(SWE))
  if (add_WaterBalVars) then
    filevardef(cdf_file_out,"DelSWE",typeof(DelSWE),getvardims(DelSWE))
  end if
  filevardef(cdf_file_out,"RoofSurfT",typeof(RoofSurfT),getvardims(RoofSurfT))
  filevardef(cdf_file_out,"WallSurfT",typeof(WallSurfT),getvardims(WallSurfT))
  filevardef(cdf_file_out,"RoadSurfT",typeof(RoadSurfT),getvardims(RoadSurfT))
  filevardef(cdf_file_out,"TairCanyon",typeof(TairCanyon),getvardims(TairCanyon))
  filevardef(cdf_file_out,"TairBuilding",typeof(TairBuilding),getvardims(TairBuilding))
  filevardef(cdf_file_out,"SnowT",typeof(SnowT),getvardims(SnowT))
  filevardef(cdf_file_out,"BaresoilT",typeof(BaresoilT),getvardims(BaresoilT))
  filevardef(cdf_file_out,"SoilMoist",typeof(SoilMoist),getvardims(SoilMoist))
  if (add_WaterBalVars) then
    filevardef(cdf_file_out,"DelSoilMoist",typeof(DelSoilMoist),getvardims(DelSoilMoist))
    filevardef(cdf_file_out,"DelIntercept",typeof(DelIntercept),getvardims(DelIntercept))
  end if
  filevardef(cdf_file_out,"SoilTemp",typeof(SoilTemp),getvardims(SoilTemp))

  filevarattdef(cdf_file_out,"latitude",lat)
  filevarattdef(cdf_file_out,"longitude",lon)
  filevarattdef(cdf_file_out,"y",y)
  filevarattdef(cdf_file_out,"x",x)
  filevarattdef(cdf_file_out,"time",time)
  filevarattdef(cdf_file_out,"area",area)
  filevarattdef(cdf_file_out,"landfrac",landfrac)
  filevarattdef(cdf_file_out,"landmask",landmask)
  filevarattdef(cdf_file_out,"pftmask",pftmask)
  filevarattdef(cdf_file_out,"levgrnd",levgrnd)
  filevarattdef(cdf_file_out,"levdcmp",levdcmp)
  filevarattdef(cdf_file_out,"nbedrock",nbedrock)
  filevarattdef(cdf_file_out,"soil_layer",soil_layer)
  filevarattdef(cdf_file_out,"SWnet",SWnet)
  filevarattdef(cdf_file_out,"LWnet",LWnet)
  filevarattdef(cdf_file_out,"SWdown",SWdown)
  filevarattdef(cdf_file_out,"LWdown",LWdown)
  filevarattdef(cdf_file_out,"SWup",SWup)
  filevarattdef(cdf_file_out,"LWup",LWup)
  filevarattdef(cdf_file_out,"Tair",Tair)
  filevarattdef(cdf_file_out,"Qair",Qair)
  filevarattdef(cdf_file_out,"PSurf",PSurf)
  filevarattdef(cdf_file_out,"Wind",Wind)
  filevarattdef(cdf_file_out,"Qle",Qle)
  filevarattdef(cdf_file_out,"Qh",Qh)
  filevarattdef(cdf_file_out,"Qstor",Qstor)
  filevarattdef(cdf_file_out,"Qanth",Qanth)
  filevarattdef(cdf_file_out,"Qanth_Qh",Qanth_Qh)
  filevarattdef(cdf_file_out,"Qanth_Qle",Qanth_Qle)
  filevarattdef(cdf_file_out,"AvgSurfT",AvgSurfT)
  filevarattdef(cdf_file_out,"Albedo",Albedo)
  filevarattdef(cdf_file_out,"SAlbedo",SAlbedo)
  filevarattdef(cdf_file_out,"SnowFrac",SnowFrac)
  filevarattdef(cdf_file_out,"Rainf",Rainf)
  filevarattdef(cdf_file_out,"Snowf",Snowf)
  filevarattdef(cdf_file_out,"Evap",Evap)
  filevarattdef(cdf_file_out,"Qs",Qs)
  filevarattdef(cdf_file_out,"Qsb",Qsb)
  filevarattdef(cdf_file_out,"Qsm",Qsm)
  filevarattdef(cdf_file_out,"Qtau",Qtau)
  filevarattdef(cdf_file_out,"SWE",SWE)
  if (add_WaterBalVars) then
    filevarattdef(cdf_file_out,"DelSWE",DelSWE)
  end if
  filevarattdef(cdf_file_out,"RoofSurfT",RoofSurfT)
  filevarattdef(cdf_file_out,"WallSurfT",WallSurfT)
  filevarattdef(cdf_file_out,"RoadSurfT",RoadSurfT)
  filevarattdef(cdf_file_out,"TairCanyon",TairCanyon)
  filevarattdef(cdf_file_out,"TairBuilding",TairBuilding)
  filevarattdef(cdf_file_out,"SnowT",SnowT)
  filevarattdef(cdf_file_out,"BaresoilT",BaresoilT)
  filevarattdef(cdf_file_out,"SoilMoist",SoilMoist)
  if (add_WaterBalVars) then
    filevarattdef(cdf_file_out,"DelSoilMoist",DelSoilMoist)
    filevarattdef(cdf_file_out,"DelIntercept",DelIntercept)
  end if
  filevarattdef(cdf_file_out,"SoilTemp",SoilTemp)

  setfileoption(cdf_file_out,"DefineMode",False)

  cdf_file_out->latitude = lat
  cdf_file_out->longitude = lon
  cdf_file_out->y = y
  cdf_file_out->x = x
  cdf_file_out->time = time
  cdf_file_out->area = area
  cdf_file_out->landfrac = landfrac
  cdf_file_out->landmask = landmask
  cdf_file_out->pftmask = pftmask
  cdf_file_out->levgrnd = levgrnd
  cdf_file_out->levdcmp = levdcmp
  cdf_file_out->nbedrock = nbedrock
  cdf_file_out->soil_layer = (/soil_layer/)
  cdf_file_out->SWnet = (/SWnet/) 
  cdf_file_out->LWnet = (/LWnet/) 
  cdf_file_out->SWdown = (/SWdown/) 
  cdf_file_out->LWdown = (/LWdown/) 
  cdf_file_out->SWup = (/SWup/) 
  cdf_file_out->LWup = (/LWup/) 
  cdf_file_out->Tair = (/Tair/) 
  cdf_file_out->Qair = (/Qair/) 
  cdf_file_out->PSurf = (/PSurf/) 
  cdf_file_out->Wind = (/Wind/) 
  cdf_file_out->Qle = (/Qle/)
  cdf_file_out->Qh = (/Qh/)
  cdf_file_out->Qstor = (/Qstor/) 
  cdf_file_out->Qanth = (/Qanth/) 
  cdf_file_out->Qanth_Qh = (/Qanth_Qh/) 
  cdf_file_out->Qanth_Qle = (/Qanth_Qle/) 
  cdf_file_out->AvgSurfT = (/AvgSurfT/) 
  cdf_file_out->Albedo = (/Albedo/)
  cdf_file_out->SAlbedo = (/SAlbedo/)
  cdf_file_out->SnowFrac = (/SnowFrac/)
  cdf_file_out->Rainf = (/Rainf/)
  cdf_file_out->Snowf = (/Snowf/)
  cdf_file_out->Evap = (/Evap/)
  cdf_file_out->Qs = (/Qs/)
  cdf_file_out->Qsb = (/Qsb/)
  cdf_file_out->Qsm = (/Qsm/)
  cdf_file_out->Qtau = (/Qtau/)
  cdf_file_out->SWE = (/SWE/)
  if (add_WaterBalVars) then
    cdf_file_out->DelSWE = (/DelSWE/)
  end if
  cdf_file_out->RoofSurfT = (/RoofSurfT/)
  cdf_file_out->WallSurfT = (/WallSurfT/)
  cdf_file_out->RoadSurfT = (/RoadSurfT/)
  cdf_file_out->TairCanyon = (/TairCanyon/)
  cdf_file_out->TairBuilding = (/TairBuilding/)
  cdf_file_out->SnowT = (/SnowT/)
  cdf_file_out->BaresoilT = (/BaresoilT/)
  cdf_file_out->SoilMoist = (/SoilMoist/)
  if (add_WaterBalVars) then
    cdf_file_out->DelSoilMoist = (/DelSoilMoist/)
    cdf_file_out->DelIntercept = (/DelIntercept/)
  end if
  cdf_file_out->SoilTemp = (/SoilTemp/)

  end if

  if (do_csv) then 

    if (add_WaterBalVars) then
      header = (/seconds_since+","+"SWnet,LWnet,SWdown,LWdown,SWup,LWup,Tair,Qair,PSurf,Wind,Qle,Qh,Qstor,Qanth,Qanth_Qh,Qanth_Qle,AvgSurfT,Albedo,SAlbedo,SnowFrac,Rainf,SnowF,Evap,Qs,Qsb,Qsm,Qtau,SWE,DelSWE,RoofSurfT,WallSurfT,RoadSurfT,TairCanyon,TairBuilding,SnowT,BaresoilT,SoilMoist0,SoilMoist1,SoilMoist2,SoilMoist3,SoilMoist4,SoilMoist5,SoilMoist6,SoilMoist7,SoilMoist8,SoilMoist9,SoilMoist10,SoilMoist11,SoilMoist12,SoilMoist13,SoilMoist14,SoilMoist15,SoilMoist16,SoilMoist17,SoilMoist18,SoilMoist19,DelSoilMoist,DelIntercept,SoilTemp0,SoilTemp1,SoilTemp2,SoilTemp3,SoilTemp4,SoilTemp5,SoilTemp6,SoilTemp7,SoilTemp8,SoilTemp9,SoilTemp10,SoilTemp11,SoilTemp12,SoilTemp13,SoilTemp14,SoilTemp15,SoilTemp16,SoilTemp17,SoilTemp18,SoilTemp19"/)
      hlist = [/header/]
      write_table(f_out_dir+f_out_csv, "w", hlist, "%s") 
      alist = [/time,SWnet,LWnet,SWdown,LWdown,SWup,LWup,Tair,Qair,PSurf,Wind,Qle,Qh,Qstor,Qanth,Qanth_Qh,Qanth_Qle,AvgSurfT,Albedo,SAlbedo,SnowFrac,Rainf,Snowf,Evap,Qs,Qsb,Qsm,Qtau,SWE,DelSWE,RoofSurfT,WallSurfT,RoadSurfT,TairCanyon,TairBuilding,SnowT,BaresoilT,SoilMoist(:,0,0,0),SoilMoist(:,1,0,0),SoilMoist(:,2,0,0),SoilMoist(:,3,0,0),SoilMoist(:,4,0,0),SoilMoist(:,5,0,0),SoilMoist(:,6,0,0),SoilMoist(:,7,0,0),SoilMoist(:,8,0,0),SoilMoist(:,9,0,0),SoilMoist(:,10,0,0),SoilMoist(:,11,0,0),SoilMoist(:,12,0,0),SoilMoist(:,13,0,0),SoilMoist(:,14,0,0),SoilMoist(:,15,0,0),SoilMoist(:,16,0,0),SoilMoist(:,17,0,0),SoilMoist(:,18,0,0),SoilMoist(:,19,0,0),DelSoilMoist(:,0,0),DelIntercept(:,0,0),SoilTemp(:,0,0,0),SoilTemp(:,1,0,0),SoilTemp(:,2,0,0),SoilTemp(:,3,0,0),SoilTemp(:,4,0,0),SoilTemp(:,5,0,0),SoilTemp(:,6,0,0),SoilTemp(:,7,0,0),SoilTemp(:,8,0,0),SoilTemp(:,9,0,0),SoilTemp(:,10,0,0),SoilTemp(:,11,0,0),SoilTemp(:,12,0,0),SoilTemp(:,13,0,0),SoilTemp(:,14,0,0),SoilTemp(:,15,0,0),SoilTemp(:,16,0,0),SoilTemp(:,17,0,0),SoilTemp(:,18,0,0),SoilTemp(:,19,0,0)/]
      printVarSummary(alist)
      write_table(f_out_dir+f_out_csv, "a", alist, "%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e")
    else
      header = (/seconds_since+","+"SWnet,LWnet,SWdown,LWdown,SWup,LWup,Tair,Qair,PSurf,Wind,Qle,Qh,Qstor,Qanth,Qanth_Qh,Qanth_Qle,AvgSurfT,Albedo,SAlbedo,SnowFrac,Rainf,SnowF,Evap,Qs,Qsb,Qsm,Qtau,SWE,RoofSurfT,WallSurfT,RoadSurfT,TairCanyon,TairBuilding,SnowT,BaresoilT,SoilMoist0,SoilMoist1,SoilMoist2,SoilMoist3,SoilMoist4,SoilMoist5,SoilMoist6,SoilMoist7,SoilMoist8,SoilMoist9,SoilMoist10,SoilMoist11,SoilMoist12,SoilMoist13,SoilMoist14,SoilMoist15,SoilMoist16,SoilMoist17,SoilMoist18,SoilMoist19,SoilTemp0,SoilTemp1,SoilTemp2,SoilTemp3,SoilTemp4,SoilTemp5,SoilTemp6,SoilTemp7,SoilTemp8,SoilTemp9,SoilTemp10,SoilTemp11,SoilTemp12,SoilTemp13,SoilTemp14,SoilTemp15,SoilTemp16,SoilTemp17,SoilTemp18,SoilTemp19"/)
      hlist = [/header/]
      write_table(f_out_dir+f_out_csv, "w", hlist, "%s") 
      alist = [/time,SWnet,LWnet,SWdown,LWdown,SWup,LWup,Tair,Qair,PSurf,Wind,Qle,Qh,Qstor,Qanth,Qanth_Qh,Qanth_Qle,AvgSurfT,Albedo,SAlbedo,SnowFrac,Rainf,Snowf,Evap,Qs,Qsb,Qsm,Qtau,SWE,RoofSurfT,WallSurfT,RoadSurfT,TairCanyon,TairBuilding,SnowT,BaresoilT,SoilMoist(:,0,0,0),SoilMoist(:,1,0,0),SoilMoist(:,2,0,0),SoilMoist(:,3,0,0),SoilMoist(:,4,0,0),SoilMoist(:,5,0,0),SoilMoist(:,6,0,0),SoilMoist(:,7,0,0),SoilMoist(:,8,0,0),SoilMoist(:,9,0,0),SoilMoist(:,10,0,0),SoilMoist(:,11,0,0),SoilMoist(:,12,0,0),SoilMoist(:,13,0,0),SoilMoist(:,14,0,0),SoilMoist(:,15,0,0),SoilMoist(:,16,0,0),SoilMoist(:,17,0,0),SoilMoist(:,18,0,0),SoilMoist(:,19,0,0),SoilTemp(:,0,0,0),SoilTemp(:,1,0,0),SoilTemp(:,2,0,0),SoilTemp(:,3,0,0),SoilTemp(:,4,0,0),SoilTemp(:,5,0,0),SoilTemp(:,6,0,0),SoilTemp(:,7,0,0),SoilTemp(:,8,0,0),SoilTemp(:,9,0,0),SoilTemp(:,10,0,0),SoilTemp(:,11,0,0),SoilTemp(:,12,0,0),SoilTemp(:,13,0,0),SoilTemp(:,14,0,0),SoilTemp(:,15,0,0),SoilTemp(:,16,0,0),SoilTemp(:,17,0,0),SoilTemp(:,18,0,0),SoilTemp(:,19,0,0)/]
      printVarSummary(alist)
      write_table(f_out_dir+f_out_csv, "a", alist, "%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e,%11.5e")
    end if
  end if

  ; Energy Balance
  bal = SWnet + LWnet + Qanth - Qle - Qh - Qstor
  print((/"Max energy imbalance: "+max((bal))/))
  print((/"Min energy imbalance: "+min((bal))/))
  swbal = SWnet - (SWdown - SWup)
  print((/swbal(0:1,0,0)/))
  print((/SWnet(0:1,0,0)/))
  print((/SWdown(0:1,0,0)/))
  print((/SWup(0:1,0,0)/))
  print((/Albedo(0:1,0,0)/))
  print((/"Max SW imbalance: "+max((swbal))/))
  print((/"Min SW imbalance: "+min((swbal))/))
  lwbal = -LWnet - (LWup - LWdown)
  print((/"Max LW imbalance: "+max((lwbal))/))
  print((/"Min LW imbalance: "+min((lwbal))/))
  print((/"Max Albedo: "+max((Albedo))/))
  print((/"Min Albedo: "+min((Albedo))/))
  swup_derived = SWdown - SWnet
  print((/"Max swup_derived: "+max((swup_derived))/))
  print((/"Min swup_derived: "+min((swup_derived))/))

  ; Water Balance
; balw = (Rainf + Snowf - Evap - Qs - Qsb)  ; mm/s
; print((/"Max water imbalance (all time steps, mm/s): "+max((balw))/))
; print((/"Min water imbalance (all time steps, mm/s): "+min((balw))/))

  if (lu_soilmoist) then

  balw = new(ntim,float)
  do t = 0,ntim-1
    ; We don't have SoilMoist at the beginning of the model run, so just check balance under the assumption that
    ; soil moisture didn't change during the first time step.  In other words, just ignore the balance in the first
    ; time step
    if (t .eq. 0) then
      balw(t) = Rainf(t,0,0) + Snowf(t,0,0) - Evap(t,0,0) - Qs(t,0,0) - Qsb(t,0,0) - (dim_sum_n(SoilMoist(t,:,0,0),0)-dim_sum_n(SoilMoist(t,:,0,0),0))/(time(1)-time(0))
    else
      balw(t) = Rainf(t,0,0) + Snowf(t,0,0) - Evap(t,0,0) - Qs(t,0,0) - Qsb(t,0,0) - (dim_sum_n(SoilMoist(t,:,0,0),0)-dim_sum_n(SoilMoist(t-1,:,0,0),0))/(time(1)-time(0))
    end if
    print((/"Rainf: "+Rainf(t,0,0)/))
    print((/"Snowf: "+Snowf(t,0,0)/))
    print((/"Evap: "+Evap(t,0,0)/))
    print((/"Qs: "+Qs(t,0,0)/))
    print((/"Qsb: "+Qsb(t,0,0)/))
    if (t .eq. 0) then
      print((/"SoilMoistDelta: "+((dim_sum_n(SoilMoist(t,:,0,0),0)-dim_sum_n(SoilMoist(t,:,0,0),0))/(time(1)-time(0)))/))
    else
      print((/"SoilMoistDelta: "+((dim_sum_n(SoilMoist(t,:,0,0),0)-dim_sum_n(SoilMoist(t-1,:,0,0),0))/(time(1)-time(0)))/))
      print((/"SoilMoist(t): "+dim_sum_n(SoilMoist(t,:,0,0),0)/))
      print((/"SoilMoist(t-1): "+dim_sum_n(SoilMoist(t-1,:,0,0),0)/))
    end if
    print((/"deltaflux: "+(Rainf(t,0,0) + Snowf(t,0,0) - Evap(t,0,0) - Qs(t,0,0) - Qsb(t,0,0))/))
    print((/"Water imbalance: "+balw(t)/))
  end do

  print((/"Mean water imbalance: "+avg(balw)/))

  end if

  if (add_WaterBalVars) then

  balw = new(ntim,float)
  do t = 0,ntim-1
    balw(t) = Rainf(t,0,0) + Snowf(t,0,0) - Evap(t,0,0) - Qs(t,0,0) - Qsb(t,0,0) - (DelSoilMoist(t,0,0)/(time(1)-time(0)) + DelSWE(t,0,0)/(time(1)-time(0)) + DelIntercept(t,0,0)/(time(1)-time(0)))
    print((/"Rainf: "+Rainf(t,0,0)/))
    print((/"Snowf: "+Snowf(t,0,0)/))
    print((/"Evap: "+Evap(t,0,0)/))
    print((/"Qs: "+Qs(t,0,0)/))
    print((/"Qsb: "+Qsb(t,0,0)/))
    print((/"DelSoilMoist: "+DelSoilMoist(t,0,0)/(time(1)-time(0))/))
    print((/"DelSWE: "+DelSWE(t,0,0)/(time(1)-time(0))/))
    print((/"DelIntercept: "+DelIntercept(t,0,0)/(time(1)-time(0))/))
    print((/"deltaflux: "+(Rainf(t,0,0) + Snowf(t,0,0) - Evap(t,0,0) - Qs(t,0,0) - Qsb(t,0,0))/))
    print((/"Water imbalance: "+balw(t)/))
  end do

  ; Latent heat of vaporization
; lamda_vap = 2.501e6     ; Our value at 0C
  lamda_vap = 2.264705e6  ; Harro's value at 20C which is wrong
  lamda_vap = 2.452718e6  ; Harro's value at 20C from his equation which is correct
  ; For this site we want the sums for 2012-01-01 to the end of the simulation (1826 days)
  if (Psites(sel) .eq. "CA-Sunset") then
    sttime = ind(mcdate .eq. 20120101 .and. mcsec .eq. 1800)
    print(sttime)
    Evap_analy_period = Evap(sttime:,0,0)
    inddew = ind(Evap_analy_period(:) .lt. 0.)
    dew = Evap_analy_period(inddew)
    print((/"Mean water imbalance: "+avg(balw)/))
    print((/"Sum(Rainf): "+dim_sum(Rainf(sttime:,0,0)*(time(1)-time(0)))/))
    print((/"Sum(Snowf): "+dim_sum(Snowf(sttime:,0,0)*(time(1)-time(0)))/))
    print((/"Sum(Evap): "+dim_sum(Evap(sttime:,0,0)*(time(1)-time(0)))/))
    print((/"Sum(Qle/lamda_vap): "+dim_sum(Qle(sttime:,0,0)/lamda_vap)*(time(1)-time(0))/))
    print((/"Sum(Dew): "+dim_sum(dew)*(time(1)-time(0))/)) 
    print((/"Sum(Qs): "+dim_sum(Qs(sttime:,0,0)*(time(1)-time(0)))/))
    print((/"Sum(Qsb): "+dim_sum(Qsb(sttime:,0,0)*(time(1)-time(0)))/))
    print((/"Sum(DelSWE): "+dim_sum(DelSWE(sttime:,0,0))/))
    print((/"Sum(DelIntercept): "+dim_sum(DelIntercept(sttime:,0,0))/))
    print((/"Sum(DelSoilMoist): "+dim_sum(DelSoilMoist(sttime:,0,0))/))
    print((/"Soil Storage: "+(dim_sum(SoilMoist(ntim-1,:,0,0))-dim_sum(SoilMoist(sttime,:,0,0)))/))
  end if

  end if

  print ("=========================================")
  print ("Finish Time: "+systemfunc("date") )
  print ("=========================================")

  end
