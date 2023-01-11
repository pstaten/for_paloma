;*************************************************
; calc_div_aqua.ncl
;*************************************************
begin
;*************************************************
; open file and read in data
;*************************************************
  f1 = addfile (infile, "r")
  if isfilevar(f1, "U") then
    u    = f1->U
    v    = f1->V
    lat  = f1->lat
    lon  = f1->lon
    xofyear = f1->xofyear
    psfc = f1->ps
    psfc0 = f1->ps/100.
    pfull = f1->pfull
    if (.not.isatt(u,"_FillValue")) then
      u@_FillValue = -1e20
      v@_FillValue = u@_FillValue
    end if
  else
    u    = f1->ucomp
    v    = f1->vcomp
    lat  = f1->lat
    lon  = f1->lon
    xofyear = f1->xofyear
    psfc = f1->ps
    pfull = f1->pfull
    u@_FillValue = -1e20
    v@_FillValue = u@_FillValue
    u@units = "ms-1"
    v@units = u@units
    replace_ieeenan(u, u@_FillValue, 0)
    replace_ieeenan(v, v@_FillValue, 0)
    ; fill in missing values via interpolation
    ; for CMIP5 data
    cyclic=True
    guess_Type=1
    nscan = 10
    epsx = 0.1
    relc = 0.5
    opt = 0
    poisson_grid_fill(u, cyclic, guess_Type, nscan, epsx, relc, opt)
    poisson_grid_fill(v, cyclic, guess_Type, nscan, epsx, relc, opt)
  end if
  nlat = dimsizes(lat)
  nlon = dimsizes(lon)
  nlev = dimsizes(pfull)
  ntim = dimsizes(xofyear)

;*************************************************
; calculate divergence on a gaussian grid
;*************************************************
  printVarSummary(u)
  div_all = uv2dv_cfd(u,v,u&lat,u&lon,3)                ; u,v ==> divergence
  copy_VarCoords(u,div_all)

;*************************************************
; calculate divergent wind component 
;*************************************************    
  print ("NCL calculating divergent wind component")
  ud    = new ( dimsizes(u), typeof(u), u@_FillValue )
  vd    = new ( dimsizes(v), typeof(v), v@_FillValue )
  dv2uvg(div_all,ud,vd) ; div_all  ==> divergent  wind components
  copy_VarCoords(u, ud ) 
  copy_VarCoords(u, vd ) 
  ud@long_name  = "Zonal Divergent Wind"
  ud@units      = u@units
  vd@long_name  = "Meridional Divergent Wind"
  vd@units      = v@units

;*************************************************
; calculate meridional and zonal divergence on a gaussian grid
;*************************************************
  div_merid = uv2dv_cfd(u*0,vd,u&lat,u&lon,3)                ; u,v ==> divergence
  copy_VarCoords(div_all,div_merid)
  div_zonal = uv2dv_cfd(ud,v*0,u&lat,u&lon,3)                ; u,v ==> divergence
  copy_VarCoords(div_all,div_zonal)

;************************************************************
; Compute delta omega at each hybrid level via kinematic method
; .  see eqn 8.20 in Wallace and Hobbs
;************************************************************

  print ("NCL calculating omega via kinematic method")
  dpi = dpres_plevel(pfull, psfc, 0, 0)

  omega_all = div_all
  omega_all = div_all*dpi
  omega_merid = div_merid
  omega_merid = div_merid*dpi
  omega_zonal = div_zonal
  omega_zonal = div_zonal*dpi
  dimo  = dimsizes(omega_all)
  klev  = dimo(1)
  ;printVarSummary(div_all)

;************************************************************
; integrate [sum] the delta-omega [bottom->top]
; .  store back into the omega variable array
; .  assume that omega at sfc of atmosphere is 0.0
; .  the funky omega subscript is because omega goes top->bottom
; .  while temp goes from sfc to top [after reordering]
;************************************************************
  print ("NCL calculating omega, step 2")

  temp_all = omega_all(xofyear|:,lat|:,lon|:,pfull|::-1)      ; reshape array for function
  temp_merid = omega_merid(xofyear|:,lat|:,lon|:,pfull|::-1)      ; reshape array for function
  temp_zonal = omega_zonal(xofyear|:,lat|:,lon|:,pfull|::-1)      ; reshape array for function
  dimt  = dimsizes(temp_all)
  klev  = dimt(3)
  oBot  = 0.0                                     ; explicit
  
  do kl=0,klev-1                             
     omega_all(:,klev-1-kl,:,:) = oBot - dim_sum(temp_all(:,:,:,0:kl))   
     omega_merid(:,klev-1-kl,:,:) = oBot - dim_sum(temp_merid(:,:,:,0:kl))   
     omega_zonal(:,klev-1-kl,:,:) = oBot - dim_sum(temp_zonal(:,:,:,0:kl))   
  end do

  omega_all@units     = "Pa/s"
  omega_merid@units     = "Pa/s"
  omega_zonal@units     = "Pa/s"
  omega_all@long_name = "OMEGA_ALL"
  omega_all@long_name = "OMEGA_ALL"
  omega_merid@long_name = "OMEGA_MERID"
  omega_zonal@long_name = "OMEGA_ZONAL"
  omega_all@method    = "derived via kinematic method: eqn 8.20 Wallace and Hobbs"
  omega_all@method    = "derived via kinematic method: eqn 8.20 Wallace and Hobbs"
  omega_merid@method    = "derived via kinematic method: eqn 8.20 Wallace and Hobbs"
  omega_zonal@method    = "derived via kinematic method: eqn 8.20 Wallace and Hobbs"



;************************************************************
; Define the netCDF file explicitly, rather than just saving
; the variable and then hoping that NCL will call the
; appropriate subroutines.
; Upgrading from 6.4 (on Karst) to 6.6.2 (on Carbonate)
; Caused a problem with that approach on 8/29/2019.
;************************************************************

print ("NCL saving output to netcdf")
setfileoption("nc","format","netcdf4")
system("/bin/rm -f " + divfile)   ; remove any pre-existing file
ncdf = addfile(divfile ,"c")  ; open output netCDF file
setfileoption(ncdf,"DefineMode",True)

; create global attributes of the file (optional)
fAtt               = True            ; assign file attributes
fAtt@title         = "Divergent wind field and kinematic vertical velocities"
fAtt@source_file   = infile
fAtt@Conventions   = "None"
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )            ; copy file attributes

dimNames = (/"xofyear", "pfull", "lat", "lon"/)
dimSizes = (/ntim, nlev, nlat, nlon /)
dimUnlim = (/ True, False, False, False/)
filedimdef(ncdf,dimNames, dimSizes, dimUnlim)

filevardef(ncdf, "xofyear", typeof(xofyear), getvardims(xofyear))
filevarattdef(ncdf,"xofyear", xofyear) ; copy time attributes
ncdf->xofyear   = (/xofyear/)

filevardef(ncdf, "pfull", typeof(pfull), getvardims(pfull) )
filevarattdef(ncdf,"pfull", pfull) ; copy lev attributes
ncdf->pfull    = (/pfull/)

filevardef(ncdf, "lat", typeof(lat), getvardims(lat))
filevarattdef(ncdf,"lat", lat) ; copy lat attributes
ncdf->lat    = (/lat/)

filevardef(ncdf, "lon", typeof(lon), getvardims(lon))
filevarattdef(ncdf,"lon", lon) ; copy lon attributes
ncdf->lon    = (/lon/)

filevardef(ncdf, "div_all", typeof(div_all), getvardims(div_all))
filevarattdef(ncdf,"div_all", div_all) ; copy T attributes
ncdf->div_all      = (/div_all/)

filevardef(ncdf, "div_merid", typeof(div_merid), getvardims(div_merid))
filevarattdef(ncdf,"div_merid", div_merid) ; copy T attributes
ncdf->div_merid      = (/div_merid/)

filevardef(ncdf, "div_zonal", typeof(div_zonal), getvardims(div_zonal))
filevarattdef(ncdf,"div_zonal", div_zonal) ; copy T attributes
ncdf->div_zonal      = (/div_zonal/)

filevardef(ncdf, "ud", typeof(ud), getvardims(ud))
filevarattdef(ncdf,"ud", ud) ; copy T attributes
ncdf->ud      = (/ud/)

filevardef(ncdf, "vd", typeof(vd), getvardims(vd))
filevarattdef(ncdf,"vd", vd) ; copy T attributes
ncdf->vd      = (/vd/)

filevardef(ncdf, "omega_all", typeof(omega_all), getvardims(omega_all))
filevarattdef(ncdf,"omega_all", omega_all) ; copy T attributes
ncdf->omega_all      = (/omega_all/)

filevardef(ncdf, "omega_merid", typeof(omega_merid), getvardims(omega_merid))
filevarattdef(ncdf,"omega_merid", omega_merid) ; copy T attributes
ncdf->omega_merid      = (/omega_merid/)

filevardef(ncdf, "omega_zonal", typeof(omega_zonal), getvardims(omega_zonal))
filevarattdef(ncdf,"omega_zonal", omega_zonal) ; copy T attributes
ncdf->omega_zonal      = (/omega_zonal/)
end

