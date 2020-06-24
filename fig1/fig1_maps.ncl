
; ------------------------------------------------------------------------------
; fig1_maps.m
;
; Description: script to create the maps of Fig. 1 in Deser et al., 2020,
; https://doi.org/10.1038/s41558-020-0731-2
;
; The full data is too large to share on github, but is available via:
; - Large Ensembles: http://www.cesm.ucar.edu/projects/community-projects/MMLEA/
;
; Required data:
; - the MMLEA archive
;   (currently, the folder structure is as on UCAR's Cheyenne supercomputer,
;   where the archive is stored --> path is set under "pathin")
; - Observations
;   (here, Berkeley Earth for temperature: http://berkeleyearth.lbl.gov/auto/Global/Gridded/Land_and_Ocean_LatLong1.nc)
;
; Author: Flavio Lehner, May 2020, flehner@ucar.edu
;
; ------------------------------------------------------------------------------

begin
  vars        = (/"tas","pr","psl"/)    ; tas or pr
  seasons     = (/"DJF","JJA","annual"/)    ; DJF JJA annual
  comp        = "Amon"
  domain      = "NA"    ; global or NA
  pathin      = "/project/yampa03/clivar_wg_le/"
  pathin_obs  = "/project/cas/flehner/observations/temperature/best/"
  pathout     = "/project/cas/flehner/"
  start_trend = 1951
  ende_trend  = 2010
  trend_time  = fspan(start_trend,ende_trend,ende_trend-start_trend+1)
  clivar_tag     = (/"cesm","canesm2","csiro_mk36","gfdl_cm3",\
                     "gfdl_esm2m","ec_earth","mpi","ipsl_cm5a"/)
  model_name     = (/"CESM1-CAM5","CanESM2","CSIRO-Mk3-6-0","GFDL-CM3",\
                     "GFDL-ESM2M","EC-EARTH","MPI-ESM","IPSL-CM5A-LR"/)
  model_nametags2= (/"CESM1","CanESM2","CSIRO-Mk3-6-0","GFDL-CM3",\
                    "GFDL-ESM2M","EC-EARTH","MPI","IPSL-CM5A-LR"/)
  model_nametags = (/"NCAR CESM1 (40)","CCCma CanESM2 (50)","CSIRO Mk36 (30)","GFDL CM3 (20)",\
                     "GFDL ESM2M (30)","EC-EARTH (16)","MPI ESM-MR (100)","IPSL CM5A-LR (30)"/)
  start0         = (/1920,1950,1850,1920,1950,1860,1850,1940/)
  ende0          = (/2100,2100,2100,2100,2100,2100,2099,2100/)
  ensmem0        = (/40,50,30,20,30,16,100,30/)

  ; -- recalc trends and SN?
  recalc1     = "no" ; cesm_lens
  recalc2     = "no" ; canesm2_lens
  recalc3     = "no" ; csiro_mk36_lens
  recalc4     = "no" ; gfdl_cm3_lens
  recalc5     = "no" ; gfdl_esm2m_lens
  recalc6     = "no" ; ec_earth_lens
  recalc7     = "no" ; mpi_lens

  ; -- recalc trends in OLENS?
  recalc10    = "no"

  ; -- recalc ensemble mean?
  recalc_em   = "no"

  ; -- recalc climatology and trend significance?
  recalc_clim = "no"

  ; -- regrid OLENS to MPI and add MPI EM
  recalc_olens_mpi = "no"

  ; -- choose bookends based on (1) trend or (2) pattern correlation?
  be = 1

  ; -- centered or uncentered pattern correlation?
  cent = 0; 0=centered, 1=uncentered
;-----------------------------------------------------------

do v = 0,0 ; full: 0,2
  var = vars(v)
  do s = 2,2 ; full: 0,2
    seas = seasons(s)
    print_clock(""+var)
    print_clock(""+seas)

  ; -- cesm_lens --
  print_clock("cesm_lens")
  if recalc1 .eq. "yes"
    start   = 1920
    ende    = 2100
    ensmem  = 40
    a = addfile(pathin+"cesm_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(0)+"_historical_rcp85_r2i1p1_192001-210012.nc","r")
    dummy   = a->$var$
    dims    = dimsizes(dummy)
    ntime   = dims(0)
    nlat    = dims(1)
    nlon    = dims(2)
    cesm_trends  = new((/ensmem,nlat,nlon/),float)
    do e = 1,ensmem
      print_clock("  ensmem = "+e+"/"+ensmem)
      if e .eq. 1
        a = addfile(pathin+"cesm_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(0)+"_historical_rcp85_r"+e+"i1p1_185001-210012.nc","r")
        tmp0     = a->$var$
        tmp00    = tmp0(840::,:,:)
        delete(tmp0)
        tmp0     = tmp00
        delete(tmp00)
      else
        a = addfile(pathin+"cesm_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(0)+"_historical_rcp85_r"+e+"i1p1_192001-210012.nc","r")
        tmp0     = a->$var$
      end if
      if seas .eq. "annual"
        tmp1 = month_to_annual(tmp0,1)
      else
        tmp1 = month_to_season(tmp0,seas)
      end if
      cesm_trends(e-1,:,:)  = regCoef_n(trend_time, tmp1(start_trend-start:ende_trend-start,:,:), 0, 0) * (ende_trend-start_trend+1)
    end do
    copy_VarCoords(dummy(0,:,:),cesm_trends(0,:,:))
    cesmEM = dim_avg_n_Wrap(cesm_trends,0)
    cesmSD = dim_stddev_n_Wrap(cesm_trends,0)
    cesmSN = (/ abs(cesmEM)/cesmSD /)
    copy_VarCoords(dummy(0,:,:),cesmSN)
    ; -- save data --
    print("save data...")
    ofile1 = pathout+"cesm_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(0)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
    system("mkdir -p "+pathout+"cesm_lens/"+comp+"/"+var)
    system("rm -rf "+ofile1)
    o = addfile(ofile1,"c")
    o->cesm_trends         = cesm_trends
    o->cesmEM         = cesmEM
    o->cesmSD         = cesmSD
    o->cesmSN         = cesmSN
    delete([/dummy,dims,ntime,nlat,nlon,tmp0,tmp1/])
  else
    ; -- load data --
    print("load data...")
    ifile1 = pathout+"cesm_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(0)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
    i = addfile(ifile1,"r")
    cesm_trends    = i->cesm_trends
    cesmEM         = i->cesmEM
    cesmSD         = i->cesmSD
    cesmSN         = i->cesmSN
  end if


  ; -- canesm2_lens --
  print_clock("canesm2_lens")
  if recalc2 .eq. "yes"
    start   = 1950
    ende    = 2100
    ensmem  = 50
    a = addfile(pathin+"canesm2_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(1)+"_historical_rcp85_r2i1p1_195001-210012.nc","r")
    dummy   = a->$var$
    dims    = dimsizes(dummy)
    ntime   = dims(0)
    nlat    = dims(1)
    nlon    = dims(2)
    canesm2_trends  = new((/ensmem,nlat,nlon/),float)
    do e = 1,ensmem
      print_clock("  ensmem = "+e+"/"+ensmem)
      a = addfile(pathin+"canesm2_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(1)+"_historical_rcp85_r"+e+"i1p1_195001-210012.nc","r")
      tmp0     = a->$var$
      if seas .eq. "annual"
        tmp1 = month_to_annual(tmp0,1)
      else
        tmp1 = month_to_season(tmp0,seas)
      end if
      canesm2_trends(e-1,:,:)  = regCoef_n(trend_time, tmp1(start_trend-start:ende_trend-start,:,:), 0, 0) * (ende_trend-start_trend+1)
    end do
    copy_VarCoords(dummy(0,:,:),canesm2_trends(0,:,:))
    canesm2EM = dim_avg_n_Wrap(canesm2_trends,0)
    canesm2SD = dim_stddev_n_Wrap(canesm2_trends,0)
    canesm2SN = (/ abs(canesm2EM)/canesm2SD /)
    copy_VarCoords(dummy(0,:,:),canesm2SN)
    ; -- save data --
    print("save data...")
    ofile1 = pathout+"canesm2_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(1)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
    system("mkdir -p "+pathout+"canesm2_lens/"+comp+"/"+var)
    system("rm -rf "+ofile1)
    o = addfile(ofile1,"c")
    o->canesm2_trends         = canesm2_trends
    o->canesm2EM         = canesm2EM
    o->canesm2SD         = canesm2SD
    o->canesm2SN         = canesm2SN
    delete([/dummy,dims,ntime,nlat,nlon,tmp0,tmp1/])
  else
    ; -- load data --
    print("load data...")
    ifile1 = pathout+"canesm2_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(1)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
    i = addfile(ifile1,"r")
    canesm2_trends    = i->canesm2_trends
    canesm2EM         = i->canesm2EM
    canesm2SD         = i->canesm2SD
    canesm2SN         = i->canesm2SN
  end if


  ; -- csiro_mk36_lens --
  print_clock("csiro_mk36_lens")
  if recalc3 .eq. "yes"
    start   = 1850
    ende    = 2100
    ensmem  = 30
    a = addfile(pathin+"csiro_mk36_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(2)+"_historical_rcp85_r2i1p1_185001-210012.nc","r")
    dummy   = a->$var$
    dims    = dimsizes(dummy)
    ntime   = dims(0)
    nlat    = dims(1)
    nlon    = dims(2)
    csiro_mk36_trends  = new((/ensmem,nlat,nlon/),float)
    do e = 1,ensmem
      print_clock("  ensmem = "+e+"/"+ensmem)
      a = addfile(pathin+"csiro_mk36_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(2)+"_historical_rcp85_r"+e+"i1p1_185001-210012.nc","r")
      tmp0     = a->$var$
      if seas .eq. "annual"
        tmp1 = month_to_annual(tmp0,1)
      else
        tmp1 = month_to_season(tmp0,seas)
      end if
      csiro_mk36_trends(e-1,:,:)  = regCoef_n(trend_time, tmp1(start_trend-start:ende_trend-start,:,:), 0, 0) * (ende_trend-start_trend+1)
    end do
    copy_VarCoords(dummy(0,:,:),csiro_mk36_trends(0,:,:))
    csiro_mk36EM = dim_avg_n_Wrap(csiro_mk36_trends,0)
    csiro_mk36SD = dim_stddev_n_Wrap(csiro_mk36_trends,0)
    csiro_mk36SN = (/ abs(csiro_mk36EM)/csiro_mk36SD /)
    copy_VarCoords(dummy(0,:,:),csiro_mk36SN)
    ; -- save data --
    print("save data...")
    ofile1 = pathout+"csiro_mk36_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(2)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
    system("mkdir -p "+pathout+"csiro_mk36_lens/"+comp+"/"+var)
    system("rm -rf "+ofile1)
    o = addfile(ofile1,"c")
    o->csiro_mk36_trends         = csiro_mk36_trends
    o->csiro_mk36EM         = csiro_mk36EM
    o->csiro_mk36SD         = csiro_mk36SD
    o->csiro_mk36SN         = csiro_mk36SN
    delete([/dummy,dims,ntime,nlat,nlon,tmp0,tmp1/])
  else
    ; -- load data --
    print("load data...")
    ifile1 = pathout+"csiro_mk36_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(2)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
    i = addfile(ifile1,"r")
    csiro_mk36_trends    = i->csiro_mk36_trends
    csiro_mk36EM         = i->csiro_mk36EM
    csiro_mk36SD         = i->csiro_mk36SD
    csiro_mk36SN         = i->csiro_mk36SN
  end if


  ; -- gfdl_cm3_lens --
  print_clock("gfdl_cm3_lens")
  if recalc4 .eq. "yes"
    start   = 1920
    ende    = 2100
    ensmem  = 20
    a = addfile(pathin+"gfdl_cm3_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(3)+"_historical_rcp85_r2i1p1_192001-210012.nc","r")
    dummy   = a->$var$
    dims    = dimsizes(dummy)
    ntime   = dims(0)
    nlat    = dims(1)
    nlon    = dims(2)
    gfdl_cm3_trends  = new((/ensmem,nlat,nlon/),float)
    do e = 1,ensmem
      print_clock("  ensmem = "+e+"/"+ensmem)
      a = addfile(pathin+"gfdl_cm3_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(3)+"_historical_rcp85_r"+e+"i1p1_192001-210012.nc","r")
      tmp0     = a->$var$
      if seas .eq. "annual"
        tmp1 = month_to_annual(tmp0,1)
      else
        tmp1 = month_to_season(tmp0,seas)
      end if
      gfdl_cm3_trends(e-1,:,:)  = regCoef_n(trend_time, tmp1(start_trend-start:ende_trend-start,:,:), 0, 0) * (ende_trend-start_trend+1)
    end do
    copy_VarCoords(dummy(0,:,:),gfdl_cm3_trends(0,:,:))
    gfdl_cm3EM = dim_avg_n_Wrap(gfdl_cm3_trends,0)
    gfdl_cm3SD = dim_stddev_n_Wrap(gfdl_cm3_trends,0)
    gfdl_cm3SN = (/ abs(gfdl_cm3EM)/gfdl_cm3SD /)
    copy_VarCoords(dummy(0,:,:),gfdl_cm3SN)
    ; -- save data --
    print("save data...")
    ofile1 = pathout+"gfdl_cm3_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(3)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
    system("mkdir -p "+pathout+"gfdl_cm3_lens/"+comp+"/"+var)
    system("rm -rf "+ofile1)
    o = addfile(ofile1,"c")
    o->gfdl_cm3_trends         = gfdl_cm3_trends
    o->gfdl_cm3EM         = gfdl_cm3EM
    o->gfdl_cm3SD         = gfdl_cm3SD
    o->gfdl_cm3SN         = gfdl_cm3SN
    delete([/dummy,dims,ntime,nlat,nlon,tmp0,tmp1/])
  else
    ; -- load data --
    print("load data...")
    ifile1 = pathout+"gfdl_cm3_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(3)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
    i = addfile(ifile1,"r")
    gfdl_cm3_trends    = i->gfdl_cm3_trends
    gfdl_cm3EM         = i->gfdl_cm3EM
    gfdl_cm3SD         = i->gfdl_cm3SD
    gfdl_cm3SN         = i->gfdl_cm3SN
  end if


  ; -- gfdl_esm2m_lens --
  print_clock("gfdl_esm2m_lens")
  if recalc5 .eq. "yes"
    start   = 1950
    ende    = 2100
    ensmem  = 30
    a = addfile(pathin+"gfdl_esm2m_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(4)+"_historical_rcp85_r2i1p1_195001-210012.nc","r")
    dummy   = a->$var$
    dims    = dimsizes(dummy)
    ntime   = dims(0)
    nlat    = dims(1)
    nlon    = dims(2)
    gfdl_esm2m_trends  = new((/ensmem,nlat,nlon/),float)
    do e = 1,ensmem
      print_clock("  ensmem = "+e+"/"+ensmem)
      a = addfile(pathin+"gfdl_esm2m_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(4)+"_historical_rcp85_r"+e+"i1p1_195001-210012.nc","r")
      tmp0     = a->$var$
      if seas .eq. "annual"
        tmp1 = month_to_annual(tmp0,1)
      else
        tmp1 = month_to_season(tmp0,seas)
      end if
      gfdl_esm2m_trends(e-1,:,:)  = regCoef_n(trend_time, tmp1(start_trend-start:ende_trend-start,:,:), 0, 0) * (ende_trend-start_trend+1)
    end do
    copy_VarCoords(dummy(0,:,:),gfdl_esm2m_trends(0,:,:))
    gfdl_esm2mEM = dim_avg_n_Wrap(gfdl_esm2m_trends,0)
    gfdl_esm2mSD = dim_stddev_n_Wrap(gfdl_esm2m_trends,0)
    gfdl_esm2mSN = (/ abs(gfdl_esm2mEM)/gfdl_esm2mSD /)
    copy_VarCoords(dummy(0,:,:),gfdl_esm2mSN)
    ; -- save data --
    print("save data...")
    ofile1 = pathout+"gfdl_esm2m_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(4)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
    system("mkdir -p "+pathout+"gfdl_esm2m_lens/"+comp+"/"+var)
    system("rm -rf "+ofile1)
    o = addfile(ofile1,"c")
    o->gfdl_esm2m_trends         = gfdl_esm2m_trends
    o->gfdl_esm2mEM         = gfdl_esm2mEM
    o->gfdl_esm2mSD         = gfdl_esm2mSD
    o->gfdl_esm2mSN         = gfdl_esm2mSN
    delete([/dummy,dims,ntime,nlat,nlon,tmp0,tmp1/])
  else
    ; -- load data --
    print("load data...")
    ifile1 = pathout+"gfdl_esm2m_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(4)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
    i = addfile(ifile1,"r")
    gfdl_esm2m_trends    = i->gfdl_esm2m_trends
    gfdl_esm2mEM         = i->gfdl_esm2mEM
    gfdl_esm2mSD         = i->gfdl_esm2mSD
    gfdl_esm2mSN         = i->gfdl_esm2mSN
  end if


  ; -- ec_earth_lens --
  print_clock("ec_earth_lens")
  if recalc6 .eq. "yes"
    start   = 1860
    ende    = 2100
    ensmem  = 16
    a = addfile(pathin+"ec_earth_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(5)+"_historical_rcp85_r2i1p1_186001-210012.nc","r")
    dummy   = a->$var$
    dims    = dimsizes(dummy)
    ntime   = dims(0)
    nlat    = dims(1)
    nlon    = dims(2)
    ec_earth_trends  = new((/ensmem,nlat,nlon/),float)
    do e = 1,ensmem
      print_clock("  ensmem = "+e+"/"+ensmem)
      a = addfile(pathin+"ec_earth_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(5)+"_historical_rcp85_r"+e+"i1p1_186001-210012.nc","r")
      tmp0     = a->$var$
      if seas .eq. "annual"
        tmp1 = month_to_annual(tmp0,1)
      else
        tmp1 = month_to_season(tmp0,seas)
      end if
      ec_earth_trends(e-1,:,:)  = regCoef_n(trend_time, tmp1(start_trend-start:ende_trend-start,:,:), 0, 0) * (ende_trend-start_trend+1)
    end do
    copy_VarCoords(dummy(0,:,:),ec_earth_trends(0,:,:))
    ec_earthEM = dim_avg_n_Wrap(ec_earth_trends,0)
    ec_earthSD = dim_stddev_n_Wrap(ec_earth_trends,0)
    ec_earthSN = (/ abs(ec_earthEM)/ec_earthSD /)
    copy_VarCoords(dummy(0,:,:),ec_earthSN)
    ; -- save data --
    print("save data...")
    ofile1 = pathout+"ec_earth_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(5)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
    system("mkdir -p "+pathout+"ec_earth_lens/"+comp+"/"+var)
    system("rm -rf "+ofile1)
    o = addfile(ofile1,"c")
    o->ec_earth_trends         = ec_earth_trends
    o->ec_earthEM         = ec_earthEM
    o->ec_earthSD         = ec_earthSD
    o->ec_earthSN         = ec_earthSN
    delete([/dummy,dims,ntime,nlat,nlon,tmp0,tmp1/])
  else
    ; -- load data --
    print("load data...")
    ifile1 = pathout+"ec_earth_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(5)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
    i = addfile(ifile1,"r")
    ec_earth_trends    = i->ec_earth_trends
    ec_earthEM         = i->ec_earthEM
    ec_earthSD         = i->ec_earthSD
    ec_earthSN         = i->ec_earthSN
  end if


  ; -- mpi_lens --
  print_clock("mpi_lens")
  if recalc7 .eq. "yes"
    start   = 1850
    ende    = 2099
    ensmem  = 100
    a = addfile(pathin+"mpi_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(6)+"_historical_rcp85_r2i1p1_185001-209912.nc","r")
    dummy   = a->$var$
    dims    = dimsizes(dummy)
    ntime   = dims(0)
    nlat    = dims(1)
    nlon    = dims(2)
    mpi_trends  = new((/ensmem,nlat,nlon/),float)
    do e = 1,ensmem
      print_clock("  ensmem = "+e+"/"+ensmem)
      a = addfile(pathin+"mpi_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(6)+"_historical_rcp85_r"+e+"i1p1_185001-209912.nc","r")
      tmp0     = a->$var$
      if seas .eq. "annual"
        tmp1 = month_to_annual(tmp0,1)
      else
        tmp1 = month_to_season(tmp0,seas)
      end if
      mpi_trends(e-1,:,:)  = regCoef_n(trend_time, tmp1(start_trend-start:ende_trend-start,:,:), 0, 0) * (ende_trend-start_trend+1)
    end do
    copy_VarCoords(dummy(0,:,:),mpi_trends(0,:,:))
    mpiEM = dim_avg_n_Wrap(mpi_trends,0)
    mpiSD = dim_stddev_n_Wrap(mpi_trends,0)
    mpiSN = (/ abs(mpiEM)/mpiSD /)
    copy_VarCoords(dummy(0,:,:),mpiSN)
    ; -- save data --
    print("save data...")
    ofile1 = pathout+"mpi_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(6)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
    system("mkdir -p "+pathout+"mpi_lens/"+comp+"/"+var)
    system("rm -rf "+ofile1)
    o = addfile(ofile1,"c")
    o->mpi_trends        = mpi_trends
    o->mpiEM         = mpiEM
    o->mpiSD         = mpiSD
    o->mpiSN         = mpiSN
    delete([/dummy,dims,ntime,nlat,nlon,tmp0,tmp1/])
  else
    ; -- load data --
    print("load data...")
    ifile1 = pathout+"mpi_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(6)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
    i = addfile(ifile1,"r")
    mpi_trends    = i->mpi_trends
    mpiEM         = i->mpiEM
    mpiSD         = i->mpiSD
    mpiSN         = i->mpiSN
  end if




  ; -- climatology (over same period as trend period) --------------------------
  m = 0
  if recalc_em .eq. "yes"
    print_clock("calc ensemble mean...")
    system("ncea -O "+pathin+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_*192001-210012*.nc "+pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_ensmean.nc")
  end if
  if recalc_clim .eq. "yes"
    print_clock("calc climatology...")
    ifile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_ensmean.nc"
    i = addfile(ifile1,"r")
    tmp0         = i->$var$
    if seas .eq. "annual"
      tmp1 = month_to_annual(tmp0,1)
    else
      tmp1 = month_to_season(tmp0,seas)
    end if
    cesm_clim   = dim_avg_n_Wrap(tmp1(start_trend-start0(m):ende_trend-start0(m)+1,:,:),0)
    tmp2 = trend_manken(tmp1(start_trend-start0(m):ende_trend-start0(m)+1,:,:),False,0)
    cesm_trend_signi = tmp2(0,:,:)
    copy_VarCoords(tmp0(0,:,:),cesm_trend_signi)
    ofile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
    o = addfile(ofile1,"w")
    o->cesm_clim         = cesm_clim
    o->cesm_trend_signi  = cesm_trend_signi
    delete([/tmp0,tmp1,tmp2/])
  end if
  ofile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
  o = addfile(ofile1,"r")
  cesm_clim         = o->cesm_clim
  cesm_trend_signi  = o->cesm_trend_signi
  files = systemfunc("ls -1 /project/cmip5/piControl/fx/sftlf/"+model_name(m)+"/r0i0p0/*.nc")
  ifile1 = files(0)
  i = addfile(ifile1,"r")
  tmp           = i->sftlf
  cesm_landfrac  = where(tmp .gt. 50,1,tmp@_FillValue)
  delete(tmp)

  m = 1
  if recalc_em .eq. "yes"
    print_clock("calc ensemble mean...")
    system("ncea -O "+pathin+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_*.nc "+pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_ensmean.nc")
  end if
  if recalc_clim .eq. "yes"
    print_clock("calc climatology...")
    ifile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_ensmean.nc"
    i = addfile(ifile1,"r")
    tmp0         = i->$var$
    if seas .eq. "annual"
      tmp1 = month_to_annual(tmp0,1)
    else
      tmp1 = month_to_season(tmp0,seas)
    end if
    canesm2_clim   = dim_avg_n_Wrap(tmp1(start_trend-start0(m):ende_trend-start0(m)+1,:,:),0)
    tmp2 = trend_manken(tmp1(start_trend-start0(m):ende_trend-start0(m)+1,:,:),False,0)
    canesm2_trend_signi = tmp2(0,:,:)
    copy_VarCoords(tmp0(0,:,:),canesm2_trend_signi)
    ofile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
    o = addfile(ofile1,"w")
    o->canesm2_clim         = canesm2_clim
    o->canesm2_trend_signi  = canesm2_trend_signi
    delete([/tmp0,tmp1,tmp2/])
  end if
  ofile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
  o = addfile(ofile1,"r")
  canesm2_clim         = o->canesm2_clim
  canesm2_trend_signi  = o->canesm2_trend_signi
  files = systemfunc("ls -1 /project/cmip5/piControl/fx/sftlf/"+model_name(m)+"/r0i0p0/*.nc")
  ifile1 = files(0)
  i = addfile(ifile1,"r")
  tmp           = i->sftlf
  canesm2_landfrac  = where(tmp .gt. 50,1,tmp@_FillValue)
  delete(tmp)

  m = 2
  if recalc_em .eq. "yes"
    print_clock("calc ensemble mean...")
    system("ncea -O "+pathin+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_*.nc "+pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_ensmean.nc")
  end if
  if recalc_clim .eq. "yes"
    print_clock("calc climatology...")
    ifile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_ensmean.nc"
    i = addfile(ifile1,"r")
    tmp0         = i->$var$
    if seas .eq. "annual"
      tmp1 = month_to_annual(tmp0,1)
    else
      tmp1 = month_to_season(tmp0,seas)
    end if
    csiro_mk36_clim   = dim_avg_n_Wrap(tmp1(start_trend-start0(m):ende_trend-start0(m)+1,:,:),0)
    tmp2 = trend_manken(tmp1(start_trend-start0(m):ende_trend-start0(m)+1,:,:),False,0)
    csiro_mk36_trend_signi = tmp2(0,:,:)
    copy_VarCoords(tmp0(0,:,:),csiro_mk36_trend_signi)
    ofile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
    o = addfile(ofile1,"w")
    o->csiro_mk36_clim         = csiro_mk36_clim
    o->csiro_mk36_trend_signi  = csiro_mk36_trend_signi
    delete([/tmp0,tmp1,tmp2/])
  end if
  ofile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
  o = addfile(ofile1,"r")
  csiro_mk36_clim         = o->csiro_mk36_clim
  csiro_mk36_trend_signi  = o->csiro_mk36_trend_signi
  files = systemfunc("ls -1 /project/cmip5/piControl/fx/sftlf/"+model_name(m)+"/r0i0p0/*.nc")
  ifile1 = files(0)
  i = addfile(ifile1,"r")
  tmp           = i->sftlf
  csiro_mk36_landfrac  = where(tmp .gt. 50,1,tmp@_FillValue)
  delete(tmp)

  m = 3
  if recalc_em .eq. "yes"
    print_clock("calc ensemble mean...")
    system("ncea -O "+pathin+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_*.nc "+pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_ensmean.nc")
  end if
  if recalc_clim .eq. "yes"
    print_clock("calc climatology...")
    ifile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_ensmean.nc"
    i = addfile(ifile1,"r")
    tmp0         = i->$var$
    if seas .eq. "annual"
      tmp1 = month_to_annual(tmp0,1)
    else
      tmp1 = month_to_season(tmp0,seas)
    end if
    gfdl_cm3_clim   = dim_avg_n_Wrap(tmp1(start_trend-start0(m):ende_trend-start0(m)+1,:,:),0)
    tmp2 = trend_manken(tmp1(start_trend-start0(m):ende_trend-start0(m)+1,:,:),False,0)
    gfdl_cm3_trend_signi = tmp2(0,:,:)
    copy_VarCoords(tmp0(0,:,:),gfdl_cm3_trend_signi)
    ofile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
    o = addfile(ofile1,"w")
    o->gfdl_cm3_clim         = gfdl_cm3_clim
    o->gfdl_cm3_trend_signi  = gfdl_cm3_trend_signi
    delete([/tmp0,tmp1,tmp2/])
  end if
  ofile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
  o = addfile(ofile1,"r")
  gfdl_cm3_clim         = o->gfdl_cm3_clim
  gfdl_cm3_trend_signi  = o->gfdl_cm3_trend_signi
  files = systemfunc("ls -1 /project/cmip5/piControl/fx/sftlf/"+model_name(m)+"/r0i0p0/*.nc")
  ifile1 = files(0)
  i = addfile(ifile1,"r")
  tmp           = i->sftlf
  gfdl_cm3_landfrac  = where(tmp .gt. 50,1,tmp@_FillValue)
  delete(tmp)

  m = 4
  if recalc_em .eq. "yes"
    print_clock("calc ensemble mean...")
    system("ncea -O "+pathin+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_*.nc "+pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_ensmean.nc")
  end if
  if recalc_clim .eq. "yes"
    print_clock("calc climatology...")
    ifile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_ensmean.nc"
    i = addfile(ifile1,"r")
    tmp0         = i->$var$
    if seas .eq. "annual"
      tmp1 = month_to_annual(tmp0,1)
    else
      tmp1 = month_to_season(tmp0,seas)
    end if
    gfdl_esm2m_clim   = dim_avg_n_Wrap(tmp1(start_trend-start0(m):ende_trend-start0(m)+1,:,:),0)
    tmp2 = trend_manken(tmp1(start_trend-start0(m):ende_trend-start0(m)+1,:,:),False,0)
    gfdl_esm2m_trend_signi = tmp2(0,:,:)
    copy_VarCoords(tmp0(0,:,:),gfdl_esm2m_trend_signi)
    ofile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
    o = addfile(ofile1,"w")
    o->gfdl_esm2m_clim         = gfdl_esm2m_clim
    o->gfdl_esm2m_trend_signi  = gfdl_esm2m_trend_signi
    delete([/tmp0,tmp1,tmp2/])
  end if
  ofile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
  o = addfile(ofile1,"r")
  gfdl_esm2m_clim         = o->gfdl_esm2m_clim
  gfdl_esm2m_trend_signi  = o->gfdl_esm2m_trend_signi
  files = systemfunc("ls -1 /project/cmip5/piControl/fx/sftlf/"+model_name(m)+"/r0i0p0/*.nc")
  ifile1 = files(0)
  i = addfile(ifile1,"r")
  tmp           = i->sftlf
  gfdl_esm2m_landfrac  = where(tmp .gt. 50,1,tmp@_FillValue)
  delete(tmp)

  m = 5
  if recalc_em .eq. "yes"
    print_clock("calc ensemble mean...")
    system("ncea -O "+pathin+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_*.nc "+pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_ensmean.nc")
  end if
  if recalc_clim .eq. "yes"
    print_clock("calc climatology...")
    ifile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_ensmean.nc"
    i = addfile(ifile1,"r")
    tmp0         = i->$var$
    if seas .eq. "annual"
      tmp1 = month_to_annual(tmp0,1)
    else
      tmp1 = month_to_season(tmp0,seas)
    end if
    ec_earth_clim   = dim_avg_n_Wrap(tmp1(start_trend-start0(m):ende_trend-start0(m)+1,:,:),0)
    tmp2 = trend_manken(tmp1(start_trend-start0(m):ende_trend-start0(m)+1,:,:),False,0)
    ec_earth_trend_signi = tmp2(0,:,:)
    copy_VarCoords(tmp0(0,:,:),ec_earth_trend_signi)
    ofile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
    o = addfile(ofile1,"w")
    o->ec_earth_clim         = ec_earth_clim
    o->ec_earth_trend_signi  = ec_earth_trend_signi
    delete([/tmp0,tmp1,tmp2/])
  end if
  ofile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
  o = addfile(ofile1,"r")
  ec_earth_clim         = o->ec_earth_clim
  ec_earth_trend_signi  = o->ec_earth_trend_signi
  files = systemfunc("ls -1 /project/cmip5/piControl/fx/sftlf/"+model_name(m)+"/r0i0p0/*.nc")
  ifile1 = files(0)
  i = addfile(ifile1,"r")
  tmp           = i->sftlf
  ec_earth_landfrac  = where(tmp .gt. 50,1,tmp@_FillValue)
  delete(tmp)

  m = 6
  if recalc_em .eq. "yes"
    print_clock("calc ensemble mean...")
    system("ncea -O "+pathin+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_*.nc "+pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_ensmean.nc")
  end if
  if recalc_clim .eq. "yes"
    print_clock("calc climatology...")
    ifile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_ensmean.nc"
    i = addfile(ifile1,"r")
    tmp0         = i->$var$
    if seas .eq. "annual"
      tmp1 = month_to_annual(tmp0,1)
    else
      tmp1 = month_to_season(tmp0,seas)
    end if
    mpi_clim   = dim_avg_n_Wrap(tmp1(start_trend-start0(m):ende_trend-start0(m)+1,:,:),0)
    tmp2 = trend_manken(tmp1(start_trend-start0(m):ende_trend-start0(m)+1,:,:),False,0)
    mpi_trend_signi = tmp2(0,:,:)
    copy_VarCoords(tmp0(0,:,:),mpi_trend_signi)
    ofile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
    o = addfile(ofile1,"w")
    o->mpi_clim         = mpi_clim
    o->mpi_trend_signi  = mpi_trend_signi
    delete([/tmp0,tmp1,tmp2/])
  end if
  ofile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_"+seas+"_trends_"+start_trend+"-"+ende_trend+".nc"
  o = addfile(ofile1,"r")
  mpi_clim         = o->mpi_clim
  mpi_trend_signi  = o->mpi_trend_signi
  files = systemfunc("ls -1 /project/cmip5/piControl/fx/sftlf/MPI-ESM-MR/r0i0p0/*.nc")
  ifile1 = files(0)
  i = addfile(ifile1,"r")
  tmp           = i->sftlf
  mpi_landfrac  = where(tmp .gt. 50,1,tmp@_FillValue)
  delete(tmp)


  ; -- Unit conversions --
  if var .eq. "pr" .or. var .eq. "psl"
    if var .eq. "pr"
      f = 86400; unit conversion factor: mm/s to mm
    end if
    if var .eq. "psl"
      f = 0.01; unit conversion factor: Pa to hPa
    end if
    cesm_trends       = cesm_trends*f
    cesmSD            = cesmSD*f
    cesm_clim         = cesm_clim*f
    canesm2_trends    = canesm2_trends*f
    canesm2SD         = canesm2SD*f
    canesm2_clim      = canesm2_clim*f
    csiro_mk36_trends = csiro_mk36_trends*f
    csiro_mk36SD      = csiro_mk36SD*f
    csiro_mk36_clim   = csiro_mk36_clim*f
    gfdl_cm3_trends   = gfdl_cm3_trends*f
    gfdl_cm3SD        = gfdl_cm3SD*f
    gfdl_cm3_clim     = gfdl_cm3_clim*f
    gfdl_esm2m_trends = gfdl_esm2m_trends*f
    gfdl_esm2mSD      = gfdl_esm2mSD*f
    gfdl_esm2m_clim   = gfdl_esm2m_clim*f
    ec_earth_trends   = ec_earth_trends*f
    ec_earthSD        = ec_earthSD*f
    ec_earth_clim     = ec_earth_clim*f
    mpi_trends        = mpi_trends*f
    mpiSD             = mpiSD*f
    mpi_clim          = mpi_clim*f
    delete(f)
  end if


  ; -- domain --
  if domain .eq. "NA"
    minlon = 200.
    maxlon = 298.
    minlat = 24.
    maxlat = 72.
  else
    minlon = 0
    maxlon = 360
    minlat = -90
    maxlat = 90
  end if
  rad      = 4.0*atan(1.0)/180.0

  ; -- function to find bookend ensemble members for each model --
  function max_min_trends(infield,m) ; needs dimensions [ensmem,lat,lon]
  begin
    tmp0        = infield
    ; -- land fraction file:
    if m .eq. 6
      files = systemfunc("ls -1 /project/cmip5/piControl/fx/sftlf/MPI-ESM-MR/r0i0p0/*.nc")
    else
      files = systemfunc("ls -1 /project/cmip5/piControl/fx/sftlf/"+model_name(m)+"/r0i0p0/*.nc")
    end if
    ifile1 = files(0)
    i = addfile(ifile1,"r")
    tmp       = i->sftlf
    if max(max(tmp)) .lt. 2
      landfrac0  = where(tmp .gt. .5,1,tmp0@_FillValue)
    else
      landfrac0  = where(tmp .gt. 50,1,tmp0@_FillValue)
    end if
    landfrac = conform_dims(dimsizes(tmp0),landfrac0,(/1,2/))
    if var .eq. "tas" .or. var .eq. "pr"
      tmp1 = where(landfrac .eq. 1,tmp0,tmp0@_FillValue)
      copy_VarCoords(tmp0,tmp1)
      delete(tmp0)
      tmp0 = tmp1
      delete(tmp1)
    end if
    ; -- weights:
    lat       = tmp0&lat
    clat      = cos(lat*rad)
    clat!0    = "lat"
    clat&lat  = lat
    ; -- area-average:
    tmp0!0       = "trends"
    tmp0&trends  = fspan(1,dimsizes(tmp0(:,0,0)),dimsizes(tmp0(:,0,0)))
    tmp1         = wgt_areaave_Wrap(tmp0({trends|:},{lat|minlat:maxlat},{lon|minlon:maxlon}),clat({lat|minlat:maxlat}),1.0,0)
    max_min_idx  = new((/5/),float)
    max_min_idx(0) = maxind(tmp1)
    max_min_idx(1) = minind(tmp1)
    max_min_idx(2) = max(tmp1)
    max_min_idx(3) = min(tmp1)
    max_min_idx(4) = avg(tmp1)
    delete([/tmp,tmp0,tmp1/])
    return(max_min_idx)
  end

  cesm_idx        = max_min_trends(cesm_trends,0)
  canesm2_idx     = max_min_trends(canesm2_trends,1)
  csiro_mk36_idx  = max_min_trends(csiro_mk36_trends,2)
  gfdl_cm3_idx    = max_min_trends(gfdl_cm3_trends,3)
  gfdl_esm2m_idx  = max_min_trends(gfdl_esm2m_trends,4)
  ec_earth_idx    = max_min_trends(ec_earth_trends,5)
  mpi_idx         = max_min_trends(mpi_trends,6)

  function trend_values(infield,m) ; infield needs dimensions [ensmem,lat,lon]
  begin
    tmp0        = infield
    ; -- land fraction file:
    if m .eq. 6
      files = systemfunc("ls -1 /project/cmip5/piControl/fx/sftlf/MPI-ESM-MR/r0i0p0/*.nc")
    else
      files = systemfunc("ls -1 /project/cmip5/piControl/fx/sftlf/"+model_name(m)+"/r0i0p0/*.nc")
    end if
    ifile1 = files(0)
    i = addfile(ifile1,"r")
    tmp       = i->sftlf
    if max(max(tmp)) .lt. 2
      landfrac0  = where(tmp .gt. .5,1,tmp0@_FillValue)
    else
      landfrac0  = where(tmp .gt. 50,1,tmp0@_FillValue)
    end if
    landfrac = conform_dims(dimsizes(tmp0),landfrac0,(/1,2/))
    if var .eq. "tas" .or. var .eq. "pr"
      tmp1 = where(landfrac .eq. 1,tmp0,tmp0@_FillValue)
      copy_VarCoords(tmp0,tmp1)
      delete(tmp0)
      tmp0 = tmp1
      delete(tmp1)
    end if
    ; -- weights:
    lat       = tmp0&lat
    clat      = cos(lat*rad)
    clat!0    = "lat"
    clat&lat  = lat
    ; -- area-average:
    tmp0!0       = "trends"
    tmp0&trends  = fspan(1,ensmem0(m),ensmem0(m))
    tmp1         = wgt_areaave_Wrap(tmp0({trends|:},{lat|minlat:maxlat},{lon|minlon:maxlon}),clat({lat|minlat:maxlat}),1.0,0)
    delete([/tmp,tmp0/])
    return(tmp1)
  end

  cesm_trend_values        = trend_values(cesm_trends,0)
  delete(tmp1)
  canesm2_trend_values     = trend_values(canesm2_trends,1)
  delete(tmp1)
  csiro_mk36_trend_values  = trend_values(csiro_mk36_trends,2)
  delete(tmp1)
  gfdl_cm3_trend_values    = trend_values(gfdl_cm3_trends,3)
  delete(tmp1)
  gfdl_esm2m_trend_values  = trend_values(gfdl_esm2m_trends,4)
  delete(tmp1)
  ec_earth_trend_values    = trend_values(ec_earth_trends,5)
  delete(tmp1)
  mpi_trend_values         = trend_values(mpi_trends,6)
  delete(tmp1)
  ; -- save data:
  m = 0
  ofile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_"+seas+"_trend_values_"+start_trend+"-"+ende_trend+".nc"
  system("rm -rf "+ofile1)
  o = addfile(ofile1,"c")
  o->trend_values         = cesm_trend_values
  m = 1
  ofile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_"+seas+"_trend_values_"+start_trend+"-"+ende_trend+".nc"
  system("rm -rf "+ofile1)
  o = addfile(ofile1,"c")
  o->trend_values         = canesm2_trend_values
  m = 2
  ofile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_"+seas+"_trend_values_"+start_trend+"-"+ende_trend+".nc"
  system("rm -rf "+ofile1)
  o = addfile(ofile1,"c")
  o->trend_values         = csiro_mk36_trend_values
  m = 3
  ofile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_"+seas+"_trend_values_"+start_trend+"-"+ende_trend+".nc"
  system("rm -rf "+ofile1)
  o = addfile(ofile1,"c")
  o->trend_values         = gfdl_cm3_trend_values
  m = 4
  ofile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_"+seas+"_trend_values_"+start_trend+"-"+ende_trend+".nc"
  system("rm -rf "+ofile1)
  o = addfile(ofile1,"c")
  o->trend_values         = gfdl_esm2m_trend_values
  m = 5
  ofile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_"+seas+"_trend_values_"+start_trend+"-"+ende_trend+".nc"
  system("rm -rf "+ofile1)
  o = addfile(ofile1,"c")
  o->trend_values         = ec_earth_trend_values
  m = 6
  ofile1 = pathout+clivar_tag(m)+"_lens/"+comp+"/"+var+"/"+var+"_"+comp+"_"+model_name(m)+"_"+seas+"_trend_values_"+start_trend+"-"+ende_trend+".nc"
  system("rm -rf "+ofile1)
  o = addfile(ofile1,"c")
  o->trend_values         = mpi_trend_values



  ; -- load obs --
  if var .eq. "pr"
    obs_name = "GPCC"
    start   = 1901
    ende    = 2017
    pathin_obs = "/project/cas/flehner/observations/precipitation/gpcc/"
    ifile1  = pathin_obs+"precip.mon.combined.total.v7.190101-201712.nc"
    i = addfile(ifile1,"r")
    tmp0    = i->precip
    mv      = tmp0@_FillValue
    tmp0    = tmp0(:,::-1,:)
    ; tmp0    = lonFlip(tmp0)
    lat_o   = tmp0&lat
    lon_o   = tmp0&lon
    if seas .eq. "annual"
      tmp1 = month_to_annual(tmp0,1)
    else
      tmp1 = month_to_season(tmp0,seas)
    end if
    tmp1 = tmp1/30.4
    obs_trend_raw  = regCoef_n(trend_time, tmp1(start_trend-start:ende_trend-start,:,:), 0, 0) * (ende_trend-start_trend+1)
    copy_VarCoords(tmp0(0,:,:),obs_trend_raw)
  else if var .eq. "tas"
    obs_name = "BEST"
    start   = 1950
    ende    = 2018
    ifile1  = pathin_obs+"Complete_TAVG_LatLong1_195001-201812.nc"
    i = addfile(ifile1,"r")
    tmp0    = i->temperature
    mv      = tofloat(-999)
    tmp0@_FillValue = mv
    tmp0@missing_value = mv
    ; tmp0    = tmp0(:,::-1,:)
    tmp0    = lonFlip(tmp0)
    lat_o   = tmp0&latitude
    lon_o   = tmp0&longitude
    if seas .eq. "annual"
      tmp1 = month_to_annual(tmp0,1)
    else
      tmp1 = month_to_season(tmp0,seas)
    end if
    obs_trend_raw  = regCoef_n(trend_time, tmp1(start_trend-start:ende_trend-start,:,:), 0, 0) * (ende_trend-start_trend+1)
    tmp        = new((/180,360/),float)
    tmp(:,:)   = obs_trend_raw(:,:)
    tmp!0 = "lat"
    tmp!1 = "lon"
    lati = fspan(-89.5,89.5,180)
    loni = fspan(0.5,359.5,360)
    lati@units = "degrees_north"
    loni@units = "degrees_east"
    tmp&lat    = lati
    tmp&lon    = loni
    delete(obs_trend_raw)
    obs_trend_raw = tmp
    delete(tmp)
    replace_ieeenan (obs_trend_raw, mv, 0)
    obs_trend_raw@_FillValue = mv
    obs_trend_raw@missing_value = mv
  end if
  end if
  ; -- weights:
  lat       = obs_trend_raw&lat
  clat      = cos(lat*rad)
  clat!0    = "lat"
  clat&lat  = lat
  ; -- area-average:
  obs_trend_values         = wgt_areaave_Wrap(obs_trend_raw({lat|minlat:maxlat},{lon|minlon:maxlon}),clat({lat|minlat:maxlat}),1.0,0)
  ; -- save data:
  ofile1 = pathin_obs+var+"_"+seas+"_trend_values_"+start_trend+"-"+ende_trend+".nc"
  system("rm -rf "+ofile1)
  o = addfile(ofile1,"c")
  o->trend_values         = obs_trend_values


  if be .eq. 2
    ; -- calculate pattern correlation with observations --
    print_clock("pattern correlation...")
    ; -- regrid to individual model and olens grids --
    ; -- olens
    lat_m = olens_trends&lat
    lon_m = olens_trends&lon
    ; obs_trend   = area_conserve_remap_Wrap (lon_o,lat_o,obs_trend_raw,lon_m,lat_m,False)
    obs_trend   = linint2_Wrap (lon_o,lat_o,obs_trend_raw,False,lon_m,lat_m,0)
    msk         = where(obs_trend .lt. -900,mv,1)
    olens_trends_masked_obs = olens_trends * conform_dims(dimsizes(olens_trends),msk,(/1,2/))
    replace_ieeenan (obs_trend, mv, 0)
    obs_trend@_FillValue = mv
    obs_trend@missing_value = mv
    olens_pc     = new((/1000/),float)
    do e = 0,999
      tmp = tofloat(olens_trends_masked_obs(e,:,:))
      replace_ieeenan (tmp, mv, 0)
      tmp@_FillValue = mv
      tmp@missing_value = mv
      tmp!0 = "lat"
      tmp!1 = "lon"
      tmp&lat   = olens_trends&lat
      tmp&lon   = olens_trends&lon
      clat = tofloat(cos(0.01745329*olens_trends&lat))
      clat!0 = "lat"
      clat&lat = olens_trends&lat
      olens_pc(e)     = pattern_cor(tmp({lat|minlat:maxlat},{lon|minlon:maxlon}),obs_trend({lat|minlat:maxlat},{lon|minlon:maxlon}),clat({lat|minlat:maxlat}),cent)
    end do
    olens_idx = new((/4/),float)
    olens_idx(0) = maxind(olens_pc)
    olens_idx(1) = minind(olens_pc)
    olens_idx(2) = olens_pc(maxind(olens_pc))
    olens_idx(3) = olens_pc(minind(olens_pc))
    delete([/lat_m,lon_m,tmp,msk,obs_trend,clat/])
    ; -- cesm
    m = 0
    lat_m = cesm_trends&lat
    lon_m = cesm_trends&lon
    obs_trend   = linint2_Wrap (lon_o,lat_o,obs_trend_raw,False,lon_m,lat_m,0)
    msk = where(obs_trend .ne. obs_trend@_FillValue,1,obs_trend@_FillValue)
    cesm_trends_masked_obs = cesm_trends * conform_dims(dimsizes(cesm_trends),msk,(/1,2/))
    cesm_pc     = new((/ensmem0(m)/),float)
    replace_ieeenan (obs_trend, mv, 0)
    obs_trend@_FillValue = mv
    obs_trend@missing_value = mv
    do e = 0,ensmem0(m)-1
      tmp = tofloat(cesm_trends_masked_obs(e,:,:))
      replace_ieeenan (tmp, mv, 0)
      tmp@_FillValue = mv
      tmp@missing_value = mv
      tmp!0 = "lat"
      tmp!1 = "lon"
      tmp&lat = cesm_trends&lat
      tmp&lon = cesm_trends&lon
      clat = tofloat(cos(0.01745329*obs_trend&lat))
      clat!0 = "lat"
      clat&lat = obs_trend&lat
      cesm_pc(e)     = pattern_cor(tmp({lat|minlat:maxlat},{lon|minlon:maxlon}),obs_trend({lat|minlat:maxlat},{lon|minlon:maxlon}),clat({lat|minlat:maxlat}),cent)
    end do
    cesm_idx(0) = maxind(cesm_pc)
    cesm_idx(1) = minind(cesm_pc)
    cesm_idx(2) = cesm_pc(maxind(cesm_pc))
    cesm_idx(3) = cesm_pc(minind(cesm_pc))
    delete([/lat_m,lon_m,tmp,msk,obs_trend,clat/])
    ; -- canesm2
    m = 1
    lat_m = canesm2_trends&lat
    lon_m = canesm2_trends&lon
    obs_trend   = linint2_Wrap (lon_o,lat_o,obs_trend_raw,False,lon_m,lat_m,0)
    msk = where(obs_trend .ne. obs_trend@_FillValue,1,obs_trend@_FillValue)
    canesm2_trends_masked_obs = canesm2_trends * conform_dims(dimsizes(canesm2_trends),msk,(/1,2/))
    canesm2_pc     = new((/ensmem0(m)/),float)
    replace_ieeenan (obs_trend, mv, 0)
    obs_trend@_FillValue = mv
    obs_trend@missing_value = mv
    do e = 0,ensmem0(m)-1
      tmp = tofloat(canesm2_trends_masked_obs(e,:,:))
      replace_ieeenan (tmp, mv, 0)
      tmp@_FillValue = mv
      tmp@missing_value = mv
      tmp!0 = "lat"
      tmp!1 = "lon"
      tmp&lat = canesm2_trends&lat
      tmp&lon = canesm2_trends&lon
      clat = tofloat(cos(0.01745329*obs_trend&lat))
      clat!0 = "lat"
      clat&lat = obs_trend&lat
      canesm2_pc(e)     = pattern_cor(tmp({lat|minlat:maxlat},{lon|minlon:maxlon}),obs_trend({lat|minlat:maxlat},{lon|minlon:maxlon}),clat({lat|minlat:maxlat}),cent)
    end do
    canesm2_idx(0) = maxind(canesm2_pc)
    canesm2_idx(1) = minind(canesm2_pc)
    canesm2_idx(2) = canesm2_pc(maxind(canesm2_pc))
    canesm2_idx(3) = canesm2_pc(minind(canesm2_pc))
    delete([/lat_m,lon_m,tmp,msk,obs_trend,clat/])
    ; -- csiro_mk36
    m = 2
    lat_m = csiro_mk36_trends&lat
    lon_m = csiro_mk36_trends&lon
    obs_trend   = linint2_Wrap (lon_o,lat_o,obs_trend_raw,False,lon_m,lat_m,0)
    msk = where(obs_trend .ne. obs_trend@_FillValue,1,obs_trend@_FillValue)
    csiro_mk36_trends_masked_obs = csiro_mk36_trends * conform_dims(dimsizes(csiro_mk36_trends),msk,(/1,2/))
    csiro_mk36_pc     = new((/ensmem0(m)/),float)
    replace_ieeenan (obs_trend, mv, 0)
    obs_trend@_FillValue = mv
    obs_trend@missing_value = mv
    do e = 0,ensmem0(m)-1
      tmp = tofloat(csiro_mk36_trends_masked_obs(e,:,:))
      replace_ieeenan (tmp, mv, 0)
      tmp@_FillValue = mv
      tmp@missing_value = mv
      tmp!0 = "lat"
      tmp!1 = "lon"
      tmp&lat = csiro_mk36_trends&lat
      tmp&lon = csiro_mk36_trends&lon
      clat = tofloat(cos(0.01745329*obs_trend&lat))
      clat!0 = "lat"
      clat&lat = obs_trend&lat
      csiro_mk36_pc(e)     = pattern_cor(tmp({lat|minlat:maxlat},{lon|minlon:maxlon}),obs_trend({lat|minlat:maxlat},{lon|minlon:maxlon}),clat({lat|minlat:maxlat}),cent)
    end do
    csiro_mk36_idx(0) = maxind(csiro_mk36_pc)
    csiro_mk36_idx(1) = minind(csiro_mk36_pc)
    csiro_mk36_idx(2) = csiro_mk36_pc(maxind(csiro_mk36_pc))
    csiro_mk36_idx(3) = csiro_mk36_pc(minind(csiro_mk36_pc))
    delete([/lat_m,lon_m,tmp,msk,obs_trend,clat/])
    ; -- gfdl_cm3
    m = 3
    lat_m = gfdl_cm3_trends&lat
    lon_m = gfdl_cm3_trends&lon
    obs_trend   = linint2_Wrap (lon_o,lat_o,obs_trend_raw,False,lon_m,lat_m,0)
    msk = where(obs_trend .ne. obs_trend@_FillValue,1,obs_trend@_FillValue)
    gfdl_cm3_trends_masked_obs = gfdl_cm3_trends * conform_dims(dimsizes(gfdl_cm3_trends),msk,(/1,2/))
    gfdl_cm3_pc     = new((/ensmem0(m)/),float)
    replace_ieeenan (obs_trend, mv, 0)
    obs_trend@_FillValue = mv
    obs_trend@missing_value = mv
    do e = 0,ensmem0(m)-1
      tmp = tofloat(gfdl_cm3_trends_masked_obs(e,:,:))
      replace_ieeenan (tmp, mv, 0)
      tmp@_FillValue = mv
      tmp@missing_value = mv
      tmp!0 = "lat"
      tmp!1 = "lon"
      tmp&lat = gfdl_cm3_trends&lat
      tmp&lon = gfdl_cm3_trends&lon
      clat = tofloat(cos(0.01745329*obs_trend&lat))
      clat!0 = "lat"
      clat&lat = obs_trend&lat
      gfdl_cm3_pc(e)     = pattern_cor(tmp({lat|minlat:maxlat},{lon|minlon:maxlon}),obs_trend({lat|minlat:maxlat},{lon|minlon:maxlon}),clat({lat|minlat:maxlat}),cent)
    end do
    gfdl_cm3_idx(0) = maxind(gfdl_cm3_pc)
    gfdl_cm3_idx(1) = minind(gfdl_cm3_pc)
    gfdl_cm3_idx(2) = gfdl_cm3_pc(maxind(gfdl_cm3_pc))
    gfdl_cm3_idx(3) = gfdl_cm3_pc(minind(gfdl_cm3_pc))
    delete([/lat_m,lon_m,tmp,msk,obs_trend,clat/])
    ; -- gfdl_esm2m
    m = 4
    lat_m = gfdl_esm2m_trends&lat
    lon_m = gfdl_esm2m_trends&lon
    obs_trend   = linint2_Wrap (lon_o,lat_o,obs_trend_raw,False,lon_m,lat_m,0)
    msk = where(obs_trend .ne. obs_trend@_FillValue,1,obs_trend@_FillValue)
    gfdl_esm2m_trends_masked_obs = gfdl_esm2m_trends * conform_dims(dimsizes(gfdl_esm2m_trends),msk,(/1,2/))
    gfdl_esm2m_pc     = new((/ensmem0(m)/),float)
    replace_ieeenan (obs_trend, mv, 0)
    obs_trend@_FillValue = mv
    obs_trend@missing_value = mv
    do e = 0,ensmem0(m)-1
      tmp = tofloat(gfdl_esm2m_trends_masked_obs(e,:,:))
      replace_ieeenan (tmp, mv, 0)
      tmp@_FillValue = mv
      tmp@missing_value = mv
      tmp!0 = "lat"
      tmp!1 = "lon"
      tmp&lat = gfdl_esm2m_trends&lat
      tmp&lon = gfdl_esm2m_trends&lon
      clat = tofloat(cos(0.01745329*obs_trend&lat))
      clat!0 = "lat"
      clat&lat = obs_trend&lat
      gfdl_esm2m_pc(e)     = pattern_cor(tmp({lat|minlat:maxlat},{lon|minlon:maxlon}),obs_trend({lat|minlat:maxlat},{lon|minlon:maxlon}),clat({lat|minlat:maxlat}),cent)
    end do
    gfdl_esm2m_idx(0) = maxind(gfdl_esm2m_pc)
    gfdl_esm2m_idx(1) = minind(gfdl_esm2m_pc)
    gfdl_esm2m_idx(2) = gfdl_esm2m_pc(maxind(gfdl_esm2m_pc))
    gfdl_esm2m_idx(3) = gfdl_esm2m_pc(minind(gfdl_esm2m_pc))
    delete([/lat_m,lon_m,tmp,msk,obs_trend,clat/])
    ; -- ec_earth
    m = 5
    ec_earth_trends = ec_earth_trends(:,::-1,:)
    lat_m = ec_earth_trends&lat
    lon_m = ec_earth_trends&lon
    obs_trend   = linint2_Wrap (lon_o,lat_o,obs_trend_raw,False,lon_m,lat_m,0)
    msk = where(obs_trend .ne. obs_trend@_FillValue,1,obs_trend@_FillValue)
    ec_earth_trends_masked_obs = ec_earth_trends * conform_dims(dimsizes(ec_earth_trends),msk,(/1,2/))
    ec_earth_pc     = new((/ensmem0(m)/),float)
    replace_ieeenan (obs_trend, mv, 0)
    obs_trend@_FillValue = mv
    obs_trend@missing_value = mv
    do e = 0,ensmem0(m)-1
      tmp = tofloat(ec_earth_trends_masked_obs(e,:,:))
      replace_ieeenan (tmp, mv, 0)
      tmp@_FillValue = mv
      tmp@missing_value = mv
      tmp!0 = "lat"
      tmp!1 = "lon"
      tmp&lat = ec_earth_trends&lat
      tmp&lon = ec_earth_trends&lon
      clat = tofloat(cos(0.01745329*obs_trend&lat))
      clat!0 = "lat"
      clat&lat = obs_trend&lat
      ec_earth_pc(e)     = pattern_cor(tmp({lat|minlat:maxlat},{lon|minlon:maxlon}),obs_trend({lat|minlat:maxlat},{lon|minlon:maxlon}),clat({lat|minlat:maxlat}),cent)
    end do
    ec_earth_idx(0) = maxind(ec_earth_pc)
    ec_earth_idx(1) = minind(ec_earth_pc)
    ec_earth_idx(2) = ec_earth_pc(maxind(ec_earth_pc))
    ec_earth_idx(3) = ec_earth_pc(minind(ec_earth_pc))
    delete([/lat_m,lon_m,tmp,msk,obs_trend,clat/])
    ; -- mpi
    m = 6
    lat_m = mpi_trends&lat
    lon_m = mpi_trends&lon
    obs_trend   = linint2_Wrap (lon_o,lat_o,obs_trend_raw,False,lon_m,lat_m,0)
    msk = where(obs_trend .ne. obs_trend@_FillValue,1,obs_trend@_FillValue)
    mpi_trends_masked_obs = mpi_trends * conform_dims(dimsizes(mpi_trends),msk,(/1,2/))
    mpi_pc     = new((/ensmem0(m)/),float)
    replace_ieeenan (obs_trend, mv, 0)
    obs_trend@_FillValue = mv
    obs_trend@missing_value = mv
    do e = 0,ensmem0(m)-1
      tmp = tofloat(mpi_trends_masked_obs(e,:,:))
      replace_ieeenan (tmp, mv, 0)
      tmp@_FillValue = mv
      tmp@missing_value = mv
      tmp!0 = "lat"
      tmp!1 = "lon"
      tmp&lat = mpi_trends&lat
      tmp&lon = mpi_trends&lon
      clat = tofloat(cos(0.01745329*obs_trend&lat))
      clat!0 = "lat"
      clat&lat = obs_trend&lat
      mpi_pc(e)     = pattern_cor(tmp({lat|minlat:maxlat},{lon|minlon:maxlon}),obs_trend({lat|minlat:maxlat},{lon|minlon:maxlon}),clat({lat|minlat:maxlat}),cent)
    end do
    mpi_idx(0) = maxind(mpi_pc)
    mpi_idx(1) = minind(mpi_pc)
    mpi_idx(2) = mpi_pc(maxind(mpi_pc))
    mpi_idx(3) = mpi_pc(minind(mpi_pc))
    delete([/lat_m,lon_m,tmp,msk,obs_trend,clat/])

    print_clock("save...")
    ofile1 = pathout+"olens_mckinnon/"+var+"_"+seas+"_trend_"+start_trend+"-"+ende_trend+"_pattern_cor_"+domain+".nc"
    system("rm -rf "+ofile1)
    o = addfile(ofile1,"c")
    o->pc         = olens_pc

    ofile1 = pathout+"cesm_lens/"+comp+"/"+var+"/"+var+"_"+seas+"_trend_"+start_trend+"-"+ende_trend+"_pattern_cor_"+domain+".nc"
    system("rm -rf "+ofile1)
    o = addfile(ofile1,"c")
    o->pc         = cesm_pc

    ofile1 = pathout+"canesm2_lens/"+comp+"/"+var+"/"+var+"_"+seas+"_trend_"+start_trend+"-"+ende_trend+"_pattern_cor_"+domain+".nc"
    system("rm -rf "+ofile1)
    o = addfile(ofile1,"c")
    o->pc         = canesm2_pc

    ofile1 = pathout+"csiro_mk36_lens/"+comp+"/"+var+"/"+var+"_"+seas+"_trend_"+start_trend+"-"+ende_trend+"_pattern_cor_"+domain+".nc"
    system("rm -rf "+ofile1)
    o = addfile(ofile1,"c")
    o->pc         = csiro_mk36_pc

    ofile1 = pathout+"gfdl_cm3_lens/"+comp+"/"+var+"/"+var+"_"+seas+"_trend_"+start_trend+"-"+ende_trend+"_pattern_cor_"+domain+".nc"
    system("rm -rf "+ofile1)
    o = addfile(ofile1,"c")
    o->pc         = gfdl_cm3_pc

    ofile1 = pathout+"gfdl_esm2m_lens/"+comp+"/"+var+"/"+var+"_"+seas+"_trend_"+start_trend+"-"+ende_trend+"_pattern_cor_"+domain+".nc"
    system("rm -rf "+ofile1)
    o = addfile(ofile1,"c")
    o->pc         = gfdl_esm2m_pc

    ofile1 = pathout+"ec_earth_lens/"+comp+"/"+var+"/"+var+"_"+seas+"_trend_"+start_trend+"-"+ende_trend+"_pattern_cor_"+domain+".nc"
    system("rm -rf "+ofile1)
    o = addfile(ofile1,"c")
    o->pc         = ec_earth_pc

    ofile1 = pathout+"mpi_lens/"+comp+"/"+var+"/"+var+"_"+seas+"_trend_"+start_trend+"-"+ende_trend+"_pattern_cor_"+domain+".nc"
    system("rm -rf "+ofile1)
    o = addfile(ofile1,"c")
    o->pc         = mpi_pc
  end if ; end of "be" if clause







; exit






;========================================================================================
  print_clock("plotting...")


  wks_type = "eps" ; eps png
  if (wks_type.eq."png") then
     wks_type@wkWidth = 1500
     wks_type@wkHeight = 1500
  end if
  if (domain.eq."NA") then
     ofile = "/home/flehner/publication/clivar19_perspective/"+get_script_prefix_name()+"_"+var+"_"+seas+"_"+start_trend+"-"+ende_trend
  else
     ofile = "/home/flehner/publication/clivar19_perspective/"+get_script_prefix_name()+"_"+var+"_"+seas+"_"+start_trend+"-"+ende_trend
  end if
  wks = gsn_open_wks(wks_type,ofile)
  plot = new(9,graphic)

  panres = True

  res = True
  res@mpProjection = "Robinson"
  res@mpFillOn = True
  ; res@gsnAddCyclic = True
  if var .eq. "psl"
    res@mpFillColors = (/-1,-1,-1,-1/)
  else
    res@mpFillColors = (/-1,0,-1,-1/)
  end if
  res@mpGeophysicalLineColor = "gray37"
  res@mpUSStateLineColor = res@mpGeophysicalLineColor
  res@mpNationalLineColor = "gray20"
  res@mpPerimOn    = False
  res@mpGridLatSpacingF =  90            ; change latitude  line spacing
  res@mpGridLonSpacingF = 180.           ; change longitude line spacing
  res@mpGridLineColor   = "transparent"  ; trick ncl into drawing perimeter
  res@mpGridAndLimbOn   = True           ; turn on lat/lon lines
  res@mpOutlineOn = True
  res@mpDataSetName = "Earth..4"
  res@mpGeophysicalLineThicknessF = 1.5
  res@mpPerimDrawOrder = "PostDraw"
  res@mpOutlineDrawOrder = "Draw"
  res@mpFillDrawOrder = "PostDraw"
  res@cnFillDrawOrder = "PreDraw"

  if (domain.eq."NA") then
    res@mpDataBaseVersion     = "LowRes";"MediumRes"
    res@mpOutlineBoundarySets = "National"
    ; res@mpOutlineBoundarySets = "GeophysicalAndUSStates"
    ; res@mpOutlineSpecifiers = (/"Canada : Provinces","Mexico : States"/)
    res@mpLimitMode = "LatLon"
    res@mpMinLonF = minlon
    res@mpMaxLonF = maxlon
    res@mpCenterLonF = (res@mpMinLonF+res@mpMaxLonF)/2.0
    res@mpMinLatF = minlat
    res@mpMaxLatF = maxlat
    res@mpPerimOn    = True
    res@gsnLeftStringOrthogonalPosF = -.96
    res@gsnLeftStringParallelPosF = 0.025
    res@gsnLeftStringFontHeightF = 0.022
    res@gsnLeftStringFontColor = "gray25"
    tw = "gsnLeftString"
    res@gsnRightStringOrthogonalPosF = -.96
    res@gsnRightStringParallelPosF = 0.97
    res@gsnRightStringFontHeightF = 0.022
    res@gsnRightStringFontColor = "gray25"
  else
    res@mpCenterLonF = 210. ; 0
    res@mpFillOn = False
    res@gsnCenterStringFontHeightF = 0.022
    res@gsnCenterStringFontColor = "gray25"
    tw = "gsnCenterString"
  end if

  res@gsnFrame     =  False
  res@gsnDraw = False
  res@cnLineLabelsOn = False
  res@cnFillOn        = True
  res@cnLinesOn       = False
  res@gsnRightString = ""
  res@gsnLeftString = ""
  res@tiMainOn = False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@lbLabelBarOn = False
  res@gsnLeftString = ""
  res@gsnRightString = ""
  ; res@cnLabelBarEndStyle = "IncludeMinMaxLabels"   ; force a label at the end of both boxes


  ; ----------------------------------------------------------------------------
  ; -- Trends --
  delete([/res@cnFillColors,res@cnLevels/])
  if (var.eq."tas" .or. var.eq."psl") then
    gsn_define_colormap(wks,"amwg_blueyellowred")
    cmap = gsn_retrieve_colormap(wks)
    hsv = rgbhsv(cmap)
    hsv(10,1) = .47
    hsv(8,1) = .45    ; 1c+1d
    hsv(9,1) = .30    ; 1d
    hsv(4,:) = (/ (hsv(3,:)+hsv(5,:))/2.0 /)
    cmap2 = hsvrgb(hsv)
    gsn_define_colormap(wks,cmap2)
    res@cnFillColors = (/3,4,5,6,7,8,9,10,11,12,13,14,15,16/)
    ; -- ‘levs’ has to be contour levels (example below) —
    ; levs                =  (/-4,-3,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,3,4/)
    levs                =  (/-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3/)
    res@cnLevels 	   		= levs
  end if
  if (var.eq."pr") then
    gsn_define_colormap(wks,"precip_diff_12lev")
    res@cnFillColors = (/2,3,4,5,6,7,8,9,10,11,12,13,14/)
    ; res@cnFillColors = (/2,3,4,5,6,8,10,11,12,13,14/)
    ; -- ‘levs’ has to be contour levels (example below) —
    if domain .eq. "NA"
      levs                =  (/-.5,-.3,-.2,-.15,-.1,-.05,0.05,.1,.15,.2,.3,.5/)
    else
      levs                =  (/-2,-1.5,-1,-0.5,-0.2,-0.1,0.1,0.2,0.5,1,1.5,2/)
    end if
    res@cnLevels 	   		= levs
  end if
  ; -- consider an expantion of the array at both ends by the first and last increment, respectively --
  dl1 = dimsizes(levs)-1
  dl2 = dl1+2
  levs1 = new((/dimsizes(levs)+2/),float)
  levs1(0) = levs(0)-levs(1) + levs(0)
  levs1(dl2) = levs(dl1)-levs(dl1-1) + levs(dl1)
  levs1(1:dl2-1) = levs
  ; -- move levs so that first value is zero --
  levs2 = levs1+abs(levs1(0))
  ; -- fill in fractions --
  levs1_3 = levs2
  levs1_3(1:dl2) = levs2(1:dl2)/levs2(dl2)

  ; -- contour resources
  resc                    = res
  resc@cnFillOn           = False
  resc@cnLinesOn          = True         ; Turn lines off
  resc@cnMonoLineThickness= True
  resc@cnLineColor        = "black"
  delete(resc@cnLevels)
  resc@cnLevels = (/5/)
  resc@cnLineThicknessF   = 2
  resc@gsnContourPosLineDashPattern = 0
  resc@cnInfoLabelOn      = False
  ; resc@cnLevels           = fspan(-1,1,3)
  resc@$tw$                 = ""

  ; -- stippling resources
  ress                    = res
  ress@cnInfoLabelOn      = False
  delete(ress@cnLevels)
  ; ress@cnFillPalette        = "default"
  ress@cnMonoFillColor      = True
  ress@cnFillColor          = "gray26"
  ress@cnLevels             = (/.1/)
  ress@cnMonoFillPattern    = False            ; want multiple patterns
  ress@cnFillPatterns       = (/-1,17/) ; the patterns
  ress@cnMonoFillScale      = False            ; want different densities
  ress@cnFillDotSizeF       = 0.002
  if domain .eq. "global"
    ress@cnFillScales         = (/.4,.4/) ; the densities
  else
    ress@cnFillScales         = (/.8,.8/) ; the densities
  end if

  ; -- hatching resources
  resh                    = ress
  resh@cnInfoLabelOn      = False
  delete(resh@cnLevels)
  ; resh@cnFillPalette        = "default"
  resh@cnMonoFillColor      = True
  resh@cnFillColor          = "gray26"
  resh@cnLevels             = (/.1/)
  resh@cnMonoFillPattern    = False            ; want multiple patterns
  resh@cnFillPatterns       = (/-1,8/) ; the patterns
  resh@cnMonoFillScale      = False            ; want different densities
  resh@cnFillDotSizeF       = 0.002
  if domain .eq. "global"
    resh@cnFillScales         = (/.4,.4/) ; the densities
  else
    resh@cnFillScales         = (/.8,.8/) ; the densities
  end if

  if (var.eq."tas") then
     panres@lbTitleString = "Temperature (~S~o~N~C "+(ende_trend-start_trend+1)+"yr~S~-1~N~)"
  end if
  if (var.eq."pr") then
     panres@lbTitleString = "Precipitation (mm day~S~-1~N~ "+(ende_trend-start_trend+1)+"yr~S~-1~N~)"
  end if
  if (var.eq."psl") then
     panres@lbTitleString = "Sea level pressure (hPa~N~ "+(ende_trend-start_trend+1)+"yr~S~-1~N~)"
  end if

  res@$tw$ = "Observations ("+obs_name+")"
  res@gsnRightString = decimalPlaces(obs_trend_values,2,True)
  plot(2) = gsn_csm_contour_map(wks,obs_trend_raw,res)

  m = 0
  res@$tw$ = model_nametags2(m)+" EM ("+ensmem0(m)+")"
  res@gsnRightString = decimalPlaces(avg(cesm_idx(4)),2,True)
  plot(0) = gsn_csm_contour_map(wks,dim_avg_n_Wrap(cesm_trends,0),res)
  ; -- stippling/hatching where obs trend is outside all model trends
  obs_trend_regridded   = linint2_Wrap (lon_o,lat_o,obs_trend_raw,False,cesm_trends&lon,cesm_trends&lat,0)
  too_low  = where(obs_trend_regridded .gt. dim_max_n(cesm_trends,0),1,0)
  too_high   = where(obs_trend_regridded .lt. dim_min_n(cesm_trends,0),1,0)
  copy_VarCoords(cesm_trends(0,:,:),too_high)
  copy_VarCoords(cesm_trends(0,:,:),too_low)
  ; ; -- plot hatching/stippling:
  ; plot_st             = gsn_csm_contour(wks,too_high,resh)
  ; overlay(plot(0),plot_st)
  ; plot_st             = gsn_csm_contour(wks,too_low,ress)
  ; overlay(plot(0),plot_st)
  delete([/too_low,too_high,obs_trend_regridded/])
  ; ---
  res@$tw$ = model_nametags2(m)+" #"+toint(cesm_idx(0)+1)
  res@gsnRightString = decimalPlaces(cesm_idx(0+2),2,True)
  plot(3) = gsn_csm_contour_map(wks,cesm_trends(toint(cesm_idx(0)),:,:),res)
  res@$tw$ = model_nametags2(m)+" #"+toint(cesm_idx(1)+1)
  res@gsnRightString = decimalPlaces(cesm_idx(1+2),2,True)
  plot(6) = gsn_csm_contour_map(wks,cesm_trends(toint(cesm_idx(1)),:,:),res)

  m = 6
  res@$tw$ = model_nametags2(m)+" EM ("+ensmem0(m)+")"
  res@gsnRightString = decimalPlaces(avg(mpi_idx(4)),2,True)
  plot(1) = gsn_csm_contour_map(wks,dim_avg_n_Wrap(mpi_trends,0),res)
  ; -- stippling/hatching where obs trend is outside all model trends
  obs_trend_regridded   = linint2_Wrap (lon_o,lat_o,obs_trend_raw,False,mpi_trends&lon,mpi_trends&lat,0)
  too_low  = where(obs_trend_regridded .gt. dim_max_n(mpi_trends,0),1,0)
  too_high   = where(obs_trend_regridded .lt. dim_min_n(mpi_trends,0),1,0)
  copy_VarCoords(mpi_trends(0,:,:),too_high)
  copy_VarCoords(mpi_trends(0,:,:),too_low)
  delete([/too_low,too_high,obs_trend_regridded/])
  ; ---
  res@$tw$ = model_nametags2(m)+" #"+toint(mpi_idx(0)+1)
  res@gsnRightString = decimalPlaces(mpi_idx(0+2),2,True)
  plot(4) = gsn_csm_contour_map(wks,mpi_trends(toint(mpi_idx(0)),:,:),res)
  res@$tw$ = model_nametags2(m)+" #"+toint(mpi_idx(1)+1)
  res@gsnRightString = decimalPlaces(mpi_idx(1+2),2,True)
  plot(7) = gsn_csm_contour_map(wks,mpi_trends(toint(mpi_idx(1)),:,:),res)

  ; res@$tw$ = "OLENS #"+toint(olens_idx(0)+1)+" + CESM forced"
  ; res@gsnRightString = decimalPlaces(olens_idx(0+2),2,True)
  ; plot(4) = gsn_csm_contour_map(wks,olens_trends(toint(olens_idx(0)),:,:),res)
  ; res@$tw$ = "OLENS #"+toint(olens_idx(1)+1)+" + CESM forced"
  ; res@gsnRightString = decimalPlaces(olens_idx(1+2),2,True)
  ; plot(5) = gsn_csm_contour_map(wks,olens_trends(toint(olens_idx(1)),:,:),res)
  ;
  ; res@$tw$ = "Observations ("+obs_name+")"
  ; res@gsnRightString = decimalPlaces(obs_trend_values,2,True)
  ; plot(9) = gsn_csm_contour_map(wks,obs_trend_raw,res)
  ; res@$tw$ = "OLENS #"+toint(olens_idx(0)+1)+" + MPI forced"
  ; res@gsnRightString = decimalPlaces(olens_mpi_idx(0+2),2,True)
  ; plot(10) = gsn_csm_contour_map(wks,olens_trends_mpi(toint(olens_mpi_idx(0)),:,:),res)
  ; res@$tw$ = "OLENS #"+toint(olens_idx(1)+1)+" + MPI forced"
  ; res@gsnRightString = decimalPlaces(olens_mpi_idx(1+2),2,True)
  ; plot(11) = gsn_csm_contour_map(wks,olens_trends_mpi(toint(olens_mpi_idx(1)),:,:),res)

  ;-- polymarker resources
  pmres                        =  True
  pmres@gsMarkerColor          = "black"       ;-- marker color
  pmres@gsMarkerIndex          = 4          ; marker style, 4=circle
  pmres@gsMarkerSizeF          =  10.        ;-- set size of marker
  pmres@gsLineThicknessF       =  3.          ;-- marker line thickness
  dum = new(1,graphic)
  dum = gsn_add_polymarker(wks, plot(2), -96.8+360, 32.7, pmres)
  ; txres                        = True
  ; txres@txFontHeightF          = .022             ; Set the font height
  ; text = gsn_add_text(wks,plot(2),"Dallas, TX",-84+360, 36 ,txres)


  panres@gsnMaximize = False
  panres@gsnPaperOrientation = "portrait"
  panres@gsnPanelLabelBar = True
  panres@gsnPanelBottom = .10
  panres@pmLabelBarWidthF = 0.62
  panres@pmLabelBarHeightF = 0.062
  panres@lbBoxEndCapStyle = "TriangleBothEnds"
  panres@cnLabelBarEndStyle = "IncludeMinMaxLabels"   ; force a label at the end of both boxes
  panres@lbBoxLineColor = "gray30"
  panres@lbTitlePosition = "bottom"
  panres@lbTitleOn = True

  panres@lbLabelFontHeightF = 0.0145
  panres@lbTitleFontHeightF = panres@lbLabelFontHeightF
  panres@pmLabelBarOrthogonalPosF = -0.013
;  panres@amJust   = "TopLeft"
;  panres@amParallelPosF = -0.475
;  panres@amOrthogonalPosF = 0.33
;  panres@gsnPanelFigureStringsFontHeightF = 0.0145
;  panres@gsnPanelFigureStringsPerimOn = False
;  panres@gsnPanelFigureStringsBackgroundFillColor = "transparent"  ;"white" (/1,1,1,0.8/)
  panres@gsnPanelMainFontHeightF = 0.016
  panres@gsnPanelYWhiteSpacePercent = 10.;2.0
  panres@gsnPanelXWhiteSpacePercent = 5.
  panres@gsnPanelMainString = changeCaseChar(var)+" Trends "+seas+" "+start_trend+"-"+ende_trend
  ; panres@lbBoxEndCapStyle   = "TriangleHighEnd"
  panres@cnLabelBarEndStyle = "IncludeMinMaxLabels"   ; force a label at the end of both boxes
  panres@cnExplicitLabelBarLabelsOn = True
  panres@lbBoxSizing        = "ExplicitSizing"
  panres@lbBoxFractions     = levs1_3
  panres@gsnPanelRowSpec    = True
  gsn_panel(wks,plot,(/3,3,3/),panres)
  delete([/plot/])
  delete([/levs,levs1,levs2,levs1_3,res@cnLevels,res@cnFillColors,panres@lbTitleString,panres@lbBoxFractions/])

  ; system("convert -density 144 -trim +repage -border 8 -bordercolor white -flatten "+ofile+".000001.png "+ofile+"_"+tword(gg-1)+".png")
  delete(wks)
  ; system("convert -density 150x150 -crop 0x0 -trim +repage "+ofile+".eps "+ofile+".jpg")
  system("convert -density 250x250 -crop 0x0 -trim +repage "+ofile+".eps "+ofile+".jpg")

  ; system("convert "+ofile+".eps "+ofile+".jpg")
  ; system("rm -rf "+ofile+"*.00000*.png")

  ; system("gv "+ofile+".eps &")

  print_clock("...done")

end do
end do

end
