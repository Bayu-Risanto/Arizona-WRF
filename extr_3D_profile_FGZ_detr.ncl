;*******************************************************************
;This script is to get 3 dimensional QVAPOR, T, PH, PHB, P, PB variables from FLAGSTAFF
;WRF simulation DETERMINISTIC. --- STORE ALL ENS in ONE FILE OUT 
;Created by C. Bayu Risanto, S.J. (9 July 2024)
;*******************************************************************
time = 0
opt            = True
opt@latlon     = True
opt@linecoords = True    ; returns lat/lon coordinates on line as attributes "lats","lons"

dir_det = "20210815_06Z_40mem_0035_6cyc_vertisundef_assim"
files = systemfunc("ls " + "/st6/bayu/ARIZONA_PRO/DART_deterministic/" + dir_det + "/" + "wrfout_d02_*")

fout_pre = dir_det

WRF = "/st6/bayu/ARIZONA_PRO/DART_deterministic/20210815_06Z_40mem_06Z_0035_6cyc_assim/wrfout_d02_2021-08-15_18:00:00/"
dummy = addfile(WRF,"r")

nf = dimsizes(files)
dirout = "/br1/castrogroup/bayu/UA-WRF/ARIZONA_PRO/WRF_EXTR/WRF-RAOB/"

;************************************************************************
;; start here
;************************************************************************

ii = 0
fgz_PRESS = new((/19,43/),float)
fgz_TEMPS = new((/19,43/),float)
fgz_TDEWS = new((/19,43/),float)
fgz_QVAPS = new((/19,43/),float)
fgz_ZALTS = new((/19,43/),float)

do while(ii.le.nf-1)

    wcday = systemfunc("date")
    a = addfile(files(ii),"r") 

    print(""+files(ii))

    opt@file_handle = a
    
    lat2d = a->XLAT(0,:,:)
    lon2d = a->XLONG(0,:,:)
    temp    = wrf_user_getvar(a, "tc",0)              ; temperature C
    t2      = wrf_user_getvar(a, "T2",0) - 273.15     ; temperature at 2 m [c]
    tdew    = wrf_user_getvar(a, "td",0)              ; dew point temp C
    tdew2   = wrf_user_getvar(a, "td2",0)             ; dew point temp at 2 m [C]
    zalt    = wrf_user_getvar(a, "height",0)          ; grid point height
    pres    = wrf_user_getvar(a, "pres",0)            ; pressure [Pa]
    psfc    = wrf_user_getvar(a, "PSFC",0)            ; pressure [Pa]
    qvap    = wrf_user_getvar(a, "QVAPOR",0)          ; specific humidity
    q2      = wrf_user_getvar(a, "Q2",0)              ; mixing ratio [kg/kg]

    fgz = wrf_user_ll_to_xy(dummy,-111.82,35.23,opt)

    fgz_q2    = q2(fgz(1),fgz(0))
    fgz_qvap = qvap(:,fgz(1),fgz(0))
    fgz_psfc = psfc(fgz(1),fgz(0))
    fgz_pres = pres(:,fgz(1),fgz(0))
    fgz_zalt = zalt(:,fgz(1),fgz(0))
    fgz_t2   = t2(fgz(1),fgz(0))
    fgz_temp = temp(:,fgz(1),fgz(0))
    fgz_td2  = tdew2(fgz(1),fgz(0))
    fgz_td   = tdew(:,fgz(1),fgz(0))

    fgz_press = array_append_record(fgz_psfc,fgz_pres,0)
    fgz_tempt = array_append_record(fgz_t2,fgz_temp,0)
    fgz_tdew  = array_append_record(fgz_td2,fgz_td,0)
    fgz_qvaps = array_append_record(fgz_q2,fgz_qvap,0)
    fgz_zalts = array_append_record(0.0,fgz_zalt,0)

    fgz_PRESS(ii,:) = fgz_press
    fgz_TEMPS(ii,:) = fgz_tempt
    fgz_TDEWS(ii,:) = fgz_tdew
    fgz_QVAPS(ii,:) = fgz_qvaps
    fgz_ZALTS(ii,:) = fgz_zalts

    delete(fgz_q2)
    delete(fgz_qvap)
    delete(fgz_psfc)
    delete(fgz_pres)
    delete(fgz_zalt)
    delete(fgz_t2)
    delete(fgz_temp)
    delete(fgz_td)
    delete(fgz_td2)
    delete(fgz_press)
    delete(fgz_tempt)
    delete(fgz_qvaps)
    delete(fgz_zalts)

    ii = ii + 1
end do

fout1a = dirout + "/FGZ-3D-" + fout_pre + ".nc"
system("/bin/rm -f " + fout1a)   ; remove any pre-existing file

print(fout1a+"")

nc1a = addfile(fout1a ,"c")  ; open output netCDF file
nc1a->FGZ_QVAP  = fgz_QVAPS
nc1a->FGZ_PRES  = fgz_PRESS
nc1a->FGZ_TEMPT = fgz_TEMPS
nc1a->FGZ_TDEWS = fgz_TDEWS
nc1a->FGZ_ZALTS = fgz_ZALTS

wallClockElapseTime(wcday, "one day", 0)    
delete(nc1a)

  
exit

