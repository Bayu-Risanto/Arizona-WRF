;*******************************************************************
;This script is to get 3 dimensional QVAPOR, T, PH, PHB, P, PB variables from TUCSON
;WRF simulation DETERMINISTIC. --- STORE ALL ENS in ONE FILE OUT 
;Created by C. Bayu Risanto, S.J. (2 July 2024)
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
sa46_PRESS = new((/19,43/),float)
sa46_TEMPS = new((/19,43/),float)
sa46_TDEWS = new((/19,43/),float)
sa46_QVAPS = new((/19,43/),float)
sa46_ZALTS = new((/19,43/),float)

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

    sa46 = wrf_user_ll_to_xy(dummy,-110.953926,32.229839,opt)

    sa46_q2   = q2(sa46(1),sa46(0))
    sa46_qvap = qvap(:,sa46(1),sa46(0))
    sa46_psfc = psfc(sa46(1),sa46(0))
    sa46_pres = pres(:,sa46(1),sa46(0))
    sa46_zalt = zalt(:,sa46(1),sa46(0))
    sa46_t2   = t2(sa46(1),sa46(0))
    sa46_temp = temp(:,sa46(1),sa46(0))
    sa46_td2  = tdew2(sa46(1),sa46(0))
    sa46_td   = tdew(:,sa46(1),sa46(0))

    sa46_press = array_append_record(sa46_psfc,sa46_pres,0)
    sa46_tempt = array_append_record(sa46_t2,sa46_temp,0)
    sa46_tdew  = array_append_record(sa46_td2,sa46_td,0)
    sa46_qvaps = array_append_record(sa46_q2,sa46_qvap,0)
    sa46_zalts = array_append_record(0.0,sa46_zalt,0)

    sa46_PRESS(ii,:) = sa46_press
    sa46_TEMPS(ii,:) = sa46_tempt
    sa46_TDEWS(ii,:) = sa46_tdew
    sa46_QVAPS(ii,:) = sa46_qvaps
    sa46_ZALTS(ii,:) = sa46_zalts

    delete(sa46_q2)
    delete(sa46_qvap)
    delete(sa46_psfc)
    delete(sa46_pres)
    delete(sa46_zalt)
    delete(sa46_t2)
    delete(sa46_temp)
    delete(sa46_td)
    delete(sa46_td2)
    delete(sa46_press)
    delete(sa46_tempt)
    delete(sa46_qvaps)
    delete(sa46_zalts)

    ii = ii + 1
end do

fout1a = dirout + "/TUS-3D-" + fout_pre + ".nc"
system("/bin/rm -f " + fout1a)   ; remove any pre-existing file

print(fout1a+"")

nc1a = addfile(fout1a ,"c")  ; open output netCDF file
nc1a->SA46_QVAP  = sa46_QVAPS
nc1a->SA46_PRES  = sa46_PRESS
nc1a->SA46_TEMPT = sa46_TEMPS
nc1a->SA46_TDEWS = sa46_TDEWS
nc1a->SA46_ZALTS = sa46_ZALTS

wallClockElapseTime(wcday, "one day", 0)    
delete(nc1a)

  
exit

