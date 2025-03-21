; This NCL script is to regrid WRFOUT directly using the file that contains the weights
; The weight file is created by ESMF_upscale_weight_to_MRMS_V3.ncl
; this file is created because early regridding script in MATLAB does no longer work (-_-)
; created by C. Bayu Risanto, S.J. (2 July 2024)

;**************************************************************************************
;load "$NCARG/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
;---Data file containing source grid
    target = "20210901_06Z_mc08_40mem_0035_6cyc_noassim"
    dirin = "/st6/bayu/ARIZONA_PRO/DART_deterministic/" + target + "/"
    files = systemfunc("ls -d " + dirin + "wrfout_d02* | cut -f7 -d/")
    ns = dimsizes(files)

    regridded = new((/ns,651,612/),float)
    times = new((/ns/),string)
    print(dimsizes(times)) 
    
    do i = 0,ns-1
      print(dirin + files(i))
      src_file = dirin + files(i)    
      sfile    = addfile(src_file,"r")
      rain = sfile->RAINNC(0,:,:)
   
;---Regrid here
      rainnc_reg2 = ESMF_regrid_with_weights(rain,"upscale_weight_WRF2MRMS.AZ+.nc",False)
      ;print(dimsizes(rainnc_reg2))

;---stack here
      regridded(i,:,:) = rainnc_reg2 

;---get temporal var
      times(i) = sfile->Times(0,:)
      print(times(i))

    end do

;--- safe file
    dirout = "/br1/castrogroup/bayu/UA-WRF/ARIZONA_PRO/WRF_EXTR/regRAINNC/"+ target +"/"
    new_file_name = "RAINNC_d02_reg2MRMS.AZ+.nc"
    fout = dirout + new_file_name
    system("/bin/rm -f " + fout)   ; remove any pre-existing file
    
    o = addfile(fout ,"c")
    o->RAINNC_reg2 = regridded
    o->TIMES = stringtochar(times)

end
