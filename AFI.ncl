begin

models = (/"ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5", "CMCC-ESM2", \
           "CNRM-CM6-1", "CNRM-ESM2-1", "EC-Earth3", "EC-Earth3-Veg-LR", "FGOALS-g3", \
           "GFDL-ESM4", "GISS-E2-1-G", "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", \
           "KACE-1-0-G", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", \
           "MRI-ESM2-0", "NorESM2-LM", "NorESM2-MM", "TaiESM1", "UKESM1-0-LL"/)

variants = (/"r1i1p1f1_gn", "r1i1p1f1_gn", "r1i1p1f1_gn", "r1i1p1f2_gn", "r1i1p1f1_gn", \
             "r1i1p1f2_gr", "r1i1p1f2_gr", "r1i1p1f1_gr", "r1i1p1f1_gr", "r1i1p1f1_gn", \
             "r1i1p1f1_gr", "r1i1p1f2_gn", "r1i1p1f1_gr", "r1i1p1f1_gr", "r1i1p1f1_gr", \
             "r1i1p1f1_gr", "r1i1p1f2_gn", "r1i1p1f1_gn", "r1i1p1f1_gn", "r1i1p1f1_gn", \
             "r1i1p1f1_gr", "r1i1p1f1_gn", "r1i1p1f1_gn", "r1i1p1f1_gn", "r1i1p1f2_gn"/)

scenarios = (/"ssp126", "ssp245", "ssp370", "ssp585"/)

f_ref = addfile("tas_day_"+models(0)+"_"+scenarios(0)+"_"+variants(0)+"_2015.nc", "r")
lat = f_ref->lat
lon = f_ref->lon
delete(f_ref)

nmodels = dimsizes(models)
nlat = dimsizes(lat)
nlon = dimsizes(lon)

do s = 0, dimsizes(scenarios)-1
    scenario = scenarios(s)
    
    output_dir = "./"+str_upper(scenario)+"/"
    system("mkdir -p " + output_dir)
    
    do y = 1950, 2100
        
        ensemble_sum = new((/nlat,nlon/), float, "No_FillValue")
        ensemble_sum = 0.0
        model_count = 0
        
        do m = 0, nmodels-1
            model = models(m)
            variant = variants(m)
            
            filename = "tas_day_"+model+"_"+scenario+"_"+variant+"_"+y+".nc"
            
            FI_3 = new((/nlat,nlon/), float, "No_FillValue")
            FI_3 = 0.0
            
            a = addfile(filename, "r")
            time = a->time
            time1 = cd_calendar(time, 0)
            IND = ind(time1(:,1).eq.7)
            
            tas1 = (a->tas(IND(0):,:,:)) - 273.15
            tas2 = where(tas1.gt.0, tas1, 0)
            tas3 = dim_sum_n(tas2, 0)
            FI_3 = FI_3 + tas3
            
            delete(tas1)
            delete(tas2)
            delete(tas3)
            delete(time)
            delete(time1)
            delete(IND)
            
            b = addfile(filename, "r")
            time_b = b->time
            time1_b = cd_calendar(time_b, 0)
            IND_b = ind(time1_b(:,1).eq.7)
            
            if (IND_b(0) .gt. 0) then
                tas4 = (b->tas(0:IND_b(0)-1,:,:)) - 273.15
                tas5 = where(tas4.gt.0, tas4, 0)
                tas6 = dim_sum_n(tas5, 0)
                FI_3 = FI_3 + tas6
                
                delete(tas4)
                delete(tas5)
                delete(tas6)
            end if
            
            delete(time_b)
            delete(time1_b)
            delete(IND_b)
            
            FI_4 = abs(FI_3)
            FI_4!0 = "lat"
            FI_4!1 = "lon"
            FI_4&lat = lat
            FI_4&lon = lon
            
            output_filename = output_dir + "TI_" + model + "_" + y + ".nc"
            fout = addfile(output_filename, "c")
            fout->TI = FI_4
            fout->lat = lat
            fout->lon = lon
            
            ensemble_sum = ensemble_sum + FI_4
            model_count = model_count + 1
            
            delete(FI_3)
            delete(FI_4)
            delete(a)
            delete(b)
            delete(fout)
            
        end do
        
        ensemble_mean = ensemble_sum / model_count
        ensemble_mean!0 = "lat"
        ensemble_mean!1 = "lon"
        ensemble_mean&lat = lat
        ensemble_mean&lon = lon
        
        ensemble_filename = output_dir + "TI_ensemble_mean_" + y + ".nc"
        fout_ens = addfile(ensemble_filename, "c")
        fout_ens->TI = ensemble_mean
        fout_ens->lat = lat
        fout_ens->lon = lon
        
        delete(ensemble_sum)
        delete(ensemble_mean)
        delete(fout_ens)
        
    end do
    
end do

end