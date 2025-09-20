
begin

models = (/"ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5", "CMCC-ESM2", \
           "CNRM-CM6-1", "CNRM-ESM2-1", "EC-Earth3", "EC-Earth3-Veg-LR", "FGOALS-g3", \
           "GFDL-ESM4", "GISS-E2-1-G", "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", \
           "KACE-1-0-G", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", \
           "MRI-ESM2-0", "NorESM2-LM", "NorESM2-MM", "TaiESM1", "UKESM1-0-LL"/)

scenarios = (/"ssp126", "ssp245", "ssp370", "ssp585"/)

mask1 = addfile("GTOPOP_v4.nc", "r")
mask2 = mask1->elevation

do s = 0, dimsizes(scenarios)-1
    scenario = scenarios(s)
    scenario_upper = str_upper(scenario)
    
    data_dir = "./"+scenario_upper+"/"
    output_dir = "./Frozen_"+scenario_upper+"/"
    system("mkdir -p " + output_dir)
    
    do y = 1950, 2099, 1
        
        do m = 0, dimsizes(models)-1
            model = models(m)
            
            fi_filename = data_dir + "FI_" + model + "_" + y + ".nc"
            ti_filename = data_dir + "TI_" + model + "_" + y + ".nc"
            
            a = addfile(fi_filename, "r")
            FI = a->FI
            b = addfile(ti_filename, "r")
            TI = b->TI
            
            lat = a->lat
            lon = a->lon
            
            nlat = dimsizes(lat)
            nlon = dimsizes(lon)
            
            frozen_TP = sqrt(FI)/(sqrt(FI)+sqrt(TI))
            
            frozen_TP@_FillValue = 9999
            frozen_TP = where(frozen_TP.ge.0.55.and.frozen_TP.le.1, 3, frozen_TP)
            frozen_TP = where(frozen_TP.lt.0.55.and.frozen_TP.ge.0.01, 2, frozen_TP)
            frozen_TP = where(frozen_TP.lt.0.01, 0, frozen_TP)
            frozen_TP@_FillValue = -9999
            
            frozen_TP!0 = "lat"
            frozen_TP&lat = lat
            frozen_TP!1 = "lon"
            frozen_TP&lon = lon
            
            frozen_TP_outside = sqrt(FI)/(sqrt(FI)+sqrt(TI))
            
            frozen_TP_outside@_FillValue = 9999
            frozen_TP_outside = where(frozen_TP_outside.ge.0.53.and.frozen_TP_outside.le.1, 3, frozen_TP_outside)
            frozen_TP_outside = where(frozen_TP_outside.lt.0.53.and.frozen_TP_outside.ge.0.1, 2, frozen_TP_outside)
            frozen_TP_outside = where(frozen_TP_outside.lt.0.1.and.frozen_TP_outside.ge.0.01, 1, frozen_TP_outside)
            frozen_TP_outside = where(frozen_TP_outside.lt.0.01, 0, frozen_TP_outside)
            frozen_TP_outside@_FillValue = -9999
            
            frozen_TP_outside!0 = "lat"
            frozen_TP_outside&lat = lat
            frozen_TP_outside!1 = "lon"
            frozen_TP_outside&lon = lon
            
            SFD = where(mask2.eq.1, frozen_TP, frozen_TP_outside)
            SFD@_FillValue = -9999
            
            SFD!0 = "lat"
            SFD&lat = lat
            SFD!1 = "lon"
            SFD&lon = lon
            
            output_filename = output_dir + model + "_Frozen_" + scenario + "_" + y + ".nc"
            fout = addfile(output_filename, "c")
            fout->frozen = SFD
            fout->lat = lat
            fout->lon = lon
            
            delete(FI)
            delete(TI)
            delete(frozen_TP)
            delete(frozen_TP_outside)
            delete(SFD)
            delete(a)
            delete(b)
            delete(fout)
        end do
        
        fi_ensemble_filename = data_dir + "FI_ensemble_mean_" + y + ".nc"
        ti_ensemble_filename = data_dir + "TI_ensemble_mean_" + y + ".nc"
        
        a_ens = addfile(fi_ensemble_filename, "r")
        FI_ens = a_ens->FI
        b_ens = addfile(ti_ensemble_filename, "r")
        TI_ens = b_ens->TI
        
        lat = a_ens->lat
        lon = a_ens->lon
        
        frozen_TP_ens = sqrt(FI_ens)/(sqrt(FI_ens)+sqrt(TI_ens))
        
        frozen_TP_ens@_FillValue = 9999
        frozen_TP_ens = where(frozen_TP_ens.ge.0.55.and.frozen_TP_ens.le.1, 3, frozen_TP_ens)
        frozen_TP_ens = where(frozen_TP_ens.lt.0.55.and.frozen_TP_ens.ge.0.01, 2, frozen_TP_ens)
        frozen_TP_ens = where(frozen_TP_ens.lt.0.01, 0, frozen_TP_ens)
        frozen_TP_ens@_FillValue = -9999
        
        frozen_TP_ens!0 = "lat"
        frozen_TP_ens&lat = lat
        frozen_TP_ens!1 = "lon"
        frozen_TP_ens&lon = lon
        
        frozen_TP_outside_ens = sqrt(FI_ens)/(sqrt(FI_ens)+sqrt(TI_ens))
        
        frozen_TP_outside_ens@_FillValue = 9999
        frozen_TP_outside_ens = where(frozen_TP_outside_ens.ge.0.53.and.frozen_TP_outside_ens.le.1, 3, frozen_TP_outside_ens)
        frozen_TP_outside_ens = where(frozen_TP_outside_ens.lt.0.53.and.frozen_TP_outside_ens.ge.0.1, 2, frozen_TP_outside_ens)
        frozen_TP_outside_ens = where(frozen_TP_outside_ens.lt.0.1.and.frozen_TP_outside_ens.ge.0.01, 1, frozen_TP_outside_ens)
        frozen_TP_outside_ens = where(frozen_TP_outside_ens.lt.0.01, 0, frozen_TP_outside_ens)
        frozen_TP_outside_ens@_FillValue = -9999
        
        frozen_TP_outside_ens!0 = "lat"
        frozen_TP_outside_ens&lat = lat
        frozen_TP_outside_ens!1 = "lon"
        frozen_TP_outside_ens&lon = lon
        
        SFD_ens = where(mask2.eq.1, frozen_TP_ens, frozen_TP_outside_ens)
        SFD_ens@_FillValue = -9999
        
        SFD_ens!0 = "lat"
        SFD_ens&lat = lat
        SFD_ens!1 = "lon"
        SFD_ens&lon = lon
        
        ensemble_output_filename = output_dir + "Ensemble_Mean_Frozen_" + scenario + "_" + y + ".nc"
        fout_ens = addfile(ensemble_output_filename, "c")
        fout_ens->frozen = SFD_ens
        fout_ens->lat = lat
        fout_ens->lon = lon
        
        delete(FI_ens)
        delete(TI_ens)
        delete(frozen_TP_ens)
        delete(frozen_TP_outside_ens)
        delete(SFD_ens)
        delete(a_ens)
        delete(b_ens)
        delete(fout_ens)
        
    end do
    
end do

end