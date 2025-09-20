
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"

begin

; 第一部分：生成不同阈值的冻土分布文件
FillValue = -9999
f1 = 0
FI_m = new((/17,600,1440/), float, -9999)
TI_m = new((/17,600,1440/), float, -9999)

do f = 2000, 2016, 1
    a = addfile("FI_ensemble_mean_"+f+".nc", "r")
    FI = a->FI
    b = addfile("TI_ensemble_mean_"+f+".nc", "r")
    TI = b->TI
    
    FI_m(f1,:,:) = FI
    TI_m(f1,:,:) = TI
    f1 = f1 + 1
end do

lat = a->lat
lon = a->lon

nlat = dimsizes(lat)
nlon = dimsizes(lon)

FI_m1 = dim_avg_n_Wrap(FI_m, 0)
TI_m1 = dim_avg_n_Wrap(TI_m, 0)

N0 = 0.01
N = tofloat(N0)

do n = 0, 99
    frozen = sqrt(FI_m1)/(sqrt(FI_m1)+sqrt(TI_m1))
    
    frozen@_FillValue = 9999
    frozen = where(frozen.ge.N.and.frozen.le.1, 3, frozen)
    frozen = where(frozen.lt.N, 0, frozen)
    frozen@_FillValue = -9999
    
    mask1 = addfile("GTOPOP_v4.nc", "r")
    mask2 = mask1->elevation
    frozen@_FillValue = -9999
    frozen = where(mask2.eq.1, frozen, frozen@_FillValue)
    
    frozen!0 = "lat"
    frozen&lat = lat
    frozen!1 = "lon"
    frozen&lon = lon
    
    fout = addfile("Frozen_ground_N="+N+".nc", "c")
    fout->frozen = frozen
    
    N = N + 0.01
    delete(frozen)
end do

; 第二部分：Ft1分析 (Ran_permafrost_area_observe.nc)
print("开始Ft1分析...")
K1 = new(100, float, 9999)
POD1 = new(100, float, 9999)
N11_1 = new(100, float, 9999)

N1 = tofloat(0.01)
ij = 0

do n = 0, 99
    a1 = 0
    a2 = 0
    a3 = 0
    a4 = 0
    b1 = 0
    b2 = 0
    b3 = 0
    b4 = 0
    a = 0
    
    in = addfile("Frozen_ground_N="+N1+".nc", "r")
    frozen = in->frozen
    lat_temp = in->lat
    lon_temp = in->lon
    obser = addfile("Ran_permafrost_area_observe.nc", "r")
    obser_frozen = obser->frozen
    obser_lat = obser->lat
    obser_lon = obser->lon
    
    nlat = dimsizes(lat_temp)
    nlon = dimsizes(lon_temp)
    
    N_count = 0
    
    do i = 0, nlat-1, 1
        do j = 0, nlon-1, 1
            pre = frozen(i,j)
            if(.not.ismissing(pre)) then
                longitude = lon_temp(j)
                latitude = lat_temp(i)
                
                index_lat = closest_val(latitude, obser_lat)
                index_lon = closest_val(longitude, obser_lon)
                
                rain = obser_frozen(index_lat, index_lon)
                if(.not.ismissing(rain)) then
                    oo = pre
                    mm = rain
                    mm_oo = mm - oo
                    
                    if(mm_oo.eq.0)
                        a = a + 1
                    end if
                    
                    if(oo.eq.0)
                        a1 = a1 + 1
                    end if
                    
                    if(oo.eq.3)
                        a4 = a4 + 1
                    end if
                    
                    if(mm.eq.0)
                        b1 = b1 + 1
                    end if
                    
                    if(mm.eq.3)
                        b4 = b4 + 1
                    end if
                    
                    N_count = N_count + 1
                end if
            end if
        end do
    end do
    
    a_1 = tofloat(a)
    N_1 = tofloat(N_count)
    a1_1 = tofloat(a1)
    a4_1 = tofloat(a4)
    b1_1 = tofloat(b1)
    b4_1 = tofloat(b4)
    
    p0 = a_1/N_1
    pe = (a1_1*b1_1+a4_1*b4_1)/(N_1*N_1)
    
    k = (p0-pe)/(1-pe)
    
    K1(ij) = k
    POD1(ij) = p0
    N11_1(ij) = N1
    
    ij = ij + 1
    N1 = N1 + 0.01
end do

; 第三部分：Ft2分析 (frozen_2000_code_3_2_0.nc)
print("开始Ft2分析...")
K2 = new(50, float, 9999)
POD2 = new(50, float, 9999)
N11_2 = new(50, float, 9999)

N1 = tofloat(0.01)
ij = 0

do n = 0, 49
    a1 = 0
    a2 = 0
    a3 = 0
    a4 = 0
    b1 = 0
    b2 = 0
    b3 = 0
    b4 = 0
    a = 0
    
    in = addfile("Frozen_ground_N="+N1+".nc", "r")
    frozen = in->frozen
    lat_temp = in->lat
    lon_temp = in->lon
    obser = addfile("frozen_2000_code_3_2_0.nc", "r")
    obser_frozen = obser->frozen
    obser_lat = obser->lat
    obser_lon = obser->lon
    
    nlat = dimsizes(lat_temp)
    nlon = dimsizes(lon_temp)
    
    N_count = 0
    
    do i = 0, nlat-1, 1
        do j = 0, nlon-1, 1
            pre = frozen(i,j)
            if(.not.ismissing(pre)) then
                longitude = lon_temp(j)
                latitude = lat_temp(i)
                
                index_lat = closest_val(latitude, obser_lat)
                index_lon = closest_val(longitude, obser_lon)
                
                rain = obser_frozen(index_lat, index_lon)
                if(.not.ismissing(rain)) then
                    oo = pre
                    mm = rain
                    mm_oo = mm - oo
                    
                    if(mm_oo.eq.0)
                        a = a + 1
                    end if
                    
                    if(oo.eq.0)
                        a1 = a1 + 1
                    end if
                    
                    if(oo.eq.2)
                        a2 = a2 + 1
                    end if
                    
                    if(oo.eq.3)
                        a4 = a4 + 1
                    end if
                    
                    if(mm.eq.0)
                        b1 = b1 + 1
                    end if
                    
                    if(mm.eq.2)
                        b2 = b2 + 1
                    end if
                    
                    if(mm.eq.3)
                        b4 = b4 + 1
                    end if
                    
                    N_count = N_count + 1
                end if
            end if
        end do
    end do
    
    a_1 = tofloat(a)
    N_1 = tofloat(N_count)
    a1_1 = tofloat(a1)
    a2_1 = tofloat(a2)
    a4_1 = tofloat(a4)
    b1_1 = tofloat(b1)
    b2_1 = tofloat(b2)
    b4_1 = tofloat(b4)
    
    p0 = a_1/N_1
    pe = (a1_1*b1_1+a2_1*b2_1+a4_1*b4_1)/(N_1*N_1)
    
    k = (p0-pe)/(1-pe)
    
    K2(ij) = k
    POD2(ij) = p0
    N11_2(ij) = N1
    
    ij = ij + 1
    N1 = N1 + 0.01
end do

; 第四部分：Ft3分析 (frozen_2000.nc)
K3 = new(6, float, 9999)
POD3 = new(6, float, 9999)
N11_3 = new(6, float, 9999)

N1 = tofloat(0.01)
ij = 0

do n = 0, 5
    a1 = 0
    a2 = 0
    a3 = 0
    a4 = 0
    b1 = 0
    b2 = 0
    b3 = 0
    b4 = 0
    a = 0
    
    in = addfile("Frozen_ground_N="+N1+".nc", "r")
    frozen = in->frozen
    lat_temp = in->lat
    lon_temp = in->lon
    obser = addfile("frozen_2000.nc", "r")
    obser_frozen = obser->frozen
    obser_lat = obser->lat
    obser_lon = obser->lon
    
    nlat = dimsizes(lat_temp)
    nlon = dimsizes(lon_temp)
    
    N_count = 0
    
    do i = 0, nlat-1, 1
        do j = 0, nlon-1, 1
            pre = frozen(i,j)
            if(.not.ismissing(pre)) then
                longitude = lon_temp(j)
                latitude = lat_temp(i)
                
                index_lat = closest_val(latitude, obser_lat)
                index_lon = closest_val(longitude, obser_lon)
                
                rain = obser_frozen(index_lat, index_lon)
                if(.not.ismissing(rain)) then
                    oo = pre
                    mm = rain
                    mm_oo = mm - oo
                    
                    if(mm_oo.eq.0)
                        a = a + 1
                    end if
                    
                    if(oo.eq.0)
                        a1 = a1 + 1
                    end if
                    
                    if(oo.eq.1)
                        a3 = a3 + 1
                    end if
                    
                    if(oo.eq.2)
                        a2 = a2 + 1
                    end if
                    
                    if(oo.eq.3)
                        a4 = a4 + 1
                    end if
                    
                    if(mm.eq.0)
                        b1 = b1 + 1
                    end if
                    
                    if(mm.eq.1)
                        b3 = b3 + 1
                    end if
                    
                    if(mm.eq.2)
                        b2 = b2 + 1
                    end if
                    
                    if(mm.eq.3)
                        b4 = b4 + 1
                    end if
                    
                    N_count = N_count + 1
                end if
            end if
        end do
    end do
    
    a_1 = tofloat(a)
    N_1 = tofloat(N_count)
    a1_1 = tofloat(a1)
    a2_1 = tofloat(a2)
    a3_1 = tofloat(a3)
    a4_1 = tofloat(a4)
    b1_1 = tofloat(b1)
    b2_1 = tofloat(b2)
    b3_1 = tofloat(b3)
    b4_1 = tofloat(b4)
    
    p0 = a_1/N_1
    pe = (a1_1*b1_1+a2_1*b2_1+a3_1*b3_1+a4_1*b4_1)/(N_1*N_1)
    
    k = (p0-pe)/(1-pe)
    
    K3(ij) = k
    POD3(ij) = p0
    N11_3(ij) = N1
    
    ij = ij + 1
    N1 = N1 + 0.01
end do

; 第五部分：输出结果和绘图
print("Ft1 分析结果:")
print("K1 = " + K1)
print("POD1 = " + POD1)

y1_max = max(K1)
y1_ind = ind(K1.eq.y1_max)
x1_max = N11_1(y1_ind)
print("Ft1 最优阈值: F = " + x1_max + ", Kappa = " + y1_max)

print("Ft2 分析结果:")
print("K2 = " + K2)
print("POD2 = " + POD2)

y2_max = max(K2)
y2_ind = ind(K2.eq.y2_max)
x2_max = N11_2(y2_ind)
print("Ft2 最优阈值: F = " + x2_max + ", Kappa = " + y2_max)

print("Ft3 分析结果:")
print("K3 = " + K3)
print("POD3 = " + POD3)

y3_max = max(K3)
y3_ind = ind(K3.eq.y3_max)
x3_max = N11_3(y3_ind)
print("Ft3 最优阈值: F = " + x3_max + ", Kappa = " + y3_max)

; 绘制Ft1结果
wks_type = "png"
wks_type@wkWidth = 2000
wks_type@wkHeight = 1500
wks1 = gsn_open_wks(wks_type, "Kappa_Ft1")

res = True
res@gsnFrame = False
res@gsnMaximize = True
res@gsnDraw = False
res@tmXBMode = "Explicit"
res@tmXBLabelFontHeightF = 0.015
res@tmYLLabelFontHeightF = 0.015
res@tmXBLabelFont = 25
res@tmYLLabelFont = 25
res@tmYLMajorLengthF = 0.01
res@tmXTMajorLengthF = 0.01
res@tmYLMinorLengthF = 0.01
res@tmXTMinorLengthF = 0.01
res@tmXBValues = (/0.1,0.3,0.5,0.7,0.9/)
res@tmXBLabels = (/"0.1","0.3","0.5","0.7","0.9"/)
res@tiMainFontHeightF = 0.015
res@tiYAxisString = "Kappa"
res@tiXAxisString = "F"
res@tiXAxisFont = 25
res@tiXAxisFontHeightF = 0.015
res@tiYAxisFont = 25
res@tiYAxisFontHeightF = 0.015
res@trYMinF = 0
res@trYMaxF = 1
res@xyLineColors = "black"
res@xyLineThicknesses = 5.0
res@xyDashPatterns = 0
res@xyMarkLineModes = "Lines"

plot1 = gsn_xy(wks1, N11_1, K1, res)

x1_refer = (/x1_max, x1_max/)
y1_refer = (/0, y1_max/)
res@xyDashPatterns = 1
plot1_ref1 = gsn_csm_xy(wks1, x1_refer, y1_refer, res)

x2_refer = (/0, x1_max/)
y2_refer = (/y1_max, y1_max/)
plot1_ref2 = gsn_csm_xy(wks1, x2_refer, y2_refer, res)

ts = True
ts@txFontHeightF = 0.02
ts@txFont = 25
ts@txJust = "CenterLeft"
text1 = gsn_add_text(wks1, plot1, "F = "+sprintf("%4.2f", x1_max), x1_max+0.01, y1_max+0.075, ts)
text2 = gsn_add_text(wks1, plot1, "Kappa = "+sprintf("%4.2f", y1_max), x1_max+0.01, y1_max+0.035, ts)

overlay(plot1, plot1_ref1)
overlay(plot1, plot1_ref2)
draw(plot1)
frame(wks1)

; 绘制Ft2结果
wks2 = gsn_open_wks(wks_type, "Kappa_Ft2")
res@tmXBValues = (/0.1,0.3,0.5,0.7/)
res@tmXBLabels = (/"0.1","0.3","0.5","0.7"/)
res@xyDashPatterns = 0

plot2 = gsn_xy(wks2, N11_2, K2, res)

x1_refer = (/x2_max, x2_max/)
y1_refer = (/0, y2_max/)
res@xyDashPatterns = 1
plot2_ref1 = gsn_csm_xy(wks2, x1_refer, y1_refer, res)

x2_refer = (/0, x2_max/)
y2_refer = (/y2_max, y2_max/)
plot2_ref2 = gsn_csm_xy(wks2, x2_refer, y2_refer, res)

text3 = gsn_add_text(wks2, plot2, "F = "+sprintf("%4.2f", x2_max), x2_max+0.01, y2_max+0.075, ts)
text4 = gsn_add_text(wks2, plot2, "Kappa = "+sprintf("%4.2f", y2_max), x2_max+0.01, y2_max+0.035, ts)

overlay(plot2, plot2_ref1)
overlay(plot2, plot2_ref2)
draw(plot2)
frame(wks2)

; 绘制Ft3结果
wks3 = gsn_open_wks(wks_type, "Kappa_Ft3")
res@tmXBValues = (/0.01,0.03,0.05,0.07/)
res@tmXBLabels = (/"0.01","0.03","0.05","0.07"/)
res@xyDashPatterns = 0

plot3 = gsn_xy(wks3, N11_3, K3, res)

x1_refer = (/x3_max, x3_max/)
y1_refer = (/0, y3_max/)
res@xyDashPatterns = 1
plot3_ref1 = gsn_csm_xy(wks3, x1_refer, y1_refer, res)

x2_refer = (/0, x3_max/)
y2_refer = (/y3_max, y3_max/)
plot3_ref2 = gsn_csm_xy(wks3, x2_refer, y2_refer, res)

text5 = gsn_add_text(wks3, plot3, "F = "+sprintf("%4.2f", x3_max), x3_max+0.01, y3_max+0.075, ts)
text6 = gsn_add_text(wks3, plot3, "Kappa = "+sprintf("%4.2f", y3_max), x3_max+0.01, y3_max+0.035, ts)

overlay(plot3, plot3_ref1)
overlay(plot3, plot3_ref2)
draw(plot3)
frame(wks3)

end
