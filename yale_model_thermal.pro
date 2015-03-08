;+
;       Modified to do the R value calculation 03-MAR-2009
;-
;openr, lunit, 'ssm455.sink', /get_lun
;npts = 378
;data_rec = { shell: 0, Radius: 0., Mass: 0., Radius_cm: 0., Mass_g: 0., Lum_erg_s: 0., T: 0., rho: 0., $
;        P: 0.,  Cv: 0., Cp: 0., Gamma1: 0.,  mu: 0., c_s: 0., Eps: 0., $
;        X: 0.,  Z: 0., kappa: 0.,   dKdr: 0.,   dkdt: 0., Del: 0., del_ad: 0., $
;        dlnP_dlnRho: 0.,   dlnP_dlnT: 0.,     dlnEps_dlnRho: 0., dlnEps_dlnT: 0.,   Lamb2: 0.,   N2: 0., $
;        F_crit_p_2: 0.,   F_crit_m_2: 0.,   alpha: 0.,           H_P: 0.,  lambda: 0.}
;data_st = replicate( data_rec, npts)
openr, lunit, '../../solar/evolution/ssm_d3_7d_Grey.txt', /get_lun
npts = 2012
data_rec = { shell: 0, Radius: 0.d0, Mass: 0.d0, Radius_cm: 0.d0, Mass_g: 0.d0, Lum_erg_s: 0.d0, T: 0.d0, rho: 0.d0, $
        P: 0.d0,  Cv: 0.d0, Cp: 0.d0, Gamma1: 0.d0,  mu: 0.d0, c_s: 0.d0, Eps: 0.d0, $
        X: 0.d0,  Z: 0.d0, kappa: 0.d0,   dKdr: 0.d0,   dkdt: 0.d0, Del: 0.d0, del_ad: 0.d0, $
        dlnP_dlnRho: 0.d0,   dlnP_dlnT: 0.d0,     dlnEps_dlnRho: 0.d0, dlnEps_dlnT: 0.d0,   Lamb2: 0.d0,   N2: 0.d0, $
        F_crit_p_2: 0.d0,   F_crit_m_2: 0.d0}
data_st = replicate( data_rec, npts)
a = ' '
for i=0, 100 do begin
   readf, lunit, a
   if( strpos( a, 'Freq-') gt 0 ) then goto, read_vals
endfor
goto, skip_file
;
read_vals:
;
readf, lunit, a
for i=0, npts-1 do begin
   readf, lunit, data_rec
   data_st[i] = data_rec
endfor
free_lun, lunit
;
!p.multi = [0,2,2]
!x.title = '!8r!6/!8R!9!dn!n!6'
n_win = 0
!p.thick = 3.75
!x.thick = !p.thick
!y.thick = !p.thick
!p.charsize = 2.5
!p.charthick = 2.5
start_plot, n_win, 'yale_mod_basic'
   !y.title = '!8T!6 (K)'
   plot_io, data_st.radius, data_st.t, xrange=[0,1.05], /xstyle
   !y.title = '!7q!6 (gr cm!u-3!n)'
   plot_io, data_st.radius, data_st.rho, xrange=[0,1.05], /xstyle
   !y.title = '!8P!6 (dyn cm!u-2!n)'
   plot_io, data_st.radius, data_st.P, xrange=[0,1.05], /xstyle
   !y.title = '!8L!6 (erg s!u-1!n)'
   plot, data_st.radius, data_st.lum_erg_s, xrange=[0,1.05], /xstyle
end_plot
!p.multi = 0
start_plot, n_win, 'yale_mod_thermo'
   !y.title = '!6Various Thermodynamic Quantities'
   plot, data_st.radius, data_st.gamma1, xrange=[0,1.05], /xstyle, yrange=[0,5]
   oplot, data_st.radius, data_st.dlnP_dlnRho, line=1
   oplot, data_st.radius, data_st.dlnP_dlnT, line=2
   oplot, data_st.radius, data_st.del, line=4
   oplot, data_st.radius, data_st.del_ad, line=3
   oplot, data_st.radius, alog10(data_st.kappa), line=5
;
   x_fig = 0.15
   y_fig = 0.65
   del_x = 0.175
   del_y = 6*0.04
   labels = [ ['%0', '%1', '%2', '%3', '%4', '%5'], $
              ['%L' + ['!7C!6!d1!n','!7v!dq!6!n','!7v!6!dT!n','!9G!6!dad!n','!9G!6', '!6log(!7j!6)'] ] ]
   table, x_fig, y_fig, x_fig+del_x, y_fig+del_y, labels ;, line_color=0
;
end_plot
iout = 3
lval = 2
gor = 6.67d-8*data_st.Mass_g/data_st.radius_cm^2
H_p = data_st.p/data_st.rho/gor
dm1 = (data_st[1:npts-1].radius_cm^3 - data_st[0:npts-2].radius_cm^3)*!pi*4/3.*data_st.rho
dm2 = [dm1[0], (dm1[1:npts-2] + dm1[0:npts-3])/2., dm1[npts-2]/2.]
rzone = (data_st[1:npts-1].radius_cm + data_st[0:npts-2].radius_cm)/2.
v = 1./data_st[0:npts-2].rho
g3m1 = (data_st.gamma1 - data_st.dlnP_dlnrho)/data_st.dlnP_dlnT
radial_modes = castor(0, 1, 3, npts-1, data_st.radius_cm, data_st.t, v, data_st.cv, data_st.dkdr, $
        data_st.dkdt, dm1, data_st.kappa, dm2, data_st.Mass_g, data_st.N2, data_st.p, $
        data_st.gamma1, g3m1, data_st.rho, rzone, xo)
goto, skip_file
test = matset( npts-1, iout, lval, data_st.radius, data_st.t, v, data_st.cv, data_st.dkdr, $
        data_st.dkdt, dm1, data_st.kappa, dm2, data_st.Mass, data_st.N2, data_st.p, $
        data_st.gamma1, g3m1, data_st.rho, rzone)
skip_file:                     
;sigma_SB = 5.66961e-5
;thermal_time_dm = sqrt(3.*data_st.kappa*data_st.cv*data_st.T/(4.*7.56591e-15*3.e10*data_st.T^4))*dm1/(4.*!pi*data_st.radius_cm^2)
;thermal_time_P_dm = data_st.Cv*data_st.T*(dm1/data_st.Lum_erg_s)
;npts = n_elements(thermal_time_dm)
;thermal_time = fltarr(npts)
;thermal_time[npts-1] = 1.e0*thermal_time_dm[npts-1]
;for i=npts-2,0,-1 do thermal_time[i] = thermal_time[i+1] + thermal_time_dm[i]
;thermal_time_P = fltarr(npts)
;thermal_time_P[npts-1] = 1.e4*thermal_time_P_dm[npts-1]
;for i=npts-2,0,-1 do thermal_time_P[i] = thermal_time_P[i+1] + thermal_time_P_dm[i]
;
;       Critical frequency calculation and comparison
;
omega_ac_2 = data_st.Gamma1*data_st.P/data_st.rho/H_p^2/4. ;*(1.d0 + 2.d0*2.d0*H_p/data_st.radius_cm)
om_inter_2 = (omega_ac_2+data_st.Lamb2)/2.
om_crit_2_p = om_inter_2 + sqrt( om_inter_2^2 - data_st.N2*data_st.Lamb2)
om_crit_2_m = om_inter_2 - sqrt( om_inter_2^2 - data_st.N2*data_st.Lamb2)
!x.title = '!8r!6'
!y.title = '!6Frequencies'
start_plot, n_win, 'yale_mod_freqs', color=39
   plot_io, data_st.radius, data_st.Lamb2, /nodata, xrange=[0,1.05], /xstyle, yrange=[1.e-7,1.e3], color=0, back=255
   convective = where( data_st.N2 gt 0. and data_st.radius lt 0.75, n_convective)
;   oplot, data_st[convective].radius, data_st[convective].F_crit_m_2, line=3
;   oplot, data_st[convective].radius, data_st.F_crit_p_2, line=4
   oplot, data_st.radius, om_crit_2_p, color=50, psym=1
   oplot, data_st[convective].radius, om_crit_2_m[convective], color=30, psym=1
   L_new = [5, 50, 500L]
   for i=0, n_elements(L_new)-1 do begin
      label = where( data_st.radius lt 0.75, n_label)
      Lamb20 = double(L_new[i]*(L_new[i]+1))*data_st.Gamma1*data_st.P/data_st.rho/data_st.radius_cm^2
      om_inter_20 = (omega_ac_2+Lamb20)/2.
      om_crit_20_p = om_inter_20 + sqrt( om_inter_20^2 - data_st.N2*Lamb20)
      om_crit_20_m = om_inter_20 - sqrt( om_inter_20^2 - data_st.N2*Lamb20)
      oplot, data_st.radius, om_crit_20_p, color=70+50*i, psym=1
;      oplot, data_st[convective].radius, om_crit_20_m[convective], color=70+50*i, psym=1
   endfor
   oplot, data_st[convective].radius, data_st[convective].N2, line=5, color=0
   oplot, data_st.radius, data_st.Lamb2, color=0
   oplot, !x.crange, 2.65e-6*[1,1], line=1, color=0
end_plot
;!x.title = '!61 - !8r!6/!8R!9!dn!n!6'
;!y.title = '!6Thermal Times (years)'
;start_plot, n_win, 'yale_mod_thermal_times', /encap
;   plot_oo, 1.-data_st.radius, (thermal_time/!pi)^2/(365.25*86400.), yrange=[1.e-8,1.e10], /ystyle, xrange=[1,1.e-6]
;   oplot, 1.-data_st.radius, thermal_time_P/(365.25*86400.), line=5
;   oplot, 10.^!x.crange, 1.e6 * [1,1], line=1
;   oplot, 10.^!x.crange, 1.e5 * [1,1], line=1
;   oplot, 10.^!x.crange, 1.e4 * [1,1], line=1
;   oplot, 10.^!x.crange, 11 * [1,1], line=1
;   oplot, (1-0.72)*[1,1], 10.^!y.crange, line=1
;end_plot
;
end