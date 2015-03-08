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
v = 1./data_st[0:npts-1].rho
g3m1 = (data_st.gamma1 - data_st.dlnP_dlnrho)/data_st.dlnP_dlnT
;
;       Calculate the radial modes
;
ihmin = 1   ; radial fundamental
ihmax = 51
radial_modes = castor(0, ihmin, ihmax, npts-1, data_st.radius_cm, data_st.t, v, data_st.cv, data_st.dkdr, $
        data_st.dkdt, dm1, data_st.kappa, dm2, data_st.Mass_g, data_st.N2, data_st.p, $
        data_st.gamma1, g3m1, data_st.rho, rzone, xo)
goto, skip_file
test = matset( npts-1, iout, lval, data_st.radius, data_st.t, v, data_st.cv, data_st.dkdr, $
        data_st.dkdt, dm1, data_st.kappa, dm2, data_st.Mass, data_st.N2, data_st.p, $
        data_st.gamma1, g3m1, data_st.rho, rzone)
skip_file:                     
!x.title = '!61 - !8r!6/!8R!9!dn!n!6'
!y.title = '!7d!8r!6/!8r!6'
start_plot, n_win, 'yale_mod_radial_modes', /encap
   case !d. name of
;      'X': window, n_win-1, xsize=1024, ysize=1024+64
      'X': window, n_win-1, xsize=1024, ysize=700
      'Z': device, set_resolution=[1024,1024+64]
      else:
   endcase
   plot_oi, [1], [1], /nodata, yrange=[-1.05,1.05], /ystyle, xrange=[1,1.e-5]
   for ih=ihmin, ihmax do oplot, 1.d0-data_st.radius, radial_modes[ih-1].dror/max(abs(radial_modes[ih-1].dror)) ;, line=5
      oplot, !x.crange, [0,0], line=1
end_plot
!x.title = '!61 - !8r!6/!8R!9!dn!n!6'
!y.title = '!7dq!6/!7q!6'
start_plot, n_win, 'yale_mod_radial_modes', /encap
   case !d. name of
;      'X': window, n_win-1, xsize=1024, ysize=1024+64
      'X': window, n_win-1, xsize=1024, ysize=700
      'Z': device, set_resolution=[1024,1024+64]
      else:
   endcase
   plot_oi, [1], [1], /nodata, yrange=[-2.e-2,2.e-2], /ystyle, xrange=[1,1.e-5]
   for ih=ihmin, ihmax do oplot, 1.d0-data_st.radius, radial_modes[ih-1].drho/max(abs(radial_modes[ih-1].drho)) ;, line=5
      oplot, !x.crange, [0,0], line=1
end_plot
;
;       Fit an asymptotic line to the radial frequencies
safe = where(radial_modes.k gt 0)
coeffs = linfit( radial_modes[safe].k, radial_modes[safe].omega_ad )
print, 'Delta nu = ', coeffs/2/!pi
!x.title = '!6Order'
!y.title = '!6Frequency'
start_plot, n_win, 'yale_mod_radial_modes', /encap
   plot, radial_modes.k, radial_modes.omega_ad, psym=1, xrange=[0, ihmax+1], /xstyle
   oplot, !x.crange, poly( !x.crange, coeffs)
end_plot
;
end