puls_sum = { mass_flag: 0, $
             ell: 0L, k: 0L, n_p: 0L, n_g: 0L, $
              omega_ad: 0., omsq_n: 0., nu_ad: 0., period: 0., $
	      acc_ad_w: 0., acc_ad_e: 0., rbar: 0., $
         del_om_qa: complex(0.0,0.0), $
         comega: complex(0.0,0.0), $
         kappa: 0., acc_na_k: 0., date_str: '        ', time_str: '          '}
nonradial = puls_sum
buf = '                  '
openu, lunit, 'model_s_vals_test.sum', /f77_unformatted, /get_lun
i = 0L
while( not(eof(lunit)) ) do begin
   readu, lunit, nonradial
;   print, i, nonradial.ell, nonradial.k, nonradial.nu_ad*1.e6, nonradial.period, ' ', nonradial.date_str,' ', nonradial.time_str
   i = i + 1L
endwhile
free_lun, lunit
npts = i-1L
nonradial = replicate( puls_sum, npts )
openu, lunit, 'model_s_vals_test.sum', /f77_unformatted, /get_lun
for i=0L, npts-1 do begin
   readu, lunit, puls_sum
   nonradial[i] = puls_sum
endfor
free_lun, lunit
n_repeat = -1L
if( n_repeat eq 0 ) then begin
   for i=0, npts-2 do begin
      XIABS = nonradial[i].omsq_n
      for j=I+1,npts-1 do begin
         XJABS = nonradial[j].omsq_n
         IF( (nonradial[i].ell eq nonradial[j].ell) and (nonradial[i].k eq nonradial[j].k) ) then begin ;and $
;      ABS(nonradial[i].omsq_n-nonradial[j].omsq_n) LE 1.e-4*(XIABS+XJABS)/2.D0 ) THEN begin
;         print, 'Identical roots: ', i, nonradial[i].ell, nonradial[i].k, nonradial[i].omsq_n, $
;	                             j, nonradial[j].ell, nonradial[j].k, nonradial[j].omsq_n
            n_repeat = n_repeat + 1
         endif
      endfor
   endfor
   print, n_repeat, ' repeated roots, out of ', n_elements(nonradial)
endif
;
;       Determine which modes are pesent in the file
;
k_min = min(nonradial.k)
k_max = max(nonradial.k)
fmode = where( nonradial.k eq 0 and (nonradial.nu_ad gt 0.), n_fmode)
p1 = where( nonradial.k eq 1 and (nonradial.nu_ad le 20.e-3), n_p1)
pmode = where( nonradial.k gt 0 and (nonradial.nu_ad le 20.e-3), n_pmode)
p70 = where( nonradial.k eq 70, n_p70)
gmode = where( nonradial.k lt 0, n_gmode)
g1 = where( nonradial.k eq -1, n_g1)
;
;       \omega-k diagram (actually a \nu-\ell diagram)
;
n_win = 0
!p.title = '!6' + strtrim(string(n_pmode),2) + ' !8p!6-modes, ' + $
                  strtrim(string(n_fmode),2) + ' !8f!6-modes, ' + $
                  strtrim(string(n_gmode),2) + ' !8g!6-modes'
!x.title = '!13l!6'
!y.title = '!7m!8!dn,!13l!n!6 (mHz)'
!p.thick = 2
!x.thick = !p.thick
!y.thick = !p.thick
!p.charsize = 2
!p.charthick = 2.5
start_plot, n_win, 'model_s_k_ell', color=39, /encap
   case !d.name of
      'Z': device, set_resolution=[1792, 1280] ;$ ; 6144 = 32", 8192=46"
      'X': window, n_win-1, xsize=1152, ysize=800
      'WIN': window, n_win-1, xsize=1536, ysize=640
;      'PS': device, /inches, xsize=18, ysize=12, xoffset=0., yoffset=0.
      else:
   endcase
   tvlct, r, g, b, /get
   r[1] = 128
   g[1] = 128
   b[1] = 128
   tvlct, r, g, b
;   plot, nonradial.ell, nonradial.nu_ad*1.e3, /nodata, psym=1, color=0, back=255, $
;      xr=[-1,max(nonradial.ell)+1], /xstyle, yrange=[0,12], /ystyle
   plot_oo, nonradial.ell, nonradial.nu_ad*1.e3, /nodata, psym=1, color=0, back=255, $
      xr=[1,1000], /xstyle, yrange=10.^[-1,1.2], /ystyle
   for k=k_min, k_max do begin
      safe = where( nonradial.k eq k, n_safe)
;      if( n_safe gt 15 ) then oplot, nonradial[safe].ell, nonradial[safe].nu_ad*1.e3, color=1
   endfor
   if( n_pmode gt 0 ) then oplot, nonradial[pmode].ell, nonradial[pmode].nu_ad*1.e3, psym=1, color=200
   if( n_p1 gt 0 ) then oplot, nonradial[p1].ell, nonradial[p1].nu_ad*1.e3, psym=1, color=240
   if( n_p70 gt 0 ) then oplot, nonradial[p70].ell, nonradial[p70].nu_ad*1.e3, psym=1, color=140
   if( n_fmode gt 0 ) then oplot, nonradial[fmode].ell, nonradial[fmode].nu_ad*1.e3, psym=1, color=30
   if( n_gmode gt 0 ) then oplot, nonradial[gmode].ell, nonradial[gmode].nu_ad*1.e3, psym=1, color=50
   if( n_g1 gt 0 ) then oplot, nonradial[g1].ell, nonradial[g1].nu_ad*1.e3, psym=1, color=80
;
   oplot, !x.crange, 5.2575*[1,1], line=0, color=240  ; acoustic cutoff at surface
;
;   oplot, !x.crange, 2.5e-3*!x.crange, line=0, color=230  ; initial guess in scan
   oplot, [2,!x.crange[1]], 7.2d-2*float([2,!x.crange[1]])^0.4829d0, line=0, color=230  ; initial lower frequency in scan
   oplot, !x.crange, 8.79d0 + 1.7d-2*(100. + float(!x.crange)/2.+0.25d0), line=0, color=230  ; initial upper frequency in scan
;   
   oplot, 10.^!x.crange, 5.2575*[1,1], line=0, color=240  ; acoustic cutoff at surface
;   oplot, 10.^!x.crange, 2.5e-3*10.^!x.crange, line=0, color=230  ; initial guess in scan
   oplot, [2,10.^!x.crange[1]], 7.2d-2*float([2,10.^!x.crange[1]])^0.4805d0, line=0, color=230  ; initial lower frequency in scan
   oplot, 10.^!x.crange, 8.79d0 + 1.7d-2*(100. + float(10.^!x.crange)/2.+0.25d0), line=0, color=230  ; initial upper frequency in 
;
   oplot, 10.^!x.crange, 0.1104d0*float(10.^!x.crange)^0.4829d0, line=0, color=30   ; f-mode fit
   oplot, 10.^!x.crange, 0.2289d0*float(10.^!x.crange)^0.4117d0, line=0, color=240  ; p1-mode fit
;
end_plot, /png
;
;
!p.title = '!6' + strtrim(string(n_pmode),2) + ' !8p!6-modes'
!y.title = '(!7m!8!dn,!13l!n!6-!7m!6!dmodel!n)/!7m!8!dn,!13l!n!6'
x_nu = double(nonradial[pmode].k) + double(nonradial[pmode].ell)/2. + 0.25
coeffs = linfit( x_nu, nonradial[pmode].nu_ad )
print, 'p-mode: ', coeffs*1.e3
coeffs_p = linfit( alog(x_nu), alog(nonradial[pmode].nu_ad) )
print, 'p-mode (power): ', coeffs_p
x_nu_f = double(nonradial[fmode].ell)
coeffs_f = linfit( alog(x_nu_f), alog(nonradial[fmode].nu_ad) )
print, 'f-mode: ', coeffs_f ;*1.e3
start_plot, n_win, 'model_s_k_ell_diff'
   plot, nonradial[pmode].ell, (nonradial[pmode].nu_ad-poly(x_nu, coeffs))/nonradial[pmode].nu_ad, psym=1, $
      xr=[0,max(nonradial.ell)+1], /xstyle, yrange=2.*[-1,1]
;   for k=k_min, k_max do begin
;      safe = where( nonradial.k eq k, n_safe)
;      if( n_safe gt 4 ) then oplot, nonradial[safe].ell, nonradial[safe].nu_ad*1.e6
;   endfor
end_plot
!x.title = '!8x!6 = !8n!6 + !13l!6/2 + 1/4'
start_plot, n_win, 'model_s_k_ell_fit'
   case !d.name of
      'Z': device, set_resolution=[1792, 1280] ;$ ; 6144 = 32", 8192=46"
      'X': window, n_win-1, xsize=1152, ysize=800
      'WIN': window, n_win-1, xsize=1536, ysize=640
;      'PS': device, /inches, xsize=18, ysize=12, xoffset=0., yoffset=0.
      else:
   endcase
   plot, x_nu, (nonradial[pmode].nu_ad-poly(x_nu,coeffs))/nonradial[pmode].nu_ad, psym=1, xr=[0,max(x_nu)+1], /xstyle, $
      yrange=[-3,0.5], /ystyle
end_plot
start_plot, n_win, 'model_s_k_ell_power_fit'
   case !d.name of
      'Z': device, set_resolution=[1792, 1280] ;$ ; 6144 = 32", 8192=46"
      'X': window, n_win-1, xsize=1152, ysize=800
      'WIN': window, n_win-1, xsize=1536, ysize=640
;      'PS': device, /inches, xsize=18, ysize=12, xoffset=0., yoffset=0.
      else:
   endcase
   plot, x_nu, (nonradial[pmode].nu_ad-exp(poly(alog(x_nu), coeffs_p)))/nonradial[pmode].nu_ad, psym=1, xr=[0,max(x_nu)+1], /xstyle, $
      yrange=[-3,0.5], /ystyle
end_plot
!x.title = '!8x!6 = !8n!6 + !13l!6/2 + 1/4'
!y.title = '!6<!8r!6>/!8R!9!dn!n!6'
start_plot, n_win, 'model_s_rbar', color=39
   case !d.name of
      'Z': device, set_resolution=[1792, 1280] ;$ ; 6144 = 32", 8192=46"
      'X': window, n_win-1, xsize=1152, ysize=800
      'WIN': window, n_win-1, xsize=1536, ysize=640
;      'PS': device, /inches, xsize=18, ysize=12, xoffset=0., yoffset=0.
      else:
   endcase
   plot, x_nu, nonradial[pmode].rbar, /nodata, psym=1, xr=[-10,max(x_nu)+1], /xstyle, $
      yrange=[0.0,1.01], /ystyle
   if( n_pmode gt 0 ) then oplot, x_nu, nonradial[pmode].rbar, psym=1, color=200
   if( n_fmode gt 0 ) then oplot, nonradial[fmode].ell/2.+nonradial[fmode].k+0.25, nonradial[fmode].rbar, psym=1, color=50
   if( n_gmode gt 0 ) then oplot, nonradial[gmode].ell/2.+nonradial[gmode].k+0.25, nonradial[gmode].rbar, psym=4, color=80
end_plot
!p.title = '!6' + strtrim(string(n_fmode),2) + ' !8f!6-modes and ' + strtrim(string(n_p1),2) + ' !8p!6!d1!n-modes'
!x.title = '!8x!6 = !13l!6'
start_plot, n_win, 'model_s_rbar', color=39
   case !d.name of
      'Z': device, set_resolution=[1792, 1280]
      'X': window, n_win-1, xsize=1152, ysize=800
      'WIN': window, n_win-1, xsize=1536, ysize=640
;      'PS': device, /inches, xsize=18, ysize=12, xoffset=0., yoffset=0.
      else:
   endcase
   plot, x_nu, nonradial[pmode].rbar, /nodata, psym=1, xr=[0,max(nonradial[fmode].ell)+1], /xstyle, $
      yrange=[0.8,1], /ystyle
   if( n_fmode gt 0 ) then oplot, nonradial[fmode].ell, nonradial[fmode].rbar, psym=1, color=50
   if( n_p1 gt 0 ) then oplot, nonradial[p1].ell, nonradial[p1].rbar, psym=1, color=220
end_plot
end
