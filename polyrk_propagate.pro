;pltcom, 'polplt.rk
;pltcom, 'polyrk_3000_04094.plt'
;pltcom, 'polyrk_4500_04001.plt'
pltcom, 'polyrk_4000_05000.plt'
;omsq = getvec('omsq')
;nmodes = n_elements(omsq)
;k_p = findgen(nmodes) + 1
;wpes = getvec('wpes')
;wsch = getvec('wsch')
x = getvec('x')
npts = n_elements(x)
;rz = getvec('rz')
nz = n_elements(rz)
;
n_win = 1
nmodes = 0
if( nmodes gt 1 ) then begin
   !x.title = '!8k!6'
   !y.title = '!7x!6!u2!n'
   start_plot, n_win, 'omsq'
      plot, k_p, omsq
      oplot, k_p, Wpes, line=5
     oplot, k_p, wsch, line=3
   end_plot
endif
;
bv = getvec('bv')
sl2 = getvec('sl2')
gor = getvec('g')
theta = getvec('thet')
dtdr = getvec('dtdr')
posit = where( dtdr gt 0)
H_p = theta[safe]*x[safe]^2/dtdr[safe]
my_lines = [0, 5, 3, 4, 2]
n_my_lines = n_elements(my_lines)
;
;       n = 4.5
g_modes = [7.831E+01, 8.859E+01, 9.902E+01, 1.108E+02, 1.245E+02]
f_mode = 1.391E+02
p_modes = [1.560E+02, 1.745E+02, 1.950E+02, 2.174E+02, 2.423E+02, $
           2.701E+02, 2.980E+02, 3.304E+02, 3.674E+02, 4.034E+02, $
           4.398E+02, 4.835E+02]
;
;       n = 3.00
;p_modes = [20.35, 35.63, 55.29, 79.23, 107.4]
;f_mode = 10.9
;g_modes = [6.553, 3.771, 2.430, 1.694, 1.248]
!p.thick = 3
!x.thick = !p.thick
!y.thick = !p.thick
!p.charsize = 2
!p.charthick = 2
!p.title = '!6 '
!x.title = '!8x = r/R!6!d*!n'
!y.title = '!8N!6!u2!n & !8S!6!s!d2!r!u2!n'
start_plot, n_win, 'propagate', color=39
   plot_io, x, 10.^BV, /nodata, xrange=[0,1], /xstyle, yrange=[10., 1.e4], /ystyle, color=0, back=255
   oplot, x, BV, color=0
   oplot, x, 10.^sl2, color=0
   for i=0, n_elements(p_modes)-1 do begin
      safe = where( 10.^sl2 lt p_modes[i] and BV lt p_modes[i], n_safe)
      if( n_safe gt 1 ) then oplot, [min(x[safe]), max(x[safe])], p_modes[i]*[1,1], line=5, color=50
      oplot, [min(x), max(x)], p_modes[i]*[1,1], line=5, color=50
   endfor
   oplot, !x.crange, f_mode*[1,1], line=5, color=30
   for i=0, n_elements(g_modes)-1 do begin
      safe = where( 10.^sl2 gt g_modes[i] and BV gt g_modes[i], n_safe)
      if( n_safe gt 1 ) then oplot, [min(x[safe]), max(x[safe])], g_modes[i]*[1,1], line=5, color=240
   endfor
   xyouts, 0.25, 18, '!8N!6!u2!n', color=0
   xyouts, 0.225, 75, '!8S!6!s!d2!r!u2!n', color=0, /align
end_plot
;
skip_plots:
close, 1
end
