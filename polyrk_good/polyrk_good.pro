files = 'polyrk_0000_' + ['00011', '00100', '01000'] + '.plt'
n_files = n_elements( files )

pltcom, 'polplt.rk
;omsq = getvec('omsq')
;nmodes = n_elements(omsq)
nmodes = 2
k_p = findgen(nmodes) + 1
;wpes = getvec('wpes')
;wsch = getvec('wsch')
x = getvec('x')
npts = n_elements(x)
rz = getvec('rz')
nz = n_elements(rz)
;
n_win = 0
if( nmodes gt 10 ) then begin
   !x.title = '!8k!6'
   !y.title = '!7x!6!u2!n'
   start_plot, n_win, 'omsq'
      plot, k_p, omsq
      oplot, k_p, Wpes, line=5
     oplot, k_p, wsch, line=3
   end_plot
endif
;
my_lines = [0, 5, 3, 4, 2]
n_my_lines = n_elements(my_lines)
;
!x.title = '!8x = r/R!6!d*!n'
!y.title = '!7d!8r/r!6'
start_plot, n_win, 'dr_r'
   yplot = fltarr(npts,nmodes)
   for imode=0, nmodes-1 do yplot[*,imode] = getvec('dr/r')
   plot, x, yplot[*,0], yrange=[min(yplot), max(yplot)]
   for imode=1, nmodes-1 do begin
;      dror = getvec('dr/r')
      oplot, x, yplot[*,imode], line=(imode mod n_my_lines)
   endfor
end_plot
;
bv = getvec('bv')
!y.title = '!7d!8h/h'
start_plot, n_win, 'dhor'
   yplot = fltarr(nz,nmodes)
   for imode=0, nmodes-1 do yplot[*,imode] = getvec('dh/h')
   plot, x, yplot[*,0], yrange=[min(yplot), max(yplot)]
   for imode=1, nmodes-1 do begin
      oplot, x, yplot[*,imode], line=(imode mod n_my_lines)
   endfor
end_plot
;
drho = getvec('bv')
!y.title = '!7dq/q!6'
start_plot, n_win, 'drho'
   drho = fltarr(nz,nmodes)
   for imode=0, nmodes-1 do drho[*,imode] = getvec('drho')
;   print, drho
   plot, rz, drho[*,0], yrange=[min(drho), max(drho)]
   for imode=1, nmodes-1 do begin
;      drho = getvec('drho')
      oplot, rz, drho[*,imode], line=(imode mod n_my_lines)
   endfor
end_plot
;
drho = getvec('bv')
;!y.title = '!7dq/q!6 (Eulerian)'
;start_plot, n_win, 'Erho'
;   yplot = fltarr(nz,nmodes)
;   for imode=0, nmodes-1 do yplot[*,imode] = getvec('Erho')
;;   drho = getvec('Erho')
;;   print, drho
;   plot, x, yplot[*,0], yrange=[min(yplot), max(yplot)]
;   for imode=1, nmodes-1 do begin
;      oplot, x, yplot[*,imode], line=(imode mod n_my_lines)
;   endfor
;end_plot
;
;goto, skip_plots
;
drho = getvec('bv')
!y.title = '!7c!6'
start_plot, n_win, 'gamma'
   yplot = fltarr(nz,nmodes)
   for imode=0, nmodes-1 do yplot[*,imode] = getvec('gam')
   plot, x, yplot[*,0], yrange=[min(yplot), max(yplot)]
   for imode=1, nmodes-1 do begin
      oplot, x, yplot[*,imode], line=(imode mod n_my_lines)
   endfor
end_plot
;
drho = getvec('bv')
!y.title = '!6Integrated Weight Function'
start_plot, n_win, 'drho'
   drho = getvec('wint')
;   print, drho
   plot, x, drho ;, xrange=[0.95,1]
   for imode=1, nmodes-1 do begin
      drho = getvec('wint')
      oplot, x, drho, line=(imode mod n_my_lines)
   endfor
end_plot
;
drho = getvec('bv')
!y.title = '!6Kinetic Energy Amplitude'
start_plot, n_win, 'rke'
   drho = getvec('rke')
   plot, x, drho;, yrange=[-0.03, 0.03]
   for imode=1, nmodes-1 do begin
      drho = getvec('rke')
      oplot, x, drho, line=(imode mod n_my_lines)
   endfor
end_plot
;
drho = getvec('bv')
!y.title = '!6Weight Function'
for imode=0, nmodes-1 do begin
   start_plot, n_win, 'weight_' + string(format='(i2.2)', imode)
      wait = getvec('wait')
      wthr = getvec('wthr')
      wgrv = getvec('wgrv')
      wdia = getvec('wdia')
      wcrs = getvec('wcrs')
      y_max = max( [ wait, wthr, wgrv, wdia, wcrs ] )
      y_min = min( [ wait, wthr, wgrv, wdia, wcrs ] )
      plot, rz, wait, thick=1.5, xthick=1.5, ythick=1.5, yrange=[y_min,y_max]
      oplot, rz, wthr, line=5, thick=1.5
      oplot, rz, wgrv, line=1, thick=1.5
      oplot, rz, wdia, line=3, thick=1.5
      oplot, rz, wcrs, line=2, thick=1.5
   end_plot
endfor
;
drho = getvec('bv')
!y.title = '!6Weight Function'
start_plot, n_win, 'weight'
   drho = getvec('wait')
   plot, rz, drho ;, xrange=[0.95,1], yrange=[-0.03, 0.03]
   for imode=1, nmodes-1 do begin
      drho = getvec('wait')
      oplot, rz, drho, line=(imode mod n_my_lines)
   endfor
end_plot
;
drho = getvec('bv')
!y.title = '!6Thermal Weight Function'
start_plot, n_win, 'weight'
   drho = getvec('wthr')
   plot, rz, drho, yrange=[-0.03, 0.03] ;, xrange=[0.95,1]
   for imode=1, nmodes-1 do begin
      drho = getvec('wthr')
      oplot, rz, drho, line=(imode mod n_my_lines)
   endfor
end_plot
;
drho = getvec('bv')
!y.title = '!6Gravitational Weight Function'
start_plot, n_win, 'weight'
   yplot = fltarr(npts,nmodes)
   for imode=0, nmodes-1 do yplot[*,imode] = getvec('wgrv')
   plot, x, yplot[*,0], yrange=[min(yplot), max(yplot)]
   for imode=1, nmodes-1 do begin
      oplot, x, yplot[*,imode], line=(imode mod n_my_lines)
   endfor
end_plot
;
drho = getvec('bv')
!y.title = '!6Cross Weight Function'
start_plot, n_win, 'weight'
   drho = getvec('wcrs')
   plot, x, drho ;, xrange=[0.95,1], yrange=[-0.03, 0.03]
   for imode=1, nmodes-1 do begin
      drho = getvec('wcrs')
      oplot, x, drho, line=(imode mod n_my_lines)
   endfor
end_plot
;
drho = getvec('bv')
!y.title = '!6Diagonal Weight Function'
start_plot, n_win, 'weight'
   drho = getvec('wdia')
   plot, x, drho ;, xrange=[0.95,1], yrange=[-0.03, 0.03]
   for imode=1, nmodes-1 do begin
      drho = getvec('wdia')
      oplot, x, drho, line=(imode mod n_my_lines)
   endfor
end_plot
skip_plots:
close, 1
end
