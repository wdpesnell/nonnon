pltcom, 'polplt.rk
omsq = getvec('omsq')
nmodes = n_elements(omsq)
k_p = findgen(nmodes) + 1
wpes = getvec('wpes')
wsch = getvec('wsch')
x = getvec('x')
rzone = getvec('rz')
;
n_win = 0
if( n_elements(omsq) gt 1 ) then begin
   !x.title = '!8k!6'
   !y.title = '!7x!6!u2!n'
   start_plot, n_win, 'omsq'
      plot, k_p, Wpes, line=5
      oplot, k_p, omsq
;      oplot, k_p, wsch, line=3
   end_plot
endif
;
lines = [5, 3, 4, 0]
n_lines = n_elements(lines)
!x.title = '!8x = r/R!6!d*!n'
!y.title = '!7d!8r/r!6'
start_plot, n_win, 'dr_r'
   dror = getvec('dr/r')
;   print, dror
   plot, x, dror, yrange=[-5,5], xrange=[0., 1.05], /xstyle
   for imode=1, nmodes-1 do begin
      dror = getvec('dr/r')
      oplot, x, dror, line=lines[imode mod n_lines]
   endfor
end_plot
;
;
drho = getvec('bv')
!y.title = '!7dq/q!6'
start_plot, n_win, 'drho'
   drho = getvec('drho')
;   print, drho
   plot, rzone, drho, yrange=[-50,50], xrange=[0.8, 1.05], /xstyle
   for imode=1, nmodes-1 do begin
      drho = getvec('drho')
      oplot, rzone, drho, line=lines[imode mod n_lines]
   endfor
end_plot
;
drho = getvec('bv')
!y.title = '!6Kinetic Energy Amplitude'
start_plot, n_win, 'rke'
   drho = getvec('rke')
   plot, x, drho, yrange=[0., 0.5], xrange=[0., 1.05], /xstyle
   for imode=1, nmodes-1 do begin
      drho = getvec('rke')
      oplot, x, drho, line=lines[imode mod n_lines]
   endfor
end_plot
;
drho = getvec('bv')
!y.title = '!6Weight Function'
start_plot, n_win, 'weight'
   drho = getvec('wait')
   plot, x, drho, yrange=[-0.25, 0.25], xrange=[0., 1.05], /xstyle
   for imode=1, nmodes-1 do begin
      drho = getvec('wait')
      oplot, x, drho, line=lines[imode mod n_lines]
   endfor
end_plot
end