pltcom, 'p_n0_l0.rk
omsq = getvec('omsq')
nmodes = n_elements(omsq)
k_p = findgen(nmodes) + 1
wpes = getvec('wpes')
wsch = getvec('wsch')
x = getvec('x')
;
n_win = 0
;!x.title = '!8k!6'
;!y.title = '!7x!6!u2!n'
;start_plot, n_win, 'dr_r'
;   plot, k_p, omsq
;   oplot, k_p, Wpes, line=5
;   oplot, k_p, wsch, line=3
;end_plot
;
!p.thick = 2.5
!x.thick = !p.thick
!y.thick = !p.thick
!p.charsize = 1.4
!p.charthick = 1.9

!p.title = '!6 '
!x.title = '!8x = r/R!6!d*!n'
!y.title = '!7d!8r/r!6'
start_plot, n_win, 'non_fig_1', /encap
   dror = getvec('dr/r')
;   print, dror
   plot, x, dror, yrange=[-5.,5], /ystyle
   analytic = fltarr(32)
   x_anal = findgen(32)/float(32)
   analytic(*) = 1.
   oplot, x_anal, analytic, psym=1
   analytic = -(1. - 7./5.*x_anal^2)*5./2.
   oplot, x_anal, analytic, psym=1
   analytic = (1. - (18./5.)*x_anal^2 + (99./35.)*x_anal^4)*35./8.
   oplot, x_anal, analytic, psym=1
   for imode=1, nmodes-1 do begin
      dror = getvec('dr/r')
      oplot, x, dror
   endfor
end_plot
goto, skip_rest
;
;
drho = getvec('bv')
!y.title = '!7dq/q!6'
start_plot, n_win, 'drho'
   drho = getvec('drho')
;   print, drho
   plot, x, drho, yrange=[-25,10]
   for imode=1, nmodes-1 do begin
      drho = getvec('drho')
      oplot, x, drho
   endfor
end_plot
;
drho = getvec('bv')
!y.title = '!6Kinetic Energy Amplitude'
start_plot, n_win, 'rke'
   drho = getvec('rke')
   plot, x, drho, yrange=[0., 0.025], xrange=[0.95, 1.]
   for imode=1, nmodes-1 do begin
      drho = getvec('rke')
      oplot, x, drho
   endfor
end_plot
;
drho = getvec('bv')
!y.title = '!6Weight Function'
start_plot, n_win, 'weight'
   drho = getvec('wait')
   plot, x, drho, yrange=[-0.04, 0.04], xrange=[0.95, 1.]
   for imode=1, nmodes-1 do begin
      drho = getvec('wait')
      oplot, x, drho
   endfor
end_plot
skip_rest:
end
