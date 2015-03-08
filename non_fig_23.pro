pltcom, 'non_fig_23.rk
omsq = getvec('omsq')
nmodes = n_elements(omsq)
k_p = findgen(nmodes) + 1
wpes = getvec('wpes')
wsch = getvec('wsch')
x = getvec('x')
;
n_win = 0
;
npts = 32
dx = 1./float(npts)
x_Pek = dx/2. + findgen(npts)*dx
;
!p.title = '!6'
!x.title = '!8x = r/R!6!d*!n'
!y.title = '!7d!8r/r!6'
start_plot, n_win, 'non_fig_2', /encap
   dror = getvec('dr/r')
   plot, x, dror, yrange=[-1.5,3.5];, xrange=[0.95,1]
   oplot, x_Pek, 2.36*x_Pek^2 - 1.36, psym=1
   dror = getvec('dr/r')
   oplot, x, dror
   oplot, x_Pek, -2.19*x_Pek^2 + 3.19, psym=1
end_plot, /notime
;
drho = getvec('bv')
!y.title = '!7c!6'
start_plot, n_win, 'non_fig_3', /encap
   drho = getvec('gam')
   plot, x, drho ;, yrange=[-1., 1.]
   drho = getvec('gam')
   oplot, x, drho
   oplot, x_Pek, -21./14.*x_Pek^2*(x_Pek^2-1), psym=1
end_plot, /notime
;
end
