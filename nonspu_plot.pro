pltcom, 'non001.plt'
omsq = getvec('omsq')
nmodes = n_elements(omsq)
k_p = findgen(nmodes) + 1
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
   plot, x, dror, yrange=[-0.1,0.1], /ystyle, xrange=[0.,1], /xstyle
;   oplot, x_Pek, 2.36*x_Pek^2 - 1.36, psym=1
   for i=0, 5 do begin
      dror = getvec('dr/r')
      oplot, x, dror
   endfor
   dror = getvec('dh/h')
   oplot, x, dror, line=5
   oplot, !x.crange, [0,0], line=1
end_plot, /notime
;
drho = getvec('bv')
!y.title = '!7dq!6/!7q!6'
start_plot, n_win, 'non_fig_3', /encap
   drho = getvec('drho')
   plot, x, drho ;, yrange=[-1., 1.]
;end_plot, /notime
;
end
