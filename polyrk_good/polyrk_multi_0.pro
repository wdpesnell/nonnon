;files = 'polyrk_0000_' + ['00011', '00100', '01000'] + '.plt'
files = 'polyrk_3000_' + ['00200', '01000'] + '.plt'
n_files = n_elements( files )
;
npts = intarr(n_files)
for i_file=0, n_files-1 do begin
   print, files[i_file]
   pltcom, files[i_file], /noscan
   x = getvec('x')
   npts[i_file] = n_elements(x)
   close, 1
endfor
;
x = fltarr(max(npts), n_files)
dror = fltarr(max(npts), n_files)
dhor = fltarr(max(npts), n_files)
rz = fltarr(max(npts), n_files)
drho = fltarr(max(npts), n_files)
;
for i_file=0, n_files-1 do begin
   pltcom, files[i_file], /noscan
   vec = getvec('x')
   x[0:npts[i_file]-1,i_file] = vec
   vec = getvec('rz')
   rz[0:npts[i_file]-2,i_file] = vec
   vec = getvec('dr/r')
   dror[0:npts[i_file]-1,i_file] = vec
   vec = getvec('dh/h')
   dhor[0:npts[i_file]-2,i_file] = vec
   vec = getvec('drho')
   drho[0:npts[i_file]-2,i_file] = vec
   close, 1
endfor
;
n_win = 0
;
my_lines = [0, 5, 3, 4, 2]
n_my_lines = n_elements(my_lines)
;
!x.title = '!8x = r/R!6!d*!n'
!y.title = '!7d!8r/r!6'
start_plot, n_win, 'dr_r'
   plot, [1], [1], yrange=[min(dror), max(dror)], xrange=[0,1.05], /xstyle
   for i_file=0, n_files-1 do begin
      oplot, x[0:npts[i_file]-1,i_file], dror[0:npts[i_file]-1,i_file], line=my_lines[i_file mod n_my_lines]
   endfor
end_plot
;
!y.title = '!7d!8h/h'
start_plot, n_win, 'dh_h'
   plot, [1], [1], yrange=[min(dhor), max(dhor)], xrange=[0,1.05], /xstyle
   for i_file=0, n_files-1 do begin
      oplot, rz[0:npts[i_file]-2,i_file], dhor[0:npts[i_file]-2,i_file], line=my_lines[i_file mod n_my_lines]
   endfor
end_plot
;
!y.title = '!7dq/q!6'
start_plot, n_win, 'dh_h'
   plot, [1], [1], yrange=[min(drho), max(drho)], xrange=[0,1.05], /xstyle
   for i_file=0, n_files-1 do begin
      oplot, rz[0:npts[i_file]-2,i_file], drho[0:npts[i_file]-2,i_file], line=my_lines[i_file mod n_my_lines]
   endfor
end_plot
end
;
;
drho = getvec('bv')
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
end
