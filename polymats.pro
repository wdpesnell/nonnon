openr, lunit, 'polout.rk', /get_lun
a = ' '
;
;       Scan for the 'zones' line for number of points and the 'theta' line to start output read
;
for i=0, 2000 do begin
;another_line:
   readf, lunit, a
   print, a
   if( strlen(a) gt 0 ) then begin
      if( strmid(a, 7, 5) eq 'zones' ) then npts = fix(strmid(a,1,5)) $
      else if( strmid(a, 1, 5) eq 'fnpol') then fnpol = float(strmid(a,8,5)) $
      else if( strpos(a, 'theta') ge 0 ) then goto, aligned
   endif
;goto, another_line
endfor
;
aligned:
print, 'Polytropic index (n) is ', fnpol
index_string = ', !8n!6 = ' + strtrim( string(format='(f5.2)', fnpol), 2)
;
;       Read the initial model part of the file
;
;
matread = fltarr(11,npts+2)
readf, lunit, matread
matread = transpose(matread)
;
x = matread(*,4)
;
readf, lunit, a
readf, lunit, a
readf, lunit, a
;
matread = fltarr(11,npts+1)
readf, lunit, matread
;
matread = transpose(matread)
;
ag1 = matread(*,1:3)
drdm = matread(*,4:5)
ag3 = matread(*,6:7)
ag4 = matread(*,8:9)
ah4 = matread(*,10:10)
;
readf, lunit, a
readf, lunit, a
readf, lunit, a
;
matread = fltarr(10,npts)
readf, lunit, matread
matread = transpose(matread)
;
ah1 = matread(*,1:2)
ah3 = matread(*,3:3)
ap1 = matread(*,4:5)
ap3 = matread(*,6:8)
ap4 = matread(*,9:9)
;
n_win = 0
;
!x.title = '!8x = r/R!6!d*!n'
!y.title = ' '
!p.title = '!6DRDM Elements' + index_string
yrange = [min(drdm), max( drdm )]
yrange = [-1.e5, 1.e5]
start_plot, n_win, 'drdm'
   plot, x, drdm(*,0), yrange=yrange
   oplot, x, drdm(*,1), line=2
end_plot
;
!x.title = '!8x = r/R!6!d*!n'
!y.title = ' '
!p.title = '!6AG1 Elements' + index_string
yrange = [min(ag1), max( ag1 )]
start_plot, n_win, 'ag1'
   plot, x, ag1(1:npts,0), yrange=yrange
   oplot, x, ag1(0:npts,1), line=1
   oplot, x, ag1(0:npts-1,2), line=2
end_plot
;
!x.title = '!8x = r/R!6!d*!n'
!y.title = ' '
!p.title = '!6AG3 Elements' + index_string
yrange = [min(ag3), max( ag3 )]
start_plot, n_win, 'ag3'
   plot, x, ag3(1:npts,0), yrange=yrange
   oplot, x, ag3(0:npts-1,1), line=3
end_plot
;
!x.title = '!8x = r/R!6!d*!n'
!y.title = ' '
!p.title = '!6AG4 Elements' + index_string
yrange = [min(ag4), max( ag4 )]
yrange = [-1.e6, 1.e6]
start_plot, n_win, 'ag4'
   plot, x, ag4(1:npts,0), yrange=yrange
   oplot, x, ag4(0:npts-1,1), line=2
end_plot
;
!p.title = '!6AH1 Elements' + index_string
yrange = [min(ah1), max( ah1 )]
yrange = [-1.e5,1.e5]
start_plot, n_win, 'ah1'
   plot, x, ah1(*,0), yrange=yrange
   oplot, x, ah1(*,1), line=2
end_plot
;
;       Boringly constant for polytropes
;
;!p.title = '!6AH3 Elements' + index_string
;yrange = [min(ah3), max(ah3)]
;start_plot, n_win, 'ah3'
;   plot, x, ah3(*,0), yrange=yrange
;end_plot
;
!x.title = '!8x = r/R!6!d*!n'
!y.title = ' '
!p.title = '!6AH4 Elements' + index_string
yrange = [min(ah4), max( ah4 )]
yrange = [0, 1.e5]
start_plot, n_win, 'ah4'
   plot, x, ah4(*), yrange=yrange
end_plot
;
!p.title = '!6AP1 Elements' + index_string
yrange = [min(ap1), max( ap1 )]
start_plot, n_win, 'ap1'
   plot, x, ap1(*,0), yrange=yrange
   oplot, x, ap1(*,1), line=1
end_plot
!p.title = '!6AP3 Elements' + index_string
yrange = [min(ap3), max(ap3)]
start_plot, n_win, 'ap3'
   plot, x, ap3(1:npts-1,0), yrange=yrange
   oplot, x, ap3(*,1), line=2
   oplot, x, ap3(0:npts-2,2), line=3
end_plot
!p.title = '!6AP4 Elements' + index_string
yrange = [min(ap4), max(ap4)]
start_plot, n_win, 'ap4'
   plot, x, ap4(*,0), yrange=yrange
end_plot
;
free_lun, lunit
end