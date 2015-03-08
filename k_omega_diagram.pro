;       K_OMEGA_DIAGRAM.PRO
;
;       Double precision is necessary to resolve the g-mode branch
;
;openr, lunit, 'polfreqs_n3p25.txt', /get_lun
filename = 'polyrk_4500_06000.txt'
openr, lunit, filename, /get_lun
aread = ' '
i = 0
while( not eof(lunit) ) do begin
   readf, lunit, aread
   i = i+1
endwhile
free_lun, lunit
n_modes = i-1
rec = {l: 0, n: 0, np: 0, ng: 0, omsq: 0.}
freqs = replicate( rec, n_modes)
openr, lunit, filename, /get_lun
readf, lunit, fnpol                   ; polytropic index is first line
print, fnpol
for i=0, n_modes-1 do begin
   readf, lunit, aread
   print, aread
   vals = strsplit( aread, ' ', /extract)
   freqs[i].L = fix(vals[0])
   freqs[i].n = fix(vals[1])
   freqs[i].np = fix(vals[2])
   freqs[i].ng = fix(vals[3])
   freqs[i].omsq = float(vals[4])
endfor
free_lun, lunit
safe = where( freqs.l eq 0, n_safe)
if( n_safe gt 0 ) then begin
   freqs[safe].omsq = freqs[safe].omsq*4./3.
   freqs[safe].n = freqs[safe].n + 1
endif
Gamma_1 = 5.d0/3.d0
n_win = 1
start_plot, n_win, 'perk_omsq'
   !p.thick = 2.
   !x.thick = !p.thick
   !y.thick = !p.thick
   !p.charthick = 1.5
   !p.charsize = 1.5
   !p.title = '!6Polytrope Model, !8n!6 = ' + string(fnpol,format='(f5.3)')
   !x.title = '!12l!6'
   !y.title = '!7x!8!dn,!12l!6!n'
   plot, [1], [1], /nodata, xrange=[-1,22], /xstyle, yrange=[5, 50], /ystyle
   for i=0, n_modes-1 do begin
      oplot, freqs.l, sqrt(freqs.omsq), psym=1, symsize=2
   endfor
   n_max = max(freqs.n)
   n_min = min(freqs.n)
   for i=n_min, n_max do begin
      safe = where( freqs.n eq i, n_safe)
      if( n_safe gt 1 ) then begin
         oplot, freqs[safe].l, sqrt(freqs[safe].omsq), psym=-4
         xyouts, 20.5, max(sqrt(freqs[safe].omsq)), strtrim(string(i),2)
      endif
   endfor
end_plot
end
