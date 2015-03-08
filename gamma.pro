function tridag, a, x
npts = n_elements(x)
beta = dblarr(npts)
gamm = dblarr(npts)
beta(0) = a(0,1)
gamm(0) = x(0)/beta(0)
for i=1, npts-1 do begin
   beta(i) = a(i,1) - a(i,0)*a(i-1,2)/beta(i-1)
   gamm(i) = (x(i) - a(i,0)*gamm(i-1))/beta(i)
endfor
;
v = dblarr(npts)
v(npts-1) = gamm(npts-1)
for i=npts-2, 0, -1 do v(i) = gamm(i) - a(i,2)*v(i+1)/beta(i)
return, v
end

n_win = 0
npts = 320
dr = 1.d0/double(npts+1)
r = dr + dindgen(npts+1)*dr
drint = dblarr(npts)
drint(0:npts-1) = r(1:npts)-r(0:npts-1)
rz = (r(1:npts)+r(0:npts-1))/2.d0
drz = dblarr(npts)
drz(0:npts-2) = rz(1:npts-1)-rz(0:npts-2)
drz(npts-1) = drz(npts-2)
l = 2
rl = double(l)
;dr = 1.d0
;
;       Power-law index for the source
;
k = 0
;
ap3 = dblarr(npts,3)
for i=1, npts-2 do begin
   ap3i1 = r(i)^(2*l+2)/(drint(i)*drz(i))/rz(i)^(2*l+2)
   ap3i3 = r(i+1)^(2*l+2)/(drint(i)*drz(i+1))/rz(i)^(2*l+2)
   ap3(i,0) = ap3i1
   ap3(i,1) =-(ap3i1 + ap3i3)
   ap3(i,2) = ap3i3
endfor
;
ap311 = r(0)^(2*l+2)/(drint(0)*drz(0))/(rz(0))^(2*l+2)
ap313 = r(1)^(2*l+2)/(drint(0)*drz(1))/(rz(0))^(2*l+2)
ap3(0,0) = 0.
ap3(0,1) =-(ap311 + ap313)
ap3(0,2) = ap313 + ap311
;
ap3n1 = r(npts-1)^(2*l+2)/(drint(npts-1)*drz(npts-1))/rz(npts-1)^(2*l+2)
ap3n3 = r(npts)^(2*l+2)/(drint(npts-1)*drz(npts-1))/rz(npts-1)^(2*l+2)
ap3(npts-1,0) = ap3n1
fp3 = (1.d0 - (rl+0.5)*(drz(npts-1)/r(npts)))/ $
      (1.d0 + (rl+0.5)*(drz(npts-1)/r(npts)))
ap3(npts-1,1) =-(ap3n1 + (1.d0-fp3)*ap3n3)
ap3(npts-1,2) = 0.d0
;
ap3(*,*) = dr^2*ap3(*,*)
;
;       The important thing in this problem is the surface boundary condition
;
source = dblarr(npts)
alpha = 3.d0
beta = 0.d0
if( k eq 0 ) then source(*) = alpha $
else source = alpha*rz^k
source(npts-1) = source(npts-1) + beta*ap3n3*(drz(npts-1)/r(npts))/ $
      (1.d0 + (rl+0.5)*(drz(npts-1)/r(npts)))
source(*) = source(*)*dr^2
;
k_size = size(k)
if( k_size(1) ge 4 ) then k_str = ', !12k!6 = ' + strtrim( string(format='(f5.1)', k),2) $
else                   k_str = ', !12k!6 = ' + strtrim( string(format='(i4)', k),2)
title_str = '!6' + strtrim( string(format='(i6)', npts),2) + $
            ' points, !7a!6 = ' + strtrim( string(format='(f7.2)', alpha),2) + $
            ', !7b!6 = ' + strtrim( string(format='(f7.2)', beta),2) + $
            ', !12l!6 = ' + strtrim( string(format='(i4)', l),2) + k_str
;
start_plot, n_win, 'gamm_tst'
   !p.title = title_str
   !x.title = '!8x = r/R!6!d*!n'
   !y.title = '!6AP3 Elements'
   yrange = [min(ap3), max(ap3)]
   plot, rz, ap3(*,1), yrange=yrange
   oplot, rz, ap3(*,0), line=5
   oplot, rz, ap3(*,2), line=3
end_plot
;
gamma = trisol( ap3(*,0), ap3(*,1), ap3(*,2), source, /double )
;gamm_2 = tridag( ap3, source )
start_plot, n_win, 'gamm_sol'
   !p.title = title_str
   !y.title = '!8f!6'
   if( abs(alpha) gt 0. ) then begin
      C_surf = ((2.*rl+k+3)/(2.*rl+1))*(beta/alpha + 1./(k+2.))
      exact = alpha*(rz(*)^(k+2)/(k+2.) - C_surf)/(2.*rl+k+3.)
   endif else begin
      exact = dblarr(npts)
      exact(*) = -beta/(2.*rl+1.)
   endelse
   yrange = [min([gamma,exact]), max([gamma,exact])]
   plot, rz, gamma, yrange=yrange
   oplot, rz, exact, line=5
end_plot
print, gamma(0), gamma(npts-1)
print, exact(0), exact(npts-1)
gamm_fit = poly_fit( rz, gamma, 2)
print, transpose(gamm_fit)
exact_fit = poly_fit( rz, exact, 2)
print, [transpose(exact_fit), gamm_fit(0)/exact_fit(0), gamm_fit(2)/exact_fit(2)]
;
start_plot, n_win, 'gamm_rat'
   !p.title = title_str
   !y.title = '!6(!8f!6 - Analytic)/!8f!6'
   ratio = (gamma-exact)/gamma
   plot, rz, ratio
end_plot
;
;       Test the current version of the matrices
;
ddr11 = dblarr(npts)
ddr11(0:npts-2) =-drz(1:npts-1)/(drz(0:npts-2)*(drz(0:npts-2)+drz(1:npts-1)))
ddr11(npts-1) = ddr11(npts-2)
ddr12 = dblarr(npts)
ddr12(0:npts-2) = (drz(1:npts-1)-drz(0:npts-2))/(drz(0:npts-2)*drz(1:npts-1))
ddr12(npts-1) = ddr12(npts-2)
ddr13 = dblarr(npts)
ddr13(0:npts-2) = drz(0:npts-2)/(drz(1:npts-1)*(drz(0:npts-2)+drz(1:npts-1)))
ddr13(npts-1) = ddr13(npts-2)
;
ap3o = dblarr(npts,3)
;
for i=1, npts-2 do begin
   ap3i1 = r(i)^(2)/(drint(i)*drz(i))/rz(i)^(2)
   ap3i3 = r(i+1)^(2)/(drint(i)*drz(i+1))/rz(i)^(2)
   ap3o(i,0) = ap3i1 +           ddr11(i)*2.d0*rl/rz(i)
   ap3o(i,1) =-(ap3i1 + ap3i3) + ddr12(i)*2.d0*rl/rz(i)
   ap3o(i,2) = ap3i3 +           ddr13(i)*2.d0*rl/rz(i)
endfor
;
ap311 = r(0)^(2)/(drint(0)*drz(0))/rz(0)^(2)
ap313 = r(1)^(2)/(drint(0)*drz(1))/rz(0)^(2)
ap3o(0,0) = 0.
ap3o(0,1) =-(ap311 + ap313) + ddr12(0)*2.d0*rl/rz(0)
ap3o(0,2) = ap313 + ap311   + ddr11(0)*2.d0*rl/rz(0) + ddr13(0)*2.d0*rl/rz(0)
;
ap3n1 = r(npts-1)^(2)/(drint(npts-1)*drz(npts-1))/rz(npts-1)^(2)
ap3n3 = r(npts)^(2)/(drint(npts-1)*drz(npts-1))/rz(npts-1)^(2)
ap3o(npts-1,0) = ap3n1 + ddr11(npts-1)*2.d0*rl/rz(npts-1)
fp3 = (1.d0 - (rl+0.5)*(drz(npts-1)/r(npts)))/ $
      (1.d0 + (rl+0.5)*(drz(npts-1)/r(npts)))
ap3o(npts-1,1) =-(ap3n1 + (1.d0-fp3)*ap3n3)  + ddr12(npts-1)*2.d0*rl/rz(npts-1) $
                + fp3*ddr13(npts-1)*2.d0*rl/rz(npts-1)
ap3o(npts-1,2) = 0.d0
;
ap3o(*,*) = dr^2*ap3o(*,*)
;
start_plot, n_win, 'gamm_tst'
   !p.title = title_str
   !x.title = '!8x = r/R!6!d*!n'
   !y.title = '!6Current AP3 Elements'
   yrange = [min(ap3o), max(ap3o)]
   plot, rz, ap3o(*,1), yrange=yrange
   oplot, rz, ap3o(*,0), line=5
   oplot, rz, ap3o(*,2), line=3
end_plot
;
if( k eq 0 ) then source(*) = alpha $
else source = alpha*rz^k
;
source(npts-1) = source(npts-1) + $
   beta*(ap3n3 + ddr13(npts-1)*2.d0*rl/rz(npts-1))*(drz(npts-1)/r(npts))/ $
   (1.d0 + (rl+0.5)*(drz(npts-1)/r(npts)))
source(*) = source(*)*dr^2
;
gamma_old = trisol( ap3o(*,0), ap3o(*,1), ap3o(*,2), source, /double )
;
start_plot, n_win, 'gamm_old'
   !p.title = title_str
   !y.title = '!8f!6'
   yrange = [min([gamma,exact]), max([gamma,exact,gamma_old])]
   plot, rz, gamma, yrange=yrange
   oplot, rz, gamma_old, line=3
   oplot, rz, exact, line=5
end_plot
print, gamma_old(0), gamma_old(npts-1)
skip_current:
end