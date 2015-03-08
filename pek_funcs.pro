;       PEK_FUNCS.PRO
;
;       Double precision is necessary to resolve the g-mode branch
;
k = 1
l = 2

Gamma_1 = 5.d0/3.d0
D_k = Gamma_1*k*(k+l+0.5d0) - 2.d0
disc = D_k^2 + double(l*(l+1))
Omega_2_p = D_k + sqrt(disc)
Omega_2_g = D_k - sqrt(disc)
print, l, Omega_2_g, Omega_2_p
;
B = -6.d0 - 4.d0*l + 2.d0/gamma_1*( 4.d0 + Omega_2_p - double(l*(l+1))/Omega_2_p)
B = double(2*(k-1)*(2*(k-1)+5+2*l))
cj = dblarr(abs(k))
cj[0] = 1.d0
dcj = dblarr(abs(k))
dcj[0] = double(l)
for j=1, abs(k)-1 do begin
   cj[j] = cj[j-1]*double(2*(j-1)*(2*(j-1)+5+2*l)-B)/double((2*(j))*(2*(j-1)+3+2*l))
   dcj[j] = cj[j]*double(2*(j)+l)
endfor
print, cj
npts = 120
dx = 1./double(npts)
x = dx/2.d0 + dindgen(npts)*dx
alfa = x^l*poly( x^2, cj )
dalfa_dx = x^(l-1)*poly( x^2, dcj )
dror = ( (x^2-1.d0)*(dalfa_dx + double(l*(l+1))*alfa/Omega_2_p/x) + alfa*2.d0*x)/(B+4.d0*l+6.d0)/x
print, transpose(POLY_FIT(x^2, dror/dror(npts-1), 1))
n_win = 0
!p.thick = 1.
!p.title = '!6Pekeris Model Solutions'
!x.title = '!8x!6 = !8r/R!d*!n!6'
!y.title = '!7a!6 (' + strtrim(string(format='(i4)',l),2) + $
                      ',' + strtrim(string(format='(i4)',k),2) + ')'
start_plot, n_win, 'pek_alfa'
   plot, x, alfa/alfa(npts-1)
   oplot, x, -5./2.*(3. - 7.*x^2)/10., psym=1
end_plot
!p.title = '!6Pekeris Model Solutions'
!x.title = '!8x!6 = !8r/R!d*!n!6'
!y.title = '!7d!8r/r!6 (' + strtrim(string(format='(i4)',l),2) + $
                      ',' + strtrim(string(format='(i4)',k),2) + ')'
start_plot, n_win, 'pek_dror'
   plot, x, 35./8.*(1. - 18./5.*x^2 + 99./35.*x^4), psym=1
;   plot, x, -5./2.*(1. - 7./5.*x^2), psym=1
   oplot, x, dror/dror(npts-1)
end_plot
end