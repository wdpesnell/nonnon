      implicit real*8(a-h,o-z)
c
c          comparison of omsq with the eigenvalues given by robe.
c
      dimension omsq1(21),omsq2(21),omsq3(21),omsq33(21),
     $         omsq35(21),omsq4(21),omsq15(21), omsq0(21)
      data omsq1/-1.583D-2,-1.898D-2,-2.319D-2,-2.899D-2,-3.731D-2,
     $      -4.989D-2,-7.029D-2,-0.1070D0,-.1844D0,-.4039D0,1.997D0,
     $      12.41D0,31.32D0,57.11D0,89.44D0,128.2D0,173.2D0,224.5D0,
     $      282.2D0,346.d0,416.d0/
      data omsq2/4.061D-2,4.840D-2,5.868D-2,7.265D-2,9.238D-2,.1215D0,
     $      .1672D0,0.2452D0,0.3957D0,0.7510D0,4.151D0,15.41D0,32.10D0,
     $    54.17D0,81.51D0,114.0D0,151.6D0,194.3D0,242.d0,295.d0,353.d0/
      data omsq3/ 0.429D0,0.510D0,0.616D0,0.7588D0,0.9584D0,1.248D0,
     $      1.694D0,2.430D0,3.771D0,6.553D0,10.90D0,20.35D0,35.63D0,
     $  55.29D0,79.23D0,107.4D0,139.8D0,176.5D0,217.d0,262.d0,312.d0/
      data omsq35/1.40D0,1.66D0,2.d0,2.464D0,3.099D0,4.023D0,5.417D0,
     $      7.655D0,11.39D0,16.13D0,21.55D0,27.91D0,39.68D0,57.88D0,
     $      80.39D0,106.8D0,137.2D0,171.d0,210.d0,252.d0,297.d0/
      data omsq4/ 6.2D0,7.34D0,8.82D0,10.78D0,13.44D0,17.02D0,20.48D0,
     $      23.99D0,30.67D0,36.79D0,45.77D0,56.18D0,67.75D0,83.82D0,
     $      102.2D0,116.6D0,140.3D0,171.d0,206.d0,245.d0,288.d0/
c
      lval = 2
      rl = dfloat(lval)
      rl1 = rl*(rl+1.d0)
c
c          pekeris model
c
c
c          The f mode
c
      omsq0(11) = 2.d0*rl*(rl-1.d0)/(2.d0*rl+1.d0)
      do k=1, 10
c
c          The g mode branch has an altered node count.
c
            dn =-2.d0 + dfloat(k)*(rl + 0.5d0 + dfloat(k))*(5.d0/3.d0)
            omsq0(11-k) = dn - dsqrt(dn**2 + rl1)
      enddo
c
c          The p mode branch
c
      do k=1, 10
         dn =-2.d0 + dfloat(k)*(rl + 0.5d0 + dfloat(k))*(5.d0/3.d0)
         omsq0(k+11) = dn + dsqrt(dn**2 + rl1)
      enddo
c
      open(16, file='robe_eig.tex', status='unknown', 
     $    carriagecontrol='list')
      write(16,1000)
      do k=1, 11
         write(16,1500) k-11, omsq0(k), omsq1(k), omsq2(k), omsq3(k),
     $      omsq35(k), omsq4(k)
      enddo
      do k=12, 15
         write(16,1510) k-11, omsq0(k), omsq1(k), omsq2(k), omsq3(k),
     $      omsq35(k), omsq4(k)
      enddo
      do k=16, 21
         write(16,1520) k-11, omsq0(k), omsq1(k), omsq2(k), omsq3(k),
     $      omsq35(k), omsq4(k)
      enddo
      write(16,1750)
      write(16,2000)
      stop
c
 1000 format(
     $ '\documentclass[12pt]{article}', /
     $ '\begin{document}',/,'\newcounter{T1}',/
     $ '\begin{center}',/
     $ '\addtocounter{T1}{1}', /
     $ '\begin{tabular}[hbt]{rllllll}',/
     $ '\multicolumn{7}{c}{Table~\arabic{T1}: Eigenvalues of ',
     $     'Polytropes (Robe, 1968)} \\',/
     $ '\hline \hline',/
     $ '\multicolumn{1}{c}{$k$} & ',/
     $ '\multicolumn{1}{c}{$\omega_0^2$} & ',/
     $ '\multicolumn{1}{c}{$\omega_1^2$} & ',/
     $ '\multicolumn{1}{c}{$\omega_2^2$} & ',/
     $ '\multicolumn{1}{c}{$\omega_3^2$} & ',/
     $ '\multicolumn{1}{c}{$\omega_{3.5}^2$} & ',/
     $ '\multicolumn{1}{c}{$\omega_4^2$} \\ \hline')
 1500 format(1x,i3,' & ',5(f10.4,' & '),f8.2, ' \\' )
 1510 format(1x,i3,' & '5(f8.2,' & '),f8.2, ' \\' )
 1520 format(1x,i3,' & '5(f8.1,' & '),f8.1, ' \\' )
 1750 format('\hline',/,'\end{tabular}', /, '\end{center}')
 2000 format('\end{document}' )
      end
