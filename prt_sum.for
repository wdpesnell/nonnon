      program prt_sum
      implicit real*8(a-h,o-z)
      real*4 mass
      integer*2 i
      character*8 date_buf
      character*10 time_buf
      parameter ( pi=3.141592654d0, pi2=2.d0*pi, secpday=86400.d0 )
c
      type puls_sum
         integer*2 :: mass_flag
         integer*4 :: ell, k, n_p, n_g
         real*4 :: omega_ad, omsq_n, nu_ad, period, acc_ad_w, acc_ad_e
         complex*8 :: del_om_qa
         complex*8 :: comega
         real*4 :: kappa, acc_na_k
      end type
      type(puls_sum):: info
c
      open (unit=41,file='Model_S_vals.sum',status='old',
     $      action='read',
     $      form='unformatted')
      open(unit=11,file='summary.lis',status='unknown')
      write(11,1000)
      do i=1,1000
         read(41,end=20) info, date_buf, time_buf
         mass = float(info%mass_flag)/10.
         period_na = pi2/real(info%comega)/secpday
         write(11,1050) mass, info%k, 
     $         info%period/secpday, 
     $         info%acc_ad_w, info%acc_ad_e, info%period,
     $         info%nu_ad, info%acc_na_k
         write(*,1050) mass, info%k, 
     $         info%period/secpday, 
     $         info%acc_ad_w, info%acc_ad_e, info%period,
     $         info%nu_ad*1.e3, info%acc_na_k
         write(*,1100) date_buf, time_buf
      enddo
  20  continue
c      do i=1,1000
c         read(41,end=20) info, date_buf, time_buf
c         mass = float(info%mass_flag)/10.
c         write(11,1010) mass, info%ell, info%k, info%n_p,
c     $         info%n_g, info%omega_ad, info%omsq_n, info%period, 
c     $         info%acc_ad_w, info%acc_ad_e, date_buf, time_buf
c      enddo
c  20  continue
c      rewind(41)
c
c      write(11,2000)
c      do i=1,1000
c         read(41,end=50) info, date_buf, time_buf
c         mass = float(info%mass_flag)/10.
c         write(11,2010) mass, info%ell, info%k, info%n_p,
c     $         info%n_g, info%del_om_qa, info%comega, 
c     $         info%kappa, info%acc_na_k
c      enddo
c  50  continue
      close(11)
      close(41)
      stop
c
 1000 format(1x,'Adiabatic frequency information',/,
     $  1x,'M/Ms   l    k   n_p  n_g      omega          Omsq',5x,
     $       2x,'    Period          err_w          err_e')
 1010 format(1x,f5.1,4(1x,i4),1p5e15.7,2(1x,a))
 1050 format(1x,f5.1,' & ',i4,' & ',1pg12.5,2(' & ',e8.1),' & ',
     $          g12.5,' & ',e10.3,' & ',e8.1,' \\')
 1100 format(1x,'Created on ',a8,' at ',a10)
c
 2000 format(//,1x,'Nonadiabatic information',/,
     $ 1x,'M/Ms   l    k   n_p  n_g   Real(omega_QA)'
     $ ' Imag(omega_QA) Real(omega_NA) Imag(omega_NA)     Kappa'
     $ '         error_k')
 2010 format(1x,f5.1,4(1x,i4),1p6e15.7)
c
      end
