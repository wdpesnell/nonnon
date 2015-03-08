      subroutine summary( info )
      implicit real*8(a-h,o-z)
      character*8 date_buf
      character*10 time_buf
c
      type puls_sum
         integer*2 :: mass_flag
         integer*4 :: ell, k, n_p, n_g
         real*4 :: omega_ad, omsq_n, nu_ad, period, acc_ad_w, acc_ad_e
	 real*4 :: rbar
         complex*8 :: del_om_qa
         complex*8 :: comega
         real*4 :: kappa, acc_na_k
      end type
      type(puls_sum):: info
c
      open (unit=41,file='Model_S_vals_test.sum',status='unknown',
     $      access='append',
     $      form='unformatted')
      call date_and_time( date=date_buf, time=time_buf )
      write(41) info, date_buf, time_buf
      close(41)
c      write(*,1000) date_buf, time_buf, info
      return
c
 1000 format(2a12,/,5(1x,i4),1p12e15.7)
      end
