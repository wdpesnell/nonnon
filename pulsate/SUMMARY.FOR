      subroutine summary( info )
      implicit real*8(a-h,o-z)
      character*8 dat_buf, time_buf
c
      structure /puls_sum/
         integer*2 mass_flag
         integer*4 ell, k, n_p, n_g
         real*4 omega_ad, omsq_n, period, acc_ad_w, acc_ad_e
         complex*8 del_om_qa
         complex*8 comega
         real*4 kappa, acc_na_k
      end structure
      record /puls_sum/ info
c
      open (unit=41,file='summary.zms',status='unknown',access='append',
     $      form='unformatted')
      call date( dat_buf )
      call time( tim_buf )
      write(41) info, dat_buf, tim_buf
      close(41)
      return
c
 1000 format(5(1x,i4),1p11e15.7)
      end
