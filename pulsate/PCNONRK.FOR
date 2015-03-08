c
c lanl modifications: xeq nonrk
c*byp3
c*rcft i=%me, x=%mex
c*destroy b%me m%me l%me
c      program nonrk(tty,tape5=tty,tape6=tty,cjhinp,tape28=cjhinp,
c     $   cjhout,tape11=cjhout,cjhplt,tape21=cjhplt,sumout,tape1=sumout)
      program nonrk
      implicit real*8(a-h,o-z)
      character*4 ititl
      logical switch, iray, nconv
c
c     program to solve full fourth order fluid system via runge kutta
c
c
c      setting verg and eps within program
c
      parameter ( verg=1.00e-5, eps=1.00e-6 )
      parameter ( nmax=600, nmax10=10*nmax )
      common g(nmax),rho(nmax),x(nmax),yliq(4,nmax),weight(4,nmax)
      common/misc/   grav,pi,pi4,p43,eigt,y3t
      common/pulsst/ rl1,l,lindex,nsurf
      common/ray/    iray
      common/rs/     r(nmax)
      common/modect/ y1m(nmax10),y2m(nmax10),nfine,nodes1,nodes2,modep
      common/savel/  xi(nmax),c1(4,nmax),y1(nmax),c2(4,nmax),
     $ y2(nmax),c3(4,nmax),qy3(nmax),c4(4,nmax),y4(nmax),c5(4,nmax),
     $ y5(nmax),c6(4,nmax),m
      dimension rint(4,nmax),ititl(13)
      dimension stroot(100),sty3(100),stwait(nmax)
c
c          files in this program.
c  unit   routine    name                purpose
c    1      main    sumout              summary output file
c    5      main    sys$input           terminal input file
c    6      main    sys$output          terminal output file
c   11      main    cjh_out.dat         primary output file
c
c   28     readin   cjhinp              contains the model quantities
c                                          for use by the stability
c                                          analyses as dumped by the
c                                          program evmdplt or modelk.
c   21     pltint   cjhplt              has vectors of pulsational
c                                          quantities in a form suitable
c                                          for pltcmd program.
c
c       call filerep
c
      open (unit=1,file='sumout',status='unknown')
      open (unit=11,file='cjhout',status='unknown')
c
c          A useful hint, when IRAY is true the weight functions
c       are calculated in the routine RFKLIQ.
c
      pi = 3.1415926535898d0
      pi4 = 4.d0*pi
      p43 = pi4/3.d0
      grav = 6.6723d-8
c
      call readin
      nsurf = nsurf - 1
      write(1,2000) nsurf
   2  continue
      write(6,2002)
      read(5,*,end=900,err=900) l
      if( l .le. 0 ) stop
      write(11,2003) l,verg,eps
      rl = float(l)
      rl1 = rl*(rl+1.d0)
      lindex = 2 - l
c
      write(6,2004)
      read(5,*,end=900,err=900) info
   3  continue
      if( info .gt. 0 ) then
c
c          calculate surface discriminant only
c       (in cowling approximation)
c
         write(6,2050)
         read(5,*,end=900) permin,permax,nper
         delper = (permax-permin)/nper
         nper1 = nper+1
         nroot = 0
         iray = .false.
         disc2 = zero
         do 210 idisc=1,nper1
            period = permin + delper*float(idisc-1)
            eig = (2.d0*pi/period)**2
            eigt = eig
            call bump(disc)
            disc = disc/(r(nsurf)**lindex)
            write(11,2103) period,disc
            if( idisc .eq. 1 ) then
               disc2 = disc
            elseif( disc*disc2 .lt. 0. ) then
               nroot = nroot + 1
               stroot(nroot) = period - delper*disc/(disc-disc2)
               sty3(nroot) =-1.e-2
               disc2 = disc
            endif
 210     continue
         write(6,2100)
         write(6,2101) nroot,permin,permax,l
         write(11,2100)
         write(11,2101) nroot,permin,permax,l
         if( nroot .le. 0 ) goto 2
         write(6,2102) (stroot(i),i=1,nroot)
         write(11,2102) (stroot(i),i=1,nroot)
c
c          end of discriminant calculation
c
c          given a first guess at the eigenvalue and the inner value
c       of y3, attempt to converge on the period using a two-dimensional
c       newton-rapheson method.
c
      else
         write(6,1001)
         read(5,*,end=900) period,y3
         if( period .le. 0.d0 ) goto 2
         write(11,1002)  period
         stroot(1) = period
         sty3(1) = y3
         nroot = 1
      endif
      do 60 iroot=1,nroot
         y3 = sty3(iroot)
         period = stroot(iroot)
c
         omsq = (2.d0*pi/period)**2
         nconv = .false.
         iray = .false.
         switch = .false.
         nserch = 20
         do 200 ntry=1,nserch
            eigt = omsq
            y3t = y3
            call grind(b1,b2)
            eigt = (1.d0+eps)*omsq
            call grind(dum1,dum2)
            deigb1 = dum1-b1
            deigb2 = dum2-b2
            eigt = omsq
            y3t = (1.d0+eps)*y3
            call grind(dum1,dum2)
            dy3b1 = dum1-b1
            dy3b2 = dum2-b2
            det = (deigb1*dy3b2-deigb2*dy3b1)
            dume = b2*dy3b1-b1*dy3b2
            dumy = b1*deigb2-b2*deigb1
            deig = eps*omsq*dume/det
            dy3 = eps*y3*dumy/det
            write(11,1003) ntry,omsq,deig,y3,dy3
            write(11,1020) b1,b2,yliq(3,nsurf),yliq(1,nsurf)
            write(6,1003) ntry,omsq,deig,y3,dy3
            write(6,1020) b1,b2,yliq(3,nsurf),yliq(1,nsurf)
            adeig = abs(deig/omsq)
            ady3 = abs(dy3/y3)
            if( adeig.lt.verg .and. ady3.lt.verg ) nconv = .true.
            epseig = deig/omsq
            epsy3 = dy3/y3
            if( abs(epseig) .gt. 0.1e0 ) deig = 0.1e0*sign(omsq,deig)
c
c          The value of dy3 must be limited so as to allow y3 to
c       change sign.
c
            if( ady3 .gt. 1.2e0 ) then
               if( .not. switch .and. dy3*y3 .lt. 0.d0 ) then
                  dy3 = 2.d0*sign(y3,dy3)
                  switch = .true.
               else
                  dy3 = 0.1e0*sign(y3,dy3)
                  switch = .false.
               endif
            elseif( ady3 .gt. 0.1e0 ) then
               dy3 = 0.1e0*sign(y3,dy3)
               switch = .false.
            endif
            omsq = omsq + deig
            y3 = y3 + dy3
            if( nconv ) goto 201
 200     continue
         write(11,1004)
         write(6,1004)
         goto 60
c
c          converged to an eigenvalue, calculate and print the mode
c       identification, weight functions and eigenvectors.
c
 201     continue
         eigt = omsq
         y3t = y3
         iray = .true.
         call fingrn(b1,b2)
         period = (2.d0*pi)/sqrt(abs(omsq))
         write(11,*) '************************************************'
         write(11,1005) period,omsq,y3
         write(6,1005) period,omsq,y3
         write(11,1007) b1,b2
         write(1,*) '************************************************'
         write(1,1005) period,omsq,y3
         write(1,1007) b1,b2
c
c            establish mode number
c
         call modeid
         write(1,1010) modep,nodes1,nodes2
         write(11,1010) modep,nodes1,nodes2
         write(6,1010) modep,nodes1,nodes2
c
         rke = weight(1,nsurf)
         eigtry = (weight(2,nsurf)+weight(3,nsurf)+weight(4,nsurf))/rke
         pertry = (2.d0*pi)/sqrt(abs(eigtry))
         qch = abs((eigtry-omsq)/omsq)
         qrke = p43*weight(1,nsurf)*omsq/2.d0
         robe = omsq/(pi*y1(nsurf)*p43)
         write(1,1011) eigtry,pertry,qch,qrke,robe
         write(11,1011) eigtry,pertry,qch,qrke,robe
         write(6,1011) eigtry,pertry,qch,qrke,robe
c
c          copy integrated t,c,n,g to array rint(4,nmax). thus, the
c       array rint will contain the integrated weight functions. after
c       call to eigenf, weight will contain the kernel of the weight
c       function as a function of position.
c
         do 10 iint=1,nsurf
            do 10 jint=1,4
               rint(jint,iint) = weight(jint,iint)/(omsq*rke)
  10     continue
         call eigenf(stwait)
         write(11,1030)
         do 30 i=1,nsurf
            wint = rint(2,i)+rint(3,i)+rint(4,i)
            write(11,1021) i,r(i),yliq(1,i),yliq(2,i),rint(1,i),
     $                     wint
  30     continue
         write(11,1022)
         do 20 i=1,nsurf
            write(11,1023) i,yliq(3,i),yliq(4,i),weight(1,i),
     $                     weight(2,i),weight(3,i),weight(4,i)
  20     continue
         write(11,1024)
         write(6,1024)
c
c          write the eigenvectors to the plot file.
c
         ititl(1) = 'dr/r'
         ititl(2) = 'y2  '
         ititl(3) = 'y3  '
         ititl(4) = 'y4  '
         ititl(5) = 'ke  '
         ititl(6) = 'cyr '
         ititl(7) = 'nyr '
         ititl(8) = 'gyr '
         ititl(9) = 'kei '
         ititl(10) = 'cyri'
         ititl(11) = 'nyri'
         ititl(12) = 'gyri'
         ititl(13) = 'wait'
         do 90 j=1,4
            write(21,9000) nsurf,ititl(j),(yliq(j,i),i=1,nsurf)
  90     continue
         do 95 j=1,4
            write(21,9000) nsurf,ititl(j+4),(weight(j,i),i=1,nsurf)
  95     continue
         do 96 j=1,4
            write(21,9000) nsurf,ititl(j+8),(rint(j,i),i=1,nsurf)
  96     continue
         write(21,9000) nsurf,ititl(13),(stwait(i),i=1,nsurf)
         do 97 i=1,nsurf
            yliq(2,i) = yliq(2,i)*g(i)/(omsq*r(i))
            yliq(3,i) = yliq(3,i)*g(i)*r(i)
  97     continue
         ititl(2) = 'dh/h'
         ititl(3) = 'gam '
         do 98 j=2,3
            write(21,9000) nsurf,ititl(j),(yliq(j,i),i=1,nsurf)
  98     continue
         y3 =-y3
  60  continue
      goto 3
c
c          termination due to end-of-file on any read.
c
 900  continue
      stop
c
 1001 format(1x,'Enter period (secs.) and central value of y3',/1x,
     $          '(enter a period of 0 or <crtl-z> to stop) '  )
 1002 format(1x,15hguessed period=,1pe14.6,4h sec)
 1003 format(6h ntry=,i3,5h eig=,1pe12.5,6h deig=,e12.5,4h y3=,e12.5,
     $  5h dy3=,e12.5)
 1004 format(1x,14h not converged)
 1005 format(1x,18h final values are:,/,1x,8hperiod =,1pe16.9,5h eig=,
     $   e16.9,5h y3 =,e11.4)
 1007 format(1x,28h boundary conditions are b1=,1pe11.3,4h b2=,e11.3)
 1010 format(1x,21h  phase diagram mode ,i4,2x,i4,14h nodes in y1 ,,
     $   i4,13h nodes in y2 )
 1011 format(1x,19hintegrated sig**2 =,1pe12.4,19h integrated period=,
     $   e12.4,/,2x,18h error in sig**2 =,e12.4,/,
     $   1x,21hkinetic energy amp. =,e12.4,/,
     $   1x,18hnormalized omsq = ,e12.4)
 1020 format(6x,8h and b1=,1pe12.4,4h b2=,e12.4,4h y3=,e12.4,
     $ 4h y1= ,e12.4)
 1021 format(1x,i4,f10.6,2(1pe12.4),3x,2(1pe12.4))
 1022 format(/1x,37hy3  , y4  , t(i) , c(i) , n(i) , g(i))
 1023 format(1x,i4,6(1pe12.4))
 1024 format(//)
 1030 format(1x,37hn, r/r* ,  y1  ,  y2  ,   tint , wint)
c
 2000 format(1h1,/1x,47hOutput from CJH runge-kutta pulsation analysis.,
     $           /,1x,i5,12h zone model.)
 2002 format(1x,'Enter l: '  )
 2003 format(1x,3h l=,i4,6h verg=,1pe10.2,5h eps=,e10.2)
 2004 format(1x,'Calculate discriminant or guess periods?'/
     $       3x,'Enter 1 for discriminant, -1 for guessing: '  )
c
 2100 format(1x,26hdiscriminant on file tape6)
 2101 format(1x,i4,21h roots found between ,1pe10.3,5h and ,e10.3,
     $      7h for l=,i4)
 2102 format(1x,1p,4e12.5)
 2103 format(1x,1p,2e14.6)
c
 2050 format(1x,52henter min. period, max. period, and number of points)
 9000 format(i4,10x,a4,/,(1p,6e12.4) )
      end
      subroutine readin
      implicit real*8(a-h,o-z)
c
      parameter ( nmax=600 )
      common g(nmax),rho(nmax),x(nmax),yliq(4,nmax),weight(4,nmax)
      common/pulsst/ rl1,l,lindex,nsurf
      common/rs/     r(nmax)
      common/savel/  xi(nmax),c1(4,nmax),y1(nmax),c2(4,nmax),
     $ y2(nmax),c3(4,nmax),y3(nmax),c4(4,nmax),y4(nmax),c5(4,nmax),
     $ y5(nmax),c6(4,nmax),m
      dimension rnorm(nmax)
c
c permanently assign itapec to tape28
c
      itapec=28
c
c      open (unit=itapec,file='cjhinp',status='old')
c
      read(itapec,1000) blum,nmod
      read(itapec,1010) amass
      write(11,1001) amass,nmod,blum
      write(6,1001) amass,nmod,blum
      read(itapec,1011) nsurf
      if( nsurf .gt. nmax ) goto 700
      read(itapec,1007) (x(i),r(i),g(i),rho(i),i=1,nsurf)
      read(itapec,5) m
      do 100 i=1,m
         read(itapec,6) xi(i),y1(i),y2(i),y3(i)
         read(itapec,7) y4(i),y5(i)
  100 continue
      close (unit=itapec)
      call consl(xi,y1,m,c1)
      call consl(xi,y2,m,c2)
      call consl(xi,y3,m,c3)
      call consl(xi,y4,m,c4)
      call consl(xi,y5,m,c5)
      call consl(xi,r ,m,c6)
      call pltint(amass,blum,nmod,nsurf,r,y3)
      return
 700  continue
      write(11,7000) nsurf
      write(6,7000) nsurf
      stop
c
 1000 format(e12.5,i5)
 1010 format(e12.4)
 1001 format(1x,f7.2,' solar mass, no. points',i4,7h log L=,f8.4)
 1011 format(i4)
 1007 format(4e20.12)
   5  format(i5)
   6  format(4e20.12)
   7  format(2e20.12)
 7000 format(1x,37hfrom readin...too many points, nsurf=,i4)
      end
      subroutine bump(b2)
      implicit real*8(a-h,o-z)
c
c          find the outer boundary condition satisfaction as a
c       function of frequency**2 in the cowling approximation.
c
      parameter ( nmax=600 )
      common g(nmax),rho(nmax),x(nmax),yliq(4,nmax),weight(4,nmax)
      common/misc/   grav,pi,pi4,p43,eigt,y3t
      common/pulsst/ rl1,l,lindex,nsurf
      common/rs/     r(nmax)
      common/splinq/ vq(6)
      dimension y(2),work(15),iwork(5)
      external rkfcow
c
      rl = float(l)
      call splntl( x(1), vq )
      yliq(1,1) = 1.d0
      yliq(2,1) = eigt/(rl*vq(1))
      y(1) = yliq(1,1)
      y(2) = yliq(2,1)
      neqn = 2
      nsurf1=nsurf-1
      do 10 j=1,nsurf1
         k=j+1
         xstart = x(j)
         xfin = x(k)
         relerr = 1.e-8
         abserr = dmin1(abs(y(1)),abs(y(2)))
         abserr = 1.e-8*abserr
         abserr = dmax1(abserr,1.d-20)
         iflag = 1
      call rkf(rkfcow,neqn,y,xstart,xfin,relerr,abserr,iflag,work,iwork)
      if(iflag.eq.3.or.iflag.eq.4.or.iflag.eq.5)write(11,1000)iflag,xfin
         if( iflag .ne.2 ) write(6,1000) iflag,xfin
         if( iflag .eq. -10 ) write(6,1010) j,xfin
         do 3 i=1,2
            yliq(i,k) = y(i)
   3     continue
   10 continue
      rt = x(nsurf)
      call splntl(rt,vq)
c
c      temp = 1.d0/vq(5)-1.d0
c      b2 = yliq(1,nsurf)*((4.d0+r(nsurf)*eigt/g(nsurf))/temp-1.d0)
c     $ +(1.d0-rl1*g(nsurf)/r(nsurf)/eigt/temp)*yliq(2,nsurf)
c
c      temp = 1.d0/y5(nsurf)-1.d0
c      b2 = ((4.d0+eigt/y1(nsurf))/temp-1.d0)*yliq(1,nsurf) +
c     $      (1.d0-rl1*y1(nsurf)/eigt/temp)  *yliq(2,nsurf)
c
      temp = 1.d0/vq(5)-1.d0
      om2nrm = eigt/vq(1)
      b2 = ( (4.d0 + om2nrm)/temp - 1.d0)*yliq(1,nsurf) +
     $     ( 1.d0-rl1/om2nrm/temp )      *yliq(2,nsurf)
      return
 1000 format(1x,17hwatch out, iflag=,i4,6h xfin=,1pe10.2)
 1010 format(1x,'More than 80 steps in RKF, j=',i4,6h xfin=,1pe11.3)
      end
      subroutine rkfcow(rt,y,yp)
      implicit real*8(a-h,o-z)
c
c          returns the right hand side of the first order
c       differential equations in the cowling approximation.
c
      common/misc/   grav,pi,pi4,p43,eigt,y3t
      common/pulsst/ rl1,l,lindex,nsurf
      dimension y(2),yp(2)
      common/splinq/ vq(6)
c
      rl = float(l)
      rlm2 = float(lindex)
      call splntl(rt,vq)
      yp(1) = y(1)*(vq(2)-3.d0+rlm2)+y(2)*(rl1*vq(1)/eigt-vq(2))
      yp(2) = y(1)*(eigt/vq(1)-vq(3))+y(2)*(1.d0-vq(4)+vq(3)+rlm2)
      do 10 i=1,2
         yp(i) = yp(i)*vq(5)
  10  continue
      return
      end
      subroutine grind(b1,b2)
      implicit real*8(a-h,o-z)
      logical iray
c
      parameter ( nmax=600 )
      common g(nmax),rho(nmax),x(nmax),yliq(4,nmax),weight(4,nmax)
      common/misc/   grav,pi,pi4,p43,eigt,y3t
      common/pulsst/ rl1,l,lindex,nsurf
      common/ray/    iray
      common/rs/     r(nmax)
      common/splinq/ vq(6)
      dimension y(8),work(51),iwork(5)
      external rkfliq
c
      rl = float(l)
      call splntl( x(1), vq )
      yliq(1,1) = 1.d0
      yliq(2,1) = eigt/(rl*vq(1))
      y(1) = yliq(1,1)
      y(2) = yliq(2,1)
      yliq(3,1) = y3t
      yliq(4,1) = rl*yliq(3,1)
      y(3) = yliq(3,1)
      y(4) = yliq(4,1)
      neqn = 4
      if( iray ) then
         do 6 j=1,4
            y(j+4) = 0.d0
            weight(j,1) = 0.d0
   6     continue
         neqn = 8
      endif
      nsurf1=nsurf-1
      do 1 j=1,nsurf1
         k=j+1
         xstart=x(j)
         xfin=x(k)
         relerr = 1.e-8
         abserr = abs(y(1))
         do 2 i=2,neqn
            t = abs(y(i))
           abserr = dmin1(abserr,t)
   2    continue
         abserr = 1.e-8*abserr
         abserr = dmax1(abserr,1.d-20)
         iflag=1
      call rkf(rkfliq,neqn,y,xstart,xfin,relerr,abserr,iflag,work,iwork)
         if( iflag.eq.3 .or. iflag.eq.4 .or. iflag.eq.5)
     $               write(11,1000) iflag,xfin
         if( iflag .ne. 2 ) write(6,1000) iflag,xfin
         if( iflag .eq. -10 ) write(6,1010) j,xfin
         do 3 i=1,4
            yliq(i,k)=y(i)
   3     continue
         if( iray ) then
            do 5 i=1,4
               weight(i,k) = y(i+4)
   5        continue
         endif
   1  continue
      rt=x(nsurf)
      call splntl(rt,vq)
c      temp = 1.d0/vq(5) - 1.d0
c      b2=yliq(1,nsurf)*((4.+r(nsurf)*eigt/g(nsurf))/temp-1.)
c     $ +(1.-rl1*g(nsurf)/r(nsurf)/eigt/temp)*yliq(2,nsurf)
c     $ + ((l+1.)/temp-1.)*yliq(3,nsurf)
c
c these are saio and cox bc's, changed
c from osaki and hansen's
c
c      temp = 1.d0/y5(nsurf) - 1.d0
c      b1 = yliq(1,nsurf)*y4(nsurf)+yliq(3,nsurf)*(rl+1.d0)+yliq(4,nsurf)
c      b2 = ((4.d0+eigt/y1(nsurf))/temp-1.d0)*yliq(1,nsurf) +
c     $      (1.d0-rl1*y1(nsurf)/eigt/temp)  *yliq(2,nsurf) +
c     $     ((rl+1.d0)/temp-1.d0)            *yliq(3,nsurf)
c
      temp = 1.d0/vq(5) - 1.d0
      om2nrm = eigt/vq(1)
      b1 = yliq(1,nsurf)*vq(4)+yliq(3,nsurf)*(rl+1.d0)+yliq(4,nsurf)
      b2 = ( (4.d0 + om2nrm)/temp-1.d0 )*yliq(1,nsurf) +
     $     ( 1.d0 - rl1/om2nrm/temp )   *yliq(2,nsurf) +
     $     ( (rl+1.d0)/temp - 1.d0)     *yliq(3,nsurf)
      return
 1010 format(1x,'More than 80 steps in RKF, j=',i4,6h xfin=,1pe11.3)
 1000 format(1x,16hwatch out,iflag=,i4,6h xfin=,1pe10.2)
      end
      subroutine rkfliq(rt,y,yp)
      implicit real*8(a-h,o-z)
      logical iray
      common/misc/   grav,pi,pi4,p43,eigt,y3t
      common/pulsst/ rl1,l,lindex,nsurf
      common/splinq/ vq(6)
      common/ray/    iray
      dimension y(8),yp(8)
      dimension z(4)
c
c  Order x, g/r, V/G1, -rA, U, fac, r
c
c          the vq vector has components:
c                (1) = local g/r
c                (2) = V/G1
c                (3) =-r A
c                (4) = U
c                (5) = fac (local transform to x = ln(r/p) )
c                (6) = r (radius)
c
      call splntl(rt,vq)
c
c          set up equations 17.50 - 17.53 in cox(1980)
c
      rl = float(l)
      rlm2 = float(lindex)
      yp(1) = y(1)*(vq(2)-3.d0+rlm2) + y(2)*(rl1*vq(1)/eigt-vq(2)) +
     $        y(3)*vq(2)
      yp(2) = y(1)*(eigt/vq(1)-vq(3)) + y(2)*(1.-vq(4)+vq(3)+rlm2) -
     $        y(3)*vq(3)
      yp(3) = y(3)*(1.d0-vq(4)+rlm2) + y(4)
      yp(4) = y(1)*vq(4)*vq(3)       + y(2)*vq(4)*vq(2) +
     $        y(3)*(rl1-vq(4)*vq(2)) + y(4)*(rlm2-vq(4))
      do 1 i=1,4
         yp(i)=vq(5)*yp(i)
   1  continue
      if( iray ) then
         rrt = 1.d0/vq(6)**lindex
         do 3 i=1,4
            z(i)=y(i)*rrt
   3     continue
         rhot = vq(4)*vq(1)/pi4/grav
         yp(5) = rhot*vq(6)**5*(z(1)**2 + rl1*z(2)**2*(vq(1)/eigt)**2)
         rrt = rhot*vq(1)*vq(6)**5
         yp(6) = rrt*vq(2)*(z(2)-z(3))**2
         yp(7) = rrt*vq(3)*z(1)**2
         yp(8) =-rrt*(z(4)+(rl+1.d0)*z(3))**2/vq(4)
         do 2 i=5,8
            yp(i) = yp(i)*vq(5)
   2     continue
      endif
      return
      end
      subroutine fingrn(b1,b2)
      implicit real*8(a-h,o-z)
      logical iray
c
      parameter ( nmax=600, nmax10=10*nmax )
      common g(nmax),rho(nmax),x(nmax),yliq(4,nmax),weight(4,nmax)
      common/misc/   grav,pi,pi4,p43,eigt,y3t
      common/pulsst/ rl1,l,lindex,nsurf
      common/ray/    iray
      common/rs/     r(nmax)
      common/modect/ y1m(nmax10),y2m(nmax10),nfine,nodes1,nodes2,modep
      common/splinq/ vq(6)
      parameter ( nzfine=10 )
      dimension y(8),work(51),iwork(5)
c
      external rkfliq
c
      rl = float(l)
      call splntl( x(1), vq )
      yliq(1,1) = 1.d0
      yliq(2,1) = eigt/(rl*vq(1))
      y(1)=yliq(1,1)
      y(2)=yliq(2,1)
c load point 1 of node count arrays
      y1m(1) = y(1)
      y2m(1) = y(2)
      nfine = 10*(nsurf-1)+1
c
      yliq(3,1) = y3t
      yliq(4,1) = rl*yliq(3,1)
      y(3) = yliq(3,1)
      y(4) = yliq(4,1)
      neqn=4
      if( iray ) then
         do 6 j=1,4
            y(j+4) = 0.d0
            weight(j,1) = 0.d0
   6     continue
         neqn=8
      endif
      nsurf1=nsurf-1
      do 1 j=1,nsurf1
         k=j+1
c
c break up interval between x(j) and x(j+1) into 10 parts for
c the integration and save this fine grid for node counting
c
         dxfine = (x(j+1)-x(j))/float(nzfine)
         do 10 jfine=1,nzfine
            xstart=x(j)+float(jfine-1)*dxfine
            xfin=xstart+dxfine
c multiplying original relerr by a factor of 10
            relerr = 1.e-7
c            relerr = relerr*10.d0
            abserr = abs(y(1))
            do 2 i=2,neqn
               t = abs(y(i))
               abserr = dmin1(abserr,t)
   2        continue
            abserr=1.e-8*abserr
            abserr = dmax1(abserr,1.d-20)
            iflag=1
            call rkf(rkfliq,neqn,y,xstart,xfin,relerr,abserr,iflag,
     $                                                  work,iwork)
            if( iflag.eq.3 .or. iflag.eq.4 .or. iflag.eq.5 ) 
     $                       write(11,1000) iflag,xfin
            if( iflag.ne. 2 ) write(6,1000) iflag,xfin
            if( iflag .eq. -10 ) write(6,1010) j,xfin
            y1m(10*(j-1)+jfine+1)=y(1)
            y2m(10*(j-1)+jfine+1)=y(2)
  10     continue
c
         do 3 i=1,4
            yliq(i,k)=y(i)
   3     continue
         if( iray ) then
            do 5 i=1,4
               weight(i,k) = y(i+4)
   5        continue
         endif
   1  continue
      rt=x(nsurf)
      call splntl(rt,vq)
c      temp=1.d0/vq(5) - 1.d0
c
c          These are Saio and Cox bc's, changed from Osaki and Hansen's
c
c      b2=yliq(1,nsurf)*((4.+r(nsurf)*eigt/g(nsurf))/temp-1.)
c     $ +(1.-rl1*g(nsurf)/r(nsurf)/eigt/temp)*yliq(2,nsurf)
c     $ + ((l+1.)/temp-1.)*yliq(3,nsurf)
c
c      temp = 1.d0/y5(nsurf) - 1.d0
c      b1 = yliq(1,nsurf)*y4(nsurf)+yliq(3,nsurf)*(rl+1.d0)+yliq(4,nsurf)
c      b2 = ((4.d0+eigt/y1(nsurf))/temp-1.d0)*yliq(1,nsurf) +
c     $      (1.d0-rl1*y1(nsurf)/eigt/temp)  *yliq(2,nsurf) +
c     $     ((rl+1.d0)/temp-1.d0)            *yliq(3,nsurf)
c
      temp = 1.d0/vq(5) - 1.d0
      om2nrm = eigt/vq(1)
      b1 = yliq(1,nsurf)*vq(4)+yliq(3,nsurf)*(rl+1.)+yliq(4,nsurf)
      b2 = ( (4.d0 + om2nrm)/temp - 1.d0)*yliq(1,nsurf) +
     $     (  1.d0 - rl1/om2nrm/temp )   *yliq(2,nsurf) +
     $     ( (rl+1.d0)/temp - 1.d0)      *yliq(3,nsurf)
      return
 1000 format(1x,16hwatch out,iflag=,i4,6h xfin=,1pe10.2)
 1010 format(1x,'More than 80 steps in RKF, j=',i4,6h xfin=,1pe11.3)
      end
      subroutine eigenf(stwait)
      implicit real*8(a-h,o-z)
c
c          find the normalized eigenvectors and the weight functions
c       as a function of position.
c
      parameter ( nmax=600 )
      common g(nmax),rho(nmax),x(nmax),yliq(4,nmax),weight(4,nmax)
      common/misc/   grav,pi,pi4,p43,eigt,y3t
      common/pulsst/ rl1,l,lindex,nsurf
      common/rs/     r(nmax)
      common/splinq/ vq(6)
      dimension stwait(nsurf)
c
      rl = float(l)
      do 10 i=1,nsurf
         t1 = (r(nsurf)/r(i))**lindex/yliq(1,nsurf)
         do 15 j=1,4
            yliq(j,i) = yliq(j,i)*t1
  15     continue
  10  continue
c
c          weight functions
c
      omsq = eigt
      rke = 0.d0
      twait = 0.d0
      do 20 i=1,nsurf-1
         rt = x(i)
         call splntl(rt,vq)
         if( i .eq. 1 ) then
            dx = x(2)-x(1)
         elseif( i .eq. nsurf ) then
            dx = x(nsurf) - x(nsurf-1)
         else
            dx = (x(i+1)-x(i-1))/2.d0
         endif
         dr = dx*vq(5)
         temp = rho(i)*vq(6)**5
         weight(1,i)=temp*(yliq(1,i)**2 +
     $                            rl1*(vq(1)/omsq)**2*yliq(2,i)**2)
         rke = rke + weight(1,i)*dr
         temp = temp*vq(1)
         weight(2,i) = temp*vq(2)*(yliq(2,i)-yliq(3,i))**2
         weight(3,i) = temp*vq(3)*yliq(1,i)**2
         weight(4,i) =-temp*(yliq(4,i)+(rl+1.d0)*yliq(3,i))**2/vq(4)
         stwait(i) = weight(2,i) + weight(3,i) + weight(4,i)
         twait = twait + stwait(i)*dr
  20  continue
      twait = twait/rke
      qch = abs((omsq-twait)/omsq)
      qrke = rke*omsq/2.d0
      write(11,2000) omsq,twait,qch,qrke
      write(6,2000) omsq,twait,qch,qrke
      write(1,2000) omsq,twait,qch,qrke
      do 30 i=1,nsurf
         weight(1,i) = weight(1,i)/(rke*r(nsurf))
         weight(2,i) = weight(2,i)/(omsq*rke)
         weight(3,i) = weight(3,i)/(omsq*rke)
         weight(4,i) = weight(4,i)/(omsq*rke)
         stwait(i) = stwait(i)/(omsq*rke)
  30  continue
      return
c
 2000 format(2x,32heigenvalue from matrix solution=,1pe12.3,
     $     /,2x,'eigenvalue from direct weight function=',e12.3,
     $       7h error=,e12.3,/,
     $       3x,21hkinetic energy amp. =,e11.4)
      end
      subroutine splntl(xint,yout)
      implicit real*8(a-h,o-z)
c
c     interpolates in the initial model mesh.
c
c      common/savngl/ k,kq
      parameter ( nmax=600 )
      common/savel/  xi(nmax),c1(4,nmax),y1(nmax),c2(4,nmax),
     $ y2(nmax),c3(4,nmax),y3(nmax),c4(4,nmax),y4(nmax),c5(4,nmax),
     $ y5(nmax),c6(4,nmax),m
      common/rs/     y6(nmax)
      dimension yout(6)
c
      mm = m-1
      if( xint.ge. xi(1) .and. xint.lt.xi(m) ) goto 40
      if( xint .lt. (1.d0+1.e-8)*xi(1) ) goto 110
      if( xint .ge. xi(m) ) goto 20
      k = 1
      goto 70
  20  if( xint .gt. (1.d0+1.e-6)*xi(m) ) goto 111
      k = mm
      goto 70
c
c          given x such that xi(1) < x < xi(m), find k such that
c       xi(k) < xint < xi(k+1). uses a binary search for ordered
c       data.
c
  40  continue
      il = 1
      ir = m
  50  k = il+((ir-il)/2)
      if( xint .ge. xi(k) ) goto 60
      ir = k
      goto 50
  60  if( xint .lt. xi(k+1) ) goto 70
      il = k
      goto 50
c
c          use the spline fits.
c
  70  continue
      x1 = xi(k+1)-xint
      xx = xint-xi(k)
      x12 = x1*x1
      xx2 = xx*xx
      yout(1) = x1*(c1(1,k)*x12 +c1(3,k))+xx*(c1(2,k)*xx2 +c1(4,k))
      yout(2) = x1*(c2(1,k)*x12 +c2(3,k))+xx*(c2(2,k)*xx2 +c2(4,k))
      yout(3) = x1*(c3(1,k)*x12 +c3(3,k))+xx*(c3(2,k)*xx2 +c3(4,k))
      yout(4) = x1*(c4(1,k)*x12 +c4(3,k))+xx*(c4(2,k)*xx2 +c4(4,k))
      yout(5) = x1*(c5(1,k)*x12 +c5(3,k))+xx*(c5(2,k)*xx2 +c5(4,k))
      yout(6) = x1*(c6(1,k)*x12 +c6(3,k))+xx*(c6(2,k)*xx2 +c6(4,k))
      return
c
  110 continue
      write(6,1100) xint,xi(1)
      yout(1) = y1(1)
      yout(2) = y2(1)
      yout(3) = y3(1)
      yout(4) = y4(1)
      yout(5) = y5(1)
      yout(6) = y6(1)
      return
  111 continue
      write(6,1101) xint,xi(m)
      yout(1) = y1(m)
      yout(2) = y2(m)
      yout(3) = y3(m)
      yout(4) = y4(m)
      yout(5) = y5(m)
      yout(6) = y6(m)
      return
c
 1100 format(1x,25hrange error in spline, x=,1pe16.8,/
     $       1x,29hpoint less than minimum x of ,e15.8)
 1101 format(1x,25hrange error in spline, x=,1pe16.8,/
     $       1x,32hpoint greater than maximum x of ,e15.8)
      end
      subroutine consl(x,y,m,c)
      implicit real*8(a-h,o-z)
c
c          set up coefficients for a cubic spline in the initial
c       model quantities.
c
      parameter ( nmax=600 )
      dimension c(4,m),x(m),y(m)
      dimension a(nmax,3),d(nmax),b(nmax),z(nmax),p(nmax)
c
      mm = m-1
      do 10 k=1,mm
         d(k) = x(k+1)-x(k)
         p(k) = d(k)/6.d0
         z(k) = (y(k+1)-y(k))/d(k)
  10  continue
      do 20 k=2,mm
         b(k) = z(k)-z(k-1)
  20  continue
      a(1,2) =-1.d0 - d(1)/d(2)
      a(1,3) = d(1)/d(2)
      a(2,3) = p(2)-p(1)*a(1,3)
      a(2,2) = 2.d0*(p(1)+p(2))-p(1)*a(1,2)
      a(2,3) = a(2,3)/a(2,2)
      b(2) = b(2)/a(2,2)
      do 30 k=3,mm
         a(k,2) = 2.d0*(p(k-1)+p(k))-p(k-1)*a(k-1,3)
         b(k) = b(k)-p(k-1)*b(k-1)
         a(k,3) = p(k)/a(k,2)
         b(k) = b(k)/a(k,2)
   30 continue
      q = d(m-2)/d(m-1)
      a(m,1) = 1.d0+q+a(m-2,3)
      a(m,2) =-q-a(m,1)*a(m-1,3)
      b(m) = b(m-2)-a(m,1)*b(m-1)
      z(m) = b(m)/a(m,2)
      mn = m-2
      do 40 i=1,mn
         k = m-i
         z(k) = b(k)-a(k,3)*z(k+1)
   40 continue
      z(1) =-a(1,2)*z(2)-a(1,3)*z(3)
      do 50 k=1,mm
         q = 1.d0/(6.d0*d(k))
         c(1,k) = z(k)*q
         c(2,k) = z(k+1)*q
         c(3,k) = y(k)/d(k)-z(k)*p(k)
         c(4,k) = y(k+1)/d(k)-z(k+1)*p(k)
  50  continue
      return
      end
      subroutine rkf(f,neqn,y,t,tout,relerr,abserr,iflag,work,
     $ iwork)
      implicit real*8(a-h,o-z)
c
      dimension y(neqn),work(1),iwork(5)
      external f
      k1m=neqn+1
      k1=k1m+1
      k2=k1+neqn
      k3=k2+neqn
      k4=k3+neqn
      k5=k4+neqn
      k6=k5+neqn
      call rkfs(f,neqn,y,t,tout,relerr,abserr,iflag,work(1),
     $ work(k1m),work(k1),work(k2),work(k3),work(k4),work(k5),
     $ work(k6),work(k6+1),iwork(1),iwork(2),iwork(3),iwork(4),
     $ iwork(5))
      return
      end
      subroutine rkfs(f,neqn,y,t,tout,relerr,abserr,iflag,yp,h,f1,f2,f3,
     $  f4,f5,savre,savae,nfe,kop,init,jflag,kflag)
      implicit real*8(a-h,o-z)
      logical hfaild,output
      parameter ( u26=2.e-13 ,remin=1.e-12 )
      parameter ( maxnfe=3000 )
c
      dimension y(neqn),yp(neqn),f1(neqn),f2(neqn),f3(neqn),f4(neqn),
     $         f5(neqn)
      external f
      if( neqn .lt. 1 ) then
         iflag = 7
         return
      endif
      if( (relerr .lt. 0.d0)  .or.  (abserr .lt. 0.d0) ) then
         iflag = 7
         return
      endif
      mflag = iabs(iflag)
      if( (mflag .lt. 1) .or. (mflag .gt. 7) ) then
         iflag = 7
         return
      endif
      if( t .eq. tout  ) then
         iflag = 7
         return
      endif
      if( mflag .eq. 1 ) goto 50
      if( mflag .ne. 2 ) then
         if( iflag .eq. 3 ) goto 40
         if( (iflag .eq. 4) .and. (abserr .gt. 0.d0) ) then
            iflag = jflag
            goto 50
         endif
c
         write(6,1000) iflag,t
         stop
      endif
      if( init .eq. 0 )  goto 45
      if( kflag .eq. 3 ) goto 40
      if( (kflag.eq.4) .and.  (abserr.eq.0.d0) ) then
         write(6,1000) iflag,t
         stop
      endif
      if( (kflag .eq.5)  .and. (relerr.le.savre) .and.
     $ (abserr.le.savae) ) then
         write(6,1000) iflag,t
         stop
      endif
      goto 50
   40 nfe=0
      if( mflag .eq. 2 ) goto 50
   45 iflag=jflag
   50 jflag=iflag
      kflag=0
      savre=relerr
      savae=abserr
      rer = dmax1(relerr,remin)
      dt=tout-t
      if( mflag .eq. 1) then
c
c        Entry for normal point (iflag = 1 in argument list.)
c
         init=0
         kop=0
         a=t
         call f(a,y,yp)
         nfe=1
         if( t .eq. tout) then
            iflag=2
            return
         endif
      elseif( init .ne. 0 ) then
         goto 80
      endif
   65 continue
      init=1
      ymax=0.
      ypn=0.
      do 70 k=1,neqn
         ypn = dmax1(abs(yp(k)),ypn)
         ymax = dmax1(abs(y(k)),ymax)
  70  continue
      etn=rer*ymax+abserr
      h = abs(dt)
      if( etn .lt. ypn*h**5 ) then
         h = dmax1((etn/ypn)**0.2,u26*dmax1(abs(t),h))
      endif
  80  h = sign(h,dt)
      if( abs(h) .ge. abs(dt) ) kop=kop+1
      if( kop .eq. 100 ) then
         iflag=6
         return
      endif
      if( abs(dt) .le. u26*abs(t)) then
c
c          Linearly interpolate the new y's and return to program.
c
         do 90 k=1,neqn
            y(k)=y(k)+dt*yp(k)
  90     continue
         a=tout
         call f(a,y,yp)
         nfe=nfe+1
         t=tout
         iflag=2
         return
      endif
c
c          Ready to integrate the equations.
c
      output=.false.
      scale=2./rer
      ae=scale*abserr
 100  continue
      hfaild=.false.
      hmin = u26*abs(t)
      dt=tout-t
c      if( abs(dt) .ge. 2.d0*abs(h)) goto 200
c      if( abs(dt) .gt. abs(h)/0.9e0 ) goto 150
c      output=.true.
c      h=dt
c      goto 200
c 150  h=0.5e0*dt
c
c          Limits on h are rather strange but can be seen as
c       three regions:
c          abs(h) < abs(dt)/2.                    h is okay
c          abs(dt)/2. < abs(h) < 0.9 abs(dt)      h = 0.5 dt
c          abs(h) > 0.9 abs(dt)                   h = dt
c
      if( abs(h) .gt. abs(dt)/2.d0) then
         if( abs(h) .gt. abs(dt)*0.9e0 ) then
            output=.true.
            h = dt
         else
            h = 0.5e0*dt
         endif
      endif
 200  continue
      if( nfe .gt. maxnfe ) then
         iflag=3
         kflag=3
         return
      endif
c
c          Call the integrator subroutine. Note that the high-order
c       estimate for y is returned in f1.
c
      call fehl(f,neqn,y,t,h,yp,f1,f2,f3,f4,f5,f1)
      nfe=nfe+5
      eeoet = 0.d0
      do 250 k=1,neqn
         et = abs(y(k)) + abs(f1(k)) + ae
         if( et .le. 0.d0 ) then
            iflag=4
            kflag=4
            return
         endif
         ee = abs((-2090.d0*yp(k)+(21970.d0*f3(k)-15048.d0*f4(k)))+
     $                       (22528.d0*f2(k)-27360.d0*f5(k)))
         eeoet = dmax1(eeoet,ee/et)
 250  continue
      esttol = abs(h)*eeoet*scale/752400.d0
      if( esttol .gt. 1. ) then
c
c          Current guess at h has resulted in an error that is too
c       large. Reduce h and try agian.
c
         hfaild=.true.
         output=.false.
         s=0.1e0
         if( esttol .lt. 59049. ) s=0.9/esttol**0.2
         h=s*h
         if( abs(h) .gt. hmin ) then
            goto 200
         else
            iflag=5
            kflag=5
            return
         endif
      endif
c
c          Accept this point and continue trying to complete
c       the integration up to tout.
c
      t = t+h
      do 270 k=1,neqn
         y(k)=f1(k)
 270  continue
      a=t
      call f(a,y,yp)
      nfe=nfe+1
      if( .not. hfaild) then
         s = 5.d0
         if( esttol .gt. 1.889568e-4 ) s=0.9/esttol**0.2
         h = sign(dmax1(s*abs(h),hmin),h)
      endif
      if( output ) then
c
c          Integration has been completed, set iflag and return
c
         t=tout
         iflag=2
c         if( nfe .gt. 80 ) iflag =-10
         return
      endif
      if( iflag .gt. 0 ) goto 100
      iflag=-2
      return
 1000 format(1x,15hrkf says iflag=,i5,3h x=,e12.4)
 2700 format(1x,1pe11.3,i4)
      end
      subroutine fehl(f,neqn,y,t,h,yp,f1,f2,f3,f4,f5,s)
      implicit real*8(a-h,o-z)
      dimension y(neqn),yp(neqn),f1(neqn),f2(neqn),f3(neqn),f4(neqn),
     $          f5(neqn),s(neqn)
      ch = 0.25e0*h
      do 221 k=1,neqn
         f5(k)=y(k)+ch*yp(k)
 221  continue
      call f(t+0.25e0*h,f5,f1)
      ch = 0.09375e0*h
      do 222 k=1,neqn
         f5(k)=y(k)+ch*(yp(k)+3.d0*f1(k))
 222  continue
      call f(t+0.375e0*h,f5,f2)
      ch = h/2197.d0
      do 223 k=1,neqn
         f5(k)=y(k)+ch*(1932.*yp(k)+(7296.*f2(k)-7200.*f1(k)))
 223  continue
      call f(t+12.d0/13.d0*h,f5,f3)
      ch = h/4104.d0
      do 224 k=1,neqn
         f5(k)=y(k)+ch*((8341.d0*yp(k)-845.d0*f3(k))+
     $                   (29440.d0*f2(k)-32832.d0*f1(k)))
 224  continue
      call f(t+h,f5,f4)
      ch = h/20520.d0
      do 225 k=1,neqn
         f1(k)=y(k)+ch*((-6080.d0*yp(k)+(9295.d0*f3(k)-5643.d0*f4(k)))+
     $                             (41040.d0*f1(k)-28352.d0*f2(k)))
 225  continue
      call f(t+0.5e0*h,f1,f5)
      ch = h/7618050.d0
c
c          s is returned with the higher order estimate of y(new)
c
      do 230 k=1,neqn
         s(k) = y(k) + ch*( (902880.*yp(k) +
     $                      (3855735.*f3(k)-1371249.*f4(k)) ) +
     $                      (3953664.*f2(k)+277020.*f5(k))  )
 230  continue
c
c          This s is returned with the lower order estimate of y(new)
c
c      a1 = 25./216.
c      b1 = 0.
c      b2 = 1408.d0/2565.d0
c      b3 = 2197.d0/4104.d0
c      b4 =-0.20e0
c      b5 = 0.d0
c      do 230 k=1,neqn
c         s(k) = y(k)+h*( a1*yp(k) + (b3*f3(k) + b4*f4(k)) + b2*f2(k) )
c 230  continue
      return
      end
      subroutine modeid
      implicit real*8(a-h,o-z)
c
      parameter ( nmax=600, nmax10=10*nmax )
      common/modect/ y1m(nmax10),y2m(nmax10),nfine,nodes1,nodes2,modep
      integer quad2,quad1
      modep = 0
c
c          first count the nodes in y1 and y2
c
      nodes1 = 0
      nodes2 = 0
      do 20 i=2,nfine
         fit1 = y1m(i)*y1m(i-1)
         fit2 = y2m(i)*y2m(i-1)
         if( fit1 .le. 0 ) nodes1 = nodes1+1
         if( fit2 .le. 0 ) nodes2 = nodes2+1
  20  continue
c
c          count crossings in the phase diagram
c
      quad1 = iquad(y1m(1),y2m(1))
      do 30 i=2,nfine
         quad2 = iquad(y1m(i),y2m(i))
         idq = quad2-quad1
         if( idq .eq. 0 ) goto 30
c
c          see if the quadrant change is a crossing in y1
c
         if( iabs(idq) .eq. 3 ) goto 100
         if( quad1.eq.3 .and. quad2.eq.2 ) goto 100
         if( quad1.eq.2 .and. quad2.eq.3 ) goto 100
c   not a y1 crossing
         quad1 = quad2
        goto 30
c
c          if crossing is clockwise, increment mode count.
c          if crossing is negative, decrement mode count.
c
 100     continue
          if( idq.eq.-3 .or. idq.eq.1 ) kdmod =-1
          if( idq.eq.3 .or. idq.eq.-1 ) kdmod = 1
          modep = modep+kdmod
          quad1 = quad2
  30  continue
      return
      end
      function iquad(y1,y2)
      implicit real*8(a-h,o-z)
c
c      if( (y1.gt.0.d0) .and. (y2.ge.0.d0) ) iquad = 1
c      if( (y1.le.0.d0) .and. (y2.gt.0.d0) ) iquad = 2
c      if( (y1.lt.0.d0) .and. (y2.le.0.d0) ) iquad = 3
c      if( (y1.ge.0.d0) .and. (y2.lt.0.d0) ) iquad = 4
      if( y1 .ge. 0. ) then
         if( y2 .ge. 0. ) then
            iquad = 1
         else
            iquad = 4
         endif
      else
         if( y2 .ge. 0. ) then
            iquad = 2
         else
            iquad = 3
         endif
      endif
      return
      end
      subroutine pltint(qmass,blum,nmod,n,r,ra)
      implicit real*8(a-h,o-z)
c
c          initialize the plot file and write the title line.
c
      parameter ( nmax=600 )
      common g(nmax),rho(nmax),x(nmax),yliq(4,nmax),weight(4,nmax)
      dimension r(nmax),ra(nmax),onemr(nmax)
c
c     open (unit=21,file='cjhplt',status='unknown')
c
      write(21,1000) qmass,blum,nmod
      call pltdmp(r,nmax,n,'r   ')
      call pltdmp(x,nmax,n,'x   ')
      call pltdmp(g,nmax,n,'g   ')
      do 10 i=1,n-1
         onemr(i) =-dlog10(1.d0 - r(i)/r(n))
  10  continue
      onemr(n) = onemr(n-1) + 0.2e0
      call pltdmp(onemr,nmax,n,'1-r ')
      call pltdmp(ra   ,nmax,n,'ra  ')
      do 30 i=1,n
         onemr(i) = dlog10(rho(i))
  30  continue
      call pltdmp(onemr,nmax,n,'rho ')
      return
c 1000 format(1x,f6.2,19hm^lh.5m11\h^lxhxm0\,8h l/lsun=,
c     $       1pe10.3,9h # zones=,i4)
 1000 format(1x,f6.2,'M/Msol, l/lsun= ',f8.4,9h # zones=,i4)
      end
      subroutine pltdmp(vec,nmax,n,ititl)
      implicit real*8(a-h,o-z)
      character*4 ititl
c
c          write the vector vec to the plot file (tape12) with appended
c       title ititl.
c
      dimension vec(nmax)
      write(21,1000) n,ititl,(vec(i),i=1,n)
      return
 1000 format(i4,10x,a4,/,(1p,6e12.4) )
      end
