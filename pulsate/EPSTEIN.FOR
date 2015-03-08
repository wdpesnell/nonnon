c
c	    Given the radial eigenvectors, calculate the Epstein
c	weight functions. See Epstein Ap. J. 112(6)1950.
c
      epstwt(1) = (thre*g1(1)-for)*g*rm(2)*dm2(2)/r(2)*dr(1)**2
      chwait = epstwt(1)
      do 104 i=2,n
        epstwt(i)=(thre*g1(i)-for)*g*rm(i+1)*dm2(i+1)/r(i+1)*dr(i+1)**2+
     $     p(i)*v(i)*g1(i)*((dr(i+1)-dr(i))/dlog(r(i+1)/r(i)))**2*dm1(i)
         chwait = chwait + epstwt(i)
 104  continue
      chwait = chwait/rke
      ech = dabs((omsq-chwait)/omsq)
      write(6,1001) omsq,qomsq,iv,ih,stwait(np1),qch,chwait,ech
      write(1,1001) omsq,qomsq,iv,ih,stwait(np1),qch,chwait,ech
      write(11,1001) omsq,qomsq,iv,ih,stwait(np1),qch,chwait,ech
      iv = nodes(dr,np1,1)
c
c          normalize the weight function per zone to omsq.
c
      do 105 i=1,n
         weight(i) = weight(i)/(omsq*rke)
         stwait(i) = stwait(i)/(omsq*rke)
         epstwt(i) = epstwt(i)/(omsq*rke)
 105  continue
      stwait(np1) = stwait(np1)/omsq
c
      call pltdmp(epstwt,nmax,n,'epwt')
 1001 format(/,19h linear radial mode,/,
     $ 8h  omsq =,1pe13.6,2X,7hqomsq =,e11.3,4H iv=,i3,4h ih=,i3,/,
     $ 2x,33heigenvalue from weight function =,e12.5,8H error =,e10.3,/,
     $ 1X,42h eigenvalue from epstein weight function =,e12.5,
     $ 7H error=,e10.3)
