openr, lunit, 'L0102800203.pul', /get_lun
title = ' '
readf, lunit, title
readf, lunit, npts,rlumgv,totmas,teff,rphoto,corlum
readf, lunit, irad,noburn,onemq0,onemq1,rlums,rlumc
      np1 = npts+1
rp = dblarr(np1)
readf, lunit, format='(4e20.13)', rp
tp = dblarr(npts)
readf, lunit, format='(4e20.13)', tp
vp = dblarr(npts)
readf, lunit, format='(4e20.13)', vp
cv = dblarr(npts)
readf, lunit, format='(4e20.13)', cv
dkdr = dblarr(npts)
readf, lunit, format='(4e20.13)', dkdr
dkdt = dblarr(npts)
readf, lunit, format='(4e20.13)', dkdt
dm1 = dblarr(npts)
readf, lunit, format='(4e20.13)', dm1
akap = dblarr(npts)
readf, lunit, format='(4e20.13)', akap
dm2 = dblarr(np1)
readf, lunit, format='(4e20.13)', dm2
rm = dblarr(np1)
readf, lunit, format='(4e20.13)', rm
p = dblarr(npts)
readf, lunit, format='(4e20.13)', p
g1 = dblarr(npts)
readf, lunit, format='(4e20.13)', g1
g3m1 = dblarr(npts)
readf, lunit, format='(4e20.13)', g3m1
frft = dblarr(np1)
readf, lunit, format='(4e20.13)', frft
sorce = dblarr(npts)
readf, lunit, format='(4e20.13)', sorce
dtsorc = dblarr(npts)
readf, lunit, format='(4e20.13)', dtsorc
dvsorc = dblarr(npts)
readf, lunit, format='(4e20.13)', dvsorc
bv = dblarr(np1)
readf, lunit, format='(4e20.13)', bv
readf, lunit, format='(4e20.13)', pc,rhoc,tc,cormas,rl0,chit0,chr0,q0, $
        g10,g3m10,cv0,cp0,opac0,dkdt0,dkdr0,sorc0,dedt0,dedv0
free_lun, lunit
totmas = max(rm)
dkdr = dkdr[*] + g3m1*dkdt
iout = 2
lin = 10
result = matset( npts, iout, lin, rp,tp,vp,cv,dkdr, $
        dkdt,dm1,gkap,dm2,rm, bv, p,g1,g3m1, rzone )
end