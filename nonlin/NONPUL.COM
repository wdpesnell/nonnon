$!
$!          Run the nonlinear pulsation code
$!
$ set default [.private.tdconv]
$ assign [-.model]evnondmp.dat lnrad
$ run nonpul
1000 10 5.e-2 .9
