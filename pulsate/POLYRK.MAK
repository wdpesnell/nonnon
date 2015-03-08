#
#          MAKE file for polyrk
#
polyrk.obj: polyrk.for
       for/g_f polyrk.for

cjhdmp.obj: cjhdmp.for
      for/g_f cjhdmp.for

pekeris.obj: pekeris.for
      for/g_f pekeris

polyrk.exe: polyrk.mak polyrk.obj cjhdmp.obj pekeris.obj
      link polyrk + cjhdmp + pekeris
