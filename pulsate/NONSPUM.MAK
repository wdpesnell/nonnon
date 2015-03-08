#
#       MAKE file for NONSPU
#
nonspum.obj: nonspum.for
      for/g_float nonspum

cmplna.obj: cmplna.for
      for/g_f cmplna

lnanon.obj: lnanon.for
      for/g_f lnanon

matset.obj: matset.for
      for/g_f  matset

matrixes.obj: matrixes.for
      for/g_f matrixes

nonmisc.obj: nonmisc.for
      for/g_f nonmisc

summary.obj: summary.for
      for/g_f summary

nonspum.exe: nonspum.obj cmplna.obj lnanon.obj matset.obj -
            matrixes.obj nonmisc.obj summary.obj
      link nonspum + cmplna + lnanon + matset + matrixes + nonmisc -
           + summary
