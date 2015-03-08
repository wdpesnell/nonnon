#
#          MAKE file for NONAD
#
nonad.obj: nonad.for
       for/g_f nonad.for

lnradp.obj: lnradp.for
      for/g_f lnradp.for

admatset.obj: admatset.for
      for/g_f admatset.for

matrixes.obj: matrixes.for
      for/g_f matrixes.for

adnonmsc.obj: adnonmsc.for
      for/g_f adnonmsc.for

invint.obj: invint.for
      for/g_f invint.for

nonad.exe: nonad.mak nonad.obj lnradp.obj admatset.obj invint.obj -
    matrixes.obj adnonmsc.obj
      link/map/cross nonad + lnradp + admatset + matrixes + adnonmsc + invint
