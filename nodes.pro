function nodes, w, lout
;c
;c finds no. of nodes(sign changes) in nz elements of array w
;c
;      dimension w(nz),inode(300)
;      common/const/  zero,one,two,thre,for,ten,ahf,qrt
;c
inode = replicate( -1, 300 )
nodcnt = 0
nz = n_elements(w)
for i=3,nz-4 do begin
   if( w[i-1]*w[i] lt 0. ) then begin
;c zero crossing found. is it general enough to be a node#
      if( (w[i-2]*w[i-1] gt 0.) and (w[i]*w[i+1] gt 0.) ) then begin
         inode[nodcnt] = i-1
         nodcnt = nodcnt + 1
      endif
   endif
endfor
;inode(nodcnt+1) = nz
nodes = nodcnt
;print, nodes, inode[0:nodcnt-1], format='(" nodes=",(30i5))'
if( lout gt 0 ) then begin
   print, nodes, format='(1x,i2," nodes found")'
   if( nodcnt gt 0 ) then print, inode[0:nodcnt-1], format='(" nodes=",(30i5))'
endif
return, nodes
end