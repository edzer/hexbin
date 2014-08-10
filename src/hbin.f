	subroutine hbin(x,y,cell,cnt,xcm,ycm, size, shape,
     *                  rx,ry, bnd, n, cellid)

C	Copyright 1991
C	Version Date:	September 16, 1994
C	Programmer:	Dan Carr
C	Indexing:	Left to right, bottom to top
C			bnd(1) rows, bnd(2) columns
C	Output:		cell ids for non empty cells, revised bnd(1)

c			optionally also return cellid(1:n)
c       Copyright (2004) Nicholas Lewin-Koh and Martin Maechler

	implicit none

	integer n, nc, cell(*), cnt(*), bnd(2), cellid(*)
c       cellid(*): length 1 or n
	double precision x(n), y(n), xcm(*),ycm(*), rx(2),ry(2), size
        double precision shape
	integer i, i1, i2, iinc
	integer j1, j2, jinc
	integer L, lmax, lat
	double precision c1, c2, con1, con2, dist1
	double precision sx, sy, xmin, ymin, xr, yr
	logical keepID

	keepID = (cellid(1) .eq. 0)
C_______Constants for scaling the data_____________________________

	xmin = rx(1)
	ymin = ry(1)
	xr = rx(2)-xmin
	yr = ry(2)-ymin
	c1 = size/xr
	c2 = size*shape/(yr*sqrt(3.))

	jinc= bnd(2)
	lat=jinc+1
	iinc= 2*jinc
	lmax=bnd(1)*bnd(2)
	con1=.25
	con2=1.0/3.0

C_______Binning loop________________________________________

	do i=1,n
	  sx = c1 * (x(i) - xmin)
	  sy = c2 * (y(i) - ymin)
	  j1 = sx+.5
	  i1 = sy+.5
	  dist1=(sx-j1)**2 + 3.*(sy-i1)**2

	  if(dist1 .lt. con1) then
	    L=i1*iinc + j1+1
	  elseif(dist1 .gt. con2) then
	    L=int(sy)*iinc + int(sx)+lat
	  else
	    j2 = sx
	    i2 = sy
	    if(dist1 .le. (sx-j2 -.5)**2 + 3.*(sy-i2 -.5)**2) then
	      L=i1*iinc+ j1+1
	    else
	      L=i2*iinc+ j2+lat
	    endif
	  endif

	  cnt(L)=cnt(L)+1
	  if (keepID) cellid(i)=L
	  xcm(L)=xcm(L)+ (x(i)-xcm(L))/cnt(L)
	  ycm(L)=ycm(L)+ (y(i)-ycm(L))/cnt(L)
	enddo

C_______Compression of output________________________________________

	nc=0
	do L=1,lmax
	   if(cnt(L) .gt. 0) then
	      nc=nc+1
	      cell(nc)=L
	      cnt(nc)=cnt(L)
	      xcm(nc)=xcm(L)
	      ycm(nc)=ycm(L)
	   endif
	enddo
	n=nc
	bnd(1)=(cell(nc)-1)/bnd(2)+1
	return
	end
