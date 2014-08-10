	subroutine hcell(x,y,cell,n,size,shape,rx,ry,bnd)
C	Copyright 1991
C	Version Date: 	September 16, 1994
C	Programmer: 	Dan Carr
C	Indexing:	Left to right, bottom to top
C			bnd(1) rows, bnd(2) columns
C	Output:		cell ids for none empty cells, revised bnd(1)

c	implicit none
	integer n, cell(1), bnd(2)
	double precision x(1), y(1), rx(2), ry(2), size, shape
	integer i, i1, i2, iinc
	integer j1, j2, jinc
	integer L, lat, celmax
	double precision c1, c2, con1, con2, dist1
	double precision sx, sy, xmin, ymin, xr, yr

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
	con1=.25
	con2=1./3.
	celmax=0

C_______Binning loop________________________________________

	do i=1,n
	  sx = c1 * (x(i) - xmin)
          sy = c2 * (y(i) - ymin)
	  j1 = sx+.5
	  i1 = sy+.5
	  dist1=(sx-j1)**2 + 3.*(sy-i1)**2

	  if(dist1.lt.con1)then
	    L=i1*iinc+j1+1
	  elseif(dist1.gt.con2)then
	    L=int(sy)*iinc + int(sx)+lat
	  else
	    j2 = sx
	    i2 = sy
	    if( dist1.le.(sx-j2-.5)**2 + 3. * (sy - i2 -.5)**2) then
	      L=i1*iinc+j1+1
            else
	      L=i2*iinc+j2+lat
	    endif
	  endif

	  cell(i)=L
	  celmax = max(celmax,L)
	enddo
	bnd(1)=(celmax-1)/bnd(2)+1
	return
	end
