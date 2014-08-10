C	File:	     	herode.f
C	Version date:	Jan 4, 1994
C	Programmer:	Daniel B. Carr
C
C	The vector erode returns the gray-level erosion order for hexagon cells.
C		The erosion cycle is:
C			 cycle = (erode-1)/6 + 1
C       	Many cells may be eroded in the same cycle
C		A tie break is the cell count deficit at erosion time:
C			deficit=erode - 6*cycle
C	 	The last eroded cell might be considered a bivariate median
C	
C	The algorithm:
C	Repeat until no cells are left in the list.
C	   Process list
C	      Reduce the cell counts by the a multiple of exposed sides
C	      If a cell count is zero or less after an erosion cycle
C		 let  order=order + 6
C		 report erode = order + cell count (count is <= 0)
C		 remove the cell from consideration
C		 update exposed side counts for existing neighbor cells
C		    if exposed sides was zero, temporarily store id's
C	      else
C		 compress list 
C	      endif
C	   Add temporarily stored id's to list					
C	End Repeat

	subroutine herode(cell,cnt,n,bdim,
     *                    erode,ncnt,ncell,sides,neib,exist)

C
C	
	implicit none

	integer cell(*), cnt(*) ! cell id and count
	integer n, bdim(2)	! number of cells and 2-D array bounds
	integer erode(*)	! erosion status
	integer	ncell(*),ncnt(*)  ! extracted id's and expanded counts
	integer sides(*)	! number of exposed sides
	integer neib(6,*)	! pointers to the neighbors
	logical exist(0:*)	! cell existence

	integer nrow, ncol, Lmax  ! dimensions
	integer	inc1(6), inc2(6) ! increments to get neighbors
	integer i, icell, j, k, L ! subscripts
	integer nc, nnc, nb, ninc, r, c  !more subscripts
	integer loop, order, maxcnt


C_______Zero cell ordering numbers________________________________

	order=0

C_______Load the increment arrays and constants

	nrow = bdim(1)
	ncol = bdim(2)
	Lmax = nrow * ncol
	nnc = n

C______Load increment arrays to neigbors______________
C
C order=right, up left, down left, up right, left, down right

	inc1(1)= 1
	inc1(2)= ncol-1
	inc1(3)= -ncol-1
	inc1(4)= ncol
	inc1(5)=-1
	inc1(6)= -ncol

	inc2(1)= 1
	inc2(2)= ncol
	inc2(3)= -ncol
	inc2(4)= ncol+1
	inc2(5)=-1
	inc2(6)= -ncol+1


c_______load working arrays_______________________________________________

        do i=0,Lmax
          exist(i)=.false.
        enddo

	maxcnt=0
	do i=1,n
	  icell=cell(i)
	  ncnt(icell)=cnt(i)
	  exist(icell)=.true.
	  maxcnt=max(maxcnt,cnt(i))
	enddo

C_______Store pointers to cell neighbors_________________________
C
C	A pointer of 0 means the neigbor in out of bounds 
C	Also find the max count
C		Speed:  Can avoid adding 1's to r and c
C	 		but this code is easier to follow
        do i=1,n
           L=cell(i)
	   k = L -1
	   r=k/ncol+1
	   c=mod(k,ncol)+1
	   if(mod(r,2).eq.1)then 
             do j = 1,6
	        neib(j,L) = L + inc1(j)
	     enddo

             if (c .eq. 1) then	 	
	       neib(2,L) = 0
               neib(3,L) = 0
               neib(5,L) = 0
             else if (c .eq. ncol) then
	       neib(1,L) = 0
	     endif

	     if (r .eq. 1) then 
	       neib(3,L) = 0
	       neib(6,L) = 0
	     else if(r.eq.nrow)then
	       neib(2,L) = 0
	       neib(4,L)  =0
	     endif

	   else
	     do j= 1,6
	       neib(j,L) = L + inc2(j)
	     enddo
	  
	     if (c .eq. 1) then
               neib(5,L) = 0
	     else if (c .eq. ncol) then
	       neib(1,L) = 0
	       neib(4,L) = 0
	       neib(6,L) = 0
	     endif

	     if (r .eq. nrow) then
	       neib(2,L) = 0
	       neib(4,L) = 0
	     endif

	   endif
        enddo


C_______Count exposed sides for cells in the contour_________________
	
          do i=1,n
	    icell=cell(i)
	    sides(icell)=0
            do j=1,6
              if(.not. exist( neib(j,icell) ) )then
                sides(icell)=sides(icell)+ 1
              endif
            enddo
          enddo

C________Grab surface cells___________________________________________

	nc=0
	do i=1,n
	  if(sides(cell(i)).gt.0)then
	    nc=nc+1
	    ncell(nc)=cell(i)
	  endif
	enddo
	n=nc           !n is now the number of exposed, non-empty cells

C_______The outer loop________________________________________________
C
C	temporary indices
C		nc: index for cells remaining on the list
C	        ninc:	index for newly exposed cells added to back of list 

	do while(n.gt.0)

C         Subtract exposed-side counts from the surface cell counts
C         until at least one cell is empty.

          loop=maxcnt
          do i=1,n
	    icell=ncell(i)
            loop=min( (ncnt(icell)-1)/sides(icell) , loop)
	  enddo
	  loop=loop+1     !all loop values are 1 too small

C         update the counts, rank and remove eroded cells
        
          nc=0
	  order=order+6
          ninc=n
          do i=1,n
	    icell=ncell(i)
            ncnt(icell)=ncnt(icell)-sides(icell)*loop
	    if(ncnt(icell).le.0)then

C	      Remove the empty cell and store it's order
	      exist(icell)=.false.
	      erode(icell)=order+ncnt(icell)

C	      Update the neighbors of the empty cell
	      do j=1,6
		nb=neib(j,icell)
	        if(exist(nb))then

C	          Store cells for addition to surface list
	          if(sides(nb).eq.0)then
	            ninc=ninc+1
	            ncell(ninc)=nb
	          endif

C	          Update sides for the neighbors
                  sides(nb)=sides(nb)+1
	        endif
	      enddo                     
	    else

C	      Save remaining cells
	      nc=nc+1
	      ncell(nc)=ncell(i)
	    endif
	  enddo

C	  Add new surface cells if any

	  do i=n+1,ninc,1
	    nc=nc+1
	    ncell(nc)=ncell(i)
	  enddo
	  n=nc
	enddo

C_______compress result___________________________________________


	do i=1,nnc
	  erode(i)=erode(cell(i))
	enddo
	n=nnc

	return
	end
