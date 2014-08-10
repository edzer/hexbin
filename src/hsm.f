C File:		hsm.f
C Programmer:	Daniel B. Carr
C Version Date: January 3, 1994
C
C This program is an hexagon cell smoother.  It smooths into
C neighboring cells and hence expands.  

C The kernal is a crude integer kernel.
C The boundary hexagons get weight 1, the center hexagon
C gets weight, wt, which by default is set to six.
C
C 
 
	subroutine hsm(cell,cnt,n,nmax,sm,ncol,wt)

	implicit none

	integer n, nmax, ncol
	integer cell(*), cnt(*), sm(*), wt(*)
	integer ind, ind1(6), ind2(12),ind3(6), ind4(12), loc
	integer row, cnt1, cnt2, wta, wtb, wtc
	integer i, j

C__________Constants___________________________________________

	ind1(1)=-1
	ind1(2)=ncol-1
	ind1(3)=ncol
	ind1(4)=+1
	ind1(5)=-ncol
	ind1(6)=-ncol-1

	ind2(1)=-2
	ind2(2)=ncol-2
	ind2(3)=2*ncol-1
	ind2(4)=2*ncol
	ind2(5)=2*ncol+1
	ind2(6)=ncol+1
	ind2(7)=2
	ind2(8)=-ncol+1
	ind2(9)=-2*ncol+1
	ind2(10)=-2*ncol
	ind2(11)=-2*ncol-1
	ind2(12)=-ncol-2

	ind3(1)=-1
	ind3(2)=ncol
	ind3(3)=ncol+1
	ind3(4)=+1
	ind3(5)=-ncol+1
	ind3(6)=-ncol

	ind4(1)=-2
	ind4(2)=ncol-1
	ind4(3)=2*ncol-1
	ind4(4)=2*ncol
	ind4(5)=2*ncol+1
	ind4(6)=ncol+2
	ind4(7)=2
	ind4(8)=-ncol+2
	ind4(9)=-2*ncol+1
	ind4(10)=-2*ncol
	ind4(11)=-2*ncol-1
	ind4(12)=-ncol-1
	
	wta = wt(1)
	wtb = wt(2)
	wtc = wt(3)

C_________Smoothing_____________________________________

	do i=1,n
	   sm(cell(i))=wta*cnt(i)
	enddo

	do i=1,n
	   loc=cell(i)
	   row=(loc-1)/ncol + 1
	   cnt1=wtb*cnt(i)
	   cnt2=wtc*cnt(i)

	   if(mod(row,2).eq.1)then
	      do j=1,6
	         ind=loc+ind1(j)
		 sm(ind)=sm(ind)+cnt1
	      enddo
	      do j=1,12
		 ind=loc+ind2(j)
		 sm(ind)=sm(ind)+cnt2
	      enddo
	   else
	      do j=1,6
	         ind=loc+ind3(j)
		 sm(ind)=sm(ind)+cnt1
	      enddo
	      do j=1,12
		 ind=loc+ind4(j)
		 sm(ind)=sm(ind)+cnt2
	      enddo
	   endif
	enddo

	n=0
	do i=1,nmax
	   if(sm(i).gt.0)then
	      n=n+1
	      cell(n)=i
	      cnt(n)=sm(i)
	   endif
	enddo
	return
	end


