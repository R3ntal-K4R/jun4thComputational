!~~~~~~~ NUMERICAL APPROXIMATION OF 3D MAGNETIC FIELD SURROUNDING SPACECRAFT ~~~~~~~

! This program calculates the three Cartesian coordinates of the vector potential and
! associated magnetic field permeating the space surrounding a spacecraft. The 
! boundary conditions of this magnetic field, generated around the surface of the 
! spacecraft, are variable to experiment with shape, intensity and form of the 
! resultant magnetic field.

! Program created by: Lórien MacEnulty 		August, 2018


module global
implicit none
	
	real(8), parameter :: pi = 3.141592653589793238d0
	real(8) :: m,dimr,dims,dimx,dimy,dimz,cntr,A0
	integer :: n,bcchoice
	
	contains
	subroutine input_read
		print *, 'beginning read'
		!grid parameters
		read(*,*) n
		read(*,*) m					!{decimeters} in a meter
		cntr = real(n/2 + 1)		!center of n^3 grid
		print *, 'reading  dim'
	
		!spacecraft Index dimensions
		read(*,*) dimr				!index (non-physical) spherical radius
		read(*,*) dims				!index (non-physical) cylindrical radius
		read(*,*) dimx				!index (non-physical) Cartesian dimensions x,y,z
		read(*,*) dimy
		read(*,*) dimz
	
		read(*,*) A0
		read(*,*) bcchoice
			! bcchoice ->
			!		= 1, rectangular surface
			!		= 2, rectangular volume
			!		= 3, spherical surface
			!		= 4, cylindrical surface
			!		= 5, planes at ends of cylinder
		
	end subroutine input_read

end module global

!Define shape of boundaries
!======================================================================================
function boundary(x,y,z)
!======================================================================================
use global
implicit none
integer :: x,y,z,boundary,r,s

	s = sqrt(real(abs(cntr-x)**2) + real(abs(cntr-y)**2))
	r = sqrt(real(abs(cntr-x)**2) + real(abs(cntr-y)**2) + real(abs(cntr-z)**2))
	
	boundary = 0
	!bcchoice = 5
		! bcchoice ->
		!		= 1, rectangular surface
		!		= 2, rectangular volume
		!		= 3, spherical surface
		!		= 4, cylindrical surface
		!		= 5, planes at ends of cylinder
	
	if (bcchoice.eq.1) then						!rectangular surface bc
		if ((z.EQ.(cntr+dimz).AND.x.LE.(cntr+dimx).AND.x.GE.(cntr-dimx).AND.y.LE.(cntr+dimy).AND.y.GE.(cntr-dimy)) &
			.OR.(z.EQ.(cntr-dimz).AND.x.LE.(cntr+dimx).AND.x.GE.(cntr-dimx).AND.y.LE.(cntr+dimy).AND.y.GE.(cntr-dimy)) & 
			.OR.(y.EQ.(cntr+dimy).AND.x.LE.(cntr+dimx).AND.x.GE.(cntr-dimx).AND.z.LE.(cntr+dimz).AND.z.GE.(cntr-dimz)) &		
			.OR.(y.EQ.(cntr-dimy).AND.x.LE.(cntr+dimx).AND.x.GE.(cntr-dimx).AND.z.LE.(cntr+dimz).AND.z.GE.(cntr-dimz)) &
			.OR.(x.EQ.(cntr+dimx).AND.y.LE.(cntr+dimy).AND.y.GE.(cntr-dimy).AND.z.LE.(cntr+dimz).AND.z.GE.(cntr-dimz)) &
			.OR.(x.EQ.(cntr-dimx).AND.y.LE.(cntr+dimy).AND.y.GE.(cntr-dimy).AND.z.LE.(cntr+dimz).AND.z.GE.(cntr-dimz))) then
				boundary = 1
		end if
	else if (bcchoice.eq.2) then				!rectangular volume bc
		if ((x.LE.(cntr+dimx)).AND.(x.GE.(cntr-dimx)).AND.(y.LE.(cntr+dimy)) &
		.AND.(y.GE.(cntr-dimy)).AND.(z.LE.(cntr+dimz)).AND.(z.GE.(cntr+dimz))) then
			boundary = 1
		end if
	else if (bcchoice.eq.3) then				!spherical surface bc
		if (r.EQ.dimr) then
			boundary = 1
		end if
	else if (bcchoice.eq.4) then
		!cylindrical surface oriented with z-axis
		if ((s.EQ.dimr.AND.z.LT.(cntr+dimz).AND.z.GT.(cntr-dimz)) &
			.OR.(z.EQ.(cntr+dimz).AND.s.LE.dimr).OR.(z.EQ.(cntr-dimz).AND.s.LE.dimr)) then
				boundary = 1
		end if
	else										!ends of cylinder only
		if (z.EQ.(cntr+dimz)) then
			if (s.LE.dims) then
				boundary = 2
			else
				boundary = 1
			end if
		else if (z.EQ.(cntr-dimz)) then
			if (s.LE.dims) then
				boundary = 4
			else
				boundary = 5
			end if
		end if
	end if
	
	if (x.eq.1.or.x.eq.n.or.y.eq.1.or.y.eq.n.or.z.eq.1.or.z.eq.n) then		!edge of grid
		boundary = 3
	end if

return
end function boundary

!initialize grid
!======================================================================================
function initial(x,y,z,prmtr,indx)
!======================================================================================
use global
implicit none
real(8) :: initial
real(8) :: unitx,unity,unitz,theta,rad,s,r
integer :: x,y,z,prmtr,indx,boundary

	!boundary conditions
	!clockwise from top
	
	
	!Physical distance of index from origin (cntr,cntr,cntr) in meters
	unitx = real((x-cntr)*m)
	unity = real((y-cntr)*m)	
	unitz = real((z-cntr)*m)
	s = sqrt(real(abs(cntr-x)**2) + real(abs(cntr-y)**2))					!cylindrical		
	r = sqrt(unitx**2+unity**2+unitz**2)									!spherical
	
	!Components for boundary conditions
	rad = dims*m									!currently cylindrical
	!theta = asin(unitx/sqrt(unitx**2 + unitz**2))

	!define initial boundary conditions for Afield generated in spacecraft
	if (prmtr.eq.1) then
		if (indx.eq.2) then
			initial = -A0*unity
		else if (indx.eq.1) then
			initial = -A0*rad**2/s**2*unity
		else if (indx.eq.3) then
			initial = 0.d0
		else if (indx.eq.4) then
			initial = -A0*unity
		else if (indx.eq.5) then
			initial = -A0*rad**2/s**2*unity
		end if
	else if (prmtr.eq.2) then
		if (indx.eq.2) then
			initial = A0*unitx
		else if (indx.eq.1) then
			initial = A0*rad**2/s**2*unitx
		else if (indx.eq.3) then
			initial = 0.d0
		else if (indx.eq.4) then
			initial = A0*unitx
		else if (indx.eq.5) then
			initial = A0*rad**2/s**2*unitx
		end if
	else
		if (indx.eq.2) then
			initial = 0.d0
		else if (indx.eq.1) then
			initial = 0.d0
		else if (indx.eq.3) then
			initial = 0.d0
		else
			initial = 0.d0
		end if
	end if
	
return
end function initial

!Calculate next values of grid
!======================================================================================
function next(ac,an,ae,as,aw,au,ad)
!======================================================================================
use global
implicit none
real(8) :: next,res,fopt
real(8) :: an,ae,as,aw,ac,au,ad
	
	
	fopt = 2.0d0 - 2.0d0*pi/(n-2.0d0)
	res = (an + ae + as + aw + au + ad)/6.0d0 - ac
	next = ac + fopt*res
	
return
end function next

!Create new grid
!======================================================================================
subroutine newgrid(prm,Af)
!======================================================================================
use global
implicit none
real(8) :: initial,Af(n,n,n)
integer :: i,j,k,prm,indx,boundary

	do i=1,n
		do j=1,n
			do k=1,n
				indx = boundary(i,j,k)
				Af(i,j,k) = initial(i,j,k,prm,indx)
			end do
		end do
	end do
	
return
end subroutine


!Employ relaxation method
!======================================================================================
subroutine relax(prm,Af,count)
!======================================================================================
use global
implicit none
real(8) :: Af(n,n,n),next,wbcm,diff,vpot,s,rad,initial
integer :: i,j,k,x,l,count,prm,indx,boundary

	wbcm = 0.000001d0
	rad = dims*m
	

	do x=0,1
		do k=2,n-1
			do i=2,n-1
				l = mod((k + i + x),2)
				j = 3 - l
				do while ((j.gt.1).AND.(j.lt.n))
					vpot = next(Af(i,j,k),Af(i-1,j,k),Af(i,j+1,k),Af(i+1,j,k),&
							Af(i,j-1,k),Af(i,j,k+1),Af(i,j,k-1))
					diff = abs(vpot - Af(i,j,k))
					indx = boundary(i,j,k)
					if (indx.eq.0) then
						if (diff.le.wbcm) count = count + 1
						Af(i,j,k) = vpot
					else
						Af(i,j,k) = initial(i,j,k,prm,indx)
					end if
					j = j + 2
				end do
			end do
		end do
	end do

return
end subroutine relax

!calculate magnetic field from vector potential
!======================================================================================
subroutine magnetic(Ax,Ay,Az)
!======================================================================================
use global
implicit none
real(8) :: Ax(n,n,n),Ay(n,n,n),Az(n,n,n)
real(8) :: Bfldx(n,n,n),Bfldy(n,n,n),Bfldz(n,n,n)
real(8) :: dyaz,dzay,dxaz,dzax,dxay,dyax,di,dj,dk
integer :: i,j,k

di = m	!units in meters
dj = m
dk = m

open(22,file = 'Bfield.dat')

write(22,*) n,m,dims,cntr,dimz!, NEW_LINE('A')

do i=2,n-1
	do j=2,n-1
		do k=2,n-1
			dyaz = (Az(i,j+1,k)-Az(i,j-1,k))/2.d0*dj
			dzay = (Ay(i,j,k+1)-Ay(i,j,k-1))/2.d0*dk
			dxaz = (Az(i+1,j,k)-Az(i-1,j,k))/2.d0*di
			dzax = (Ax(i,j,k+1)-Ax(i,j,k-1))/2.d0*dk
			dxay = (Ay(i+1,j,k)-Ay(i-1,j,k))/2.d0*di
			dyax = (Ax(i,j+1,k)-Ax(i,j-1,k))/2.d0*dj
			
			Bfldx(i,j,k) = dyaz - dzay
			Bfldy(i,j,k) = dzax - dxaz
			Bfldz(i,j,k) = dxay - dyax
			write(22,*) i-2,j-2,k-2,Bfldx(i,j,k),Bfldy(i,j,k),Bfldz(i,j,k)
		end do
	end do
end do

close(22)

return
end subroutine magnetic

!Compare spherical solutions to analytical model
!======================================================================================
function analytical(x,y,z)
!======================================================================================
use global
implicit none
real(8) :: rad,unitx,unity,unitz,r,theta,E0,analytical
integer :: x,y,z

	rad = dims*m
	unitx = real((x-cntr)*m) 								!value of distance of index from
	unity = real((y-cntr)*m)								!origin (cntr,cntr,cntr) in meters
	unitz = real((z-cntr)*m)
	r = sqrt(unitx**2+unity**2+unitz**2)
	theta = asin(unitx/sqrt(unitx**2 + unitz**2))
	E0 = 1.d0

	if (r.le.rad) then										!potential inside sphere
		analytical = 0.d0
	else
		analytical = -E0*(r - rad**3/r**2)*unitz/r			!potential outside sphere
	end if

return
end function analytical

!======================================================================================
program laplacian
!======================================================================================	
	use global
	implicit none
	real(8), allocatable :: A1(:,:,:),A2(:,:,:),A3(:,:,:),Bfield(:,:,:)
	real(8) :: start,finish,analytical,theory
	integer :: t,dt,i,j,l,k,count,converge,x
	print *, 'reading input'
    call input_read
    print *, 'Starting'
    call cpu_time(start)
	
	t = 1000
	converge = 0.95*(n-2)**3
	
	allocate(A1(n,n,n),A2(n,n,n),A3(n,n,n),Bfield(n,n,n))
	print *, 'allocation finished'
	do x=1,3
		if (x.eq.1) then
			print *, 'Commencing x-coordinate...'
			call newgrid(1,A1)
			do dt=1,t
				count = 0
				!print *,dt
				call relax(1,A1,count)
				if (count.ge.converge) exit
			end do
		else if (x.eq.2) then
			print *, 'Commencing y-coordinate...'
			call newgrid(2,A2)
			do dt=1,t
				!print *,dt
				count = 0
				call relax(2,A2,count)
				if (count.ge.converge) exit
			end do
		else
			print *, 'Commencing z-coordinate...'
			A3 = 0.d0
			call newgrid(3,A3)
			do dt=1,t
				count = 0
				!print *,dt
				call relax(3,A3,count)
				if (count.ge.converge) exit
			end do
		end if
	end do

	print *, 'Calculating magnetic field...'
	call magnetic(A1,A2,A3)

	!print and write grids
	print *, 'Writing files...'
!	open(21,file = 'vp.dat')
	!open(22,file = 'vp_analytical.plt')
!	do i=2,n-1
!		do j=2,n-1
!			do k=2,n-1
!				write(21,*) i-2,j-2,k-2,A1(i,j,k),A2(i,j,k),A3(i,j,k)
				!theory = analytical(i,j,k)
				!write(22,*) i,j,k,theory
!			end do
!		end do
!	end do
!	close(21)
	!close(22)

	!print run time
    call cpu_time(finish)
    if (finish-start.lt.60.00) then
    	print '("Time = ",f10.3," seconds")',(finish-start)
    else
    	print '("Time = ",f10.3," minutes")',(finish-start)/60.00
    end if
    
    print *, 'Algorithm complete. Live long and prosper.'

end program laplacian