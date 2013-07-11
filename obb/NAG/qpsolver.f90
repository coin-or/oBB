subroutine qpsolver(h,cvec,bl,bu,n,x,obj)
! Solves QP c^Tx + 1/2x^THx
! subject to l <= x <= u 
! Generate signature file: f2py -m qpsolver -h qpsolver.pyf qpsolver.f90
! Then compile with: f2py -c --fcompiler=intelem qpsolver.pyf -I$nag_dir/nag_interface_blocks qpsolver.f90 $nag_dir/lib/libnag_mkl.so

! Use Statements ..
use nag_library, only : e04nfa, e04wbf, e54nfu, nag_wp

! Implicit None Statement ..
implicit none

! Inputs & Outputs
integer, intent(in) 				:: n
real (kind=nag_wp), dimension(n,n), intent(in) 	:: h
real (kind=nag_wp), dimension(n), intent(in) 	:: cvec
real (kind=nag_wp), dimension(n), intent(in) 	:: bl
real (kind=nag_wp), dimension(n), intent(in) 	:: bu
real (kind=nag_wp), dimension(n), intent(inout) :: x 
real (kind=nag_wp), intent(out) 		:: obj

!f2py integer, intent(in) 				:: n
!f2py real (kind=nag_wp), dimension(n,n), intent(in) 	:: h
!f2py real (kind=nag_wp), dimension(n), intent(in) 	:: cvec
!f2py real (kind=nag_wp), dimension(n), intent(in) 	:: bl
!f2py real (kind=nag_wp), dimension(n), intent(in) 	:: bu
!f2py real (kind=nag_wp), dimension(n), intent(inout) :: x 
!f2py real (kind=nag_wp), intent(out) 		:: obj

! Local Parameters ..
integer, parameter              :: lcwsav = 1, liwsav = 610, llwsav = 120, lrwsav = 475, nclin = 0, sda = 1, lda = 1

! Local Scalars ..
integer                         :: ifail, iter, ldh, liwork, lwork

! Local Arrays ..
real (kind=nag_wp), dimension(lda,sda) 	:: a 
real (kind=nag_wp), dimension(1) 	:: ax
real (kind=nag_wp), dimension(n) 	:: clamda
real (kind=nag_wp), dimension(n**2 + 8*n) 	:: work			
real (kind=nag_wp), dimension(1)   	:: ruser 
real (kind=nag_wp), dimension(lrwsav)	:: rwsav
integer, dimension(n)			:: istate 
integer, dimension(2*n + 3) 		:: iwork
integer, dimension(liwsav)              :: iwsav
logical, dimension(llwsav)              :: lwsav
character (len=80), dimension(lcwsav)   :: cwsav

! Executable Statements
ldh = n 
liwork = 2*n + 3
lwork = n**2 + 8*n

! Initialise E04NFA
ifail = 0
CALL e04wbf('E04NFA',cwsav,lcwsav,lwsav,llwsav,iwsav,liwsav,rwsav,lrwsav,ifail)

! Solve the problem
ifail = -1
CALL e04nfa(n,nclin,a,lda,bl,bu,cvec,h,ldh,e54nfu,istate,x,iter,obj,ax, &
	clamda,iwork,liwork,work,lwork,iwsav,ruser,lwsav,iwsav,rwsav,ifail)

end subroutine qpsolver

