subroutine qpsolver_lincon(h,cvec,bl,bu,n,x,obj,nclin,a)
! Solves QP c^Tx + 1/2x^THx
! subject to l <= (x,Ax)^T <= u 
! Generate signature file: f2py -m qpsolver_lincon -h qpsolver_lincon.pyf qpsolver_lincon.f90
! Then Compile with: f2py -c --fcompiler=intelem qpsolver_lincon.pyf -I$nag_dir/nag_interface_blocks qpsolver_lincon.f90 $nag_dir/lib/libnag_mkl.so

! Use Statements ..
use nag_library, only : e04nfa, e04wbf, e54nfu, nag_wp

! Implicit None Statement ..
implicit none

! Inputs & Outputs
integer, intent(in) 				:: n
real (kind=nag_wp), dimension(n,n), intent(in) 	:: h
real (kind=nag_wp), dimension(n), intent(in) 	:: cvec
real (kind=nag_wp), dimension(n+nclin), intent(in) 	:: bl
real (kind=nag_wp), dimension(n+nclin), intent(in) 	:: bu
real (kind=nag_wp), dimension(n), intent(inout) :: x 
real (kind=nag_wp), intent(out) 		:: obj
integer, intent(in) 				:: nclin
real (kind=nag_wp), dimension(nclin,n), intent(in) 	:: a

!f2py integer, intent(in) 				:: n
!f2py real (kind=nag_wp), dimension(n,n), intent(in) 	:: h
!f2py real (kind=nag_wp), dimension(n), intent(in) 	:: cvec
!f2py real (kind=nag_wp), dimension(n+nclin), intent(in) 	:: bl
!f2py real (kind=nag_wp), dimension(n+nclin), intent(in) 	:: bu
!f2py real (kind=nag_wp), dimension(n), intent(inout) :: x 
!f2py real (kind=nag_wp), intent(out) 		:: obj
!f2py integer, intent(in) 				:: nclin
!f2py real (kind=nag_wp), dimension(nclin,n), intent(in) 	:: a

! Local Parameters ..
integer, parameter              :: lcwsav = 1, liwsav = 610, llwsav = 120, lrwsav = 475

! Local Scalars ..
integer                         :: ifail, iter, ldh, liwork, lwork

! Local Arrays ..
real (kind=nag_wp), dimension(nclin) 	:: ax
real (kind=nag_wp), dimension(n+nclin) 	:: clamda
real (kind=nag_wp), dimension(2*n**2 + 8*n + 5*nclin) 	:: work			
real (kind=nag_wp), dimension(1)   	:: ruser 
real (kind=nag_wp), dimension(lrwsav)	:: rwsav
integer, dimension(n+nclin)			:: istate 
integer, dimension(2*n + 3) 		:: iwork
integer, dimension(liwsav)              :: iwsav
logical, dimension(llwsav)              :: lwsav
character (len=80), dimension(lcwsav)   :: cwsav

! Executable Statements
ldh = n 
liwork = 2*n + 3
lwork = 2*n**2 + 8*n + 5*nclin

! Initialise E04NFA
ifail = 0
CALL e04wbf('E04NFA',cwsav,lcwsav,lwsav,llwsav,iwsav,liwsav,rwsav,lrwsav,ifail)

! Solve the problem
ifail = -1
CALL e04nfa(n,nclin,a,nclin,bl,bu,cvec,h,ldh,e54nfu,istate,x,iter,obj,ax, &
	clamda,iwork,liwork,work,lwork,iwsav,ruser,lwsav,iwsav,rwsav,ifail)

end subroutine qpsolver_lincon

