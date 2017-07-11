program main

use omp_lib
implicit none

integer :: div
double precision :: startTime, endTime,startTime2, endTime2, sum, sum0, sum1, sum2, gflops
character :: arg*10

! node coordinate
integer :: x1, y1, z1
integer :: x2, y2, z2

! SpMV iteration times
integer, parameter :: itermax = 256 ! number of SpMV

! Matrix
integer :: nnodes_edge =  100   ! number of nodes on a edge of cube
real(8), parameter :: mev         =  1.2d0 ! value of matrix element
real(8), parameter :: vev         = 10.2d0 ! value of vector element

integer(8) :: ncols, nnz
integer(8) :: matdim
integer(8) :: nne, nne2

integer(8),    pointer :: irow(:)
integer(8),    pointer :: jcol(:)
real(8),    pointer :: val(:)

! Vector

real(8),    pointer :: X(:) ! given vector

real(8),    pointer :: y(:)     ! result (Ab=y)

! CRS format
integer(8) :: row
integer ti1, ti2, t_rate, t_max, diff

real(kind=8) :: START_TIME, END_TIME, Tcomm
integer(kind=8) :: i, j, jS, jE, in
real(kind=8) :: YV1, YV2, YV3, XV1, XV2, XV3

integer(kind=8) :: N, NP
integer(kind=8), pointer :: indexL(:), itemL(:), indexU(:), itemU(:)
real(kind=8), pointer :: AL(:), AU(:), D(:)

integer,       parameter :: coef_block = 9     ! BCSR parameter : number of element in a block (fix to 9)
integer,       parameter :: vec_coef_block = 3 ! BCSR parameter : number of element in a block (fix to 9)
integer, parameter :: numOfBlockPerThread = 100
logical, save :: isFirst = .true.
integer, save :: numOfThread = 1
integer, save, allocatable :: startPos(:), endPos(:)
integer(kind=8), save :: sectorCacheSize0, sectorCacheSize1
integer(kind=8) :: threadNum, blockNum, numOfBlock
integer(kind=8) :: numOfElement, elementCount, blockIndex
real(kind=8) :: numOfElementPerBlock

! elap time
real(4) :: etime, tarray(2), t1, t2, tsum

!
integer(8) :: k,l,iter,i2,i3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! create matrix A
!
!

if(iargc()>0)then 
   call getarg( 1, arg )
   read(arg,*) nnodes_edge
end if

nne  = nnodes_edge
nne2 = nnodes_edge * nnodes_edge

matdim = nnodes_edge * nnodes_edge * nnodes_edge

N = matdim
ncols = matdim

write(*,"(a30,i13)") 'Num of nodes on edge: nne',nne !DEBUG
write(*,"(a30,i13)") 'nne**2',nne2 !DEBUG
write(*,"(a30,i13)") 'Num of rows',N   !DEBUG
write(*,"(a30,i13)") 'Num of columns',ncols !DEBUG




! count number of non zero elements
nnz=0

t1 = etime(tarray) ! TIMER

do z1 = 1, nne
  do y1 = 1, nne
    do x1 = 1, nne
      i = nne2*(z1 - 1) + nne*(y1 - 1) + x1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! z=-1
x2=x1-1
y2=y1-1
z2=z1-1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1
y2=y1-1
z2=z1-1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1+1
y2=y1-1
z2=z1-1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1-1
y2=y1
z2=z1-1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1
y2=y1
z2=z1-1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1+1
y2=y1
z2=z1-1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1-1
y2=y1+1
z2=z1-1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1
y2=y1+1
z2=z1-1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1+1
y2=y1+1
z2=z1-1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! z=0
x2=x1-1
y2=y1-1
z2=z1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1
y2=y1-1
z2=z1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1+1
y2=y1-1
z2=z1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1-1
y2=y1
z2=z1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1
y2=y1
z2=z1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1+1
y2=y1
z2=z1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1-1
y2=y1+1
z2=z1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1
y2=y1+1
z2=z1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1+1
y2=y1+1
z2=z1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! z=+1
x2=x1-1
y2=y1-1
z2=z1+1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1
y2=y1-1
z2=z1+1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1+1
y2=y1-1
z2=z1+1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1-1
y2=y1
z2=z1+1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1
y2=y1
z2=z1+1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1+1
y2=y1
z2=z1+1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1-1
y2=y1+1
z2=z1+1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1
y2=y1+1
z2=z1+1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1+1
y2=y1+1
z2=z1+1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  nnz = nnz+1
end if

    end do
  end do
end do


write(*,"(a30,i13)") 'nnz',nnz !DEBUG

t2 = etime(tarray) ! TIMER
!write(*,*) 'user time for count number of non-zero element (sec)', t2-t1 !TIMER
write(*,"(a30,f13.4)") 'Count non-zero elmement ',t2-t1 !DEBUG

allocate(irow(nnz), jcol(nnz))

 t1 = etime(tarray) ! TIMER
 
 ! set irow jcol matrix
 
 k=1
 
 do z1 = 1, nne
   do y1 = 1, nne
     do x1 = 1, nne
       i = nne2*(z1 - 1) + nne*(y1 - 1) + x1
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! z=-1
x2=x1-1
y2=y1-1
z2=z1-1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1
y2=y1-1
z2=z1-1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1+1
y2=y1-1
z2=z1-1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1-1
y2=y1
z2=z1-1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1
y2=y1
z2=z1-1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1+1
y2=y1
z2=z1-1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1-1
y2=y1+1
z2=z1-1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1
y2=y1+1
z2=z1-1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1+1
y2=y1+1
z2=z1-1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! z=0
x2=x1-1
y2=y1-1
z2=z1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1
y2=y1-1
z2=z1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1+1
y2=y1-1
z2=z1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1-1
y2=y1
z2=z1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1
y2=y1
z2=z1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1+1
y2=y1
z2=z1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1-1
y2=y1+1
z2=z1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1
y2=y1+1
z2=z1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1+1
y2=y1+1
z2=z1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! z=+1
x2=x1-1
y2=y1-1
z2=z1+1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1
y2=y1-1
z2=z1+1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1+1
y2=y1-1
z2=z1+1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1-1
y2=y1
z2=z1+1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1
y2=y1
z2=z1+1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1+1
y2=y1
z2=z1+1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1-1
y2=y1+1
z2=z1+1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1
y2=y1+1
z2=z1+1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x2=x1+1
y2=y1+1
z2=z1+1
j=nne2*(z2-1)+nne*(y2-1)+x2
if ((x2 .ge. 1) .and. (y2 .ge. 1) .and. (z2 .ge. 1) .and. (x2 .le. nne) .and. (y2 .le. nne) .and. (z2 .le. nne)) then
  irow(k) = i
  jcol(k) = j
  k = k + 1
end if
     end do
   end do
 end do

t2 = etime(tarray) ! TIMER
write(*,"(a30,f13.4,a)") 'create irow jcol matrix ',t2-t1," sec" !DEBUG


  t1 = etime(tarray) ! TIMER
  allocate(indexL(0:N)) 
  allocate(itemL(nnz))
  allocate(AL(nnz*coef_block))

  indexL(0) = 0

  row = 1 
  do i=1,nnz
  ! matrix format error check
    if (row > irow(i)) then
      write(*,*) 'ERROR: row is not sorted' 
      stop
    end if
    if (row < irow(i)-1) then
      write(*,*) 'ERROR: vacant row' 
      stop
    end if

  ! new row?
    if (row == irow(i)-1) then 
      indexL(row) = i-1
    end if

    row = irow(i)
    itemL(i) = jcol(i) 

  ! set block csr val
    AL(9*i-8)  = mev*i+1
    AL(9*i-7)  = mev*i+2
    AL(9*i-6)  = mev*i+3
    AL(9*i-5)  = mev*i+4
    AL(9*i-4)  = mev*i+5
    AL(9*i-3)  = mev*i+6
    AL(9*i-2)  = mev*i+7
    AL(9*i-1)  = mev*i+8
    AL(9*i  )  = mev*i+9
  end do
      indexL(N)=nnz 


t2 = etime(tarray) ! TIMER
!write(*,*) 'user time for create CSR matrix (sec)', t2-t1 !TIMER
write(*,"(a30,f13.4,a)") 'Create CSR matrix ',t2-t1, " sec" !DEBUG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! set vector b for Ab=x
!

t1 = etime(tarray) ! TIMER
allocate(X(N*vec_coef_block))

do i=1, N*vec_coef_block
  X(i) = vev
end do
!write(*,*) 'X ',X !DEBUG

t2 = etime(tarray) ! TIMER
write(*,"(a30,f13.4,a)") 'Create vector b ', t2-t1, " sec"

write(*,"(a30,i13,a)") 'Length of vector ', ncols
write(*,"(a30,i13,a)") 'Size of vector ', vec_coef_block*ncols*8/1024**2 , " MB"

write(*,"(a30,i13)") 'Num of non-zero elem. ', nnz * coef_block
write(*,"(a30,i13,a)") 'Size of non-zero elem. ', nnz * coef_block * 8/1024**2, " MB"

write(*,"(a30,i13)") 'Num of SpMV iteration ', itermax

allocate(y(N*vec_coef_block))
tsum=0

      
      




! added for turning >>>
if (.not. isFirst) then
  numOfBlock = numOfThread * numOfBlockPerThread
  if (endPos(numOfBlock-1) .ne. N-1) then
    deallocate(startPos, endPos)
    isFirst = .true.
  endif
endif
if (isFirst) then
  !$ numOfThread = omp_get_max_threads()
  numOfBlock = numOfThread * numOfBlockPerThread
  allocate (startPos(0 : numOfBlock - 1), endPos(0 : numOfBlock - 1))
  numOfElement = N + indexL(N)
  numOfElementPerBlock = dble(numOfElement) / numOfBlock
  blockNum = 0
  elementCount = 0
  startPos(blockNum) = 1
  do i= 1, N
    elementCount = elementCount + 1
    elementCount = elementCount + (indexL(i) - indexL(i-1))
    if (elementCount > (blockNum + 1) * numOfElementPerBlock) then
      endPos(blockNum) = i
      ! write(9000+hecMESH%my_rank,*) mod(blockNum, numOfThread), &
      !      startPos(blockNum), endPos(blockNum)
      blockNum = blockNum + 1
      startPos(blockNum) = i + 1
      if (blockNum == (numOfBlock - 1)) exit
    endif
  enddo
  endPos(blockNum) = N

  do i= blockNum+1, numOfBlock-1
    startPos(i) = N
    endPos(i) = N-1
  end do

  isFirst = .false.
endif

div=1

sum=0.0d0
call system_clock(ti1) 

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP & PRIVATE(i,XV1,XV2,XV3,YV1,YV2,YV3,jS,jE,j,in,threadNum,blockNum,blockIndex,sum,endTime,& 
!$OMP & startTime,endTime2,startTime2,sum1,sum2) &
!$OMP & SHARED(D,AL,AU,indexL,itemL,indexU,itemU,X,Y,startPos,endPos,numOfThread,N,div,tsum)
      threadNum = 0
      sum=0.0d0
      !$ threadNum = omp_get_thread_num()
  do blockNum = 0 , numOfBlockPerThread - 1
      blockIndex = blockNum * numOfThread  + threadNum


    do iter=1, itermax
!$omp barrier
      startTime2 = omp_get_wtime()

      !t1 = etime(tarray)
      
      do i = startPos(blockIndex), endPos(blockIndex)
        !Y = 0.0d0
        X(3*i-2) = vev * iter* i
        X(3*i-1) = vev * iter*0.5/i
        X(3*i-0) = -vev * iter *0.334*i
      end do
      endTime2 = omp_get_wtime()
      sum2 = sum2 + (endTime2 - startTime2)
!$omp barrier
      startTime = omp_get_wtime()

      YV1 = 0.0d0
      YV2 = 0.0d0
      YV3 = 0.0d0
      do i = startPos(blockIndex), endPos(blockIndex)
        XV1= X(3*i-2)
        XV2= X(3*i-1)
        XV3= X(3*i  )

        jS= indexL(i-1) + 1
        jE= indexL(i  )
        do j= jS, jE
          in  = itemL(j)
          XV1= X(3*in-2)
          XV2= X(3*in-1)
          XV3= X(3*in  )
          YV1= YV1 + AL(9*j-8)*XV1 + AL(9*j-7)*XV2 + AL(9*j-6)*XV3
          YV2= YV2 + AL(9*j-5)*XV1 + AL(9*j-4)*XV2 + AL(9*j-3)*XV3
          YV3= YV3 + AL(9*j-2)*XV1 + AL(9*j-1)*XV2 + AL(9*j  )*XV3
        enddo
        Y(3*i-2)= YV1
        Y(3*i-1)= YV2
        Y(3*i  )= YV3
      enddo
      endTime = omp_get_wtime()
      sum = sum + (endTime - startTime)
!$omp barrier
      
    end do
    
  enddo
    if (omp_get_thread_num()==0) then 
      tsum=sum
      write(*,"(a30,f13.4,a)") 'Time for Clearing Vect', sum2, " sec"

    end if 
!$OMP END PARALLEL

call system_clock(ti2, t_rate)

!      write(*,*) 'y' ! DEBUG
!      write(*,*) y ! DEBUG

      write(*,"(a30,f13.4,a)") 'Time for SpMV iteration.', ((ti2-ti1)/dble(t_rate)), " sec"
      write(*,"(a30,f13.4,a)") 'Time for only part of SpMV ', tsum, " sec"

      gflops = (((nnz*(2.0*coef_block)*itermax)/div)/tsum)/(10.0**9.0)
      write(*,"(a30,f13.4,a)") 'Performance:   ', gflops, ' GFLOPS'

      gflops = (((nnz*(8*coef_block+8+8*2*3)*itermax)/div)/tsum)/(10.0**9.0)
      write(*,"(a30,f13.4,a)") 'Memory(Cache) throughput:   ', gflops, ' GB/s'
      
      gflops = (((nnz*(8*coef_block+8)*itermax)/div)/tsum)/(10.0**9.0)
      write(*,"(a30,f13.4,a)") 'If X,Y Vector cached then ', gflops, ' GB/s'
    !write(*,*) 'Came here 1'
    !write(*,*) y
    !write(*,*) 'Came here 2'

      stop
      end

