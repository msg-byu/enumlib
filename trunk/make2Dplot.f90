! This program reads in a file in the format of struct_enum.out and
! makes a postscript output of the 2D structures
PROGRAM make2Dplot
use num_types
use vector_matrix_utilities
implicit none
character(80) fname, title, labeling, dummy
integer ioerr,  i,j,  is,js, iD, nD
integer k, strN, hnfN, sizeN, nAt, diag(3), a,b,c,d,e,f,  fullL(3,3), pgOps,&
     & label, L(2,2), temp(2), Nspots, rows, cols, indx
real(sp) :: p(3,3),  scale, spotsize
real(sp) :: xorig, yorig, parLat(2,2),point(2),yoff,xoff
real(dp), allocatable :: dvec(:,:)

if(iargc()/=1) stop "2Dplot.x requires as an argument &
     & the filename to read from"
call getarg(1,dummy)
read(dummy,*) fname

open(11,file=fname,status='old',iostat=ioerr)
if(ioerr/=0)then; write(*,'("Input file doesn''t exist:",a80)') trim(fname);endif
! Read in the title from the struct_enum.out file
read(11,'(a80)') title; title = trim(title)
read(11,*) dummy ! Read in surf/bulk
 
! Read in the parent lattice vectors
do i = 1,3; read(11,*) p(:,i); enddo
parLat = p(2:3,2:3)

! Read in the basis vectors for the d-set
read(11,*) nD
allocate(dvec(3,nD))
do i = 1,nD
   read(11,*) dvec(:,i)
enddo

! Read in the number of labels, i.e., binary, ternary, etc.
read(11,'(i2)') k
do i = 1,4;read(11,*);enddo ! Skip 4 lines

scale = .19
Nspots = 10
rows = 16
cols = 10
spotsize = .35
call init(scale,-10.0,-2.0)
call init_colors

xorig = 0.; yorig = 0
! Need to add a loop over each structure
outer: do js = 1,rows
   do is = 1,cols
      xoff = (is-1)*Nspots
      yoff = 28/scale - (js-1)*Nspots
      !if ((js-1)*cols+is>160) exit outer ! Ternary
      !if ((js-1)*cols+is>155) exit outer ! binary
      if (yoff < 0) exit outer
      ! Read in the info for the given structure
      read(11,*,iostat=ioerr) strN, hnfN,sizeN, nAt, pgOps, diag, a,b,c,d,e,f,&
      fullL,labeling
      if (ioerr/=0) exit
      L = reshape((/fullL(2,2),fullL(2,3),fullL(3,2),fullL(3,3)/),(/2,2/))
      do i = 1,2
         write(*,'(2i2)') L(i,:)
      enddo
      indx = diag(2)*diag(3)
      do i = 0,Nspots-2
         do j = 0,Nspots-2
            do iD = 1, nD
            ! point is a 2-vector, Cartesian coordinates of the point
            point = i*parLat(:,1)+j*parLat(:,2)+(/xoff,yoff/)+dvec(2:3,iD)
            ! i, j are the parent lattice direct coordinates of a point 
            ! (not including the d-vector index). In other words, (i,j)=A^-1.r
            ! So temp is (LA^-1)r which maps the vector into the group
            temp = modulo(matmul(L,((/i,j/))),&
                 ((/diag(2),diag(3)/)))
            !write(*,'("i,j: ",2(i2,1x))') i,j
            ! temp will now contain the "group member index" of the point (i,j)

            if (any(temp<0)) stop "less than zero"
            ! the rhs of label= converts the "group member index" to a single index in the labeling
            !string. It also accounts for the d-vector index
            ! The coordinates in "point" tell us where to draw a circle; "label" tells us which
            ! color belongs to that site.
            label = (iD-1)*indx+(temp(2)+temp(1)*diag(3))+1  
            write(*,'("i,j:",2(1x,i2),5x,"g coords:",2(1x,i1),5x,"label:",i1)') i,j,temp,label
            if (labeling(label:label)=='1') then
                call color('blue')
                call bullet(point(1),point(2),spotsize)
            else if (labeling(label:label)=='2') then
                  call color('green')
                  call bullet(point(1),point(2),spotsize)
            else
               if (diag(2)==1) then; call color('red');
                  else; call color('red'); endif
                  call bullet(point(1),point(2),spotsize)
               endif
               !print *, point, labeling(label:label)
            enddo
            enddo
         enddo

enddo
enddo outer
call draw
close(11)

END PROGRAM make2Dplot
