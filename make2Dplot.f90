! This program reads in a file in the format of struct_enum.out and
! makes a postscript output of the 2D structures
PROGRAM make2Dplot
use num_types
use vector_matrix_utilities
implicit none
character(80) fname, title, labeling
integer ioerr,  i,j,  is,js
integer k, strN, sizeN, nAt, diag(3), a,b,c,d,e,f,  fullL(3,3), pgOps,&
     & label, L(2,2), temp(2), Nspots, rows, cols, indx
real(sp) :: p(3,3),  scale
real(sp) :: xorig, yorig, parLat(2,2),point(2),yoff,xoff
 
read(*,*) fname
open(11,file=fname,status='old',iostat=ioerr)
if(ioerr/=0)then; write(*,'("Input file doesn''t exist:",a80)') trim(fname);endif
! Read in the title from the struct_enum.out file
read(11,'(a80)') title; title = trim(title)

! Read in the parent lattice vectors
do i = 1,3; read(11,*) p(:,i); enddo
parLat = p(2:3,2:3)

! Read in the number of labels, i.e., binary, ternary, etc.
read(11,'(i2)') k
read(11,*) ! Skip the line with the header labels

scale = .16
Nspots = 10
rows = 16
cols = 10
call init(scale,-10.0,-2.0)
call init_colors

xorig = 0.; yorig = 0.
! Need to add a loop over each structure
outer: do js = 1,rows
   do is = 1,cols
      xoff = (is-1)*Nspots
      yoff = 25.3/scale - (js-1)*Nspots
      !if ((js-1)*cols+is>160) exit outer ! Ternary
      !if ((js-1)*cols+is>155) exit outer ! binary
      if (yoff < 0) exit outer
      ! Read in the info for the given structure
      read(11,*,iostat=ioerr) strN, sizeN, nAt, pgOps, diag, a,b,c,d,e,f,&
      fullL,labeling
      if (ioerr/=0) exit
      L = reshape((/fullL(2,2),fullL(2,3),fullL(3,2),fullL(3,3)/),(/2,2/))
      !do i = 1,2
      !   write(*,'(2i2)') L(i,:)
      !enddo
      indx = diag(2)*diag(3)
      do i = 0,Nspots-2
         do j = 0,Nspots-2
            point = i*parLat(:,1)+j*parLat(:,2)+(/xoff,yoff/)
            temp = modulo(matmul(L,((/i,j-1/))),&
                 ((/diag(2),diag(3)/)))
            !if (temp(1) < 0) temp(1) = temp(1) + diag(2)
            !if (temp(2) < 0) temp(2) = temp(2) + diag(3)
            if (any(temp<0)) stop "less than zero"
            label = (temp(2)+temp(1)*diag(3))+1  
            if (labeling(label:label)=='1') then
                call color('black')
                call bullet(point(1),point(2),.45)
            else if (labeling(label:label)=='2') then
                  call color('blue')
                  call bullet(point(1),point(2),.35)
            else
               if (diag(2)==1) then; call color('red');
                  else; call color('green'); endif
                  call bullet(point(1),point(2),.35)
               endif
               !print *, point, labeling(label:label)
            enddo
         enddo

enddo
enddo outer
call draw
close(11)

END PROGRAM make2Dplot
