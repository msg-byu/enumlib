C 2-dimensional plotting routine

c initialization

	subroutine init(scale_in,xorigin,yorigin)

	real origin(2)
	logical ex

	common scale,origin
	common /common_psfile/filename
	common /common_bbox/bbxmin,bbxmax,bbymin,bbymax

* open list file
	open(30,file='splot.lis')
* copy input parameters to list file
	write(30,*)'init'
	write(30,*)'scale = ',scale_in
	write(30,*)'origin = ',xorigin,yorigin
	write(30,*)
* open plot file
        inquire(31,opened=ex)
        if(.not.ex)then
          open(31,file='splot.ps')
	  write(31,'(a)')'%!PS-Adobe-2.0'
        endif
* origin
	origin(1)=xorigin
	origin(2)=yorigin
* scale
	scale=scale_in*72/2.54
* line width
	widthline=0.08
	write(31,*)widthline,' setlinewidth'
* black
	write(31,*)'0 setgray'
* bounding box
	bbxmin=1e6
	bbxmax=-1e6
	bbymin=1e6
	bbymax=-1e6
	end

* end plot
	subroutine draw
	common scale,origin
	common /common_bbox/bbxmin,bbxmax,bbymin,bbymax
	write(30,*)'end of plot'
	write(31,*)'showpage'
	write(6,*)bbxmin,bbxmax,bbymin,bbymax,scale
	end

c error routine

	subroutine splot_err(x)
	character*(*) x
	character format_string*9
	format_string='(1x,a   )'
	length=len(x)
	write(format_string(6:8),'(i3)')length
	write(6,*)
	write(6,*)'**** Fatal Error ****'
	write(6,*)
	write(6,format_string)x
	write(6,*)
	stop 'This program has bombed!'
	end

c open specific output ps file

	subroutine psfile(filename)
	character(*) :: filename
        open(31,file=filename)
        end

* translate origin

	subroutine translate(x,y)
	call point(x,y,xp,yp)
	write(30,*)'translate',x,y
	write(31,*)xp,yp,' translate'
	end

c line

	subroutine line(n,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,
	2     x8,y8,x9,y9,x10,y10,x11,y11,x12,y12,x13,y13,x14,y14,x15,
	2     y15,x16,y16,x17,y17,x18,y18,x19,y19,x20,y20)

	dimension x(20),y(20)

	write(30,*)'line: ',n,' points'
	if(n.lt.2.or.n.gt.20)
	2     call splot_err('1st argument of LINE is invalid')

	goto (2,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),n

20	x(20)=x20
	y(20)=y20
19	x(19)=x19
	y(19)=y19
18	x(18)=x18
	y(18)=y18
17	x(17)=x17
	y(17)=y17
16	x(16)=x16
	y(16)=y16
15	x(15)=x15
	y(15)=y15
14	x(14)=x14
	y(14)=y14
13	x(13)=x13
	y(13)=y13
12	x(12)=x12
	y(12)=y12
11	x(11)=x11
	y(11)=y11
10	x(10)=x10
	y(10)=y10
9	x(9)=x9
	y(9)=y9
8	x(8)=x8
	y(8)=y8
7	x(7)=x7
	y(7)=y7
6	x(6)=x6
	y(6)=y6
5	x(5)=x5
	y(5)=y5
4	x(4)=x4
	y(4)=y4
3	x(3)=x3
	y(3)=y3
2	x(2)=x2
	y(2)=y2
	x(1)=x1
	y(1)=y1

	do i=1,n
	  write(30,*)x(i),y(i)
	enddo
	write(30,*)

	write(31,*)'newpath'
	call point(x(1),y(1),xp,yp)
	call bbox_limits(xp,yp)
	write(31,*)xp,yp,' moveto'
	do i=2,n
	  call point(x(i),y(i),xp,yp)
	  call bbox_limits(xp,yp)
	  write(31,*)xp,yp,' lineto'
	enddo
	write(31,*)'stroke'

	end

c line

	subroutine line_white(n,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,
	1    x7,y7,
	2     x8,y8,x9,y9,x10,y10,x11,y11,x12,y12,x13,y13,x14,y14,x15,
	2     y15,x16,y16,x17,y17,x18,y18,x19,y19,x20,y20)

	dimension x(20),y(20)

	write(30,*)'line: ',n,' points'
	if(n.lt.2.or.n.gt.20)
	2     call splot_err('1st argument of LINE is invalid')

	goto (2,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),n

20	x(20)=x20
	y(20)=y20
19	x(19)=x19
	y(19)=y19
18	x(18)=x18
	y(18)=y18
17	x(17)=x17
	y(17)=y17
16	x(16)=x16
	y(16)=y16
15	x(15)=x15
	y(15)=y15
14	x(14)=x14
	y(14)=y14
13	x(13)=x13
	y(13)=y13
12	x(12)=x12
	y(12)=y12
11	x(11)=x11
	y(11)=y11
10	x(10)=x10
	y(10)=y10
9	x(9)=x9
	y(9)=y9
8	x(8)=x8
	y(8)=y8
7	x(7)=x7
	y(7)=y7
6	x(6)=x6
	y(6)=y6
5	x(5)=x5
	y(5)=y5
4	x(4)=x4
	y(4)=y4
3	x(3)=x3
	y(3)=y3
2	x(2)=x2
	y(2)=y2
	x(1)=x1
	y(1)=y1

	do i=1,n
	  write(30,*)x(i),y(i)
	enddo
	write(30,*)

	write(31,*)'1 setgray'
	write(31,*)'newpath'
	call point(x(1),y(1),xp,yp)
	write(31,*)xp,yp,' moveto'
	do i=2,n
	  call point(x(i),y(i),xp,yp)
	  write(31,*)xp,yp,' lineto'
	enddo
	write(31,*)'stroke'
	write(31,*)'0 setgray'

	end

c closed line

	subroutine line_closed(n,
	2     x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,
	2     x8,y8,x9,y9,x10,y10,x11,y11,x12,y12,x13,y13,x14,y14,x15,
	2     y15,x16,y16,x17,y17,x18,y18,x19,y19,x20,y20)

	dimension x(20),y(20)

	write(30,*)'line: ',n,' points'
	if(n.lt.2.or.n.gt.20)
	2     call splot_err('1st argument of LINE is invalid')

	goto (2,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),n

20	x(20)=x20
	y(20)=y20
19	x(19)=x19
	y(19)=y19
18	x(18)=x18
	y(18)=y18
17	x(17)=x17
	y(17)=y17
16	x(16)=x16
	y(16)=y16
15	x(15)=x15
	y(15)=y15
14	x(14)=x14
	y(14)=y14
13	x(13)=x13
	y(13)=y13
12	x(12)=x12
	y(12)=y12
11	x(11)=x11
	y(11)=y11
10	x(10)=x10
	y(10)=y10
9	x(9)=x9
	y(9)=y9
8	x(8)=x8
	y(8)=y8
7	x(7)=x7
	y(7)=y7
6	x(6)=x6
	y(6)=y6
5	x(5)=x5
	y(5)=y5
4	x(4)=x4
	y(4)=y4
3	x(3)=x3
	y(3)=y3
2	x(2)=x2
	y(2)=y2
	x(1)=x1
	y(1)=y1

	do i=1,n
	  write(30,*)x(i),y(i)
	enddo
	write(30,*)

	write(31,*)'newpath'
	call point(x(1),y(1),xp,yp)
	call bbox_limits(xp,yp)
	write(31,*)xp,yp,' moveto'
	do i=2,n
	  call point(x(i),y(i),xp,yp)
	  call bbox_limits(xp,yp)
	  write(31,*)xp,yp,' lineto'
	enddo
	write(31,*)'closepath'
	write(31,*)'stroke'

	end

c closed line

	subroutine line_fill(n,
	2     x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,
	2     x8,y8,x9,y9,x10,y10,x11,y11,x12,y12,x13,y13,x14,y14,x15,
	2     y15,x16,y16,x17,y17,x18,y18,x19,y19,x20,y20)

	dimension x(20),y(20)

	write(30,*)'line: ',n,' points'
	if(n.lt.2.or.n.gt.20)
	2     call splot_err('1st argument of LINE is invalid')

	goto (2,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),n

20	x(20)=x20
	y(20)=y20
19	x(19)=x19
	y(19)=y19
18	x(18)=x18
	y(18)=y18
17	x(17)=x17
	y(17)=y17
16	x(16)=x16
	y(16)=y16
15	x(15)=x15
	y(15)=y15
14	x(14)=x14
	y(14)=y14
13	x(13)=x13
	y(13)=y13
12	x(12)=x12
	y(12)=y12
11	x(11)=x11
	y(11)=y11
10	x(10)=x10
	y(10)=y10
9	x(9)=x9
	y(9)=y9
8	x(8)=x8
	y(8)=y8
7	x(7)=x7
	y(7)=y7
6	x(6)=x6
	y(6)=y6
5	x(5)=x5
	y(5)=y5
4	x(4)=x4
	y(4)=y4
3	x(3)=x3
	y(3)=y3
2	x(2)=x2
	y(2)=y2
	x(1)=x1
	y(1)=y1

	do i=1,n
	  write(30,*)x(i),y(i)
	enddo
	write(30,*)

	write(31,*)'newpath'
	call point(x(1),y(1),xp,yp)
	call bbox_limits(xp,yp)
	write(31,*)xp,yp,' moveto'
	do i=2,n
	  call point(x(i),y(i),xp,yp)
	  call bbox_limits(xp,yp)
	  write(31,*)xp,yp,' lineto'
	enddo
	write(31,*)'closepath'
	write(31,*)'0.7 setgray'
	write(31,*)'fill'
	write(31,*)'0 setgray'

	end

c closed line

	subroutine line_fill2(n,gray,
	2     x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,
	2     x8,y8,x9,y9,x10,y10,x11,y11,x12,y12,x13,y13,x14,y14,x15,
	2     y15,x16,y16,x17,y17,x18,y18,x19,y19,x20,y20)

	dimension x(20),y(20)

	write(30,*)'line: ',n,' points'
	if(n.lt.2.or.n.gt.20)
	2     call splot_err('1st argument of LINE is invalid')

	goto (2,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),n

20	x(20)=x20
	y(20)=y20
19	x(19)=x19
	y(19)=y19
18	x(18)=x18
	y(18)=y18
17	x(17)=x17
	y(17)=y17
16	x(16)=x16
	y(16)=y16
15	x(15)=x15
	y(15)=y15
14	x(14)=x14
	y(14)=y14
13	x(13)=x13
	y(13)=y13
12	x(12)=x12
	y(12)=y12
11	x(11)=x11
	y(11)=y11
10	x(10)=x10
	y(10)=y10
9	x(9)=x9
	y(9)=y9
8	x(8)=x8
	y(8)=y8
7	x(7)=x7
	y(7)=y7
6	x(6)=x6
	y(6)=y6
5	x(5)=x5
	y(5)=y5
4	x(4)=x4
	y(4)=y4
3	x(3)=x3
	y(3)=y3
2	x(2)=x2
	y(2)=y2
	x(1)=x1
	y(1)=y1

	do i=1,n
	  write(30,*)x(i),y(i)
	enddo
	write(30,*)

	write(31,*)'newpath'
	call point(x(1),y(1),xp,yp)
	call bbox_limits(xp,yp)
	write(31,*)xp,yp,' moveto'
	do i=2,n
	  call point(x(i),y(i),xp,yp)
	  call bbox_limits(xp,yp)
	  write(31,*)xp,yp,' lineto'
	enddo
	write(31,*)'closepath'
	write(31,*)gray,' setgray'
	write(31,*)'fill'
	write(31,*)'0 setgray'

	end

c closed line

	subroutine line_fill3(n,
	2     x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,
	2     x8,y8,x9,y9,x10,y10,x11,y11,x12,y12,x13,y13,x14,y14,x15,
	2     y15,x16,y16,x17,y17,x18,y18,x19,y19,x20,y20)

	dimension x(20),y(20)

	write(30,*)'line: ',n,' points'
	if(n.lt.2.or.n.gt.20)
	2     call splot_err('1st argument of LINE is invalid')

	goto (2,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),n

20	x(20)=x20
	y(20)=y20
19	x(19)=x19
	y(19)=y19
18	x(18)=x18
	y(18)=y18
17	x(17)=x17
	y(17)=y17
16	x(16)=x16
	y(16)=y16
15	x(15)=x15
	y(15)=y15
14	x(14)=x14
	y(14)=y14
13	x(13)=x13
	y(13)=y13
12	x(12)=x12
	y(12)=y12
11	x(11)=x11
	y(11)=y11
10	x(10)=x10
	y(10)=y10
9	x(9)=x9
	y(9)=y9
8	x(8)=x8
	y(8)=y8
7	x(7)=x7
	y(7)=y7
6	x(6)=x6
	y(6)=y6
5	x(5)=x5
	y(5)=y5
4	x(4)=x4
	y(4)=y4
3	x(3)=x3
	y(3)=y3
2	x(2)=x2
	y(2)=y2
	x(1)=x1
	y(1)=y1

	do i=1,n
	  write(30,*)x(i),y(i)
	enddo
	write(30,*)

	write(31,*)'newpath'
	call point(x(1),y(1),xp,yp)
	call bbox_limits(xp,yp)
	write(31,*)xp,yp,' moveto'
	do i=2,n
	  call point(x(i),y(i),xp,yp)
	  call bbox_limits(xp,yp)
	  write(31,*)xp,yp,' lineto'
	enddo
	write(31,*)'closepath'
	write(31,*)'fill'

	end

c line using vectors

	subroutine line2(n,x,y)

	dimension x(100000),y(100000)

	write(30,*)'line: ',n,' points'
	if(n.lt.2.or.n.gt.100000)
	2     call splot_err('1st argument of LINE2 is invalid')

	do i=1,n
	  write(30,*)x(i),y(i)
	enddo
	write(30,*)

	write(31,*)'newpath'
	call point(x(1),y(1),xp,yp)
	call bbox_limits(xp,yp)
	write(31,*)xp,yp,' moveto'
	do i=2,n
	  call point(x(i),y(i),xp,yp)
	  call bbox_limits(xp,yp)
	  write(31,*)xp,yp,' lineto'
	enddo
	write(31,*)'stroke'

	return
	end

c line using vectors, closed path

	subroutine line2_closed(n,x,y)

	dimension x(100000),y(100000)

	write(30,*)'line: ',n,' points'
	if(n.lt.2.or.n.gt.100000)
	2     call splot_err('1st argument of LINE2 is invalid')

	do i=1,n
	  write(30,*)x(i),y(i)
	enddo
	write(30,*)

	write(31,*)'newpath'
	call point(x(1),y(1),xp,yp)
	call bbox_limits(xp,yp)
	write(31,*)xp,yp,' moveto'
	do i=2,n
	  call point(x(i),y(i),xp,yp)
	  call bbox_limits(xp,yp)
	  write(31,*)xp,yp,' lineto'
	enddo
	write(31,*)'closepath'
	write(31,*)'stroke'

	return
	end

c line using vectors, closed path filled with gray

	subroutine line2_fill(n,x,y)

	dimension x(100000),y(100000)

	write(30,*)'line: ',n,' points'
	if(n.lt.2.or.n.gt.100000)
	2     call splot_err('1st argument of LINE2 is invalid')

	do i=1,n
	  write(30,*)x(i),y(i)
	enddo
	write(30,*)

	write(31,*)'newpath'
	call point(x(1),y(1),xp,yp)
	call bbox_limits(xp,yp)
	write(31,*)xp,yp,' moveto'
	do i=2,n
	  call point(x(i),y(i),xp,yp)
	  call bbox_limits(xp,yp)
	  write(31,*)xp,yp,' lineto'
	enddo
	write(31,*)'closepath'
	write(31,*)'0.7 setgray'
	write(31,*)'fill'
	write(31,*)'0 setgray'

	return
	end

c line using vectors, closed path filled with gray

	subroutine line2_fill2(n,gray,x,y)

	dimension x(100000),y(100000)

	write(30,*)'line: ',n,' points'
	if(n.lt.2.or.n.gt.100000)
	2     call splot_err('1st argument of LINE2 is invalid')

	do i=1,n
	  write(30,*)x(i),y(i)
	enddo
	write(30,*)

	write(31,*)'newpath'
	call point(x(1),y(1),xp,yp)
	call bbox_limits(xp,yp)
	write(31,*)xp,yp,' moveto'
	do i=2,n
	  call point(x(i),y(i),xp,yp)
	  write(31,*)xp,yp,' lineto'
	  call bbox_limits(xp,yp)
	enddo
	write(31,*)'closepath'
	write(31,*)gray,' setgray'
	write(31,*)'fill'
	write(31,*)'0 setgray'

	return
	end

c dashed line

	subroutine dash(n,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,
	2     x8,y8,x9,y9,x10,y10,x11,y11,x12,y12,x13,y13,x14,y14,x15,
	2     y15,x16,y16,x17,y17,x18,y18,x19,y19,x20,y20)

	dimension x(20),y(20)

	write(30,*)'dashed line: ',n,' points'
	if(n.lt.2.or.n.gt.20)
	2     call splot_err('1st argument of DASH is invalid')

	goto (2,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),n

20	x(20)=x20
	y(20)=y20
19	x(19)=x19
	y(19)=y19
18	x(18)=x18
	y(18)=y18
17	x(17)=x17
	y(17)=y17
16	x(16)=x16
	y(16)=y16
15	x(15)=x15
	y(15)=y15
14	x(14)=x14
	y(14)=y14
13	x(13)=x13
	y(13)=y13
12	x(12)=x12
	y(12)=y12
11	x(11)=x11
	y(11)=y11
10	x(10)=x10
	y(10)=y10
9	x(9)=x9
	y(9)=y9
8	x(8)=x8
	y(8)=y8
7	x(7)=x7
	y(7)=y7
6	x(6)=x6
	y(6)=y6
5	x(5)=x5
	y(5)=y5
4	x(4)=x4
	y(4)=y4
3	x(3)=x3
	y(3)=y3
2	x(2)=x2
	y(2)=y2
	x(1)=x1
	y(1)=y1

	do i=1,n
	  write(30,*)x(i),y(i)
	enddo
	write(30,*)

	write(31,*)'gsave'
	write(31,*)'[5 5] 0 setdash'
	write(31,*)'newpath'
	call point(x(1),y(1),xp,yp)
	call bbox_limits(xp,yp)
	write(31,*)xp,yp,' moveto'
	do i=2,n
	  call point(x(i),y(i),xp,yp)
	  call bbox_limits(xp,yp)
	  write(31,*)xp,yp,' lineto'
	enddo
	write(31,*)'stroke'
	write(31,*)'grestore'

	return
	end

c dashed line with small dashes

	subroutine dash_size(n,size1,size2,
	2     x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,
	2     x8,y8,x9,y9,x10,y10,x11,y11,x12,y12,x13,y13,x14,y14,x15,
	2     y15,x16,y16,x17,y17,x18,y18,x19,y19,x20,y20)

	dimension x(20),y(20)

	write(30,*)'dashed line: ',n,' points'
	write(30,*)'length of dash: ',size1
	write(30,*)'distance between dashes: ',size2
	if(n.lt.2.or.n.gt.20)
	2     call splot_err('1st argument of DASH is invalid')

	goto (2,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),n

20	x(20)=x20
	y(20)=y20
19	x(19)=x19
	y(19)=y19
18	x(18)=x18
	y(18)=y18
17	x(17)=x17
	y(17)=y17
16	x(16)=x16
	y(16)=y16
15	x(15)=x15
	y(15)=y15
14	x(14)=x14
	y(14)=y14
13	x(13)=x13
	y(13)=y13
12	x(12)=x12
	y(12)=y12
11	x(11)=x11
	y(11)=y11
10	x(10)=x10
	y(10)=y10
9	x(9)=x9
	y(9)=y9
8	x(8)=x8
	y(8)=y8
7	x(7)=x7
	y(7)=y7
6	x(6)=x6
	y(6)=y6
5	x(5)=x5
	y(5)=y5
4	x(4)=x4
	y(4)=y4
3	x(3)=x3
	y(3)=y3
2	x(2)=x2
	y(2)=y2
	x(1)=x1
	y(1)=y1

	do i=1,n
	  write(30,*)x(i),y(i)
	enddo
	write(30,*)

	write(31,*)'gsave'
	write(31,*)'[ ',size1,size2,' ] 0 setdash'
	write(31,*)'newpath'
	call point(x(1),y(1),xp,yp)
	call bbox_limits(xp,yp)
	write(31,*)xp,yp,' moveto'
	do i=2,n
	  call point(x(i),y(i),xp,yp)
	  call bbox_limits(xp,yp)
	  write(31,*)xp,yp,' lineto'
	enddo
	write(31,*)'stroke'
	write(31,*)'grestore'

	return
	end

c dashed line using vectors

	subroutine dash2(n,x,y)

	dimension x(100000),y(100000)

	write(30,*)'dashed line: ',n,' points'
	if(n.lt.2.or.n.gt.100000)
	2     call splot_err('1st argument of DASH2 is invalid')

	do i=1,n
	  write(30,*)x(i),y(i)
	enddo
	write(30,*)

	write(31,*)'gsave'
	write(31,*)'[5 5] 0 setdash'
	write(31,*)'newpath'
	call point(x(1),y(1),xp,yp)
	call bbox_limits(xp,yp)
	write(31,*)xp,yp,' moveto'
	do i=2,n
	  call point(x(i),y(i),xp,yp)
	  call bbox_limits(xp,yp)
	  write(31,*)xp,yp,' lineto'
	enddo
	write(31,*)'stroke'
	write(31,*)'grestore'

	end

c dashed line using vectors

	subroutine dash2_size(n,size1,size2,x,y)

	dimension x(100000),y(100000)

	write(30,*)'dashed line: ',n,' points'
	write(30,*)'length of dash: ',size1
	write(30,*)'distance between dashes: ',size2
	if(n.lt.2.or.n.gt.100000)
	2     call splot_err('1st argument of DASH2 is invalid')

	do i=1,n
	  write(30,*)x(i),y(i)
	enddo
	write(30,*)

	write(31,*)'gsave'
	write(31,*)'[ ',size1,size2,' ] 0 setdash'
	write(31,*)'newpath'
	call point(x(1),y(1),xp,yp)
	call bbox_limits(xp,yp)
	write(31,*)xp,yp,' moveto'
	do i=2,n
	  call point(x(i),y(i),xp,yp)
	  call bbox_limits(xp,yp)
	  write(31,*)xp,yp,' lineto'
	enddo
	write(31,*)'stroke'
	write(31,*)'grestore'

	end

c arrow

	subroutine arrow(xcenter,ycenter,xpoint,ypoint,arrowlength)

	common scale

	write(30,*)'arrow'
	write(30,*)'point at ',xcenter,ycenter
	write(30,*)'direction: ',xpoint,ypoint
	write(30,*)'length: ',arrowlength
	write(30,*)

	call point(xcenter,ycenter,x,y)
	AMAG=SQRT(xpoint**2+ypoint**2)
	if(amag.eq.0.)call splot_err('Invalid direction in ARROW')
	xpointunit=xpoint/AMAG
	ypointunit=ypoint/AMAG
	xpointperp=ypointunit
	ypointperp=-xpointunit

	write(31,*)'newpath'
	write(31,*)x,y,' moveto'
	write(31,*)x-scale*arrowlength*(xpointunit+xpointperp/8),
	2     y-scale*arrowlength*(ypointunit+ypointperp/8),
	2     ' lineto'
	write(31,*)x-scale*arrowlength*(xpointunit-xpointperp/8),
	2     y-scale*arrowlength*(ypointunit-ypointperp/8),
	2     ' lineto'
        write(31,*)'closepath'
	write(31,*)'fill'

	end

c arrow

	subroutine widearrow(xcenter,ycenter,xpoint,ypoint,arrowlength)

	common scale,origin

	write(30,*)'wide arrow'
	write(30,*)'point at ',xcenter,ycenter
	write(30,*)'direction: ',xpoint,ypoint
	write(30,*)'length: ',arrowlength
	write(30,*)

	call point(xcenter,ycenter,x,y)
	AMAG=SQRT(xpoint**2+ypoint**2)
	if(amag.eq.0.)call splot_err('Invalid direction in ARROW')
	xpointunit=xpoint/AMAG
	ypointunit=ypoint/AMAG
	xpointperp=ypointunit
	ypointperp=-xpointunit

	write(31,*)'newpath'
	write(31,*)x,y,' moveto'
	write(31,*)x-scale*arrowlength*(xpointunit+xpointperp/4),
	2     y-scale*arrowlength*(ypointunit+ypointperp/4),
	2     ' lineto'
	write(31,*)x-scale*arrowlength*(xpointunit-xpointperp/4),
	2     y-scale*arrowlength*(ypointunit-ypointperp/4),
	2     ' lineto'
        write(31,*)'closepath'
	write(31,*)'fill'

	return
	end

c arrow

	subroutine arrow2(xcenter,ycenter,xpoint,ypoint,
	1    arrowlength,arrowwidth)

	common scale,origin

	write(30,*)'arrow2'
	write(30,*)'point at ',xcenter,ycenter
	write(30,*)'direction: ',xpoint,ypoint
	write(30,*)'length: ',arrowlength
	write(30,*)'width: ',arrowwidth
	write(30,*)

	call point(xcenter,ycenter,x,y)
	AMAG=SQRT(xpoint**2+ypoint**2)
	if(amag.eq.0.)call splot_err('Invalid direction in ARROW')
	xpointunit=xpoint/AMAG
	ypointunit=ypoint/AMAG
	xpointperp=ypointunit
	ypointperp=-xpointunit

	write(31,*)'newpath'
	write(31,*)x,y,' moveto'
	write(31,*)x-scale*(arrowlength*xpointunit
	1    +arrowwidth/2*xpointperp),
	2     y-scale*(arrowlength*ypointunit
	2    +arrowwidth/2*ypointperp),
	2     ' lineto'
	write(31,*)x-scale*(arrowlength*xpointunit
	3    -arrowwidth/2*xpointperp),
	2     y-scale*(arrowlength*ypointunit
	4    -arrowwidth/2*ypointperp),
	2     ' lineto'
        write(31,*)'closepath'
	write(31,*)'fill'

	return
	end

c bullet

	subroutine bullet(xcenter,ycenter,radius)

	common scale

	write(30,*)'bullet'
	write(30,*)'center at ',xcenter,ycenter
	write(30,*)'radius: ',radius
	write(30,*)

	if(radius.le.0.)call splot_err('Invalid radius in BULLET')

	call point(xcenter,ycenter,x,y)
	write(31,*)'newpath'
	write(31,*)x,y,radius*scale,' 0 360 arc fill'

	end

c pie

	subroutine pie(xcenter,ycenter,radius,anglebegin,angleend)

	common scale

	write(30,*)'pie'
	write(30,*)'center at ',xcenter,ycenter
	write(30,*)'begin at angle: ',anglebegin
	write(30,*)'end at angle:   ',angleend
	write(30,*)'radius: ',radius
	write(30,*)

	if(radius.le.0.)call splot_err('Invalid radius in BULLET')

	call point(xcenter,ycenter,x,y)
	write(31,*)'newpath'
	write(31,*)x,y,radius*scale,anglebegin,angleend,' arc'
	write(31,*)x,y,' lineto'
	write(31,*)'closepath'
	write(31,*)'fill'

	end

c circle

	subroutine circle(xcenter,ycenter,radius)

	common scale

	write(30,*)'circle:'
	write(30,*)'center at ',xcenter,ycenter
	write(30,*)'radius: ',radius
	write(30,*)

	if(radius.le.0.)call splot_err('Invalid radius in CIRCLE')

	call point(xcenter,ycenter,x,y)
	write(31,*)'newpath'
	write(31,*)x,y,radius*scale,' 0 360 arc stroke'
	rp=radius*scale
	call bbox_limits(x+rp,y)
	call bbox_limits(x-rp,y)
	call bbox_limits(x,y+rp)
	call bbox_limits(x,y-rp)

	end

c circle

	subroutine circle_fill(xcenter,ycenter,radius,gray)

	common scale

	write(30,*)'filled circle:'
	write(30,*)'center at ',xcenter,ycenter
	write(30,*)'radius: ',radius
	write(30,*)

	if(radius.le.0.)call splot_err('Invalid radius in CIRCLE')

	call point(xcenter,ycenter,x,y)
	write(31,*)'newpath'
	write(31,*)x,y,radius*scale,' 0 360 arc'
	write(31,*)gray,' setgray'
	write(31,*)'fill'
	write(31,*)'0 setgray'
	rp=radius*scale
	call bbox_limits(x+rp,y)
	call bbox_limits(x-rp,y)
	call bbox_limits(x,y+rp)
	call bbox_limits(x,y-rp)

	end

c arc

	subroutine arc(xcenter,ycenter,radius,anglebegin,angleend)

	common scale

	write(30,*)'arc:'
	write(30,*)'center at ',xcenter,ycenter
	write(30,*)'radius: ',radius
	write(30,*)'begin at angle: ',anglebegin
	write(30,*)'end at angle:   ',angleend
	write(30,*)

	if(radius.le.0.)call splot_err('Invalid radius in ARC')
	if(angleend.lt.anglebegin)call splot_err('Invalid angle in ARC')
	if(angleend-anglebegin.gt.360.)
	2     call splot_err('Invalid angle in ARC')

	call point(xcenter,ycenter,x,y)
	write(31,*)'newpath'
	write(31,*)x,y,radius*scale,anglebegin,angleend,' arc stroke'
	rp=radius*scale
	call bbox_limits(x+rp*cosd(anglebegin),y+rp*sind(anglebegin))
	call bbox_limits(x+rp*cosd(angleend),y+rp*sind(angleend))
	do i=-360,360,90
	  angle=i
	  if(angle.gt.anglebegin.and.angle.lt.angleend)
	1      call bbox_limits(x+rp*cosd(angle),y+rp*sind(angle))
	enddo
	end

c dashed arc

	subroutine arc_dash(xcenter,ycenter,radius,anglebegin,angleend)

	common scale

	write(30,*)'dashed arc:'
	write(30,*)'center at ',xcenter,ycenter
	write(30,*)'radius: ',radius
	write(30,*)'begin at angle: ',anglebegin
	write(30,*)'end at angle:   ',angleend
	write(30,*)

	if(radius.le.0.)call splot_err('Invalid radius in ARC')
	if(angleend.lt.anglebegin)call splot_err('Invalid angle in ARC')
	if(angleend-anglebegin.gt.360.)
	2     call splot_err('Invalid angle in ARC')

	call point(xcenter,ycenter,x,y)
	write(31,*)'gsave'
	write(31,*)'[5 5] 0 setdash'
	write(31,*)'newpath'
	write(31,*)x,y,radius*scale,anglebegin,angleend,' arc stroke'
	write(31,*)'grestore'

	end

c triangle, pointed side up

	subroutine triangleup(xcenter,ycenter,radius)

	common scale

	write(30,*)'triangle:'
	write(30,*)'center at ',xcenter,ycenter
	write(30,*)'radius: ',radius
	write(30,*)

	if(radius.le.0.)call splot_err('Invalid radius in CIRCLE')

	a=radius/0.643
	write(31,*)'newpath'
	call point(xcenter,ycenter+a,x,y)
	write(31,*)x,y,' moveto'
	call point(xcenter+a*sqrt(3.0)/2,ycenter-a/2,x,y)
	write(31,*)x,y,' lineto'
	call point(xcenter-a*sqrt(3.0)/2,ycenter-a/2,x,y)
	write(31,*)x,y,' lineto'
	write(31,*)'closepath'
	write(31,*)'fill'

	end

c open triangle, pointed side up

	subroutine triangleupopen(xcenter,ycenter,radius)

	common scale

	write(30,*)'triangle:'
	write(30,*)'center at ',xcenter,ycenter
	write(30,*)'radius: ',radius
	write(30,*)

	if(radius.le.0.)call splot_err('Invalid radius in CIRCLE')

	a=radius/0.643
	write(31,*)'newpath'
	call point(xcenter,ycenter+a,x,y)
	write(31,*)x,y,' moveto'
	call point(xcenter+a*sqrt(3.0)/2,ycenter-a/2,x,y)
	write(31,*)x,y,' lineto'
	call point(xcenter-a*sqrt(3.0)/2,ycenter-a/2,x,y)
	write(31,*)x,y,' lineto'
	write(31,*)'closepath'
	write(31,*)'stroke'

	end

* box with tick marks

	subroutine box(xmin,xmax,ymin,ymax,nx,ny,ticksize)

	write(30,*)'box'
	write(30,*)'minimum value of x: ',xmin
	write(30,*)'maximum value of x: ',xmax
	write(30,*)'minimum value of y: ',ymin
	write(30,*)'maximum value of y: ',ymax
	write(30,*)'number of divisions on x axis: ',nx
	write(30,*)'number of divisions on y axis: ',ny
	write(30,*)'length of tick mark: ',ticksize

* box
	write(31,*)'newpath'
	call point(xmin,ymin,x,y)
	write(31,*)x,y,' moveto'
	call point(xmax,ymin,x,y)
	write(31,*)x,y,' lineto'
	call point(xmax,ymax,x,y)
	write(31,*)x,y,' lineto'
	call point(xmin,ymax,x,y)
	write(31,*)x,y,' lineto'
	write(31,*)'closepath'
	write(31,*)'gsave'
	write(31,*)'stroke'
* clip everything outside box
	write(31,*)'grestore'
	write(31,*)'clip'
* ticks
	if(nx.gt.1)then
	  dx=(xmax-xmin)/nx
	  a=xmin
	  do i=1,nx-1
	    a=a+dx
	    write(31,*)'newpath'
	    call point(a,ymin,x,y)
	    write(31,*)x,y,' moveto'
	    call point(a,ymin+ticksize,x,y)
	    write(31,*)x,y,' lineto'
	    write(31,*)'stroke'
	    write(31,*)'newpath'
	    call point(a,ymax,x,y)
	    write(31,*)x,y,' moveto'
	    call point(a,ymax-ticksize,x,y)
	    write(31,*)x,y,' lineto'
	    write(31,*)'stroke'
	  enddo
	endif
	if(ny.gt.1)then
	  dy=(ymax-ymin)/ny
	  a=ymin
	  do i=1,ny-1
	    a=a+dy
	    write(31,*)'newpath'
	    call point(xmin,a,x,y)
	    write(31,*)x,y,' moveto'
	    call point(xmin+ticksize,a,x,y)
	    write(31,*)x,y,' lineto'
	    write(31,*)'stroke'
	    write(31,*)'newpath'
	    call point(xmax,a,x,y)
	    write(31,*)x,y,' moveto'
	    call point(xmax-ticksize,a,x,y)
	    write(31,*)x,y,' lineto'
	    write(31,*)'stroke'
	  enddo
	endif

	end

* box containing data points

	subroutine box_points(xmin,xmax,ymin,ymax,n,xdata,ydata)
	
	dimension xdata(100000),ydata(100000)
	common/common_getx/scalex,offsetx
	common/common_gety/scaley,offsety

	write(30,*)'box containing data points'
	write(30,*)'minimum value of x: ',xmin
	write(30,*)'maximum value of x: ',xmax
	write(30,*)'minimum value of y: ',ymin
	write(30,*)'maximum value of y: ',ymax
	write(30,*)n,' data points:'
	
* find minimum and maximum values of data
	x1=xdata(1)
	x2=xdata(1)
	y1=ydata(1)
	y2=ydata(1)
	do i=2,n
	  if(xdata(i).lt.x1)x1=xdata(i)
	  if(xdata(i).gt.x2)x2=xdata(i)
	  if(ydata(i).lt.y1)y1=ydata(i)
	  if(ydata(i).gt.y2)y2=ydata(i)
	enddo
* find the scale and offset of the x coordinate
	dx=(x2-x1)/20
	x1=x1-dx
	x2=x2+dx
	scalex=(xmax-xmin)/(x2-x1)
	offsetx=xmin/scalex-x1
* find the scale and offset of the y coordinate
	dy=(y2-y1)/20
	y1=y1-dy
	y2=y2+dy
	scaley=(ymax-ymin)/(y2-y1)
	offsety=ymin/scaley-y1

* box
	write(31,*)'newpath'
	call point(xmin,ymin,x,y)
	write(31,*)x,y,' moveto'
	call point(xmax,ymin,x,y)
	write(31,*)x,y,' lineto'
	call point(xmax,ymax,x,y)
	write(31,*)x,y,' lineto'
	call point(xmin,ymax,x,y)
	write(31,*)x,y,' lineto'
	write(31,*)'closepath'
	write(31,*)'gsave'
	write(31,*)'stroke'
* clip everything outside box
	write(31,*)'grestore'
	write(31,*)'clip'

	end

* clip figure

	subroutine clip(xmin,xmax,ymin,ymax)

	write(30,*)'clip'
	write(30,*)'minimum value of x: ',xmin
	write(30,*)'maximum value of x: ',xmax
	write(30,*)'minimum value of y: ',ymin
	write(30,*)'maximum value of y: ',ymax

* box defining clip
	write(31,*)'newpath'
	call point(xmin,ymin,x,y)
	write(31,*)x,y,' moveto'
	call point(xmax,ymin,x,y)
	write(31,*)x,y,' lineto'
	call point(xmax,ymax,x,y)
	write(31,*)x,y,' lineto'
	call point(xmin,ymax,x,y)
	write(31,*)x,y,' lineto'
	write(31,*)'closepath'
	write(31,*)'clip'

	end

* scale points

	function getx(x)
	common/common_getx/scalex,offsetx
	getx=(x+offsetx)*scalex
	end

	function gety(y)
	common/common_gety/scaley,offsety
	gety=(y+offsety)*scaley
	end

c pen width

	subroutine pen(x)

	write(30,*)'pen: ',x
	write(30,*)

	widthline=x*0.24
	write(31,*)widthline,' setlinewidth'

	return
	end

c text

	subroutine text(xtext,ytext,word,n)

	character*(*) word

	write(30,*)'text: ',word
	write(30,*)'x,y: ',xtext,ytext
	write(30,*)'height: ',n,' points'
	write(30,*)
	call point(xtext,ytext,x,y)

	write(31,*)'/Helvetica findfont'
	write(31,*)n,' scalefont'
	write(31,*)'setfont'
	write(31,*)x,y,' moveto'
	write(31,*)'(',word,') show'

 	return
	end


c textup

	subroutine textup(xtext,ytext,word,n)

	character*(*) word

	write(30,*)'textup: ',word
	write(30,*)'x,y: ',xtext,ytext
	write(30,*)'height: ',n,' points'
	write(30,*)
	call point(xtext,ytext,x,y)

	write(31,*)'/Helvetica findfont'
	write(31,*)n,' scalefont'
	write(31,*)'setfont'
	write(31,*)'gsave'
	write(31,*)x,y,' moveto'
	write(31,*)'90 rotate'
	write(31,*)'(',word,') show'
	write(31,*)'grestore'

	return
	end


c text_rotate

	subroutine text_rotate(xtext,ytext,angle,word,n)

	character*(*) word

	write(30,*)'text_rotate: ',word
	write(30,*)'x,y: ',xtext,ytext
	write(30,*)'angle: ',angle
	write(30,*)'height: ',n,' points'
	write(30,*)
	call point(xtext,ytext,x,y)

	write(31,*)'/Helvetica findfont'
	write(31,*)n,' scalefont'
	write(31,*)'setfont'
	write(31,*)'gsave'
	write(31,*)x,y,' moveto'
	write(31,*)angle,' rotate'
	write(31,*)'(',word,') show'
	write(31,*)'grestore'

	return
	end


c text with different font

	subroutine textfont(xtext,ytext,word,n,font)

	character*(*) word,font

	write(30,*)'text: ',word
	write(30,*)'font: ',font
	write(30,*)'x,y: ',xtext,ytext
	write(30,*)'height: ',n,' points'
	write(30,*)
	call point(xtext,ytext,x,y)

	write(31,*)'/',font,' findfont'
	write(31,*)n,' scalefont'
	write(31,*)'setfont'
	write(31,*)x,y,' moveto'
	write(31,*)'(',word,') show'

 	return
	end


c resistor

	subroutine resistor(x1,y1,x2,y2)

	dimension dirup(2),dirdown(2),x(8),y(8)

	write(30,*)'resistor:'
	write(30,*)x1,y1
	write(30,*)x2,y2
	write(30,*)

	segmentlength=sqrt((y2-y1)**2+(x2-x1)**2)/6
	angle=atan2d(y2-y1,x2-x1)
	angleup=angle+60
	angledown=angle-60
	dirup(1)=segmentlength*cosd(angleup)
	dirup(2)=segmentlength*sind(angleup)
	dirdown(1)=segmentlength*cosd(angledown)
	dirdown(2)=segmentlength*sind(angledown)
	x(1)=x1
	y(1)=y1
	x(2)=x(1)+dirup(1)
	y(2)=y(1)+dirup(2)
	x(3)=x(2)+2*dirdown(1)
	y(3)=y(2)+2*dirdown(2)
	x(4)=x(3)+2*dirup(1)
	y(4)=y(3)+2*dirup(2)
	x(5)=x(4)+2*dirdown(1)
	y(5)=y(4)+2*dirdown(2)
	x(6)=x(5)+2*dirup(1)
	y(6)=y(5)+2*dirup(2)
	x(7)=x(6)+2*dirdown(1)
	y(7)=y(6)+2*dirdown(2)
	x(8)=x2
	y(8)=y2
	call line2(8,x,y)

	return
	end

c ground

	subroutine ground(xcenter,ycenter,groundwidth)

	dimension x(2),y(2)
	logical init,listfile
*	external pl2ca
	common/splot_common/xmin,xmax,ymin,ymax,init,listfile

	write(30,*)'ground:'
	write(30,*)'center: ',xcenter,ycenter
	write(30,*)'width:  ',groundwidth
	write(30,*)

	n=2

	x(1)=xcenter-groundwidth/2
	y(1)=ycenter
	x(2)=xcenter+groundwidth/2
	y(2)=ycenter
	call line2(n,x,y)

	x(1)=xcenter-groundwidth/3
	y(1)=ycenter-groundwidth/4
	x(2)=xcenter+groundwidth/3
	y(2)=ycenter-groundwidth/4
	call line2(n,x,y)

	x(1)=xcenter-groundwidth/10
	y(1)=ycenter-groundwidth/2
	x(2)=xcenter+groundwidth/10
	y(2)=ycenter-groundwidth/2
	call line2(n,x,y)

	return
	end

c battery

	subroutine battery(x1,y1,x2,y2)

	dimension center(2),x(2),y(2)
	logical init,listfile
*	external pl2ca
	common/splot_common/xmin,xmax,ymin,ymax,init,listfile

	write(30,*)'battery:'
	write(30,*)x1,y1
	write(30,*)x2,y2
	write(30,*)

	center(1)=(x1+x2)/2
	center(2)=(y1+y2)/2
	batterylength=sqrt((y2-y1)**2+(x2-x1)**2)
	cornerlength1=batterylength/2/cosd(79.)
	cornerlength2=batterylength/2/cosd(63.)
	angle=atan2d(y1-y2,x1-x2)

	x(1)=center(1)+cornerlength2*cosd(angle+63)
	y(1)=center(2)+cornerlength2*sind(angle+63)
	x(2)=center(1)+cornerlength2*cosd(angle-63)
	y(2)=center(2)+cornerlength2*sind(angle-63)
	n=2
	call line2(n,x,y)

	x(1)=center(1)-cornerlength1*cosd(angle+79)
	y(1)=center(2)-cornerlength1*sind(angle+79)
	x(2)=center(1)-cornerlength1*cosd(angle-79)
	y(2)=center(2)-cornerlength1*sind(angle-79)
	n=2
	call line2(n,x,y)

	return
	end

* ac source

	subroutine acsource(x,y,r)

	dimension x1(0:40),y1(0:40)

	write(30,*)'ac source'
	write(30,*)x,y
	write(30,*)'radius = ',r

	call circle(x,y,r)
	do i=0,40
	  x1(i)=x+((i-20)*0.0375)*r
	  y1(i)=y+r/2*sind(9.0*i)
	enddo
	call line2(41,x1,y1)
	end

c capacitor

	subroutine capacitor(x1,y1,x2,y2)

	dimension center(2),x(2),y(2)
	logical init,listfile
*	external pl2ca
	common/splot_common/xmin,xmax,ymin,ymax,init,listfile

	write(30,*)'capacitor:'
	write(30,*)x1,y1
	write(30,*)x2,y2
	write(30,*)

	center(1)=(x1+x2)/2
	center(2)=(y1+y2)/2
	capacitorlength=sqrt((y2-y1)**2+(x2-x1)**2)
	cornerlength=capacitorlength/2/cosd(76.)
	angle=atan2d(y1-y2,x1-x2)

	x(1)=center(1)+cornerlength*cosd(angle+76)
	y(1)=center(2)+cornerlength*sind(angle+76)
	x(2)=center(1)+cornerlength*cosd(angle-76)
	y(2)=center(2)+cornerlength*sind(angle-76)
	n=2
	call line2(n,x,y)

	x(1)=center(1)-cornerlength*cosd(angle+76)
	y(1)=center(2)-cornerlength*sind(angle+76)
	x(2)=center(1)-cornerlength*cosd(angle-76)
	y(2)=center(2)-cornerlength*sind(angle-76)
	n=2
	call line2(n,x,y)

	return
	end

c diode

	subroutine diode(x1,y1,x2,y2)

	dimension uppoint(2),downpoint(2),x(8),y(8)
	logical init,listfile
*	external pl2ca
	common/splot_common/xmin,xmax,ymin,ymax,init,listfile

	write(30,*)'diode:'
	write(30,*)x1,y1
	write(30,*)x2,y2
	write(30,*)

	diodelength=sqrt((y2-y1)**2+(x2-x1)**2)
	sidelength=diodelength/cosd(30.)
	angle=atan2d(y1-y2,x1-x2)
	uppoint(1)=x2+sidelength*cosd(angle+30)
	uppoint(2)=y2+sidelength*sind(angle+30)
	downpoint(1)=x2+sidelength*cosd(angle-30)
	downpoint(2)=y2+sidelength*sind(angle-30)

	x(1)=x2
	y(1)=y2
	x(2)=uppoint(1)
	y(2)=uppoint(2)
	x(3)=downpoint(1)
	y(3)=downpoint(2)
	x(4)=x2
	y(4)=y2
	n=4
	call line2(n,x,y)

	do i=1,15
	  fdown=i/16.
	  fup=1-fdown
	  x(2)=fup*uppoint(1)+fdown*downpoint(1)
	  y(2)=fup*uppoint(2)+fdown*downpoint(2)
	  n=2
	  call line2(n,x,y)
	enddo

	x(1)=uppoint(1)+x2-x1
	y(1)=uppoint(2)+y2-y1
	x(2)=downpoint(1)+x2-x1
	y(2)=downpoint(2)+y2-y1
	n=2
	call line2(n,x,y)

	return
	end

c zener diode

	subroutine zener(x1,y1,x2,y2)

	dimension uppoint(2),downpoint(2),x(8),y(8)
	logical init,listfile
*	external pl2ca
	common/splot_common/xmin,xmax,ymin,ymax,init,listfile

	write(30,*)'diode:'
	write(30,*)x1,y1
	write(30,*)x2,y2
	write(30,*)

	diodelength=sqrt((y2-y1)**2+(x2-x1)**2)
	sidelength=diodelength/cosd(30.)
	angle=atan2d(y1-y2,x1-x2)
	uppoint(1)=x2+sidelength*cosd(angle+30)
	uppoint(2)=y2+sidelength*sind(angle+30)
	downpoint(1)=x2+sidelength*cosd(angle-30)
	downpoint(2)=y2+sidelength*sind(angle-30)

	x(1)=x2
	y(1)=y2
	x(2)=uppoint(1)
	y(2)=uppoint(2)
	x(3)=downpoint(1)
	y(3)=downpoint(2)
	x(4)=x2
	y(4)=y2
	n=4
	call line2(n,x,y)

	do i=1,15
	  fdown=i/16.
	  fup=1-fdown
	  x(2)=fup*uppoint(1)+fdown*downpoint(1)
	  y(2)=fup*uppoint(2)+fdown*downpoint(2)
	  n=2
	  call line2(n,x,y)
	enddo

	zlength=diodelength/2.5
	x(1)=uppoint(1)+x2-x1+zlength*cosd(angle+150)
	y(1)=uppoint(2)+y2-y1+zlength*sind(angle+150)
	x(2)=uppoint(1)+x2-x1
	y(2)=uppoint(2)+y2-y1
	x(3)=downpoint(1)+x2-x1
	y(3)=downpoint(2)+y2-y1
	x(4)=downpoint(1)+x2-x1+zlength*cosd(angle-30)
	y(4)=downpoint(2)+y2-y1+zlength*sind(angle-30)
	n=4
	call line2(n,x,y)

	return
	end

c inductor

	subroutine inductor(x1,y1,x2,y2)

	dimension x(0:70),y(0:70)
	logical init,listfile
*	external pl2ca
	common/splot_common/xmin,xmax,ymin,ymax,init,listfile

	write(30,*)'inductor:'
	write(30,*)x1,y1
	write(30,*)'height: ',height
	write(30,*)

	n=70
	d=sqrt((x2-x1)**2+(y2-y1)**2)
	ry=d/6
	rx=d/10
	angle=atan2d(y2-y1,x2-x1)
	a=-180
	da=1260.0/n
	dd=(d-2*rx)/n
	d0=rx
	do i=0,n
	  x3=d0+rx*cosd(a)
	  y3=ry*sind(a)
	  x(i)=x1+x3*cosd(angle)-y3*sind(angle)
	  y(i)=y1+x3*sind(angle)+y3*cosd(angle)
	  a=a+da
	  d0=d0+dd
	enddo
	call line2(n+1,x,y)

	return
	end

c transformer

	subroutine transformer(x1,y1,height)

	dimension circlepoints(0:20,2),x(50),y(50)
	logical init,listfile
*	external pl2ca
	common/splot_common/xmin,xmax,ymin,ymax,init,listfile

	write(30,*)'transformer:'
	write(30,*)x1,y1
	write(30,*)'height: ',height
	write(30,*)

	scale=height/9
	r=1.125*scale
	do i=0,20
	  angle=90-18*i
	  circlepoints(i,1)=r*cosd(angle)
	  circlepoints(i,2)=r*sind(angle)
	enddo

	x(1)=x1
	y(1)=y1
	n=1
	do i=1,4
	  xcenter=x1
	  ycenter=y1+r-2*i*r
	  do j=1,10
	    n=n+1
	    x(n)=xcenter+circlepoints(j,1)
	    y(n)=ycenter+circlepoints(j,2)
	  enddo
	enddo
	call line2(n,x,y)

	x(1)=x1+7*scale
	y(1)=y1
	n=1
	do i=1,4
	  xcenter=x1+7*scale
	  ycenter=y1+r-2*i*r
	  do j=19,10,-1
	    n=n+1
	    x(n)=xcenter+circlepoints(j,1)
	    y(n)=ycenter+circlepoints(j,2)
	  enddo
	enddo
	call line2(n,x,y)

	x(1)=x1+3
	y(1)=y1
	x(2)=x1+3
	y(2)=y1-height
	n=2
	call line2(n,x,y)

	x(1)=x1+4
	y(1)=y1
	x(2)=x1+4
	y(2)=y1-height
	n=2
	call line2(n,x,y)

	return
	end

c npn transistor

	subroutine npn(x1,y1,x2,y2,nturn)

	dimension x(2),y(2)
	logical init,listfile
*	external pl2ca
	common/splot_common/xmin,xmax,ymin,ymax,init,listfile

	write(30,*)'npn transistor:'
	write(30,*)x1,y1
	write(30,*)x2,y2
	if(nturn.gt.0)write(30,*)'arrow in positive direction'
	if(nturn.lt.0)write(30,*)'arrow in negative direction'
	write(30,*)

	twidth=sqrt((y2-y1)**2+(x2-x1)**2)
	theight=twidth*8/4.5
	tbar=twidth*3/4.5
	tbarwidth=twidth/9
	tdiag=theight/2-twidth*tand(30.)
	diaglength=twidth/cosd(30.)
	angle=atan2d(y2-y1,x2-x1)

	do i=0,2
	  x(1)=x1+(x2-x1)/twidth*tbarwidth*i/2+tbar*cosd(angle+90)
	  y(1)=y1+(y2-y1)/twidth*tbarwidth*i/2+tbar*sind(angle+90)
	  x(2)=x1+(x2-x1)/twidth*tbarwidth*i/2+tbar*cosd(angle-90)
	  y(2)=y1+(y2-y1)/twidth*tbarwidth*i/2+tbar*sind(angle-90)
	  n=2
	  call line2(n,x,y)
	enddo

	x(1)=x1+tdiag*cosd(angle+90)
	y(1)=y1+tdiag*sind(angle+90)
	x(2)=x(1)+diaglength*cosd(angle+30)
	y(2)=y(1)+diaglength*sind(angle+30)
	n=2
	call line2(n,x,y)

	if(nturn.gt.0)then
	  listfile=.false.
	  call widearrow(x(2),y(2),x(2)-x(1),y(2)-y(1),0.4*twidth)
	  listfile=.true.
	endif

	x(1)=x1+tdiag*cosd(angle-90)
	y(1)=y1+tdiag*sind(angle-90)
	x(2)=x(1)+diaglength*cosd(angle-30)
	y(2)=y(1)+diaglength*sind(angle-30)
	n=2
	call line2(n,x,y)

	if(nturn.lt.0)then
	  listfile=.false.
	  call widearrow(x(2),y(2),x(2)-x(1),y(2)-y(1),0.4*twidth)
	  listfile=.true.
	endif

	return
	end

c pnp transistor

	subroutine pnp(x1,y1,x2,y2,nturn)

	dimension x(2),y(2)
	logical init,listfile
*	external pl2ca
	common/splot_common/xmin,xmax,ymin,ymax,init,listfile

	write(30,*)'npn transistor:'
	write(30,*)x1,y1
	write(30,*)x2,y2
	if(nturn.gt.0)write(30,*)'arrow in positive direction'
	if(nturn.lt.0)write(30,*)'arrow in negative direction'
	write(30,*)

	twidth=sqrt((y2-y1)**2+(x2-x1)**2)
	theight=twidth*8/4.5
	tbar=twidth*3/4.5
	tbarwidth=twidth/9
	tdiag=theight/2-twidth*tand(30.)
	diaglength=twidth/cosd(30.)
	angle=atan2d(y2-y1,x2-x1)

	do i=0,2
	  x(1)=x1+(x2-x1)/twidth*tbarwidth*i/2+tbar*cosd(angle+90)
	  y(1)=y1+(y2-y1)/twidth*tbarwidth*i/2+tbar*sind(angle+90)
	  x(2)=x1+(x2-x1)/twidth*tbarwidth*i/2+tbar*cosd(angle-90)
	  y(2)=y1+(y2-y1)/twidth*tbarwidth*i/2+tbar*sind(angle-90)
	  n=2
	  call line2(n,x,y)
	enddo

	x(1)=x1+tdiag*cosd(angle+90)
	y(1)=y1+tdiag*sind(angle+90)
	x(2)=x(1)+diaglength*cosd(angle+30)
	y(2)=y(1)+diaglength*sind(angle+30)
	n=2
	call line2(n,x,y)

	if(nturn.gt.0)then
	  listfile=.false.
	  call widearrow(x(1)+0.2*diaglength*cosd(angle+30),
	2                y(1)+0.2*diaglength*sind(angle+30),
	2                x(1)-x(2),y(1)-y(2),0.4*twidth)
	  listfile=.true.
	endif

	x(1)=x1+tdiag*cosd(angle-90)
	y(1)=y1+tdiag*sind(angle-90)
	x(2)=x(1)+diaglength*cosd(angle-30)
	y(2)=y(1)+diaglength*sind(angle-30)
	n=2
	call line2(n,x,y)

	if(nturn.lt.0)then
	  listfile=.false.
	  call widearrow(x(1)+0.2*diaglength*cosd(angle+30),
	2                y(1)+0.2*diaglength*sind(angle+30),
	2                x(1)-x(2),y(1)-y(2),0.4*twidth)
	  listfile=.true.
	endif

	return
	end

* MOSFET

	subroutine mosfet(x1,y1,xmosfetwidth)

	dimension x(4),y(4)
	logical init,listfile
*	external pl2ca
	common/splot_common/xmin,xmax,ymin,ymax,init,listfile

	write(30,*)'MOSFET:'
	write(30,*)x1,y1
	write(30,*)'width: ',xmosfetwidth
	write(30,*)

	scale=xmosfetwidth/5

	x(1)=x1
	y(1)=y1
	x(2)=x1
	y(2)=y1+6*scale
	x(3)=x1-2*scale
	y(3)=y(2)
	n=3
	call line2(n,x,y)

	y(1)=y1-scale
	y(2)=y1+7*scale

	x(1)=x1+1.5*scale
	x(2)=x(1)
	n=2
	call line2(n,x,y)

	x(1)=x1+1.25*scale
	x(2)=x(1)
	n=2
	call line2(n,x,y)

	x(1)=x1+1.75*scale
	x(2)=x(1)
	n=2
	call line2(n,x,y)

	x(1)=x1+5*scale
	y(1)=y1+6*scale
	x(2)=x1+1.5*scale
	y(2)=y(1)
	n=2
	call line2(n,x,y)

	x(1)=x1+1.5*scale
	y(1)=y1
	x(2)=x1+5*scale
	y(2)=y1
	x(3)=x(2)
	y(3)=y1+3*scale
	x(4)=x(1)
	y(4)=y(3)
	n=4
	call line2(n,x,y)

	listfile=.false.
	call widearrow(x1+1.9*scale,y1+3*scale,-1.,0.,1.8*scale)
	call bullet(x1+5*scale,y1,0.5*scale)
	listfile=.true.

	return
	end

* op amp

	subroutine opamp(x1,y1,opampwidth)

	dimension x(4),y(4)
	logical init,listfile
*	external pl2ca
	common/splot_common/xmin,xmax,ymin,ymax,init,listfile

	write(30,*)'op amp:'
	write(30,*)x1,y1
	write(30,*)'width: ',opampwidth
	write(30,*)

	scale=opampwidth/15

	x(1)=x1
	y(1)=y1
	x(2)=x1-15*scale
	y(2)=y1+8.5*scale
	x(3)=x(2)
	y(3)=y1-8.5*scale
	x(4)=x1
	y(4)=y1
	n=4
	call line2(n,x,y)

	x(1)=x1-14*scale
	y(1)=y1+4*scale
	x(2)=x(1)+1.2*scale
	y(2)=y(1)
	n=2
	call line2(n,x,y)

	x(1)=x1-14*scale
	y(1)=y1-4*scale
	x(2)=x(1)+1.2*scale
	y(2)=y(1)
	n=2
	call line2(n,x,y)

	x(1)=x1-13.4*scale
	y(1)=y1-4.6*scale
	x(2)=x(1)
	y(2)=y1-3.4*scale
	n=2
	call line2(n,x,y)

	return
	end

c and gate

	subroutine andgate(x1,y1,x2,y2)

	dimension x(4),y(4)
	logical init,listfile
*	external pl2ca
	common/splot_common/xmin,xmax,ymin,ymax,init,listfile

	write(30,*)'and gate:'
	write(30,*)x1,y1
	write(30,*)x2,y2
	write(30,*)

	gatewidth=sqrt((y2-y1)**2+(x2-x1)**2)
	scale=gatewidth/8
	angle=atan2d(y2-y1,x2-x1)
	xforward=(x2-x1)/8
	yforward=(y2-y1)/8
	xupward=scale*cosd(angle+90)
	yupward=scale*sind(angle+90)

	x(1)=x1+4.5*xforward+3.5*xupward
	y(1)=y1+4.5*yforward+3.5*yupward
	x(2)=x1+3.5*xupward
	y(2)=y1+3.5*yupward
	x(3)=x1-3.5*xupward
	y(3)=y1-3.5*yupward
	x(4)=x1+4.5*xforward-3.5*xupward
	y(4)=y1+4.5*yforward-3.5*yupward
	n=4
	call line2(n,x,y)

	listfile=.false.
	call arc(x1+4.5*xforward,y1+4.5*yforward,3.5*scale,
	2     angle-90,angle+90,18.)
	listfile=.true.

	return
	end

c or gate

	subroutine orgate(x1,y1,x2,y2)

	dimension x(2),y(2)
	logical init,listfile
*	external pl2ca
	common/splot_common/xmin,xmax,ymin,ymax,init,listfile

	write(30,*)'or gate:'
	write(30,*)x1,y1
	write(30,*)x2,y2
	write(30,*)

	gatewidth=sqrt((y2-y1)**2+(x2-x1)**2)
	scale=gatewidth/9
	angle=atan2d(y2-y1,x2-x1)
	xforward=(x2-x1)/9
	yforward=(y2-y1)/9
	xupward=scale*cosd(angle+90)
	yupward=scale*sind(angle+90)

	x(1)=x1-(sqrt(77.)-sqrt(68.75))*xforward+3.5*xupward
	y(1)=y1-(sqrt(77.)-sqrt(68.75))*yforward+3.5*yupward
	x(2)=x2-3.5*sqrt(3.)*xforward+3.5*xupward
	y(2)=y2-3.5*sqrt(3.)*yforward+3.5*yupward
	n=2
	call line2(n,x,y)

	x(1)=x1-(sqrt(77.)-sqrt(68.75))*xforward-3.5*xupward
	y(1)=y1-(sqrt(77.)-sqrt(68.75))*yforward-3.5*yupward
	x(2)=x2-3.5*sqrt(3.)*xforward-3.5*xupward
	y(2)=y2-3.5*sqrt(3.)*yforward-3.5*yupward
	n=2
	call line2(n,x,y)

	listfile=.false.
	call arc(x1-sqrt(77.)*xforward,y1-sqrt(77.)*yforward,
	2     9*scale,angle-asind(3.5/9),angle+asind(3.5/9),5.)

	call arc(x2-3.5*sqrt(3.)*xforward+3.5*xupward,
	2        y2-3.5*sqrt(3.)*yforward+3.5*yupward,
	2        7*scale,angle-90,angle-30,5.)

	call arc(x2-3.5*sqrt(3.)*xforward-3.5*xupward,
	2        y2-3.5*sqrt(3.)*yforward-3.5*yupward,
	2        7*scale,angle+30,angle+90,5.)
	listfile=.true.

	return
	end

c xor gate

	subroutine xorgate(x1,y1,x2,y2)

	dimension x(2),y(2)
	logical init,listfile
*	external pl2ca
	common/splot_common/xmin,xmax,ymin,ymax,init,listfile

	write(30,*)'or gate:'
	write(30,*)x1,y1
	write(30,*)x2,y2
	write(30,*)

	gatewidth=sqrt((y2-y1)**2+(x2-x1)**2)
	scale=gatewidth/9
	angle=atan2d(y2-y1,x2-x1)
	xforward=(x2-x1)/9
	yforward=(y2-y1)/9
	xupward=scale*cosd(angle+90)
	yupward=scale*sind(angle+90)

	x(1)=x1-(sqrt(77.)-sqrt(68.75))*xforward+3.5*xupward
	y(1)=y1-(sqrt(77.)-sqrt(68.75))*yforward+3.5*yupward
	x(2)=x2-3.5*sqrt(3.)*xforward+3.5*xupward
	y(2)=y2-3.5*sqrt(3.)*yforward+3.5*yupward
	n=2
	call line2(n,x,y)

	x(1)=x1-(sqrt(77.)-sqrt(68.75))*xforward-3.5*xupward
	y(1)=y1-(sqrt(77.)-sqrt(68.75))*yforward-3.5*yupward
	x(2)=x2-3.5*sqrt(3.)*xforward-3.5*xupward
	y(2)=y2-3.5*sqrt(3.)*yforward-3.5*yupward
	n=2
	call line2(n,x,y)

	listfile=.false.
	call arc(x1-sqrt(77.)*xforward,y1-sqrt(77.)*yforward,
	2     9*scale,angle-asind(3.5/9),angle+asind(3.5/9),5.)

	call arc(x1-(2+sqrt(77.))*xforward,y1-(2+sqrt(77.))*yforward,
	2     9*scale,angle-asind(3.5/9),angle+asind(3.5/9),5.)

	call arc(x2-3.5*sqrt(3.)*xforward+3.5*xupward,
	2        y2-3.5*sqrt(3.)*yforward+3.5*yupward,
	2        7*scale,angle-90,angle-30,5.)

	call arc(x2-3.5*sqrt(3.)*xforward-3.5*xupward,
	2        y2-3.5*sqrt(3.)*yforward-3.5*yupward,
	2        7*scale,angle+30,angle+90,5.)
	listfile=.true.

	return
	end

c buffer gate

	subroutine buffergate(x1,y1,x2,y2)

	dimension x(4),y(4)
	logical init,listfile
*	external pl2ca
	common/splot_common/xmin,xmax,ymin,ymax,init,listfile

	write(30,*)'buffer gate:'
	write(30,*)x1,y1
	write(30,*)x2,y2
	write(30,*)

	gatewidth=sqrt((y2-y1)**2+(x2-x1)**2)
	scale=gatewidth/5
	angle=atan2d(y2-y1,x2-x1)
	xupward=5/sqrt(3.)*scale*cosd(angle+90)
	yupward=5/sqrt(3.)*scale*sind(angle+90)

	x(1)=x2
	y(1)=y2
	x(2)=x1+xupward
	y(2)=y1+yupward
	x(3)=x1-xupward
	y(3)=y1-yupward
	x(4)=x2
	y(4)=y2
	n=4
	call line2(n,x,y)

	return
	end

c lamp

	subroutine lamp(x1,y1,x2,y2)

	dimension x(4),y(4)
	logical init,listfile
*	external pl2ca
	common/splot_common/xmin,xmax,ymin,ymax,init,listfile

	write(30,*)'lamp:'
	write(30,*)x1,y1
	write(30,*)x2,y2
	write(30,*)

	xlampwidth=sqrt((y2-y1)**2+(x2-x1)**2)
	scale=xlampwidth/2
	angle=atan2d(y2-y1,x2-x1)
	xforward=(x2-x1)/2
	yforward=(y2-y1)/2
	xupward=scale*cosd(angle+90)
	yupward=scale*sind(angle+90)

	x(1)=x1
	y(1)=y1
	x(2)=x1+5*xupward
	y(2)=y1+5*yupward
	n=2
	call line2(n,x,y)

	x(1)=x2
	y(1)=y2
	x(2)=x2+5*xupward
	y(2)=y2+5*yupward
	n=2
	call line2(n,x,y)

	x(1)=x1-1.25*xforward+5.75*xupward
	y(1)=y1-1.25*yforward+5.75*yupward
	x(2)=x1-1.25*xforward+1.5*xupward
	y(2)=y1-1.25*yforward+1.5*yupward
	x(3)=x2+1.25*xforward+1.5*xupward
	y(3)=y2+1.25*yforward+1.5*yupward
	x(4)=x2+1.25*xforward+5.75*xupward
	y(4)=y2+1.25*yforward+5.75*yupward
	n=4
	call line2(n,x,y)
	
	listfile=.false.
	call arc(x1+xforward+5*xupward,y1+yforward+5*yupward,
	2     scale,angle,angle+180,18.)
	call arc(x1+xforward+5.75*xupward,y1+yforward+5.75*yupward,
	2     2.25*scale,angle,angle+180,18.)
	listfile=.true.

	return
	end

c thread

	subroutine thread(x1,y1,x2,y2,nteeth)

	dimension x(102),y(102),dx(2),dy(2)
	logical init,listfile
*	external pl2ca
	common/splot_common/xmin,xmax,ymin,ymax,init,listfile
	if(nteeth.le.0)
	2     call splot_err('last argument in THREAD is invalid')
	if(nteeth.gt.50)
	2     call splot_err('last argument in THREAD is too large')

	write(30,*)'thread:'
	write(30,*)x1,y1
	write(30,*)x2,y2
	write(30,*)'number of teeth: ',nteeth
	write(30,*)

	d=sqrt((y2-y1)**2+(x2-x1)**2)
	teethsize=d/nteeth/sqrt(2.)
	angle=atan2d(y2-y1,x2-x1)
	dx(1)=teethsize*cosd(angle+45)
	dy(1)=teethsize*sind(angle+45)
	dx(2)=teethsize*cosd(angle-45)
	dy(2)=teethsize*sind(angle-45)

	x(1)=x1
	y(1)=y1
	x(2)=x1+dx(2)/2
	y(2)=y1+dy(2)/2
	n=2
	do i=1,nteeth
	do j=1,2
	  n=n+1
	  x(n)=x(n-1)+dx(j)
	  y(n)=y(n-1)+dy(j)
	enddo
	enddo
	x(n)=x2
	y(n)=y2
	call line2(n,x,y)

	return
	end

***************************************************************************

* draw a dimension in the diagram

	subroutine showdistance(x1,y1,x2,y2,ndecimal,size,xtext,ytext)
	character num*11,format_string*7
	logical insidearrow
	common scale

	if(ndecimal.le.0)
	2     call splot_err('last argument in SHOWDISTANCE is invalid')

	write(30,*)'showdistance:'
	write(30,*)x1,y1
	write(30,*)x2,y2
	write(30,*)'number of decimal places: ',ndecimal
	write(30,*)'size: ',size
	write(30,*)

* 10-point text

	itextheight=10
	textheight=itextheight*0.7*size
	itextheight=itextheight*size

* length of arrows

	arrowlength=0.4*size/scale*72/2.54
	arrowline=0.2*size/scale*72/2.54
	bordersize=0.1*size*scale
	
* find distance and direction

	distance=sqrt((x2-x1)**2+(y2-y1)**2)
	angle=atan2d(y2-y1,x2-x1)

* get string

        format_string='(f11. )'
        write(format_string(6:6),'(i1)')ndecimal
	write(num,format_string)distance
	do i=1,11
	  if(num(i:i).ne.' ')then
	    if(num(i:i).eq.'.')then
	      num(i-1:i-1)='0'
	      numstart=i-1
	    else
	      numstart=i
	    endif
	    goto 1
	  endif
	enddo
1	continue

* can arrows fit inside? 

	  if(distance.gt.4*arrowlength)then
	    insidearrow=.true.
	  else
	    insidearrow=.false.
	  endif

* draw arrows

	if(insidearrow)then
	  call arrow(x1,y1,x1-x2,y1-y2,arrowlength)
	  call arrow(x2,y2,x2-x1,y2-y1,arrowlength)
	else
	  call arrow(x1,y1,x2-x1,y2-y1,arrowlength)
	  call arrow(x2,y2,x1-x2,y1-y2,arrowlength)
	endif

* draw small line on each arrow

	call line(2,x1+arrowline/2*cosd(angle+90),
	2           y1+arrowline/2*sind(angle+90),
	2           x1-arrowline/2*cosd(angle+90),
	2           y1-arrowline/2*sind(angle+90))
	call line(2,x2+arrowline/2*cosd(angle+90),
	2           y2+arrowline/2*sind(angle+90),
	2           x2-arrowline/2*cosd(angle+90),
	2           y2-arrowline/2*sind(angle+90))

* draw tail to each arrow

	if(insidearrow)then
	  call line(2,x1,y1,x2,y2)
	else
	  taillength=2*arrowlength
	  call line(2,x1,y1,x1-taillength*cosd(angle),
	2     y1-taillength*sind(angle))
	  call line(2,x2,y2,x2+taillength*cosd(angle),
	2     y2+taillength*sind(angle))
	endif

* number

	write(31,*)'/Helvetica findfont'
	write(31,*)itextheight,' scalefont'
	write(31,*)'setfont'
	
	call point((x1+x2)/2+xtext,(y1+y2)/2+ytext,x0,y0)
	
	if(insidearrow)then
	write(31,*)'newpath'
	write(31,*)'(',num(numstart:11),') stringwidth'
	write(31,*)'pop -2 div ',x0-bordersize,' add'
	write(31,*)y0-textheight/2-bordersize,' moveto'
	write(31,*)'0 ',textheight+2*bordersize,' rlineto'
	write(31,*)'(',num(numstart:11),') stringwidth'
	write(31,*)'pop ',2*bordersize,' add'
	write(31,*)'0 rlineto'
	write(31,*)'0 ',-textheight-2*bordersize,' rlineto'
	write(31,*)'closepath 1 setgray fill 0 setgray'
	endif
	
	write(31,*)'(',num(numstart:11),') stringwidth'
	write(31,*)'pop -2 div ',x0,' add'
	write(31,*)y0-textheight/2,' moveto'
	write(31,*)'(',num(numstart:11),') show'
	
	end

* convert units of coordinates of a point

	subroutine point(x,y,xp,yp)
	dimension origin(2)
	common scale,origin
	xp=(x-origin(1))*scale
	yp=(y-origin(2))*scale
	end

* limits of bounding box

	subroutine bbox_limits(xp,yp)
	common /common_bbox/bbxmin,bbxmax,bbymin,bbymax
	if(xp.gt.bbxmax)bbxmax=xp
	if(xp.lt.bbxmin)bbxmin=xp
	if(yp.gt.bbymax)bbymax=yp
	if(yp.lt.bbymin)bbymin=yp
	end

*------------------------------------------------------------------------------

* smooth curve using cubic spline interpolation

      subroutine cubic_spline(n1,x1,y1,n2,x2,y2)
      parameter(maxpoints=1000)
      dimension x1(maxpoints),y1(maxpoints),x2(maxpoints),y2(maxpoints),
     +     deriv(maxpoints)
      call spline(x1,y1,n1,1.0e30,1.0e30,deriv)
      dx=(x1(n1)-x1(1))/(n2-1)
      x0=x1(1)
      ktry=1
      do i=1,n2
        x2(i)=x0
        x0=x0+dx
	if(i.eq.n2)x2(i)=x1(n1)
        call splint(x1,y1,deriv,n1,x2(i),y2(i),ktry)
      enddo
      end

* cubic spline interpolation
* see Numerical Recipes

      subroutine spline(x,y,n,yp1,ypn,y2)
* Given arrays x and y of length n containing a tabulated function, i.e.,
* y(i)=f(x(i)), with x(1) < x(2) < ... < x(n), and given values yp1 and
* ypn for the first derivative of the interpolating function and points
* 1 and n, respectively, this routine returns an array y2 of length n
* which contains the second derivatives of the interpolating function
* at the tabulated points x(i).  If yp1 and/or ypn are equal to 1.0d30
* or larger, the routine is signaled to set the corresponding boundary
* condition for a natural spline, with zero second derivative on that
* boundary.

      implicit none
      integer maxpoints,n
      parameter(maxpoints=100000)
      real x(n),y(n),y2(n),u(maxpoints),yp1,ypn,sig,p,
     +     qn,un
      integer i,k

* the lower boundary is set either to be natural
      if(yp1.gt.0.99e30)then
        y2(1)=0
        u(1)=0
* or else to have a specified first derivative.
      else
        y2(1)=-0.5e0
        u(1)=(3/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
* This is the decomposition loop of the tridiagonal algorithm.  y2 and u
* are used for temporary storage of the decomposed factors.
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2
        y2(i)=(sig-1)/p
        u(i)=(6*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     +       /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11      continue
* The upper boundary condition is set either to be natural
      if(ypn.gt.0.99e30)then
        qn=0
        un=0
* or else to have a specified first derivative.
      else
        qn=0.5
        un=(3/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1)
* This is the backsubstitution loop of the tridiagonal algorithm.
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12      continue
      end

      subroutine splint(xa,ya,y2a,n,x,y,ktry)
* Given the arrays xa and ya of length n, which tabulate a function
* and give the array y2a, which is the output from spline, and given a
* value of x, this routine returns a cubic-spline interpolated value y.
* The argument ktry is an approximate value of the index.
      implicit none
      integer n
      real xa(n),ya(n),y2a(n),x,y,h,a,b
      integer klo,khi,ktry

* This part differs from Numerical Recipes
* Find the index starting with a trial value ktry.
      if(xa(ktry).lt.x)then
        khi=ktry
        do while(xa(khi).lt.x)
          khi=khi+1
        enddo
        if(khi.gt.n)then
          call splot_err('Error in splint: khi out of range')
        endif
        klo=khi-1
      else
        klo=ktry
        do while(xa(klo).gt.x)
          klo=klo-1
        enddo
        if(klo.lt.1)then
          call splot_err('Error in splint: klo out of range')
        endif
        khi=klo+1
      endif
* return value at beginning of interval
      ktry=klo
      h=xa(khi)-xa(klo)
      if(h.eq.0.0)then
        call splot_err('Error in splint: invalid values of xa')
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)
     +     +((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*h**2/6
      end

c----------------------------------------------------------------------

c colors

	subroutine init_colors
	dimension rgb(3,20)
	character names(20)*10
	
	data ncolors/14/
	data (names(j),(rgb(i,j),i=1,3),j=1,14)/
	2     'white',1,1,1,
	2     'black',0,0,0,
	2     'grey',0.5,0.5,0.5,
	2     'red',1,0,0,
	2     'yellow',1,1,0,
	2     'green',0,1,0,
	2     'cyan',0,1,1,
	2     'blue',0,0,1,
	2     'magenta',1,0,1,
	2     'orange',1,0.5,0,
	2     'aqua',0,1,0.5,
	2     'purple',0.5,0,1,
	2     'pink',1,0.85,0.95,
	2     'brown',0.5,0.2,1/

	write(30,*)'initialize colors'
	do i=1,ncolors
	  write(30,'(a10,3f6.3)')names(i),(rgb(j,i),j=1,3)
	enddo
	
	do i=1,ncolors
	  write(31,100)names(i),(rgb(j,i),j=1,3)
100	  format('/',a10,'{',3f6.3,' setrgbcolor} def')
        enddo
        
        end

	subroutine color(name)
	character name*(*)
	write(30,*)'set color ',name
	write(31,*)name
	end

	subroutine wavelength(x)
	write(30,'(a,f5.0)')'set wavelength ',x
	if(x.lt.400.0)
	1    call splot_err('wavelength less than 400 nm')
	if(x.gt.750.0)
	2    call splot_err('wavelength greater than 750 nm')
	call wavelength_to_rgb(x,r,g,b)
	write(31,'(3f6.3,a)')r,g,b,' setrgbcolor'
	end

      subroutine wavelength_to_rgb(wavelength,r,g,b)
      r=0
      g=0
      b=0
      red=640
      yellow=590
      green=535
      blue=485
      violet=440
      w0=600
      dw=200
* red
      if(wavelength.gt.red)then
        r=1
*red to yellow
      else if(wavelength.gt.yellow)then
        r=1
        g=sind((red-wavelength)/(red-yellow)*90)
*yellow to green
      else if(wavelength.gt.green)then
        r=sind((wavelength-green)/(yellow-green)*90)
        g=1
* green to blue
      else if(wavelength.gt.blue)then 
        g=cosd((green-wavelength)/(green-blue)*90)
        b=sind((green-wavelength)/(green-blue)*90)
* blue to violet
      else if(wavelength.gt.violet)then 
        r=(blue-wavelength)/(blue-violet)
        b=1
* violet
      else
        r=1
        b=1
      endif
* superimpose sensitivity of eye
      f=exp(-((wavelength-w0)/dw)**2)
      r=f*r
      g=f*g
      b=f*b
      end

	
