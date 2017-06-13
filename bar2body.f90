!*********************************************************
!目的：将线单元转化成实体单元                            *
!方法：在线上等距插入点然后平移得到长方体柱              *
!输入：文件输入点的坐标和及杆单元信息(端点点的编号)      *
!      屏幕输入d(分段长度)、r(平移距离即柱边长的一半)    *
!输出：点的坐标(前面是柱的端头点，后面是中间插入点)      *
!内容：主程序main()=>输入input()=>等距插值子程序insert-  *
!      point(d)=>平移变换translation(r)=>输出output()    *
!时间：2017-06-12     HIT     作者：hitiler              *
!*********************************************************

module typedef
    type transform
        real,allocatable::point(:,:)
        integer,allocatable::line(:,:)
        integer,allocatable::area(:,:)
    end type
    type(transform)::t
    integer::points,lines,nps  !nps(新点数)
    real::np(1000,3),tp(4000,3),tnp(4000,3)
    !np(插入点),tp(原端点变换后的点),tnp(插入点变换后的点)
    end module
    
!=========================================================
!=========================================================
program main
    use typedef
    implicit none
    real::d,r
    character(12)::infile,outfile
    write(infile,'(a7)') 'bar.inp'
    write(outfile,'(a8)') 'body.inp'
    open(10,file=infile,status='old')
    open(11,file=outfile,status='new')
    write(*,*)"please input distance and radius:"
    read(*,*) d,r
    call input_data()
    call insertpoint(d)
    call translation(r)
!    call get_grid()
    call output_data()
    close(10)
    close(11)
    stop
    end program     
!=========================================================
    
!读取杆件信息
subroutine input_data()
    use typedef
    implicit none
    integer::i,j
    character::c
    read(10,*)  !跳过avs信息行
    read(10,*)points,lines !读第一行，点线的个数
    allocate(t%point(points,3))
    allocate(t%line(lines,2))
    do i=1,points
        read(10,*) c,t%point(i,1:3)
    end do
    do j=1,lines
        read(10,*) (c,i=1,3),t%line(j,1:2)
    end do
    return
    end subroutine
    
!等距离插入点
subroutine insertpoint(d)  !传入距离参数d
    use typedef
    implicit none
    integer::i,j,k,n=0,nump(lines)
    real::dt=0,d
    nps=0
    do i=1,lines
        dt=(((t%point(t%line(i,2),1)-t%point(t%line(i,1),1))**2)+ &
            ((t%point(t%line(i,2),2)-t%point(t%line(i,1),2))**2)+ &
            ((t%point(t%line(i,2),3)-t%point(t%line(i,1),3))**2))**0.5
!        write(*,*)dt
        n=anint(dt/d)
!        write(*,*)n
        nump(i)=n-1   !记录杆上插入点的个数
        do k=1,n-1
            np(k+nps,:)=t%point(t%line(i,1),:)+(t%point(t%line(i,2),:)&
            -t%point(t%line(i,1),:))*k*d/dt
        end do               
        nps=nps+nump(i)
    end do
!    write(*,*)nps
    return
    end subroutine
    
!平移变换
subroutine translation(r) !传入半径参数r
    use typedef
    implicit none
    integer::i,j
    real::r
    do i=1,points
        tp(((i-1)*4+1),:)=(/t%point(i,1)-r,t%point(i,2),t%point(i,3)/)
        tp(((i-1)*4+2),:)=(/t%point(i,1)+r,t%point(i,2),t%point(i,3)/)
        tp(((i-1)*4+3),:)=(/t%point(i,1),t%point(i,2)-r,t%point(i,3)/)
        tp(((i-1)*4+4),:)=(/t%point(i,1),t%point(i,2)+r,t%point(i,3)/)
    end do
    do j=1,nps
        tnp(((j-1)*4+1),:)=(/np(j,1)-r,np(j,2),np(j,3)/)
        tnp(((j-1)*4+2),:)=(/np(j,1)+r,np(j,2),np(j,3)/)
        tnp(((j-1)*4+3),:)=(/np(j,1),np(j,2)-r,np(j,3)/)
        tnp(((j-1)*4+4),:)=(/np(j,1),np(j,2)+r,np(j,3)/)
    end do
    return
    end subroutine
    
!输出点的坐标
 subroutine output_data()
    use typedef
    implicit none
    integer::i,j,k,l,m,n
    write(11,'(5(2X,I5))') 4*(points+nps),4*(points+nps)
    do i=1,4*points
        write(11,100) i,tp(i,1),tp(i,2),tp(i,3)
    end do
    do j=1,4*nps
        write(11,100) tnp(j,1),tnp(j,2),tnp(j,3)
    end do
    do l=1,4*(points+nps)
        write(11,'(2(I5),3X,A,3X,I5)') l,1,'pt',l
    end do
100 format(I5,2X,3(3X,E15.7))
    write(11,'(2(8x,I1))')1,1
    write(11,'(A)')'not,'
!    do k=1,nps
!        write(11,100) k,np(k,1),np(k,2),np(k,3)
!    end do
    do m=1,4*(points+nps)
        write(11,'(I5,3X,E15.7)') m,0
    end do 
    write(11,'(2(8x,I1))')1,1
    write(11,'(A)')'dt,'
    do n=1,4*(points+nps)
        write(11,'(I5,5X,I5)') n,n 
    end do
    deallocate(t%point)
    deallocate(t%line)
    return
    end subroutine
    