!*********************************************************
!目的：将线单元转化成实体单元                            *
!方法：在线上等距插入点然后平移得到长方体柱              *
!输入：文件输入点的坐标和及杆单元信息(端点点的编号)      *
!      屏幕输入d(分段长度)、r(平移距离即柱边长的一半)    *
!输出：点的坐标(前面是柱的端头点，后面是中间插入点)      *
!内容：主程序main()=>输入input()=>等距插值子程序insert-  *
!      point(d)=>平移变换translation(r)=>杆上节点排序get_*
!      nodenumber()=>生成网格get_grid()=>输出output()    *
!时间：2017-06-12     HIT     作者：hitiler              *
!*********************************************************

module typedef
    type transform
        real,allocatable::point(:,:)
        integer,allocatable::line(:,:)
        integer,allocatable::area(:,:)
    end type
    type(transform)::t
    integer::points,lines,nps,sumnp,sln,slg 
    !nps(新点数),sumnp(计数),sln(节点数),slg(网格数)
    real::np(1000,3),tp(4000,4),tnp(4000,4)
    !np(插入点),tp(原端点变换后的点),tnp(插入点变换后的点)
    integer::nump(100),nub(4000,5),grid(4000,3)
    !nump(新插点个数),nub(节点),grid(网格)
    end module
    
!=========================================================
!=========================================================
program main
    use typedef
    implicit none
    real::d,r
    character(12)::infile,outpointfile,outgridfile
    write(infile,'(a7)') 'bar.inp'
    write(outpointfile,'(a9)') 'point.inp'
    write(outgridfile,'(a8)') 'grid.inp'
    open(10,file=infile,status='old')
    open(11,file=outpointfile,status='replace')
    open(12,file=outgridfile,status='replace')
    write(*,*)"please input distance and radius:"
    read(*,*) d,r
    call input_data()
    call insertpoint(d)
    call translation(r)
    call get_nodenumber()
    call get_grid()
    call output_pointdata()
    call output_griddata()
    close(10)
    close(11)
    close(12)
    stop
    end program     
!=========================================================
    
!读取杆件信息
subroutine input_data()
    use typedef
    implicit none
    integer::i,j
    character::c  !c是读取但无用的信息，相当于跳过
    read(10,*)c,points,lines !读第一行，点线的个数
    read(10,*)  !跳过avs信息行
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
    integer::i,j,k,n=0
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
    integer::i,j,k,l
    real::r
    do i=1,points
        tp(((i-1)*4+1),:)=(/(i-1)*4.0+1,t%point(i,1)-r,t%point(i,2),t%point(i,3)/)
        tp(((i-1)*4+2),:)=(/(i-1)*4.0+2,t%point(i,1),t%point(i,2)-r,t%point(i,3)/)
        tp(((i-1)*4+3),:)=(/(i-1)*4.0+3,t%point(i,1)+r,t%point(i,2),t%point(i,3)/)
        tp(((i-1)*4+4),:)=(/(i-1)*4.0+4,t%point(i,1),t%point(i,2)+r,t%point(i,3)/)
    end do
    sumnp=0
    do k=1,lines
        do j=1,nump(k)
            tnp(((j+sumnp-1)*4+1),:)=(/4*points+(j+sumnp-1)*4.0+1,np(j+sumnp,1)-r,np(j+sumnp,2),np(j+sumnp,3)/)
            tnp(((j+sumnp-1)*4+2),:)=(/4*points+(j+sumnp-1)*4.0+2,np(j+sumnp,1),np(j+sumnp,2)-r,np(j+sumnp,3)/)
            tnp(((j+sumnp-1)*4+3),:)=(/4*points+(j+sumnp-1)*4.0+3,np(j+sumnp,1)+r,np(j+sumnp,2),np(j+sumnp,3)/)
            tnp(((j+sumnp-1)*4+4),:)=(/4*points+(j+sumnp-1)*4.0+4,np(j+sumnp,1),np(j+sumnp,2)+r,np(j+sumnp,3)/)
        end do
        sumnp=sumnp+nump(k)
    end do
    return
    end subroutine
    
!获得变换后每根杆上的节点编号
subroutine get_nodenumber()
    use typedef
    implicit none
    integer i,j,k
    sumnp=0
    sln=0
    do i=1,lines
        nub(1+sln,:)=(/(t%line(i,1)-1)*4+1,(t%line(i,1)-1)*4+2, &
           (t%line(i,1)-1)*4+3,(t%line(i,1)-1)*4+4,(t%line(i,1)-1)*4+1/)
        do j=1,nump(i)
            nub(sln+j+1,:)=(/4*points+(j+sumnp-1)*4+1,4*points+(j+sumnp-1)*4+2, &
               4*points+(j+sumnp-1)*4+3,4*points+(j+sumnp-1)*4+4,4*points+(j+sumnp-1)*4+1/)
        end do
        nub(sln+nump(i)+2,:)=(/(t%line(i,2)-1)*4+1,(t%line(i,2)-1)*4+2, &
           (t%line(i,2)-1)*4+3,(t%line(i,2)-1)*4+4,(t%line(i,2)-1)*4+1/)
        sumnp=sumnp+nump(i)
        sln=sln+nump(i)+2
    end do
    do k=1,sln
        write(*,'(5(3X,I5))') nub(k,1),nub(k,2),nub(k,3),nub(k,4),nub(k,5)
    end do
    return
    end subroutine

!拓扑连接形成表面网格
subroutine get_grid()
    use typedef
    implicit none
    integer::i,j,k,sg,sbg
    slg=0
    sg=0
    sbg=0
    !上述三个变量用来计数，没有意义
    do i=1,lines
        sg=0
        do j=1,nump(i)+1
            do k=1,4
                grid(slg+sg+k,:)=(/nub(sbg+j,k),nub(sbg+j,k+1),nub(sbg+j+1,k+1)/)
                grid(slg+sg+4+k,:)=(/nub(sbg+j,k),nub(sbg+j+1,k+1),nub(sbg+j+1,k)/)
            end do
            sg=sg+8
        end do
        slg=slg+(nump(i)+1)*8
        sbg=sbg+nump(i)+2
    end do
    end subroutine
    
!输出点的坐标
subroutine output_pointdata()
    use typedef
    implicit none
    integer::i,j,k,l,m,n
    write(11,'(5(2X,I5))') 4*(points+nps),4*(points+nps),1,0,0
    do i=1,4*points
        write(11,100) int(tp(i,1)),tp(i,2),tp(i,3),tp(i,4)
    end do
    do j=1,4*nps
        write(11,100) int(tnp(j,1)),tnp(j,2),tnp(j,3),tnp(j,4)
    end do
    do l=1,4*(points+nps)
        write(11,'(2(I5),3X,A,3X,I5)') l,1,'pt',l
    end do
100 format(I5,2X,3(3X,E15.7))
    write(11,'(2(8x,I1))') 1,1
    write(11,'(A)')'not,'
!    do k=1,nps
!        write(11,100) k,np(k,1),np(k,2),np(k,3)
!    end do
    do m=1,4*(points+nps)
        write(11,'(I5,3X,E15.7)') m,0
    end do 
    write(11,'(2(8x,I1))') 1,1
    write(11,'(A)')'dt,'
    do n=1,4*(points+nps)
        write(11,'(I5,5X,I5)') n,n 
    end do
    deallocate(t%point)
    deallocate(t%line)
    return
    end subroutine
    
!输出面的拓扑信息
subroutine output_griddata()
    use typedef
    implicit none
    integer::i,j,k
    write(12,'(5(2X,I5))') 4*(points+nps),slg,1,0,0
    do i=1,4*points
        write(12,100) int(tp(i,1)),tp(i,2),tp(i,3),tp(i,4)
    end do
    do j=1,4*nps
        write(12,100) int(tnp(j,1)),tnp(j,2),tnp(j,3),tnp(j,4)
    end do
100 format(I5,2X,3(3X,E15.7))
    do k=1,slg
        write(12,200) k,'1','tri',grid(k,1),grid(k,2),grid(k,3)
    end do
200 format(I5,2(2X,A3),3(3X,I5))
    write(12,'(2(8x,I1))') 1,1
    write(12,'(A)')'temputure,'
    end subroutine
    