!+---------------------------------------------------------------------+
!| This module contains input and output subroutines.                  |
!+---------------------------------------------------------------------+
!| Writen by Jian Fang, 2020-04-16                                     |
!+---------------------------------------------------------------------+
module readwrite
  !
  use commvardefine
  use h5readwrite
  !
  implicit none
  !
  contains
  !!+------------------------------------------------------------------+
  !|this subroutine is used to read input files.                       |
  !+-------------------------------------------------------------------+
  subroutine readinput
    !
    integer :: nrank,n,i,j,k
    !
    print*,'==========================readinput=========================='
    !
    open(11,file='datin/input.3dp.dat',form='formatted',status='old')
    read(11,"(//)")
    read(11,*)flowtype
    read(11,"()")
    read(11,*)im,jm,km
    read(11,"()")
    read(11,*)lihomo,ljhomo,lkhomo
    read(11,"()")
    read(11,"()")
    read(11,"(////////)")
    read(11,*)ref_t,reynolds,mach
    read(11,"(////////)")
    read(11,*)isp,jsp,ksp
    close(11)
    print*,'>> input.3dp.dat ... done!'
    !
    open(11,file='datin/lwrite.dat',form='formatted',status='old')
    read(11,"(//////)")
    read(11,*)deltat
    close(11)
    print*,'>> lwrite.dat ... done.'
    !
    gridfile='datin/grid.h5'
    outfolder='Outdat/'
    !
    open(11,file='datin/parallel.info',form='formatted',status='old')
    read(11,"()")
    read(11,*)isize,jsize,ksize
    !
    irkm=isize-1
    jrkm=jsize-1
    krkm=ksize-1
    !
    nproc=isize*jsize*ksize
    rankmax=nproc-1
    allocate( ink(0:rankmax),jnk(0:rankmax),knk(0:rankmax),            &
              iop(0:rankmax),jop(0:rankmax),kop(0:rankmax),            &
              imp(0:rankmax),jmp(0:rankmax),kmp(0:rankmax)             )
    !
    read(11,"()")
    do n=0,rankmax
      read(11,*)nrank,ink(n),jnk(n),knk(n),imp(n),jmp(n),kmp(n),       &
                                           iop(n),jop(n),kop(n)
    end do
    close(11)
    print*,'>> parallel.info ... done!'
    !
    print*,'==========================readinput=========================='
    print*
    !
  end subroutine readinput
  !
  subroutine readinput_astrr(inputfile)
    !
    character(len=*),intent(in) :: inputfile
    !
    integer :: nrank,n,i,j,k
    !
    print*,'==========================readinput=========================='
    !
    open(11,file=inputfile,form='formatted',status='old')
    read(11,'(///////)')
    read(11,*)im,jm,km
    read(11,"(/)")
    read(11,*)lihomo,ljhomo,lkhomo
    read(11,'(//////////)')
    read(11,*)ref_t,reynolds,mach
    read(11,'(///////////////////////////)')
    read(11,'(A)')gridfile
    close(11)
    print*,'>> ',inputfile
    !
    print*,gridfile
    ! stop
    !
    outfolder='outdat/'
    ! gridfile='datin/grid.h5'
    !
    open(11,file='datin/parallel.info',form='formatted',status='old')
    read(11,"()")
    read(11,*)isize,jsize,ksize
    !
    irkm=isize-1
    jrkm=jsize-1
    krkm=ksize-1
    !
    nproc=isize*jsize*ksize
    rankmax=nproc-1
    allocate( ink(0:rankmax),jnk(0:rankmax),knk(0:rankmax),            &
              iop(0:rankmax),jop(0:rankmax),kop(0:rankmax),            &
              imp(0:rankmax),jmp(0:rankmax),kmp(0:rankmax)             )
    !
    read(11,"()")
    do n=0,rankmax
      read(11,*)nrank,ink(n),jnk(n),knk(n),imp(n),jmp(n),kmp(n),       &
                                           iop(n),jop(n),kop(n)
    end do
    close(11)
    print*,'>> parallel.info ... done!'
    !
    print*,'==========================readinput=========================='
    print*
    !
  end subroutine readinput_astrr
  !+-------------------------------------------------------------------+
  !|the end of the subroutine readinput                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to read pepin.dat                         |
  !+-------------------------------------------------------------------+
  subroutine readpepin
    !
    use commvardefine, only: im,jm,km,ref_t,Reynolds,Mach
    !
    open(11,file='datin/pepin.dat')
    read(11,"(///)")
    read(11,*)im,jm,km
    read(11,"(//)")
    read(11,*)ref_t,Reynolds,Mach
    close(11)
    print*,' >> pepin.dat ...done!'
    !
  end subroutine readpepin
  !+-------------------------------------------------------------------+
  !|the end of the subroutine Readpepin                                |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !|This subroutine is used to read grid file with HDF5                |
  !+-------------------------------------------------------------------+
  subroutine H5_ReadGrid
    !
    print*, ' >> ', trim(gridfile)
    !
    allocate(x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km))
    call H5ReadArray(x,im,jm,km,'x',trim(gridfile))
    call H5ReadArray(y,im,jm,km,'y',trim(gridfile))
    call H5ReadArray(z,im,jm,km,'z',trim(gridfile))
    !

  end subroutine H5_ReadGrid
  !+-------------------------------------------------------------------+
  !|the end of the subroutine H5_ReadGrid                              |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !|This subroutine is used to read grid file with HDF5                |
  !+-------------------------------------------------------------------+
  subroutine H5_WriteMBgrid(blkdata)
    !
    type(block) :: blkdata(:)
    !
    integer :: nblk,jblk
    character(len=3) :: blkname
    logical :: lfilext
    !
    nblk=size(blkdata)
    write(*,'(A,I0)')'  ** number of blocks= ',nblk
    !
    inquire(file='grid.h5',exist=lfilext)
    !
    if(lfilext) call system('mv  -v grid.h5 grid.bak')
    !
    call H5WriteArray(nblk,'numblk','grid.h5')
    do jblk=1,nblk
      im=blkdata(jblk)%im
      jm=blkdata(jblk)%jm
      km=blkdata(jblk)%km-1
      !
      write(blkname,"(I3.3)")jblk
      !
      call H5WriteArray(blkdata(jblk)%x,im,jm,km,'x_'//blkname,'grid.h5')
      call H5WriteArray(blkdata(jblk)%y,im,jm,km,'y_'//blkname,'grid.h5')
      call H5WriteArray(blkdata(jblk)%z,im,jm,km,'z_'//blkname,'grid.h5')
      !
    enddo
    !
    !
    ! call H5WriteArray(x,im,jm,km,'x','grid.h5')
    !
  end subroutine H5_WriteMBgrid
  !+-------------------------------------------------------------------+
  !|the end of the subroutine H5_ReadGrid                              |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !|The subroutine is used to screen print the informations of the     |
  !|post-process program.                                              |
  !+-------------------------------------------------------------------+
  subroutine informout(npost)
    !
    integer,intent(in) :: npost
    !
    character(len=53) :: whattheprogramdo
    !
    if(npost==-1) then
      whattheprogramdo='                                      grid generation'
    elseif(npost==-2) then
      whattheprogramdo='                                      flow generation'
    elseif(npost==-3) then
      whattheprogramdo='                                             test run'
    elseif(npost==-4) then
      whattheprogramdo='                                         scheme study'
    elseif(npost==-5) then
      whattheprogramdo='                                          wing design'
    elseif(npost==0) then
      whattheprogramdo='                          read nsamples from meanflow'
    elseif(npost==2) then
      whattheprogramdo='            read instant flow and output tecplot data'
    elseif(npost==201) then
      whattheprogramdo='                            output turbulent CS field'
    elseif(npost==202) then
      whattheprogramdo='                 outputing a series of 2D flow slices'
    elseif(npost==203) then
      whattheprogramdo='             get statistics of one instantaneous flow'
    elseif(npost==204) then
      whattheprogramdo='                               process monitor status'
    elseif(npost==205) then
      whattheprogramdo='                              process flow separation'
    elseif(npost==206) then
      whattheprogramdo='                                  analyze shock angle'
    elseif(npost==210) then
      whattheprogramdo='                          to view turbulent structure'
    elseif(npost==211) then
      whattheprogramdo='                              to view a 2-D x-y slice'
    elseif(npost==212) then
      whattheprogramdo='                              to view a 2-D x-z slice'
    elseif(npost==213) then
      whattheprogramdo='                          outputing a series of tecyz'
    elseif(npost==214) then
      whattheprogramdo='                     writing a fluctuatant flow field'
    elseif(npost==215) then
      whattheprogramdo='                         instant channel data process'
    elseif(npost==3) then
      whattheprogramdo='                  read raw data and output statistics'
    elseif(npost==301) then
      whattheprogramdo='                         get results of channel flows'
    elseif(npost==302) then
      whattheprogramdo='                        get results of boundary layer'
    elseif(npost==303) then
      whattheprogramdo='                        get results of diffuser flow'
    elseif(npost==304) then
      whattheprogramdo='                        get results of cylinder flow'
    elseif(npost==305) then
      whattheprogramdo='                     get results of tip leakage flow'
    elseif(npost==4) then
      whattheprogramdo='           read raw data and output budget statistics'
    elseif(npost==401) then
      whattheprogramdo='                                 Analysis budget data'
    elseif(npost==5) then
      whattheprogramdo='                  Interpolate flowfield to a new mesh'
    elseif(npost==501) then
      whattheprogramdo='                              generate a inflow plane'
    elseif(npost==502) then
      whattheprogramdo='                re-normlize a result from other solver'
    elseif(npost==6) then
      whattheprogramdo='                        calculate 2 point correlation'
    elseif(npost==7) then
      whattheprogramdo='                           calculate spatial spectrum'
    elseif(npost==701) then
      whattheprogramdo='             read spectrum data and output statistics'
    elseif(npost==702) then
      whattheprogramdo='       get the sepctra of an instantanesou flow field'
    elseif(npost==8) then
      whattheprogramdo='                                analysis stokes layer'
    elseif(npost==801) then
      whattheprogramdo='                                 analysis kslice data'
    elseif(npost==802) then
      whattheprogramdo='                                 analysis islice data'
    elseif(npost==803) then
      whattheprogramdo='                                 analysis jkslice data'
    elseif(npost==9) then
      whattheprogramdo='                   statistical data unformatted->HDF5'
    elseif(npost==901) then
      whattheprogramdo='                        islice data unformatted->HDF5'
    elseif(npost==902) then
      whattheprogramdo='                        kslice data unformatted->HDF5'
    elseif(npost==903) then
      whattheprogramdo='             meanflow & budget data unformatted->HDF5'
    elseif(npost==907) then
      whattheprogramdo='flowfield,supply,meanflow,budget data HDF5->Unformatted'
    elseif(npost==908) then
      whattheprogramdo='                          convert jslice data to hdf5'
    elseif(npost==1000) then
      whattheprogramdo='                                        user dedefine'
    elseif(npost==1001) then
      whattheprogramdo='             process flow field from the code hipstar'
    elseif(npost==1002) then
      whattheprogramdo='                process flow field from the code sbli'
    elseif(npost==1003) then
      whattheprogramdo='              generate mean flow profile as BL inflow'
    elseif(npost==1005) then
      whattheprogramdo='                                process data from UPM'
    elseif(npost==1006) then
      whattheprogramdo='            write a tecplot script file for animation'
    else
      whattheprogramdo='                                something you defined'
    end if
    !
    print*,'==========================informout=========================='
    print*,'                   im                  jm                  km'
    write(*,"(2X,3(I20))")im,jm,km
    print*,'-------------------------------------------------------------'
    print*,'             Reynolds                Mach               ref_t'
    write(*,"(2X,3(F20.5))")Reynolds,mach,ref_t
    print*,'-------------------------------------------------------------'
    print*,'    Total Blocks        I Block        J Block        K Block'
    write(*,"(2X,4(I15))")nproc,isize,jsize,ksize
    print*,'-------------------------------------------------------------'
    print*,'   npost                             what will the program do'
    write(*,"(2X,I7,A53)")npost,whattheprogramdo
    !
    print*,'==========================informout=========================='
    ! print*
    !
  end subroutine InformOut
  !+-------------------------------------------------------------------+
  !|The end of the subroutine InformOut.                               |
  !+-------------------------------------------------------------------+
  !!
  subroutine writexml
    !
    use commvardefine, only: im,jm,km
    !
    open(18,file="xdmf2d.xmf")
    write(18,*)'<?xml version="1.0" ?>'
    write(18,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(18,*)'<Xdmf Version="2.0">'
    write(18,*) " <Domain>"
    write(18,*) '   <Grid Name="mesh1" GridType="Uniform">'
    write(18,*) '     <Topology TopologyType="3DSMesh" NumberOfElements="',im+1,jm+1,1,'/>'
    write(18,*) '     <Geometry GeometryType="X_Y">'
    write(18,*) '       <DataItem Dimensions="',im+1,jm+1,1,'" NumberType="Float" Precision="8" Format="HDF">'
    write(18,*) '        datin/grid.h5:/x'
    write(18,*) '       </DataItem>'
    write(18,*) '       <DataItem Dimensions="',im+1,jm+1,1,'" NumberType="Float" Precision="8" Format="HDF">'
    write(18,*) '        datin/grid.h5:/y'
    write(18,*) '       </DataItem>'
    write(18,*) '       <DataItem Dimensions="',im+1,jm+1,1,'" NumberType="Float" Precision="8" Format="HDF">'
    write(18,*) '        datin/grid.h5:/z'
    write(18,*) '       </DataItem>'
    write(18,*) '     </Geometry>'
    write(18,*) '     <Attribute Name="Pressure" AttributeType="Scalar" Center="Cell">'
    write(18,*) '       <DataItem Dimensions="',im+1,jm+1,1,'" NumberType="Float" Precision="8" Format="HDF">'
    write(18,*) '        outdat/flowfield.h5:/p'
    write(18,*) '       </DataItem>'
    write(18,*) '     </Attribute>'
    write(18,*) '     <Attribute Name="VelocityX" AttributeType="Scalar" Center="Node">'
    write(18,*) '       <DataItem Dimensions="',im+1,jm+1,1,'" NumberType="Float" Precision="8" Format="HDF">'
    write(18,*) '        outdat/flowfield.h5:/u1'
    write(18,*) '       </DataItem>'
    write(18,*) '     </Attribute>'
    write(18,*) '   </Grid>'
    write(18,*) ' </Domain>'
    write(18,*) '</Xdmf>'
    close(18)
    print*,' << xdmf2d.xmf'
    !
  end subroutine writexml
  !!
end module readwrite
!+---------------------------------------------------------------------+
!| The end of the module readwrite.                                    |
!+---------------------------------------------------------------------+