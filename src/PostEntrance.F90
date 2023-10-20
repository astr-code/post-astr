!+---------------------------------------------------------------------+
!| This module is the main entrance for post-process.                  |
!+---------------------------------------------------------------------+
!| Writen by Jian Fang, 2020-04-16                                     |
!+---------------------------------------------------------------------+
module postentrance
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to direct the program to the specific process. |
  !+-------------------------------------------------------------------+
  subroutine gotoprocess
    !
    use controlfileread
    use fieldview
    use readwrite
    use statistic
    use boundarylayer
    use dataprocess
    use userdefine
    use gridgen
    use flowanalyse
    use preprocess
    use h5readwrite
    use numerics
    use WingDesign
    ! use cgnstest
    !$ use OMP_LIB
    !
    ! local data
    !
    integer :: npost
    integer :: numb1,numb2,numb3
    character(len=255) :: string1,string2,string3
    integer :: firstnumb,lastnumb
    character(len=255) :: inputfile,outputfile
    integer :: nthread,mthread,ncore
    !
    !+-----------+
    !| * begin * |
    !+-----------+
    !
    call readkeyboard(npost,command=string3,number1=numb1,number2=numb2,number3=numb3, &
                      filename1=string1,filename2=string2,inputfile=inputfile,         &
                      nthread=nthread)
    !
    ! nthread=128
    !
    print*,' ** npost    : ',npost
    print*,' ** command  : ',trim(string3)
    print*,' ** number1  : ',numb1
    print*,' ** number2  : ',numb2
    print*,' ** number3  : ',numb3
    print*,' ** filename1: ',trim(string1)
    print*,' ** filename2: ',trim(string2)
    print*,' ** inputfile: ',trim(inputfile)
    print*,' ** nthread  : ',nthread
    !
    print*,' ** post command: ',trim(string3)
    !
    !$ write(*,'(A,I0,A)')'  ** The program is working in OpenMP with ',nthread,' threads'
    !$ call OMP_SET_DYNAMIC(.false.)
    !$ if (OMP_GET_DYNAMIC()) print *,' **!! Warning: dynamic adjustment of threads has been set!'
    !$ call OMP_SET_NUM_THREADS(nthread)
    !$ ncore=OMP_get_num_procs()
    !$ if(nthread>ncore) print *,' **!! Warning: the thread No. is larger than processers No.!'
    !
    !$ mthread=OMP_get_thread_num()
    !
    if(npost==-5) then
      !
      ! call wing2d_design
      ! call delta_wing_design
      !
      stop
      !
    endif
    !
    if(npost==-4) then
      !
      ! call tensoropt
      call wave_propa
      ! call errana1d_multi_domain
      !
      stop
      !
    endif
    !
    if(npost==-3) then
      ! call cgnswrite_grid_str
      ! call cgnread_grid_str
      ! call tgvstate(trim(string1))
      !
      call uprofilesym
      ! call filter_test
      ! call canteratest
      ! call expdatint
      ! call mijgen
      ! call momread
      ! call h5rcw(trim(string1),trim(string2),numb1)
      !
      ! call tgverror(trim(string1))
      ! call particle_gen
      !
      stop
      !
    endif
    !
    if(npost==-2) then
      call readpepin
      !
      call refconst
      !
      ! call onedflame_process(trim(string1))
      ! call onedflame_turb_int(trim(string1))
      !
      ! call inletmixing
      ! call lamianrBL(trim(string1))
      ! call inletprofile(trim(string1))
      ! call inlet_supchan(trim(string1))
      ! call randomforce
      call LambChaplyginDipole
      !
      stop ' ** job is done **'
    endif
    !
    if(npost==-1) then
      !
      call readpepin
      !
      call refconst
      !
      call grid2dextract
      ! call doublewedge4
      !
      ! call grid_rescale
      !
      ! call gridchannelimhomo
      ! call gridcavity2
      ! call gridsandbox
      !
      ! call gridBLmvg
      !
      ! call gridcavity4
      ! call gridchannelimhomo
      ! call gridibcylinder
      ! call gridchaib
      ! call gridchamvg
      !
      ! call gridBL2
      ! call grid2h5
      ! call gridBL_simple
      ! call gridBL_extra_spang
      ! call gridswtbli
      ! call gridjet
      ! call gridswtbli
      ! call gridmixinglayer
      ! call gridmixinglayer2
      !
      ! call GridSimple(20.d0,2.d0*pi,2.d0*pi)
      ! call gridsimpcylin
      ! call gridchannel
      ! call GridComExp
      ! call GridComExp90deg
      ! call gridcyaxix
      ! call grid2dextract
      ! call gridspanext
      !
      stop ' ** job is done **'
    endif
    !
    ! 
    if(trim(inputfile)=='datin/input.3dp.dat') then
      call readinput
      !
      call refconst
    elseif(len(trim(inputfile))>5) then
      call readinput_astrr(trim(inputfile))
      !
      call refconst
    endif
    !
    call informout(npost)
    !
    select case (npost)
    !
    case (204)
      ! call psdtest
      call pwallmon
      ! call flowstat_mon_proc
      ! call plotmonpoint
      ! call flowstat_mon_resum
      ! call monitorprocess(numb1,numb2)
      !
      ! call monitordatacon(trim(string1))
      ! call monitor_points
      !
    case (206)
      call shockjump
    case (2)
      call writetec3d(trim(string1),trim(string2))
    case (210)
      if(numb1>=0 .and. numb2>=0 ) then
        call writeteccs(numb1,numb2)
      elseif((trim(string1).ne.' ') .and. (trim(string2).ne.' ')) then
        call writeteccs(trim(string1),trim(string2))
      endif
    case (211)
      if(numb1>=0 .and. numb2>=0 .and. numb3>=0) then
        call writetecxy(numb1,numb2,numb3)
      elseif((trim(string1).ne.' ') .and. (trim(string2).ne.' ')) then
        call writetecxy(trim(string1),trim(string2),numb1)
      endif
    case (212)
      if((trim(string1).ne.' ') .and. (trim(string2).ne.' ')) then
        call writetecxz(trim(string1),trim(string2),numb1)
      endif
    case (213)
      if((trim(string1).ne.' ') .and. (trim(string2).ne.' ')) then
        call writetecyz(trim(string1),trim(string2),numb1)
      endif
    case (214)
      if((trim(string1).ne.' ') .and. (trim(string2).ne.' ')) then
        call writeprofile(trim(string1),trim(string2))
      endif
    case (215)
      call techannel
    case (3)
      !
      call H5_ReadGrid
      !
      call meanflow_fav(numb1,numb2)
      ! call meanflow_rey(numb1,numb2)
      ! call betchovsta(numb1,numb2)
      ! call highorder_rey(numb1,numb2)
      !
      call dutmean
      !
      call highordermom
      !
    case (301)
      call chanmean
    case (302)
      call blmean(numb1)
    case (304)
      call cylinder
    case (306)
      call streamlinecal(numb1,numb2)
    case (307)
      call cavity(numb1)
    case (4)
      !
      call H5_ReadGrid
      !
      call budgetnorm
      !
    case (401)
      ! call TKEBudgetSta
      call TKEBudget_zm
    case (5)
      ! call flowinterp_swtbli
      ! call flowinterp_cartesian
      ! call flowinterpxy3d
      call flowinterpxy(trim(string1))
      ! call flowinterp3d
      ! call flowinterpz3d
      ! call mapislice
      ! call flowdupz(trim(string1))
    case (501)
      call inletgen
    case (6)
      call pwcorrelation(numb1,numb2)
    case (7)
      !
      if(numb1>=0 .and. numb2>=0 ) then
        call spec_homo(numb1,numb2)
      elseif((trim(string1).ne.' ') .and. (trim(string2).ne.' ')) then
        call spec_homo_onefile(trim(string1),trim(string2))
      endif
      !
    case (801)
      call ksliceview(numb1,numb2)
    case (802)
      ! call isdatabase(numb1,numb2,numb3)
      ! call inflowpro(numb1,numb2)
      ! call isliceinterp(numb1,numb2,trim(string1))
      call isliceinterp_chan(numb1,numb2,trim(string1))
      ! call isliceview(numb1,numb2,numb3)
      ! call isliceanalys(numb1,numb2,mode=trim(string1))
    case (803)
      call jsliceview(numb1,numb2,numb3)
    case (901)
      call islice2h5(numb1,numb2)
    case (902)
      call kslice2h5(numb1,numb2)
    case (908)
      call jslice2h5(numb1,numb2)
    case (903)
      call meanflow_datacon
      !
      call budget_datacon
    case (1000)
      !
      call dataextra(numb1,numb2)
      ! call flowdataclean(numb1,numb2)
      ! call momread
      ! call mongen
      ! call betchovnorm
      ! call momtana(numb1,numb2,numb3)
      ! call vortsta(numb1,numb2)
      !
      ! call cavity
      ! call shockexpcom
      ! call comexpsl
      ! call comexp3d
      ! call comexp(numb1)
      ! call dis_streamline_wall
      ! call supplymodify
      ! call datacorase
      ! call nbody
      ! call mixinglayer(numb1)
      ! call writevelofied(trim(string1),trim(string2))
      ! call skewness_cal(trim(string1))
      ! call MVG3D
      ! call ijk2intone
      ! call writexml
      !
      ! call flowrecons
      !
    case (1002)
      ! call shearcore
      ! call swtbli(numb1)
      ! call wall3d
      ! call sonicline
    case (1007)
      call particle_stat(numb1,numb2)
    case (1008)
      call incompact3d_visu(numb1,numb2,trim(string1),trim(string2))
      ! call incompact3d(numb1,numb2,trim(string1),trim(string2))
    case default
      print*,' !! post-process not define !!'
      stop
    end select
    !
    print*,' +-----------------------------+'
    print*,' |  post-process job is done.  |'
    print*,' +-----------------------------+'
    !+---------+
    !| * end * |
    !+---------+
    !
  end subroutine gotoprocess
  !+-------------------------------------------------------------------+
  !| The end of the subroutine gotoproces.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !|this subroutine calculate some reference constants.                |
  !+-------------------------------------------------------------------+
  !|writen by fang jian, 2008-09-18.                                   |
  !+-------------------------------------------------------------------+
  subroutine refconst
    !
    use commvardefine
    !
    prandtl=0.72d0
    gamma=1.4d0
    rgas=287.1d0
    !
    ref_miu=1.458d-6*ref_t**1.5d0/(ref_t+110.4d0)
    !
    tempconst=110.4d0/ref_t
    tempconst1=1.d0+tempconst
    !
    !
    const1=gamma*(gamma-1.d0)*mach**2
    const2=gamma*mach**2
    const3=(gamma-1.d0)/3.d0*prandtl*(mach**2)
    const4=(gamma-1.d0)*mach**2*reynolds*prandtl
    const5=(gamma-1.d0)*mach**2
    const6=gamma-1.d0
    const7=(gamma-1.d0)*mach**2*reynolds*prandtl
    const8=1.d0/reynolds
    !
    uinf=1.d0
    vinf=0.d0
    tinf=1.d0
    roinf=1.d0
    !
    pinf=roinf*tinf/const2
    !
    print*,' ** Reference parameters calculation done.'
    !
  end subroutine refconst
  !+-------------------------------------------------------------------+
  !|end of the subroutine refconst.                                    |
  !+-------------------------------------------------------------------+
  !!
end module postentrance
!+---------------------------------------------------------------------+
!| The end of the module postentrance.                                 |
!+---------------------------------------------------------------------+