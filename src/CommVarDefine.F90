!+---------------------------------------------------------------------+
!| This module is used to define global variables and arraies.         |
!+---------------------------------------------------------------------+
!| Writen by Jian Fang, 2020-04-16                                     |
!+---------------------------------------------------------------------+
module commvardefine
  !
  implicit none
  !
  !
  integer :: im,jm,km
  !+-------------+-----------------------------------------------------+
  !|    im,jm,km | dimension of 3-D array                              |
  !+-------------+-----------------------------------------------------+
  !
  integer :: flowtype
  !+-------------+---------------------------------------------------=-+
  !|     flowtype| The variable in ASTR used to define type of flow    |
  !+-------------+----------------------------------------------------=+
  !
  logical :: lihomo,ljhomo,lkhomo
  !+-------------+-----------------------------------------------------+
  !|      lihomo | The logical variables used to define homogeneous    |
  !|      ljhomo | directions of the computation domain.               |
  !|      lkhomo |                                                     |
  !+-------------+-----------------------------------------------------+
  !
  real(8) :: ref_t,reynolds,mach
  !+-------------+-----------------------------------------------------+
  !|        ref_t| The reference temperature of the simulation         |
  !+-------------+-----------------------------------------------------+
  !|     Reynolds| The reference Reynolds number of the simulation     |
  !+-------------+-----------------------------------------------------+
  !|         mach| The reference Mach number of the simulation         |
  !+-------------+-----------------------------------------------------+
  !
  integer :: isp,jsp,ksp
  integer :: nproc,rankmax
  integer :: isize,jsize,ksize,irkm,jrkm,krkm
  integer,allocatable,dimension(:) :: imp,jmp,kmp,iop,jop,kop,         &
                                      ink,jnk,knk
  !+-------------+-----------------------------------------------------+
  !|  isp,jsp,ksp| The effective domain size of the simulation.        |
  !+-------------+-----------------------------------------------------+
  !|        isize| The number of sub-domains in each direction         |
  !|        jsize| the parameter of domain split in parallel           |
  !|        ksize| computation                                         |
  !+-------------+-----------------------------------------------------+
  !|         irkm| irkm=isize-1                                        |
  !|         jrkm| jrkm=jsize-1                                        |
  !|         krkm| krkm=ksize-1                                        |
  !+-------------+-----------------------------------------------------+
  !|        nproc| number of processors.                               |
  !+-------------+-----------------------------------------------------+
  !|      rankmax| max number of ranks: rankmax=nproc-1                |
  !+-------------+-----------------------------------------------------+
  !|  ink,jnk,knk| serial number of MPI ranks                          |
  !+-------------+-----------------------------------------------------+
  !|  imp,jmp,kmp| the size of mesh in each sub-domain                 |
  !+-------------+-----------------------------------------------------+
  !|  iop,jop,kop| the root point in each subdomain.                   |
  !+-------------+-----------------------------------------------------+
  !
  real(8) :: deltat
  !+-------------+-----------------------------------------------------+
  !|       deltat| The time step of the simulation                     |
  !+-------------+-----------------------------------------------------+
  !
  real(8),allocatable,dimension(:,:,:) :: x,y,z,ro,u1,u2,u3,p,t
  real(8),allocatable,dimension(:,:,:) :: ro_m,u1_m,u2_m,u3_m,p_m,t_m
  real(8),allocatable,dimension(:,:)   :: ro_zm,u1_zm,u2_zm,u3_zm,     &
                                          p_zm,t_zm
  real(8),allocatable,dimension(:)     :: ro_xzm,u1_xzm,u2_xzm,u3_xzm, &
                                          p_xzm,t_xzm
  real(8),allocatable,dimension(:,:,:,:) :: du1_m,du2_m,du3_m,dt_m 
  real(8),allocatable,dimension(:,:,:) :: du1_zm,du2_zm,du3_zm,dt_zm
  real(8),allocatable,dimension(:,:) :: du1_xzm,du2_xzm,du3_xzm,dt_xzm
  !+-------------+-----------------------------------------------------+
  !|       x,y,z | 3-D spatial coordinates                             |
  !+-------------+-----------------------------------------------------+
  !|          ro | density                                             |
  !|    u1,u2,u3 | velocity vector                                     |
  !|           p | pressure                                            |
  !|           t | temperature                                         |
  !+-------------+-----------------------------------------------------+
  !|          *_m| Mean flow statistics in 3D space.                   |
  !+-------------+-----------------------------------------------------+
  !|         *_zm| Mean flow statistics in x-y space.                  |
  !+-------------+-----------------------------------------------------+
  !|         *_xm| Mean flow statistics in y-z space.                  |
  !+-------------+-----------------------------------------------------+
  !|        *_xzm| Mean flow statistics as y proilfes.                 |
  !+-------------+-----------------------------------------------------+
  !|  du1,du2,du3| Velocity gradients martix.                          |
  !+-------------+-----------------------------------------------------+
  !|           dt| Temperature gradients martix.                       |
  !+-------------+-----------------------------------------------------+
  !
  real(8) :: prandtl,gamma,rgas,ref_miu,tempconst,tempconst1,const1,   &
             const2,const3,const4,const5,const6,const7,const8,roinf,   &
             uinf,vinf,tinf,pinf
  !+-------------+-----------------------------------------------------+
  !|         rgas| Gras constant.                                      |
  !+-------------+-----------------------------------------------------+
  !|      prandtl|Prandtl number                                       |
  !+-------------+-----------------------------------------------------+
  !|        gamma|gas capticity ratio                                  |
  !+-------------+-----------------------------------------------------+
  !|      ref_miu|reference viscosity                                  |
  !+-------------+-----------------------------------------------------+
  !|   roinf,uinf| flow variales at the infinite incoming flow.        |
  !|    tinf,pinf|                                                     |
  !+-------------+-----------------------------------------------------+
  !|    tempconst| parameters used to calculate viscosity.             |
  !|   tempconst1|                                                     |
  !+-------------+-----------------------------------------------------+
  !|       const1| const1=gamma*(gamma-1.d0)*mach**2                   |
  !|       const2| const2=gamma*mach**2                                |
  !|       const3| const3=(gamma-1.d0)/3.d0*prandtl*(mach**2)          |
  !|       const4| const4=(gamma-1.d0)*mach**2*reynolds*prandtl        |
  !|       const5| const5=(gamma-1.d0)*mach**2                         |
  !|       const6| const6=gamma-1.d0                                   |
  !|       const7| const7=(gamma-1.d0)*mach**2*Reynolds*prandtl        |
  !|       const8| const8=1.d0/Reynolds                                |
  !+-------------+-----------------------------------------------------+
  !
  !+-------------+-----------------------------------------------------+
  !|          pi | the ratio of a circle's circumference to diamete    |
  !+-------------+-----------------------------------------------------+
  character(len=64) :: gridfile
  character(len=7) :: outfolder
  !
  real(8),parameter :: pi=4.d0*atan(1.0_8),                            &
                       num1d35 =1.d0/35.d0,  num1d3  =1.d0/3.d0,       &
                       num2d3  =2.d0/3.d0,   num1d24 =1.d0/24.d0,      &
                       num4d3  =4.d0/3.d0,   num1d6  =1.d0/6.d0,       &
                       num1d12 =1.d0/12.d0,  num7d12 =7.d0/12.d0,      &
                       num7d9  =7.d0/9.d0,   num1d36 =1.d0/36.d0,      &
                       num1d60 =1.d0/60.d0,  num65d3 =65.d0/3.d0,      &
                       num20d3 =20.d0/3.d0,  num1d11 =1.d0/11.d0,      &
                       num25d12=25.d0/12.d0, num11d6 =11.d0/6.d0,      &
                       num1d840=1.d0/840.d0, num13d60=13.d0/60.d0,     &
                       num1d30 =1.d0/30.d0,  num47d60=47.d0/60.d0,     &
                       num5d6  =5.d0/6.d0,   num1d18 =1.d0/18.d0,      &
                       num19d18=19.d0/18.d0, num5d9  =5.d0/9.d0,       &
                       num9d36 =9.d0/36.d0
  !
  type :: block
    integer :: im,jm,km
    real(8),allocatable :: x(:,:,:),y(:,:,:),z(:,:,:)
  end type block
  !+-------------+-----------------------------------------------------+
  !|       block | multi-block class                                   |
  !+-------------+-----------------------------------------------------+
  !!
  !
end module commvardefine
!+---------------------------------------------------------------------+
!| The end of the module commvardefine                                 |
!+---------------------------------------------------------------------+