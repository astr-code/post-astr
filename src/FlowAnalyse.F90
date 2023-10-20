!+---------------------------------------------------------------------+
!| This module contains subroutines to get statistics of different flow|
!+---------------------------------------------------------------------+
!| Writen by Jian Fang, 2020-04-16                                     |
!+---------------------------------------------------------------------+
module flowanalyse
  !
  use readwrite
  use h5readwrite
  use gradsolver
  use writetec
  use basicfunction
  !
  implicit none
  !
  contains
  !
  subroutine betchovnorm
    !
    use commvardefine,only: im,jm,km,x,y,z,Reynolds,Mach,lihomo,lkhomo

    use gradsolver
    use basicfunction, only: integration
    !
    real(8) :: m_iim_jj(0:im,0:jm,0:km),m_ijm_ji(0:im,0:jm,0:km),      &
               m_ijm_ij(0:im,0:jm,0:km),m_iim_jjm_kk(0:im,0:jm,0:km),  &
               m_ijm_jim_kk(0:im,0:jm,0:km),                           &
               m_ijm_ijm_kk(0:im,0:jm,0:km),                           &
               m_ijm_jkm_ki(0:im,0:jm,0:km),                           &
               m_ijm_kjm_ki(0:im,0:jm,0:km)
    real(8) :: d11_d11(0:im,0:jm,0:km),d11_d12(0:im,0:jm,0:km),        &
               d11_d13(0:im,0:jm,0:km),d11_d21(0:im,0:jm,0:km),        &
               d11_d22(0:im,0:jm,0:km),d11_d23(0:im,0:jm,0:km),        &
               d11_d31(0:im,0:jm,0:km),d11_d32(0:im,0:jm,0:km),        &
               d11_d33(0:im,0:jm,0:km),d12_d12(0:im,0:jm,0:km),        &
               d12_d13(0:im,0:jm,0:km),d12_d21(0:im,0:jm,0:km),        &
               d12_d22(0:im,0:jm,0:km),d12_d23(0:im,0:jm,0:km),        &
               d12_d31(0:im,0:jm,0:km),d12_d32(0:im,0:jm,0:km),        &
               d12_d33(0:im,0:jm,0:km),d13_d13(0:im,0:jm,0:km),        &
               d13_d21(0:im,0:jm,0:km),d13_d22(0:im,0:jm,0:km),        &
               d13_d23(0:im,0:jm,0:km),d13_d31(0:im,0:jm,0:km),        &
               d13_d32(0:im,0:jm,0:km),d13_d33(0:im,0:jm,0:km),        &
               d21_d21(0:im,0:jm,0:km),d21_d22(0:im,0:jm,0:km),        &
               d21_d23(0:im,0:jm,0:km),d21_d31(0:im,0:jm,0:km),        &
               d21_d32(0:im,0:jm,0:km),d21_d33(0:im,0:jm,0:km),        &
               d22_d22(0:im,0:jm,0:km),d22_d23(0:im,0:jm,0:km),        &
               d22_d31(0:im,0:jm,0:km),d22_d32(0:im,0:jm,0:km),        &
               d22_d33(0:im,0:jm,0:km),d23_d23(0:im,0:jm,0:km),        &
               d23_d31(0:im,0:jm,0:km),d23_d32(0:im,0:jm,0:km),        &
               d23_d33(0:im,0:jm,0:km),d31_d31(0:im,0:jm,0:km),        &
               d31_d32(0:im,0:jm,0:km),d31_d33(0:im,0:jm,0:km),        &
               d32_d32(0:im,0:jm,0:km),d32_d33(0:im,0:jm,0:km),        &
               d33_d33(0:im,0:jm,0:km),m11(0:im,0:jm,0:km),            &
               m12(0:im,0:jm,0:km),m13(0:im,0:jm,0:km),                &
               m21(0:im,0:jm,0:km),m22(0:im,0:jm,0:km),                &
               m23(0:im,0:jm,0:km),m31(0:im,0:jm,0:km),                &
               m32(0:im,0:jm,0:km),m33(0:im,0:jm,0:km)
    !
    real(8) :: m_iim_jj_xzm(0:jm),m_ijm_ji_xzm(0:jm),                  &
               m_ijm_ij_xzm(0:jm),m_iim_jjm_kk_xzm(0:jm),              &
               m_ijm_jim_kk_xzm(0:jm),m_ijm_ijm_kk_xzm(0:jm),          &
               m_ijm_jkm_ki_xzm(0:jm),m_ijm_kjm_ki_xzm(0:jm),          &
               d11_d11_xzm(0:jm),d11_d12_xzm(0:jm),                    &
               d11_d13_xzm(0:jm),d11_d21_xzm(0:jm),                    &
               d11_d22_xzm(0:jm),d11_d23_xzm(0:jm),                    &
               d11_d31_xzm(0:jm),d11_d32_xzm(0:jm),                    &
               d11_d33_xzm(0:jm),d12_d12_xzm(0:jm),                    &
               d12_d13_xzm(0:jm),d12_d21_xzm(0:jm),                    &
               d12_d22_xzm(0:jm),d12_d23_xzm(0:jm),                    &
               d12_d31_xzm(0:jm),d12_d32_xzm(0:jm),                    &
               d12_d33_xzm(0:jm),d13_d13_xzm(0:jm),                    &
               d13_d21_xzm(0:jm),d13_d22_xzm(0:jm),                    &
               d13_d23_xzm(0:jm),d13_d31_xzm(0:jm),                    &
               d13_d32_xzm(0:jm),d13_d33_xzm(0:jm),                    &
               d21_d21_xzm(0:jm),d21_d22_xzm(0:jm),                    &
               d21_d23_xzm(0:jm),d21_d31_xzm(0:jm),                    &
               d21_d32_xzm(0:jm),d21_d33_xzm(0:jm),                    &
               d22_d22_xzm(0:jm),d22_d23_xzm(0:jm),                    &
               d22_d31_xzm(0:jm),d22_d32_xzm(0:jm),                    &
               d22_d33_xzm(0:jm),d23_d23_xzm(0:jm),                    &
               d23_d31_xzm(0:jm),d23_d32_xzm(0:jm),                    &
               d23_d33_xzm(0:jm),d31_d31_xzm(0:jm),                    &
               d31_d32_xzm(0:jm),d31_d33_xzm(0:jm),                    &
               d32_d32_xzm(0:jm),d32_d33_xzm(0:jm),                    &
               d33_d33_xzm(0:jm),m11_xzm(0:jm),                        &
               m12_xzm(0:jm),m13_xzm(0:jm),                            &
               m21_xzm(0:jm),m22_xzm(0:jm),                            &
               m23_xzm(0:jm),m31_xzm(0:jm),                            &
               m32_xzm(0:jm),m33_xzm(0:jm)
    real(8) :: m_iim_jj_zm(0:im,0:jm),m_ijm_ji_zm(0:im,0:jm),          &
               m_ijm_ij_zm(0:im,0:jm),m_iim_jjm_kk_zm(0:im,0:jm),      &
               m_ijm_jim_kk_zm(0:im,0:jm),m_ijm_ijm_kk_zm(0:im,0:jm),  &
               m_ijm_jkm_ki_zm(0:im,0:jm),m_ijm_kjm_ki_zm(0:im,0:jm),  &
               d11_d11_zm(0:im,0:jm),d11_d12_zm(0:im,0:jm),            &
               d11_d13_zm(0:im,0:jm),d11_d21_zm(0:im,0:jm),            &
               d11_d22_zm(0:im,0:jm),d11_d23_zm(0:im,0:jm),            &
               d11_d31_zm(0:im,0:jm),d11_d32_zm(0:im,0:jm),            &
               d11_d33_zm(0:im,0:jm),d12_d12_zm(0:im,0:jm),            &
               d12_d13_zm(0:im,0:jm),d12_d21_zm(0:im,0:jm),            &
               d12_d22_zm(0:im,0:jm),d12_d23_zm(0:im,0:jm),            &
               d12_d31_zm(0:im,0:jm),d12_d32_zm(0:im,0:jm),            &
               d12_d33_zm(0:im,0:jm),d13_d13_zm(0:im,0:jm),            &
               d13_d21_zm(0:im,0:jm),d13_d22_zm(0:im,0:jm),            &
               d13_d23_zm(0:im,0:jm),d13_d31_zm(0:im,0:jm),            &
               d13_d32_zm(0:im,0:jm),d13_d33_zm(0:im,0:jm),            &
               d21_d21_zm(0:im,0:jm),d21_d22_zm(0:im,0:jm),            &
               d21_d23_zm(0:im,0:jm),d21_d31_zm(0:im,0:jm),            &
               d21_d32_zm(0:im,0:jm),d21_d33_zm(0:im,0:jm),            &
               d22_d22_zm(0:im,0:jm),d22_d23_zm(0:im,0:jm),            &
               d22_d31_zm(0:im,0:jm),d22_d32_zm(0:im,0:jm),            &
               d22_d33_zm(0:im,0:jm),d23_d23_zm(0:im,0:jm),            &
               d23_d31_zm(0:im,0:jm),d23_d32_zm(0:im,0:jm),            &
               d23_d33_zm(0:im,0:jm),d31_d31_zm(0:im,0:jm),            &
               d31_d32_zm(0:im,0:jm),d31_d33_zm(0:im,0:jm),            &
               d32_d32_zm(0:im,0:jm),d32_d33_zm(0:im,0:jm),            &
               d33_d33_zm(0:im,0:jm),m11_zm(0:im,0:jm),                &
               m12_zm(0:im,0:jm),m13_zm(0:im,0:jm),                    &
               m21_zm(0:im,0:jm),m22_zm(0:im,0:jm),                    &
               m23_zm(0:im,0:jm),m31_zm(0:im,0:jm),                    &
               m32_zm(0:im,0:jm),m33_zm(0:im,0:jm)
    real(8) :: div,d11,d12,d13,d21,d22,d23,d31,d32,d33, &
               d11d11,d11d12,d11d13,d11d21,d11d22,d11d23,d11d31,d11d32,&
               d11d33,d12d12,d12d13,d12d21,d12d22,d12d23,d12d31,d12d32,&
               d12d33,d13d13,d13d21,d13d22,d13d23,d13d31,d13d32,d13d33,&
               d21d21,d21d22,d21d23,d21d31,d21d32,d21d33,d22d22,d22d23,&
               d22d31,d22d32,d22d33,d23d23,d23d31,d23d32,d23d33,d31d31,&
               d31d32,d31d33,d32d32,d32d33,d33d33
    real(8) :: var1,var2,var3,var4
    !
    integer :: nsamples,jmm,n,i,j,k
    logical :: lfilalive
    !
    ! call H5ReadArray(u1,im,jm,km,'u1rem','Results/budget.h5')
    ! call H5ReadArray(u2,im,jm,km,'u2rem','Results/budget.h5')
    ! call H5ReadArray(u3,im,jm,km,'u3rem','Results/budget.h5')
    !
    call H5ReadArray(nsamples,'nsamples','Outdat/betchov.h5')
    call H5ReadArray(m_iim_jj,    im,jm,km,    'm_iim_jj','Outdat/betchov.h5')
    call H5ReadArray(m_ijm_ji,    im,jm,km,    'm_ijm_ji','Outdat/betchov.h5')
    call H5ReadArray(m_ijm_ij,    im,jm,km,    'm_ijm_ij','Outdat/betchov.h5')
    call H5ReadArray(m_iim_jjm_kk,im,jm,km,'m_iim_jjm_kk','Outdat/betchov.h5')
    call H5ReadArray(m_ijm_jim_kk,im,jm,km,'m_ijm_jim_kk','Outdat/betchov.h5')
    call H5ReadArray(m_ijm_ijm_kk,im,jm,km,'m_ijm_ijm_kk','Outdat/betchov.h5')
    call H5ReadArray(m_ijm_jkm_ki,im,jm,km,'m_ijm_jkm_ki','Outdat/betchov.h5')
    call H5ReadArray(m_ijm_kjm_ki,im,jm,km,'m_ijm_kjm_ki','Outdat/betchov.h5')
    !
    call H5ReadArray(d11_d11,im,jm,km,'d11_d11','Outdat/betchov.h5')
    call H5ReadArray(d11_d12,im,jm,km,'d11_d12','Outdat/betchov.h5')
    call H5ReadArray(d11_d13,im,jm,km,'d11_d13','Outdat/betchov.h5')
    call H5ReadArray(d11_d21,im,jm,km,'d11_d21','Outdat/betchov.h5')
    call H5ReadArray(d11_d22,im,jm,km,'d11_d22','Outdat/betchov.h5')
    call H5ReadArray(d11_d23,im,jm,km,'d11_d23','Outdat/betchov.h5')
    call H5ReadArray(d11_d31,im,jm,km,'d11_d31','Outdat/betchov.h5')
    call H5ReadArray(d11_d32,im,jm,km,'d11_d32','Outdat/betchov.h5')
    call H5ReadArray(d11_d33,im,jm,km,'d11_d33','Outdat/betchov.h5')
    call H5ReadArray(d12_d12,im,jm,km,'d12_d12','Outdat/betchov.h5')
    call H5ReadArray(d12_d13,im,jm,km,'d12_d13','Outdat/betchov.h5')
    call H5ReadArray(d12_d21,im,jm,km,'d12_d21','Outdat/betchov.h5')
    call H5ReadArray(d12_d22,im,jm,km,'d12_d22','Outdat/betchov.h5')
    call H5ReadArray(d12_d23,im,jm,km,'d12_d23','Outdat/betchov.h5')
    call H5ReadArray(d12_d31,im,jm,km,'d12_d31','Outdat/betchov.h5')
    call H5ReadArray(d12_d32,im,jm,km,'d12_d32','Outdat/betchov.h5')
    call H5ReadArray(d12_d33,im,jm,km,'d12_d33','Outdat/betchov.h5')
    call H5ReadArray(d13_d13,im,jm,km,'d13_d13','Outdat/betchov.h5')
    call H5ReadArray(d13_d21,im,jm,km,'d13_d21','Outdat/betchov.h5')
    call H5ReadArray(d13_d22,im,jm,km,'d13_d22','Outdat/betchov.h5')
    call H5ReadArray(d13_d23,im,jm,km,'d13_d23','Outdat/betchov.h5')
    call H5ReadArray(d13_d31,im,jm,km,'d13_d31','Outdat/betchov.h5')
    call H5ReadArray(d13_d32,im,jm,km,'d13_d32','Outdat/betchov.h5')
    call H5ReadArray(d13_d33,im,jm,km,'d13_d33','Outdat/betchov.h5')
    call H5ReadArray(d21_d21,im,jm,km,'d21_d21','Outdat/betchov.h5')
    call H5ReadArray(d21_d22,im,jm,km,'d21_d22','Outdat/betchov.h5')
    call H5ReadArray(d21_d23,im,jm,km,'d21_d23','Outdat/betchov.h5')
    call H5ReadArray(d21_d31,im,jm,km,'d21_d31','Outdat/betchov.h5')
    call H5ReadArray(d21_d32,im,jm,km,'d21_d32','Outdat/betchov.h5')
    call H5ReadArray(d21_d33,im,jm,km,'d21_d33','Outdat/betchov.h5')
    call H5ReadArray(d22_d22,im,jm,km,'d22_d22','Outdat/betchov.h5')
    call H5ReadArray(d22_d23,im,jm,km,'d22_d23','Outdat/betchov.h5')
    call H5ReadArray(d22_d31,im,jm,km,'d22_d31','Outdat/betchov.h5')
    call H5ReadArray(d22_d32,im,jm,km,'d22_d32','Outdat/betchov.h5')
    call H5ReadArray(d22_d33,im,jm,km,'d22_d33','Outdat/betchov.h5')
    call H5ReadArray(d23_d23,im,jm,km,'d23_d23','Outdat/betchov.h5')
    call H5ReadArray(d23_d31,im,jm,km,'d23_d31','Outdat/betchov.h5')
    call H5ReadArray(d23_d32,im,jm,km,'d23_d32','Outdat/betchov.h5')
    call H5ReadArray(d23_d33,im,jm,km,'d23_d33','Outdat/betchov.h5')
    call H5ReadArray(d31_d31,im,jm,km,'d31_d31','Outdat/betchov.h5')
    call H5ReadArray(d31_d32,im,jm,km,'d31_d32','Outdat/betchov.h5')
    call H5ReadArray(d31_d33,im,jm,km,'d31_d33','Outdat/betchov.h5')
    call H5ReadArray(d32_d32,im,jm,km,'d32_d32','Outdat/betchov.h5')
    call H5ReadArray(d32_d33,im,jm,km,'d32_d33','Outdat/betchov.h5')
    call H5ReadArray(d33_d33,im,jm,km,'d33_d33','Outdat/betchov.h5')
    call H5ReadArray(m11    ,im,jm,km,    'm11','Outdat/betchov.h5')     
    call H5ReadArray(m12    ,im,jm,km,    'm12','Outdat/betchov.h5')     
    call H5ReadArray(m13    ,im,jm,km,    'm13','Outdat/betchov.h5')     
    call H5ReadArray(m21    ,im,jm,km,    'm21','Outdat/betchov.h5')     
    call H5ReadArray(m22    ,im,jm,km,    'm22','Outdat/betchov.h5')     
    call H5ReadArray(m23    ,im,jm,km,    'm23','Outdat/betchov.h5')     
    call H5ReadArray(m31    ,im,jm,km,    'm31','Outdat/betchov.h5')     
    call H5ReadArray(m32    ,im,jm,km,    'm32','Outdat/betchov.h5')     
    call H5ReadArray(m33    ,im,jm,km,    'm33','Outdat/betchov.h5')
    !
    print*,' ** nsamples=',nsamples
    !
    if(lihomo .and. lkhomo) then
      !
      jmm=jm/2
      !
      call H5_ReadGrid
      !
      ! du1=grad(u1,x,y,z) 
      ! du2=grad(u2)
      ! du3=grad(u3)
      !
      do j=0,jm
        !
            m_iim_jj_xzm(j) = integration(    m_iim_jj(:,j,:))     
            m_ijm_ji_xzm(j) = integration(    m_ijm_ji(:,j,:))      
            m_ijm_ij_xzm(j) = integration(    m_ijm_ij(:,j,:))     
        m_iim_jjm_kk_xzm(j) = integration(m_iim_jjm_kk(:,j,:)) 
        m_ijm_jim_kk_xzm(j) = integration(m_ijm_jim_kk(:,j,:)) 
        m_ijm_ijm_kk_xzm(j) = integration(m_ijm_ijm_kk(:,j,:)) 
        m_ijm_jkm_ki_xzm(j) = integration(m_ijm_jkm_ki(:,j,:)) 
        m_ijm_kjm_ki_xzm(j) = integration(m_ijm_kjm_ki(:,j,:)) 
             d11_d11_xzm(j) = integration(     d11_d11(:,j,:)) 
             d11_d12_xzm(j) = integration(     d11_d12(:,j,:)) 
             d11_d13_xzm(j) = integration(     d11_d13(:,j,:)) 
             d11_d21_xzm(j) = integration(     d11_d21(:,j,:)) 
             d11_d22_xzm(j) = integration(     d11_d22(:,j,:)) 
             d11_d23_xzm(j) = integration(     d11_d23(:,j,:)) 
             d11_d31_xzm(j) = integration(     d11_d31(:,j,:)) 
             d11_d32_xzm(j) = integration(     d11_d32(:,j,:)) 
             d11_d33_xzm(j) = integration(     d11_d33(:,j,:)) 
             d12_d12_xzm(j) = integration(     d12_d12(:,j,:)) 
             d12_d13_xzm(j) = integration(     d12_d13(:,j,:)) 
             d12_d21_xzm(j) = integration(     d12_d21(:,j,:)) 
             d12_d22_xzm(j) = integration(     d12_d22(:,j,:)) 
             d12_d23_xzm(j) = integration(     d12_d23(:,j,:)) 
             d12_d31_xzm(j) = integration(     d12_d31(:,j,:)) 
             d12_d32_xzm(j) = integration(     d12_d32(:,j,:)) 
             d12_d33_xzm(j) = integration(     d12_d33(:,j,:)) 
             d13_d13_xzm(j) = integration(     d13_d13(:,j,:)) 
             d13_d21_xzm(j) = integration(     d13_d21(:,j,:)) 
             d13_d22_xzm(j) = integration(     d13_d22(:,j,:)) 
             d13_d23_xzm(j) = integration(     d13_d23(:,j,:)) 
             d13_d31_xzm(j) = integration(     d13_d31(:,j,:)) 
             d13_d32_xzm(j) = integration(     d13_d32(:,j,:)) 
             d13_d33_xzm(j) = integration(     d13_d33(:,j,:)) 
             d21_d21_xzm(j) = integration(     d21_d21(:,j,:)) 
             d21_d22_xzm(j) = integration(     d21_d22(:,j,:)) 
             d21_d23_xzm(j) = integration(     d21_d23(:,j,:)) 
             d21_d31_xzm(j) = integration(     d21_d31(:,j,:)) 
             d21_d32_xzm(j) = integration(     d21_d32(:,j,:)) 
             d21_d33_xzm(j) = integration(     d21_d33(:,j,:)) 
             d22_d22_xzm(j) = integration(     d22_d22(:,j,:)) 
             d22_d23_xzm(j) = integration(     d22_d23(:,j,:)) 
             d22_d31_xzm(j) = integration(     d22_d31(:,j,:)) 
             d22_d32_xzm(j) = integration(     d22_d32(:,j,:)) 
             d22_d33_xzm(j) = integration(     d22_d33(:,j,:)) 
             d23_d23_xzm(j) = integration(     d23_d23(:,j,:)) 
             d23_d31_xzm(j) = integration(     d23_d31(:,j,:)) 
             d23_d32_xzm(j) = integration(     d23_d32(:,j,:)) 
             d23_d33_xzm(j) = integration(     d23_d33(:,j,:)) 
             d31_d31_xzm(j) = integration(     d31_d31(:,j,:)) 
             d31_d32_xzm(j) = integration(     d31_d32(:,j,:)) 
             d31_d33_xzm(j) = integration(     d31_d33(:,j,:)) 
             d32_d32_xzm(j) = integration(     d32_d32(:,j,:)) 
             d32_d33_xzm(j) = integration(     d32_d33(:,j,:)) 
             d33_d33_xzm(j) = integration(     d33_d33(:,j,:)) 
                 m11_xzm(j) = integration(         m11(:,j,:))         
                 m12_xzm(j) = integration(         m12(:,j,:))         
                 m13_xzm(j) = integration(         m13(:,j,:))         
                 m21_xzm(j) = integration(         m21(:,j,:))         
                 m22_xzm(j) = integration(         m22(:,j,:))         
                 m23_xzm(j) = integration(         m23(:,j,:))         
                 m31_xzm(j) = integration(         m31(:,j,:))         
                 m32_xzm(j) = integration(         m32(:,j,:))         
                 m33_xzm(j) = integration(         m33(:,j,:))     
        !
        ! du1_xzm(1,j)=integration(du1(1,:,j,:))
        ! du1_xzm(2,j)=integration(du1(2,:,j,:))
        ! du1_xzm(3,j)=integration(du1(3,:,j,:))
        ! du2_xzm(1,j)=integration(du2(1,:,j,:))
        ! du2_xzm(2,j)=integration(du2(2,:,j,:))
        ! du2_xzm(3,j)=integration(du2(3,:,j,:))
        ! du3_xzm(1,j)=integration(du3(1,:,j,:))
        ! du3_xzm(2,j)=integration(du3(2,:,j,:))
        ! du3_xzm(3,j)=integration(du3(3,:,j,:))
        !
      enddo
      !    
          m_iim_jj_xzm =     m_iim_jj_xzm/dble(nsamples*im*km)     
          m_ijm_ji_xzm =     m_ijm_ji_xzm/dble(nsamples*im*km)      
          m_ijm_ij_xzm =     m_ijm_ij_xzm/dble(nsamples*im*km)     
      m_iim_jjm_kk_xzm = m_iim_jjm_kk_xzm/dble(nsamples*im*km) 
      m_ijm_jim_kk_xzm = m_ijm_jim_kk_xzm/dble(nsamples*im*km) 
      m_ijm_ijm_kk_xzm = m_ijm_ijm_kk_xzm/dble(nsamples*im*km) 
      m_ijm_jkm_ki_xzm = m_ijm_jkm_ki_xzm/dble(nsamples*im*km) 
      m_ijm_kjm_ki_xzm = m_ijm_kjm_ki_xzm/dble(nsamples*im*km) 
           d11_d11_xzm =      d11_d11_xzm/dble(nsamples*im*km)
           d11_d12_xzm =      d11_d12_xzm/dble(nsamples*im*km)
           d11_d13_xzm =      d11_d13_xzm/dble(nsamples*im*km)
           d11_d21_xzm =      d11_d21_xzm/dble(nsamples*im*km)
           d11_d22_xzm =      d11_d22_xzm/dble(nsamples*im*km)
           d11_d23_xzm =      d11_d23_xzm/dble(nsamples*im*km)
           d11_d31_xzm =      d11_d31_xzm/dble(nsamples*im*km)
           d11_d32_xzm =      d11_d32_xzm/dble(nsamples*im*km)
           d11_d33_xzm =      d11_d33_xzm/dble(nsamples*im*km)
           d12_d12_xzm =      d12_d12_xzm/dble(nsamples*im*km)
           d12_d13_xzm =      d12_d13_xzm/dble(nsamples*im*km)
           d12_d21_xzm =      d12_d21_xzm/dble(nsamples*im*km)
           d12_d22_xzm =      d12_d22_xzm/dble(nsamples*im*km)
           d12_d23_xzm =      d12_d23_xzm/dble(nsamples*im*km)
           d12_d31_xzm =      d12_d31_xzm/dble(nsamples*im*km)
           d12_d32_xzm =      d12_d32_xzm/dble(nsamples*im*km)
           d12_d33_xzm =      d12_d33_xzm/dble(nsamples*im*km)
           d13_d13_xzm =      d13_d13_xzm/dble(nsamples*im*km)
           d13_d21_xzm =      d13_d21_xzm/dble(nsamples*im*km)
           d13_d22_xzm =      d13_d22_xzm/dble(nsamples*im*km)
           d13_d23_xzm =      d13_d23_xzm/dble(nsamples*im*km)
           d13_d31_xzm =      d13_d31_xzm/dble(nsamples*im*km)
           d13_d32_xzm =      d13_d32_xzm/dble(nsamples*im*km)
           d13_d33_xzm =      d13_d33_xzm/dble(nsamples*im*km)
           d21_d21_xzm =      d21_d21_xzm/dble(nsamples*im*km)
           d21_d22_xzm =      d21_d22_xzm/dble(nsamples*im*km)
           d21_d23_xzm =      d21_d23_xzm/dble(nsamples*im*km)
           d21_d31_xzm =      d21_d31_xzm/dble(nsamples*im*km)
           d21_d32_xzm =      d21_d32_xzm/dble(nsamples*im*km)
           d21_d33_xzm =      d21_d33_xzm/dble(nsamples*im*km)
           d22_d22_xzm =      d22_d22_xzm/dble(nsamples*im*km)
           d22_d23_xzm =      d22_d23_xzm/dble(nsamples*im*km)
           d22_d31_xzm =      d22_d31_xzm/dble(nsamples*im*km)
           d22_d32_xzm =      d22_d32_xzm/dble(nsamples*im*km)
           d22_d33_xzm =      d22_d33_xzm/dble(nsamples*im*km)
           d23_d23_xzm =      d23_d23_xzm/dble(nsamples*im*km)
           d23_d31_xzm =      d23_d31_xzm/dble(nsamples*im*km)
           d23_d32_xzm =      d23_d32_xzm/dble(nsamples*im*km)
           d23_d33_xzm =      d23_d33_xzm/dble(nsamples*im*km)
           d31_d31_xzm =      d31_d31_xzm/dble(nsamples*im*km)
           d31_d32_xzm =      d31_d32_xzm/dble(nsamples*im*km)
           d31_d33_xzm =      d31_d33_xzm/dble(nsamples*im*km)
           d32_d32_xzm =      d32_d32_xzm/dble(nsamples*im*km)
           d32_d33_xzm =      d32_d33_xzm/dble(nsamples*im*km)
           d33_d33_xzm =      d33_d33_xzm/dble(nsamples*im*km)
               m11_xzm =          m11_xzm/dble(nsamples*im*km)        
               m12_xzm =          m12_xzm/dble(nsamples*im*km)        
               m13_xzm =          m13_xzm/dble(nsamples*im*km)        
               m21_xzm =          m21_xzm/dble(nsamples*im*km)        
               m22_xzm =          m22_xzm/dble(nsamples*im*km)        
               m23_xzm =          m23_xzm/dble(nsamples*im*km)        
               m31_xzm =          m31_xzm/dble(nsamples*im*km)        
               m32_xzm =          m32_xzm/dble(nsamples*im*km)        
               m33_xzm =          m33_xzm/dble(nsamples*im*km)
      !
      ! du1_xzm=du1_xzm/dble(im*km)
      ! du2_xzm=du2_xzm/dble(im*km)
      ! du3_xzm=du3_xzm/dble(im*km)
!  
      print*,' ** total number of samples:',(nsamples*im*km)
      !
      do j=0,jm
        !
        d11=m11_xzm(j);      d12=m12_xzm(j);      d13=m13_xzm(j)
        d21=m21_xzm(j);      d22=m22_xzm(j);      d23=m23_xzm(j)
        d31=m31_xzm(j);      d32=m32_xzm(j);      d33=m33_xzm(j)
        !
        d11d11=d11_d11_xzm(j)-d11*d11 
        d11d12=d11_d12_xzm(j)-d11*d12 
        d11d13=d11_d13_xzm(j)-d11*d13 
        d11d21=d11_d21_xzm(j)-d11*d21 
        d11d22=d11_d22_xzm(j)-d11*d22 
        d11d23=d11_d23_xzm(j)-d11*d23 
        d11d31=d11_d31_xzm(j)-d11*d31 
        d11d32=d11_d32_xzm(j)-d11*d32 
        d11d33=d11_d33_xzm(j)-d11*d33 
        d12d12=d12_d12_xzm(j)-d12*d12 
        d12d13=d12_d13_xzm(j)-d12*d13 
        d12d21=d12_d21_xzm(j)-d12*d21 
        d12d22=d12_d22_xzm(j)-d12*d22 
        d12d23=d12_d23_xzm(j)-d12*d23 
        d12d31=d12_d31_xzm(j)-d12*d31 
        d12d32=d12_d32_xzm(j)-d12*d32 
        d12d33=d12_d33_xzm(j)-d12*d33 
        d13d13=d13_d13_xzm(j)-d13*d13 
        d13d21=d13_d21_xzm(j)-d13*d21 
        d13d22=d13_d22_xzm(j)-d13*d22 
        d13d23=d13_d23_xzm(j)-d13*d23 
        d13d31=d13_d31_xzm(j)-d13*d31 
        d13d32=d13_d32_xzm(j)-d13*d32 
        d13d33=d13_d33_xzm(j)-d13*d33 
        d21d21=d21_d21_xzm(j)-d21*d21 
        d21d22=d21_d22_xzm(j)-d21*d22 
        d21d23=d21_d23_xzm(j)-d21*d23 
        d21d31=d21_d31_xzm(j)-d21*d31 
        d21d32=d21_d32_xzm(j)-d21*d32 
        d21d33=d21_d33_xzm(j)-d21*d33 
        d22d22=d22_d22_xzm(j)-d22*d22 
        d22d23=d22_d23_xzm(j)-d22*d23 
        d22d31=d22_d31_xzm(j)-d22*d31 
        d22d32=d22_d32_xzm(j)-d22*d32 
        d22d33=d22_d33_xzm(j)-d22*d33 
        d23d23=d23_d23_xzm(j)-d23*d23 
        d23d31=d23_d31_xzm(j)-d23*d31 
        d23d32=d23_d32_xzm(j)-d23*d32 
        d23d33=d23_d33_xzm(j)-d23*d33 
        d31d31=d31_d31_xzm(j)-d31*d31 
        d31d32=d31_d32_xzm(j)-d31*d32 
        d31d33=d31_d33_xzm(j)-d31*d33 
        d32d32=d32_d32_xzm(j)-d32*d32 
        d32d33=d32_d33_xzm(j)-d32*d33 
        d33d33=d33_d33_xzm(j)-d33*d33 
        !
        d11_d11_xzm(j)=d11d11 
        d11_d12_xzm(j)=d11d12 
        d11_d13_xzm(j)=d11d13 
        d11_d21_xzm(j)=d11d21 
        d11_d22_xzm(j)=d11d22 
        d11_d23_xzm(j)=d11d23 
        d11_d31_xzm(j)=d11d31 
        d11_d32_xzm(j)=d11d32 
        d11_d33_xzm(j)=d11d33 
        d12_d12_xzm(j)=d12d12 
        d12_d13_xzm(j)=d12d13 
        d12_d21_xzm(j)=d12d21 
        d12_d22_xzm(j)=d12d22 
        d12_d23_xzm(j)=d12d23 
        d12_d31_xzm(j)=d12d31 
        d12_d32_xzm(j)=d12d32 
        d12_d33_xzm(j)=d12d33 
        d13_d13_xzm(j)=d13d13 
        d13_d21_xzm(j)=d13d21 
        d13_d22_xzm(j)=d13d22 
        d13_d23_xzm(j)=d13d23 
        d13_d31_xzm(j)=d13d31 
        d13_d32_xzm(j)=d13d32 
        d13_d33_xzm(j)=d13d33 
        d21_d21_xzm(j)=d21d21 
        d21_d22_xzm(j)=d21d22 
        d21_d23_xzm(j)=d21d23 
        d21_d31_xzm(j)=d21d31 
        d21_d32_xzm(j)=d21d32 
        d21_d33_xzm(j)=d21d33 
        d22_d22_xzm(j)=d22d22 
        d22_d23_xzm(j)=d22d23 
        d22_d31_xzm(j)=d22d31 
        d22_d32_xzm(j)=d22d32 
        d22_d33_xzm(j)=d22d33 
        d23_d23_xzm(j)=d23d23 
        d23_d31_xzm(j)=d23d31 
        d23_d32_xzm(j)=d23d32 
        d23_d33_xzm(j)=d23d33 
        d31_d31_xzm(j)=d31d31 
        d31_d32_xzm(j)=d31d32 
        d31_d33_xzm(j)=d31d33 
        d32_d32_xzm(j)=d32d32 
        d32_d33_xzm(j)=d32d33 
        d33_d33_xzm(j)=d33d33 
        !
        div=d11+d22+d33
        !
        var1  = div*div
        m_iim_jj_xzm(j) = m_iim_jj_xzm(j) - var1
        !
        var1  = d11*d11 + d12*d21 + d13*d31 +     &
                d21*d12 + d22*d22 + d23*d32 +     &
                d31*d13 + d32*d23 + d33*d33
        m_ijm_ji_xzm(j) = m_ijm_ji_xzm(j) - var1
        !
        var1  = d11*d11 + d12*d12 + d13*d13 +     &
                d21*d21 + d22*d22 + d23*d23 +     &
                d31*d31 + d32*d32 + d33*d33
        m_ijm_ij_xzm(j) = m_ijm_ij_xzm(j) - var1
        !
        var1=div*div*div
        var2=d11d11+d11d22+d11d33 +  &
             d11d22+d22d22+d22d33 +  &
             d11d33+d22d33+d33d33  
        var3=(d11+d22+d33)*var2
        m_iim_jjm_kk_xzm(j) = m_iim_jjm_kk_xzm(j) - var1-3.d0*var3
        !
        var1=d11*d11 + d12*d21 + d13*d31 +     &
             d21*d12 + d22*d22 + d23*d32 +     &
             d31*d13 + d32*d23 + d33*d33
        var2=d11*(d11d11+d11d22+d11d33) +   &
             d12*(d11d21+d21d22+d21d33) +   &
             d13*(d11d31+d22d31+d31d33) +   &
             d21*(d11d12+d12d22+d12d33) +   &
             d22*(d11d22+d22d22+d22d33) +   &
             d23*(d11d32+d22d32+d32d33) +   &
             d31*(d11d13+d13d22+d13d33) +   &
             d32*(d11d23+d22d23+d23d33) +   &
             d33*(d11d33+d22d33+d33d33)
        var3=div*(d11d11+d12d21+d13d31  +   &
                  d12d21+d22d22+d23d32  +   &
                  d13d31+d23d32+d33d33 )
        m_ijm_jim_kk_xzm(j) = m_ijm_jim_kk_xzm(j) - var1*div - 2.d0*var2 -var3
        !
        var1  = d11*d11 + d12*d12 + d13*d13 +     &
                d21*d21 + d22*d22 + d23*d23 +     &
                d31*d31 + d32*d32 + d33*d33
        var2=d11*(d11d11+d11d22+d11d33) +   &
             d12*(d11d12+d12d22+d12d33) +   &
             d13*(d11d13+d13d22+d13d33) +   &
             d21*(d11d21+d21d22+d21d33) +   &
             d22*(d11d22+d22d22+d22d33) +   &
             d23*(d11d23+d22d23+d23d33) +   &
             d31*(d11d31+d22d31+d31d33) +   &
             d32*(d11d32+d22d32+d32d33) +   &
             d33*(d11d33+d22d33+d33d33)
        var3=div*(d11d11+d12d12+d13d13  +   &
                  d21d21+d22d22+d23d23  +   &
                  d31d31+d32d32+d33d33 )
        m_ijm_ijm_kk_xzm(j) = m_ijm_ijm_kk_xzm(j) - var1*div  - 2.d0*var2 -var3
        !
        var1= d11*d11*d11 + d11*d12*d21 + d11*d13*d31 +  &
              d12*d21*d11 + d12*d22*d21 + d12*d23*d31 +  &
              d13*d31*d11 + d13*d32*d21 + d13*d33*d31 +  &
              d21*d11*d12 + d21*d12*d22 + d21*d13*d32 +  &
              d22*d21*d12 + d22*d22*d22 + d22*d23*d32 +  &
              d23*d31*d12 + d23*d32*d22 + d23*d33*d32 +  &
              d31*d11*d13 + d31*d12*d23 + d31*d13*d33 +  &
              d32*d21*d13 + d32*d22*d23 + d32*d23*d33 +  &
              d33*d31*d13 + d33*d32*d23 + d33*d33*d33
        var2= d11*(d11d11+d12d21+d13d31) +  &
              d12*(d11d21+d21d22+d23d31) +  &
              d13*(d11d31+d21d32+d31d33) +  &
              d21*(d11d12+d12d22+d13d32) +  &
              d22*(d12d21+d22d22+d23d32) +  &
              d23*(d12d31+d22d32+d32d33) +  &
              d31*(d11d13+d12d23+d13d33) +  &
              d32*(d13d21+d22d23+d23d33) +  &
              d33*(d13d31+d23d32+d33d33)
        m_ijm_jkm_ki_xzm(j) = m_ijm_jkm_ki_xzm(j) - var1 -3.d0*var2
        !
        var1=d11*(d11*d11+d21*d21+d31*d31) + &
             d12*(d12*d11+d22*d21+d32*d31) + &
             d13*(d13*d11+d23*d21+d33*d31) + &
             d21*(d11*d12+d21*d22+d31*d32) + &
             d22*(d12*d12+d22*d22+d32*d32) + &
             d23*(d13*d12+d23*d22+d33*d32) + &
             d31*(d11*d13+d21*d23+d31*d33) + &
             d32*(d12*d13+d22*d23+d32*d33) + &
             d33*(d13*d13+d23*d23+d33*d33)
        var2=d11*(d11d11+d21d21+d31d31) + &
             d12*(d11d12+d21d22+d31d32) + &
             d13*(d11d13+d21d23+d31d33) + &
             d21*(d11d12+d21d22+d31d32) + &
             d22*(d12d12+d22d22+d32d32) + &
             d23*(d12d13+d22d23+d32d33) + &
             d31*(d11d13+d21d23+d31d33) + &
             d32*(d12d13+d22d23+d32d33) + &
             d33*(d13d13+d23d23+d33d33)
        var3=(d11*d11d11+d21*d11d21+d31*d11d31) + &
             (d12*d11d12+d22*d12d21+d32*d12d31) + &
             (d13*d11d13+d23*d13d21+d33*d13d31) + &
             (d11*d12d21+d21*d21d22+d31*d21d32) + &
             (d12*d12d22+d22*d22d22+d32*d22d32) + &
             (d13*d12d23+d23*d22d23+d33*d23d32) + &
             (d11*d13d31+d21*d23d31+d31*d31d33) + &
             (d12*d13d32+d22*d23d32+d32*d32d33) + &
             (d13*d13d33+d23*d23d33+d33*d33d33)
        var4=(d11*d11d11+d21*d11d21+d31*d11d31) + &
             (d11*d12d12+d21*d12d22+d31*d12d32) + &
             (d11*d13d13+d21*d13d23+d31*d13d33) + &
             (d12*d11d21+d22*d21d21+d32*d21d31) + &
             (d12*d12d22+d22*d22d22+d32*d22d32) + &
             (d12*d13d23+d22*d23d23+d32*d23d33) + &
             (d13*d11d31+d23*d21d31+d33*d31d31) + &
             (d13*d12d32+d23*d22d32+d33*d32d32) + &
             (d13*d13d33+d23*d23d33+d33*d33d33)
        m_ijm_kjm_ki_xzm(j) = m_ijm_kjm_ki_xzm(j) -var1-var2-var3-var4
        !
      enddo
      !
      open(18,file='m_Betchov.dat')
      write(18,"(9(1X,A15))")'y','m_iim_jj','m_ijm_ji','m_ijm_ij',       &
                          'm_iim_jjm_kk','m_ijm_jim_kk','m_ijm_ijm_kk',  &
                          'm_ijm_jkm_ki','m_ijm_kjm_ki'
      write(18,"(9(1X,E15.7E3))")(y(0,j,0),m_iim_jj_xzm(j),m_ijm_ji_xzm(j),   &
                             m_ijm_ij_xzm(j),m_iim_jjm_kk_xzm(j),        &
                             m_ijm_jim_kk_xzm(j),m_ijm_ijm_kk_xzm(j),    &
                             m_ijm_jkm_ki_xzm(j),m_ijm_kjm_ki_xzm(j),j=0,jm)
      close(18)
      print*,' << m_Betchov.dat ... done !'
      !
      inquire(file='Results/betchov.xzm.h5',exist=lfilalive)
      if(lfilalive) call system('mv  -v Results/betchov.xzm.h5 Results/betchov.xzm.bak')
      !
      call H5WriteArray(    m_iim_jj_xzm,jm,    'm_iim_jj','Results/betchov.xzm.h5')     
      call H5WriteArray(    m_ijm_ji_xzm,jm,    'm_ijm_ji','Results/betchov.xzm.h5')      
      call H5WriteArray(    m_ijm_ij_xzm,jm,    'm_ijm_ij','Results/betchov.xzm.h5')     
      call H5WriteArray(m_iim_jjm_kk_xzm,jm,'m_iim_jjm_kk','Results/betchov.xzm.h5') 
      call H5WriteArray(m_ijm_jim_kk_xzm,jm,'m_ijm_jim_kk','Results/betchov.xzm.h5') 
      call H5WriteArray(m_ijm_ijm_kk_xzm,jm,'m_ijm_ijm_kk','Results/betchov.xzm.h5') 
      call H5WriteArray(m_ijm_jkm_ki_xzm,jm,'m_ijm_jkm_ki','Results/betchov.xzm.h5') 
      call H5WriteArray(m_ijm_kjm_ki_xzm,jm,'m_ijm_kjm_ki','Results/betchov.xzm.h5') 
      call H5WriteArray(     d11_d11_xzm,jm,     'd11_d11','Results/betchov.xzm.h5')
      call H5WriteArray(     d11_d12_xzm,jm,     'd11_d12','Results/betchov.xzm.h5')
      call H5WriteArray(     d11_d13_xzm,jm,     'd11_d13','Results/betchov.xzm.h5')
      call H5WriteArray(     d11_d21_xzm,jm,     'd11_d21','Results/betchov.xzm.h5')
      call H5WriteArray(     d11_d22_xzm,jm,     'd11_d22','Results/betchov.xzm.h5')
      call H5WriteArray(     d11_d23_xzm,jm,     'd11_d23','Results/betchov.xzm.h5')
      call H5WriteArray(     d11_d31_xzm,jm,     'd11_d31','Results/betchov.xzm.h5')
      call H5WriteArray(     d11_d32_xzm,jm,     'd11_d32','Results/betchov.xzm.h5')
      call H5WriteArray(     d11_d33_xzm,jm,     'd11_d33','Results/betchov.xzm.h5')
      call H5WriteArray(     d12_d12_xzm,jm,     'd12_d12','Results/betchov.xzm.h5')
      call H5WriteArray(     d12_d13_xzm,jm,     'd12_d13','Results/betchov.xzm.h5')
      call H5WriteArray(     d12_d21_xzm,jm,     'd12_d21','Results/betchov.xzm.h5')
      call H5WriteArray(     d12_d22_xzm,jm,     'd12_d22','Results/betchov.xzm.h5')
      call H5WriteArray(     d12_d23_xzm,jm,     'd12_d23','Results/betchov.xzm.h5')
      call H5WriteArray(     d12_d31_xzm,jm,     'd12_d31','Results/betchov.xzm.h5')
      call H5WriteArray(     d12_d32_xzm,jm,     'd12_d32','Results/betchov.xzm.h5')
      call H5WriteArray(     d12_d33_xzm,jm,     'd12_d33','Results/betchov.xzm.h5')
      call H5WriteArray(     d13_d13_xzm,jm,     'd13_d13','Results/betchov.xzm.h5')
      call H5WriteArray(     d13_d21_xzm,jm,     'd13_d21','Results/betchov.xzm.h5')
      call H5WriteArray(     d13_d22_xzm,jm,     'd13_d22','Results/betchov.xzm.h5')
      call H5WriteArray(     d13_d23_xzm,jm,     'd13_d23','Results/betchov.xzm.h5')
      call H5WriteArray(     d13_d31_xzm,jm,     'd13_d31','Results/betchov.xzm.h5')
      call H5WriteArray(     d13_d32_xzm,jm,     'd13_d32','Results/betchov.xzm.h5')
      call H5WriteArray(     d13_d33_xzm,jm,     'd13_d33','Results/betchov.xzm.h5')
      call H5WriteArray(     d21_d21_xzm,jm,     'd21_d21','Results/betchov.xzm.h5')
      call H5WriteArray(     d21_d22_xzm,jm,     'd21_d22','Results/betchov.xzm.h5')
      call H5WriteArray(     d21_d23_xzm,jm,     'd21_d23','Results/betchov.xzm.h5')
      call H5WriteArray(     d21_d31_xzm,jm,     'd21_d31','Results/betchov.xzm.h5')
      call H5WriteArray(     d21_d32_xzm,jm,     'd21_d32','Results/betchov.xzm.h5')
      call H5WriteArray(     d21_d33_xzm,jm,     'd21_d33','Results/betchov.xzm.h5')
      call H5WriteArray(     d22_d22_xzm,jm,     'd22_d22','Results/betchov.xzm.h5')
      call H5WriteArray(     d22_d23_xzm,jm,     'd22_d23','Results/betchov.xzm.h5')
      call H5WriteArray(     d22_d31_xzm,jm,     'd22_d31','Results/betchov.xzm.h5')
      call H5WriteArray(     d22_d32_xzm,jm,     'd22_d32','Results/betchov.xzm.h5')
      call H5WriteArray(     d22_d33_xzm,jm,     'd22_d33','Results/betchov.xzm.h5')
      call H5WriteArray(     d23_d23_xzm,jm,     'd23_d23','Results/betchov.xzm.h5')
      call H5WriteArray(     d23_d31_xzm,jm,     'd23_d31','Results/betchov.xzm.h5')
      call H5WriteArray(     d23_d32_xzm,jm,     'd23_d32','Results/betchov.xzm.h5')
      call H5WriteArray(     d23_d33_xzm,jm,     'd23_d33','Results/betchov.xzm.h5')
      call H5WriteArray(     d31_d31_xzm,jm,     'd31_d31','Results/betchov.xzm.h5')
      call H5WriteArray(     d31_d32_xzm,jm,     'd31_d32','Results/betchov.xzm.h5')
      call H5WriteArray(     d31_d33_xzm,jm,     'd31_d33','Results/betchov.xzm.h5')
      call H5WriteArray(     d32_d32_xzm,jm,     'd32_d32','Results/betchov.xzm.h5')
      call H5WriteArray(     d32_d33_xzm,jm,     'd32_d33','Results/betchov.xzm.h5')
      call H5WriteArray(     d33_d33_xzm,jm,     'd33_d33','Results/betchov.xzm.h5')
      call H5WriteArray(         m11_xzm,jm,         'm11','Results/betchov.xzm.h5')        
      call H5WriteArray(         m12_xzm,jm,         'm12','Results/betchov.xzm.h5')        
      call H5WriteArray(         m13_xzm,jm,         'm13','Results/betchov.xzm.h5')        
      call H5WriteArray(         m21_xzm,jm,         'm21','Results/betchov.xzm.h5')        
      call H5WriteArray(         m22_xzm,jm,         'm22','Results/betchov.xzm.h5')        
      call H5WriteArray(         m23_xzm,jm,         'm23','Results/betchov.xzm.h5')        
      call H5WriteArray(         m31_xzm,jm,         'm31','Results/betchov.xzm.h5')        
      call H5WriteArray(         m32_xzm,jm,         'm32','Results/betchov.xzm.h5')        
      call H5WriteArray(         m33_xzm,jm,         'm33','Results/betchov.xzm.h5')
      !
    elseif(lkhomo) then
      !
      do i=0,im
      do j=0,jm
        !
            m_iim_jj_zm(i,j) = integration(    m_iim_jj(i,j,1:km))     
            m_ijm_ji_zm(i,j) = integration(    m_ijm_ji(i,j,1:km))      
            m_ijm_ij_zm(i,j) = integration(    m_ijm_ij(i,j,1:km))     
        m_iim_jjm_kk_zm(i,j) = integration(m_iim_jjm_kk(i,j,1:km)) 
        m_ijm_jim_kk_zm(i,j) = integration(m_ijm_jim_kk(i,j,1:km)) 
        m_ijm_ijm_kk_zm(i,j) = integration(m_ijm_ijm_kk(i,j,1:km)) 
        m_ijm_jkm_ki_zm(i,j) = integration(m_ijm_jkm_ki(i,j,1:km)) 
        m_ijm_kjm_ki_zm(i,j) = integration(m_ijm_kjm_ki(i,j,1:km)) 
             d11_d11_zm(i,j) = integration(     d11_d11(i,j,1:km)) 
             d11_d12_zm(i,j) = integration(     d11_d12(i,j,1:km)) 
             d11_d13_zm(i,j) = integration(     d11_d13(i,j,1:km)) 
             d11_d21_zm(i,j) = integration(     d11_d21(i,j,1:km)) 
             d11_d22_zm(i,j) = integration(     d11_d22(i,j,1:km)) 
             d11_d23_zm(i,j) = integration(     d11_d23(i,j,1:km)) 
             d11_d31_zm(i,j) = integration(     d11_d31(i,j,1:km)) 
             d11_d32_zm(i,j) = integration(     d11_d32(i,j,1:km)) 
             d11_d33_zm(i,j) = integration(     d11_d33(i,j,1:km)) 
             d12_d12_zm(i,j) = integration(     d12_d12(i,j,1:km)) 
             d12_d13_zm(i,j) = integration(     d12_d13(i,j,1:km)) 
             d12_d21_zm(i,j) = integration(     d12_d21(i,j,1:km)) 
             d12_d22_zm(i,j) = integration(     d12_d22(i,j,1:km)) 
             d12_d23_zm(i,j) = integration(     d12_d23(i,j,1:km)) 
             d12_d31_zm(i,j) = integration(     d12_d31(i,j,1:km)) 
             d12_d32_zm(i,j) = integration(     d12_d32(i,j,1:km)) 
             d12_d33_zm(i,j) = integration(     d12_d33(i,j,1:km)) 
             d13_d13_zm(i,j) = integration(     d13_d13(i,j,1:km)) 
             d13_d21_zm(i,j) = integration(     d13_d21(i,j,1:km)) 
             d13_d22_zm(i,j) = integration(     d13_d22(i,j,1:km)) 
             d13_d23_zm(i,j) = integration(     d13_d23(i,j,1:km)) 
             d13_d31_zm(i,j) = integration(     d13_d31(i,j,1:km)) 
             d13_d32_zm(i,j) = integration(     d13_d32(i,j,1:km)) 
             d13_d33_zm(i,j) = integration(     d13_d33(i,j,1:km)) 
             d21_d21_zm(i,j) = integration(     d21_d21(i,j,1:km)) 
             d21_d22_zm(i,j) = integration(     d21_d22(i,j,1:km)) 
             d21_d23_zm(i,j) = integration(     d21_d23(i,j,1:km)) 
             d21_d31_zm(i,j) = integration(     d21_d31(i,j,1:km)) 
             d21_d32_zm(i,j) = integration(     d21_d32(i,j,1:km)) 
             d21_d33_zm(i,j) = integration(     d21_d33(i,j,1:km)) 
             d22_d22_zm(i,j) = integration(     d22_d22(i,j,1:km)) 
             d22_d23_zm(i,j) = integration(     d22_d23(i,j,1:km)) 
             d22_d31_zm(i,j) = integration(     d22_d31(i,j,1:km)) 
             d22_d32_zm(i,j) = integration(     d22_d32(i,j,1:km)) 
             d22_d33_zm(i,j) = integration(     d22_d33(i,j,1:km)) 
             d23_d23_zm(i,j) = integration(     d23_d23(i,j,1:km)) 
             d23_d31_zm(i,j) = integration(     d23_d31(i,j,1:km)) 
             d23_d32_zm(i,j) = integration(     d23_d32(i,j,1:km)) 
             d23_d33_zm(i,j) = integration(     d23_d33(i,j,1:km)) 
             d31_d31_zm(i,j) = integration(     d31_d31(i,j,1:km)) 
             d31_d32_zm(i,j) = integration(     d31_d32(i,j,1:km)) 
             d31_d33_zm(i,j) = integration(     d31_d33(i,j,1:km)) 
             d32_d32_zm(i,j) = integration(     d32_d32(i,j,1:km)) 
             d32_d33_zm(i,j) = integration(     d32_d33(i,j,1:km)) 
             d33_d33_zm(i,j) = integration(     d33_d33(i,j,1:km)) 
                 m11_zm(i,j) = integration(         m11(i,j,1:km))         
                 m12_zm(i,j) = integration(         m12(i,j,1:km))         
                 m13_zm(i,j) = integration(         m13(i,j,1:km))         
                 m21_zm(i,j) = integration(         m21(i,j,1:km))         
                 m22_zm(i,j) = integration(         m22(i,j,1:km))         
                 m23_zm(i,j) = integration(         m23(i,j,1:km))         
                 m31_zm(i,j) = integration(         m31(i,j,1:km))         
                 m32_zm(i,j) = integration(         m32(i,j,1:km))         
                 m33_zm(i,j) = integration(         m33(i,j,1:km))     
        !
      enddo
      enddo
      !    
          m_iim_jj_zm =     m_iim_jj_zm/dble(nsamples*km)     
          m_ijm_ji_zm =     m_ijm_ji_zm/dble(nsamples*km)      
          m_ijm_ij_zm =     m_ijm_ij_zm/dble(nsamples*km)     
      m_iim_jjm_kk_zm = m_iim_jjm_kk_zm/dble(nsamples*km) 
      m_ijm_jim_kk_zm = m_ijm_jim_kk_zm/dble(nsamples*km) 
      m_ijm_ijm_kk_zm = m_ijm_ijm_kk_zm/dble(nsamples*km) 
      m_ijm_jkm_ki_zm = m_ijm_jkm_ki_zm/dble(nsamples*km) 
      m_ijm_kjm_ki_zm = m_ijm_kjm_ki_zm/dble(nsamples*km) 
           d11_d11_zm =      d11_d11_zm/dble(nsamples*km)
           d11_d12_zm =      d11_d12_zm/dble(nsamples*km)
           d11_d13_zm =      d11_d13_zm/dble(nsamples*km)
           d11_d21_zm =      d11_d21_zm/dble(nsamples*km)
           d11_d22_zm =      d11_d22_zm/dble(nsamples*km)
           d11_d23_zm =      d11_d23_zm/dble(nsamples*km)
           d11_d31_zm =      d11_d31_zm/dble(nsamples*km)
           d11_d32_zm =      d11_d32_zm/dble(nsamples*km)
           d11_d33_zm =      d11_d33_zm/dble(nsamples*km)
           d12_d12_zm =      d12_d12_zm/dble(nsamples*km)
           d12_d13_zm =      d12_d13_zm/dble(nsamples*km)
           d12_d21_zm =      d12_d21_zm/dble(nsamples*km)
           d12_d22_zm =      d12_d22_zm/dble(nsamples*km)
           d12_d23_zm =      d12_d23_zm/dble(nsamples*km)
           d12_d31_zm =      d12_d31_zm/dble(nsamples*km)
           d12_d32_zm =      d12_d32_zm/dble(nsamples*km)
           d12_d33_zm =      d12_d33_zm/dble(nsamples*km)
           d13_d13_zm =      d13_d13_zm/dble(nsamples*km)
           d13_d21_zm =      d13_d21_zm/dble(nsamples*km)
           d13_d22_zm =      d13_d22_zm/dble(nsamples*km)
           d13_d23_zm =      d13_d23_zm/dble(nsamples*km)
           d13_d31_zm =      d13_d31_zm/dble(nsamples*km)
           d13_d32_zm =      d13_d32_zm/dble(nsamples*km)
           d13_d33_zm =      d13_d33_zm/dble(nsamples*km)
           d21_d21_zm =      d21_d21_zm/dble(nsamples*km)
           d21_d22_zm =      d21_d22_zm/dble(nsamples*km)
           d21_d23_zm =      d21_d23_zm/dble(nsamples*km)
           d21_d31_zm =      d21_d31_zm/dble(nsamples*km)
           d21_d32_zm =      d21_d32_zm/dble(nsamples*km)
           d21_d33_zm =      d21_d33_zm/dble(nsamples*km)
           d22_d22_zm =      d22_d22_zm/dble(nsamples*km)
           d22_d23_zm =      d22_d23_zm/dble(nsamples*km)
           d22_d31_zm =      d22_d31_zm/dble(nsamples*km)
           d22_d32_zm =      d22_d32_zm/dble(nsamples*km)
           d22_d33_zm =      d22_d33_zm/dble(nsamples*km)
           d23_d23_zm =      d23_d23_zm/dble(nsamples*km)
           d23_d31_zm =      d23_d31_zm/dble(nsamples*km)
           d23_d32_zm =      d23_d32_zm/dble(nsamples*km)
           d23_d33_zm =      d23_d33_zm/dble(nsamples*km)
           d31_d31_zm =      d31_d31_zm/dble(nsamples*km)
           d31_d32_zm =      d31_d32_zm/dble(nsamples*km)
           d31_d33_zm =      d31_d33_zm/dble(nsamples*km)
           d32_d32_zm =      d32_d32_zm/dble(nsamples*km)
           d32_d33_zm =      d32_d33_zm/dble(nsamples*km)
           d33_d33_zm =      d33_d33_zm/dble(nsamples*km)
               m11_zm =          m11_zm/dble(nsamples*km)        
               m12_zm =          m12_zm/dble(nsamples*km)        
               m13_zm =          m13_zm/dble(nsamples*km)        
               m21_zm =          m21_zm/dble(nsamples*km)        
               m22_zm =          m22_zm/dble(nsamples*km)        
               m23_zm =          m23_zm/dble(nsamples*km)        
               m31_zm =          m31_zm/dble(nsamples*km)        
               m32_zm =          m32_zm/dble(nsamples*km)        
               m33_zm =          m33_zm/dble(nsamples*km)
      !
      print*,' ** total number of samples:',(nsamples*km)
      !
      inquire(file='Results/betchov.zm.h5',exist=lfilalive)
      if(lfilalive) call system('mv  -v Results/betchov.zm.h5 Results/betchov.zm.bak')
      !
      do i=0,im
      do j=0,jm
        !
        d11=m11_zm(i,j);      d12=m12_zm(i,j);      d13=m13_zm(i,j)
        d21=m21_zm(i,j);      d22=m22_zm(i,j);      d23=m23_zm(i,j)
        d31=m31_zm(i,j);      d32=m32_zm(i,j);      d33=m33_zm(i,j)
        !
        d11d11=d11_d11_zm(i,j)-d11*d11 
        d11d12=d11_d12_zm(i,j)-d11*d12 
        d11d13=d11_d13_zm(i,j)-d11*d13 
        d11d21=d11_d21_zm(i,j)-d11*d21 
        d11d22=d11_d22_zm(i,j)-d11*d22 
        d11d23=d11_d23_zm(i,j)-d11*d23 
        d11d31=d11_d31_zm(i,j)-d11*d31 
        d11d32=d11_d32_zm(i,j)-d11*d32 
        d11d33=d11_d33_zm(i,j)-d11*d33 
        d12d12=d12_d12_zm(i,j)-d12*d12 
        d12d13=d12_d13_zm(i,j)-d12*d13 
        d12d21=d12_d21_zm(i,j)-d12*d21 
        d12d22=d12_d22_zm(i,j)-d12*d22 
        d12d23=d12_d23_zm(i,j)-d12*d23 
        d12d31=d12_d31_zm(i,j)-d12*d31 
        d12d32=d12_d32_zm(i,j)-d12*d32 
        d12d33=d12_d33_zm(i,j)-d12*d33 
        d13d13=d13_d13_zm(i,j)-d13*d13 
        d13d21=d13_d21_zm(i,j)-d13*d21 
        d13d22=d13_d22_zm(i,j)-d13*d22 
        d13d23=d13_d23_zm(i,j)-d13*d23 
        d13d31=d13_d31_zm(i,j)-d13*d31 
        d13d32=d13_d32_zm(i,j)-d13*d32 
        d13d33=d13_d33_zm(i,j)-d13*d33 
        d21d21=d21_d21_zm(i,j)-d21*d21 
        d21d22=d21_d22_zm(i,j)-d21*d22 
        d21d23=d21_d23_zm(i,j)-d21*d23 
        d21d31=d21_d31_zm(i,j)-d21*d31 
        d21d32=d21_d32_zm(i,j)-d21*d32 
        d21d33=d21_d33_zm(i,j)-d21*d33 
        d22d22=d22_d22_zm(i,j)-d22*d22 
        d22d23=d22_d23_zm(i,j)-d22*d23 
        d22d31=d22_d31_zm(i,j)-d22*d31 
        d22d32=d22_d32_zm(i,j)-d22*d32 
        d22d33=d22_d33_zm(i,j)-d22*d33 
        d23d23=d23_d23_zm(i,j)-d23*d23 
        d23d31=d23_d31_zm(i,j)-d23*d31 
        d23d32=d23_d32_zm(i,j)-d23*d32 
        d23d33=d23_d33_zm(i,j)-d23*d33 
        d31d31=d31_d31_zm(i,j)-d31*d31 
        d31d32=d31_d32_zm(i,j)-d31*d32 
        d31d33=d31_d33_zm(i,j)-d31*d33 
        d32d32=d32_d32_zm(i,j)-d32*d32 
        d32d33=d32_d33_zm(i,j)-d32*d33 
        d33d33=d33_d33_zm(i,j)-d33*d33 
        !
        d11_d11_zm(i,j)=d11d11 
        d11_d12_zm(i,j)=d11d12 
        d11_d13_zm(i,j)=d11d13 
        d11_d21_zm(i,j)=d11d21 
        d11_d22_zm(i,j)=d11d22 
        d11_d23_zm(i,j)=d11d23 
        d11_d31_zm(i,j)=d11d31 
        d11_d32_zm(i,j)=d11d32 
        d11_d33_zm(i,j)=d11d33 
        d12_d12_zm(i,j)=d12d12 
        d12_d13_zm(i,j)=d12d13 
        d12_d21_zm(i,j)=d12d21 
        d12_d22_zm(i,j)=d12d22 
        d12_d23_zm(i,j)=d12d23 
        d12_d31_zm(i,j)=d12d31 
        d12_d32_zm(i,j)=d12d32 
        d12_d33_zm(i,j)=d12d33 
        d13_d13_zm(i,j)=d13d13 
        d13_d21_zm(i,j)=d13d21 
        d13_d22_zm(i,j)=d13d22 
        d13_d23_zm(i,j)=d13d23 
        d13_d31_zm(i,j)=d13d31 
        d13_d32_zm(i,j)=d13d32 
        d13_d33_zm(i,j)=d13d33 
        d21_d21_zm(i,j)=d21d21 
        d21_d22_zm(i,j)=d21d22 
        d21_d23_zm(i,j)=d21d23 
        d21_d31_zm(i,j)=d21d31 
        d21_d32_zm(i,j)=d21d32 
        d21_d33_zm(i,j)=d21d33 
        d22_d22_zm(i,j)=d22d22 
        d22_d23_zm(i,j)=d22d23 
        d22_d31_zm(i,j)=d22d31 
        d22_d32_zm(i,j)=d22d32 
        d22_d33_zm(i,j)=d22d33 
        d23_d23_zm(i,j)=d23d23 
        d23_d31_zm(i,j)=d23d31 
        d23_d32_zm(i,j)=d23d32 
        d23_d33_zm(i,j)=d23d33 
        d31_d31_zm(i,j)=d31d31 
        d31_d32_zm(i,j)=d31d32 
        d31_d33_zm(i,j)=d31d33 
        d32_d32_zm(i,j)=d32d32 
        d32_d33_zm(i,j)=d32d33 
        d33_d33_zm(i,j)=d33d33 
        !
        div=d11+d22+d33
        !
        var1  = div*div
        m_iim_jj_zm(i,j) = m_iim_jj_zm(i,j) - var1
        !
        var1  = d11*d11 + d12*d21 + d13*d31 +     &
                d21*d12 + d22*d22 + d23*d32 +     &
                d31*d13 + d32*d23 + d33*d33
        m_ijm_ji_zm(i,j) = m_ijm_ji_zm(i,j) - var1
        !
        var1  = d11*d11 + d12*d12 + d13*d13 +     &
                d21*d21 + d22*d22 + d23*d23 +     &
                d31*d31 + d32*d32 + d33*d33
        m_ijm_ij_zm(i,j) = m_ijm_ij_zm(i,j) - var1
        !
        var1=div*div*div
        var2=d11d11+d11d22+d11d33 +  &
             d11d22+d22d22+d22d33 +  &
             d11d33+d22d33+d33d33  
        var3=(d11+d22+d33)*var2
        m_iim_jjm_kk_zm(i,j) = m_iim_jjm_kk_zm(i,j) - var1-3.d0*var3
        !
        var1=d11*d11 + d12*d21 + d13*d31 +     &
             d21*d12 + d22*d22 + d23*d32 +     &
             d31*d13 + d32*d23 + d33*d33
        var2=d11*(d11d11+d11d22+d11d33) +   &
             d12*(d11d21+d21d22+d21d33) +   &
             d13*(d11d31+d22d31+d31d33) +   &
             d21*(d11d12+d12d22+d12d33) +   &
             d22*(d11d22+d22d22+d22d33) +   &
             d23*(d11d32+d22d32+d32d33) +   &
             d31*(d11d13+d13d22+d13d33) +   &
             d32*(d11d23+d22d23+d23d33) +   &
             d33*(d11d33+d22d33+d33d33)
        var3=div*(d11d11+d12d21+d13d31  +   &
                  d12d21+d22d22+d23d32  +   &
                  d13d31+d23d32+d33d33 )
        m_ijm_jim_kk_zm(i,j) = m_ijm_jim_kk_zm(i,j) - var1*div - 2.d0*var2 -var3
        !
        var1  = d11*d11 + d12*d12 + d13*d13 +     &
                d21*d21 + d22*d22 + d23*d23 +     &
                d31*d31 + d32*d32 + d33*d33
        var2=d11*(d11d11+d11d22+d11d33) +   &
             d12*(d11d12+d12d22+d12d33) +   &
             d13*(d11d13+d13d22+d13d33) +   &
             d21*(d11d21+d21d22+d21d33) +   &
             d22*(d11d22+d22d22+d22d33) +   &
             d23*(d11d23+d22d23+d23d33) +   &
             d31*(d11d31+d22d31+d31d33) +   &
             d32*(d11d32+d22d32+d32d33) +   &
             d33*(d11d33+d22d33+d33d33)
        var3=div*(d11d11+d12d12+d13d13  +   &
                  d21d21+d22d22+d23d23  +   &
                  d31d31+d32d32+d33d33 )
        m_ijm_ijm_kk_zm(i,j) = m_ijm_ijm_kk_zm(i,j) - var1*div  - 2.d0*var2 -var3
        !
        var1= d11*d11*d11 + d11*d12*d21 + d11*d13*d31 +  &
              d12*d21*d11 + d12*d22*d21 + d12*d23*d31 +  &
              d13*d31*d11 + d13*d32*d21 + d13*d33*d31 +  &
              d21*d11*d12 + d21*d12*d22 + d21*d13*d32 +  &
              d22*d21*d12 + d22*d22*d22 + d22*d23*d32 +  &
              d23*d31*d12 + d23*d32*d22 + d23*d33*d32 +  &
              d31*d11*d13 + d31*d12*d23 + d31*d13*d33 +  &
              d32*d21*d13 + d32*d22*d23 + d32*d23*d33 +  &
              d33*d31*d13 + d33*d32*d23 + d33*d33*d33
        var2= d11*(d11d11+d12d21+d13d31) +  &
              d12*(d11d21+d21d22+d23d31) +  &
              d13*(d11d31+d21d32+d31d33) +  &
              d21*(d11d12+d12d22+d13d32) +  &
              d22*(d12d21+d22d22+d23d32) +  &
              d23*(d12d31+d22d32+d32d33) +  &
              d31*(d11d13+d12d23+d13d33) +  &
              d32*(d13d21+d22d23+d23d33) +  &
              d33*(d13d31+d23d32+d33d33)
        m_ijm_jkm_ki_zm(i,j) = m_ijm_jkm_ki_zm(i,j) - var1 -3.d0*var2
        !
        var1=d11*(d11*d11+d21*d21+d31*d31) + &
             d12*(d12*d11+d22*d21+d32*d31) + &
             d13*(d13*d11+d23*d21+d33*d31) + &
             d21*(d11*d12+d21*d22+d31*d32) + &
             d22*(d12*d12+d22*d22+d32*d32) + &
             d23*(d13*d12+d23*d22+d33*d32) + &
             d31*(d11*d13+d21*d23+d31*d33) + &
             d32*(d12*d13+d22*d23+d32*d33) + &
             d33*(d13*d13+d23*d23+d33*d33)
        var2=d11*(d11d11+d21d21+d31d31) + &
             d12*(d11d12+d21d22+d31d32) + &
             d13*(d11d13+d21d23+d31d33) + &
             d21*(d11d12+d21d22+d31d32) + &
             d22*(d12d12+d22d22+d32d32) + &
             d23*(d12d13+d22d23+d32d33) + &
             d31*(d11d13+d21d23+d31d33) + &
             d32*(d12d13+d22d23+d32d33) + &
             d33*(d13d13+d23d23+d33d33)
        var3=(d11*d11d11+d21*d11d21+d31*d11d31) + &
             (d12*d11d12+d22*d12d21+d32*d12d31) + &
             (d13*d11d13+d23*d13d21+d33*d13d31) + &
             (d11*d12d21+d21*d21d22+d31*d21d32) + &
             (d12*d12d22+d22*d22d22+d32*d22d32) + &
             (d13*d12d23+d23*d22d23+d33*d23d32) + &
             (d11*d13d31+d21*d23d31+d31*d31d33) + &
             (d12*d13d32+d22*d23d32+d32*d32d33) + &
             (d13*d13d33+d23*d23d33+d33*d33d33)
        var4=(d11*d11d11+d21*d11d21+d31*d11d31) + &
             (d11*d12d12+d21*d12d22+d31*d12d32) + &
             (d11*d13d13+d21*d13d23+d31*d13d33) + &
             (d12*d11d21+d22*d21d21+d32*d21d31) + &
             (d12*d12d22+d22*d22d22+d32*d22d32) + &
             (d12*d13d23+d22*d23d23+d32*d23d33) + &
             (d13*d11d31+d23*d21d31+d33*d31d31) + &
             (d13*d12d32+d23*d22d32+d33*d32d32) + &
             (d13*d13d33+d23*d23d33+d33*d33d33)
        m_ijm_kjm_ki_zm(i,j) = m_ijm_kjm_ki_zm(i,j) -var1-var2-var3-var4
        !
      enddo
      enddo
      !
      call H5WriteArray(    m_iim_jj_zm,im,jm,    'm_iim_jj','Results/betchov.zm.h5')     
      call H5WriteArray(    m_ijm_ji_zm,im,jm,    'm_ijm_ji','Results/betchov.zm.h5')      
      call H5WriteArray(    m_ijm_ij_zm,im,jm,    'm_ijm_ij','Results/betchov.zm.h5')     
      call H5WriteArray(m_iim_jjm_kk_zm,im,jm,'m_iim_jjm_kk','Results/betchov.zm.h5') 
      call H5WriteArray(m_ijm_jim_kk_zm,im,jm,'m_ijm_jim_kk','Results/betchov.zm.h5') 
      call H5WriteArray(m_ijm_ijm_kk_zm,im,jm,'m_ijm_ijm_kk','Results/betchov.zm.h5') 
      call H5WriteArray(m_ijm_jkm_ki_zm,im,jm,'m_ijm_jkm_ki','Results/betchov.zm.h5') 
      call H5WriteArray(m_ijm_kjm_ki_zm,im,jm,'m_ijm_kjm_ki','Results/betchov.zm.h5') 
      call H5WriteArray(     d11_d11_zm,im,jm,     'd11_d11','Results/betchov.zm.h5')
      call H5WriteArray(     d11_d12_zm,im,jm,     'd11_d12','Results/betchov.zm.h5')
      call H5WriteArray(     d11_d13_zm,im,jm,     'd11_d13','Results/betchov.zm.h5')
      call H5WriteArray(     d11_d21_zm,im,jm,     'd11_d21','Results/betchov.zm.h5')
      call H5WriteArray(     d11_d22_zm,im,jm,     'd11_d22','Results/betchov.zm.h5')
      call H5WriteArray(     d11_d23_zm,im,jm,     'd11_d23','Results/betchov.zm.h5')
      call H5WriteArray(     d11_d31_zm,im,jm,     'd11_d31','Results/betchov.zm.h5')
      call H5WriteArray(     d11_d32_zm,im,jm,     'd11_d32','Results/betchov.zm.h5')
      call H5WriteArray(     d11_d33_zm,im,jm,     'd11_d33','Results/betchov.zm.h5')
      call H5WriteArray(     d12_d12_zm,im,jm,     'd12_d12','Results/betchov.zm.h5')
      call H5WriteArray(     d12_d13_zm,im,jm,     'd12_d13','Results/betchov.zm.h5')
      call H5WriteArray(     d12_d21_zm,im,jm,     'd12_d21','Results/betchov.zm.h5')
      call H5WriteArray(     d12_d22_zm,im,jm,     'd12_d22','Results/betchov.zm.h5')
      call H5WriteArray(     d12_d23_zm,im,jm,     'd12_d23','Results/betchov.zm.h5')
      call H5WriteArray(     d12_d31_zm,im,jm,     'd12_d31','Results/betchov.zm.h5')
      call H5WriteArray(     d12_d32_zm,im,jm,     'd12_d32','Results/betchov.zm.h5')
      call H5WriteArray(     d12_d33_zm,im,jm,     'd12_d33','Results/betchov.zm.h5')
      call H5WriteArray(     d13_d13_zm,im,jm,     'd13_d13','Results/betchov.zm.h5')
      call H5WriteArray(     d13_d21_zm,im,jm,     'd13_d21','Results/betchov.zm.h5')
      call H5WriteArray(     d13_d22_zm,im,jm,     'd13_d22','Results/betchov.zm.h5')
      call H5WriteArray(     d13_d23_zm,im,jm,     'd13_d23','Results/betchov.zm.h5')
      call H5WriteArray(     d13_d31_zm,im,jm,     'd13_d31','Results/betchov.zm.h5')
      call H5WriteArray(     d13_d32_zm,im,jm,     'd13_d32','Results/betchov.zm.h5')
      call H5WriteArray(     d13_d33_zm,im,jm,     'd13_d33','Results/betchov.zm.h5')
      call H5WriteArray(     d21_d21_zm,im,jm,     'd21_d21','Results/betchov.zm.h5')
      call H5WriteArray(     d21_d22_zm,im,jm,     'd21_d22','Results/betchov.zm.h5')
      call H5WriteArray(     d21_d23_zm,im,jm,     'd21_d23','Results/betchov.zm.h5')
      call H5WriteArray(     d21_d31_zm,im,jm,     'd21_d31','Results/betchov.zm.h5')
      call H5WriteArray(     d21_d32_zm,im,jm,     'd21_d32','Results/betchov.zm.h5')
      call H5WriteArray(     d21_d33_zm,im,jm,     'd21_d33','Results/betchov.zm.h5')
      call H5WriteArray(     d22_d22_zm,im,jm,     'd22_d22','Results/betchov.zm.h5')
      call H5WriteArray(     d22_d23_zm,im,jm,     'd22_d23','Results/betchov.zm.h5')
      call H5WriteArray(     d22_d31_zm,im,jm,     'd22_d31','Results/betchov.zm.h5')
      call H5WriteArray(     d22_d32_zm,im,jm,     'd22_d32','Results/betchov.zm.h5')
      call H5WriteArray(     d22_d33_zm,im,jm,     'd22_d33','Results/betchov.zm.h5')
      call H5WriteArray(     d23_d23_zm,im,jm,     'd23_d23','Results/betchov.zm.h5')
      call H5WriteArray(     d23_d31_zm,im,jm,     'd23_d31','Results/betchov.zm.h5')
      call H5WriteArray(     d23_d32_zm,im,jm,     'd23_d32','Results/betchov.zm.h5')
      call H5WriteArray(     d23_d33_zm,im,jm,     'd23_d33','Results/betchov.zm.h5')
      call H5WriteArray(     d31_d31_zm,im,jm,     'd31_d31','Results/betchov.zm.h5')
      call H5WriteArray(     d31_d32_zm,im,jm,     'd31_d32','Results/betchov.zm.h5')
      call H5WriteArray(     d31_d33_zm,im,jm,     'd31_d33','Results/betchov.zm.h5')
      call H5WriteArray(     d32_d32_zm,im,jm,     'd32_d32','Results/betchov.zm.h5')
      call H5WriteArray(     d32_d33_zm,im,jm,     'd32_d33','Results/betchov.zm.h5')
      call H5WriteArray(     d33_d33_zm,im,jm,     'd33_d33','Results/betchov.zm.h5')
      call H5WriteArray(         m11_zm,im,jm,         'm11','Results/betchov.zm.h5')        
      call H5WriteArray(         m12_zm,im,jm,         'm12','Results/betchov.zm.h5')        
      call H5WriteArray(         m13_zm,im,jm,         'm13','Results/betchov.zm.h5')        
      call H5WriteArray(         m21_zm,im,jm,         'm21','Results/betchov.zm.h5')        
      call H5WriteArray(         m22_zm,im,jm,         'm22','Results/betchov.zm.h5')        
      call H5WriteArray(         m23_zm,im,jm,         'm23','Results/betchov.zm.h5')        
      call H5WriteArray(         m31_zm,im,jm,         'm31','Results/betchov.zm.h5')        
      call H5WriteArray(         m32_zm,im,jm,         'm32','Results/betchov.zm.h5')        
      call H5WriteArray(         m33_zm,im,jm,         'm33','Results/betchov.zm.h5')
      !
    endif
    !
  end subroutine betchovnorm
  !
  subroutine momtana(nsta,nend,iref)
    !
    use commvardefine,only: im,jm,km,x,y,z,Reynolds,Mach
    !
    integer,intent(in) :: nsta,nend,iref
    !
    character(len=4) :: fname,fname1,fname2
    character(len=128) :: filename
    real(8),allocatable,dimension(:,:) :: rom,u1m,u2m,tm,u11,u22,u33
    real(8),allocatable,dimension(:,:,:) :: ro,u1,u2,u3
    real(8),allocatable,dimension(:,:,:,:) :: du1,du2,du3
    real(8) :: du11,du12,du13,du21,du22,du23,du31,du32,du33,          &
               s11,s12,s13,s21,s22,s23,s31,s32,s33,div,var1,var2,var3,ros
    real(8) :: m_iim_jj(0:jm),m_ijm_ji(0:jm),m_ijm_ij(0:jm),           &
               m_iim_jjm_kk(0:jm),m_ijm_jim_kk(0:jm),                  &
               m_ijm_ijm_kk(0:jm),m_ijm_jkm_ki(0:jm),                  &
               m_ijm_kjm_ki(0:jm),Marms(0:jm)
    integer :: n,i,j,k,interval
    !
    write(fname1,'(i4.4)')nsta
    write(fname2,'(i4.4)')nend
    !
    write(*,'(A,I0,4(A))')'  ** get statistics at the location i=',    &
                           iref,' flowfield from ',fname1,' to ',fname2
    !
    call H5_ReadGrid
    !
    allocate( u1m(0:im,0:jm),u2m(0:im,0:jm),tm(0:im,0:jm) )
    !
    call H5ReadArray(tm,im,jm,'t','Results/mean.fav.zm.h5')
    ! call H5ReadArray(u1m,im,jm,'u1','Results/mean.fav.zm.h5')
    ! call H5ReadArray(u2m,im,jm,'u2','Results/mean.fav.zm.h5')
    call H5ReadArray(u1m,im,jm,'u1rem','Results/budget.zm.h5')
    call H5ReadArray(u2m,im,jm,'u2rem','Results/budget.zm.h5')
    !
    allocate( u1(0:im,0:jm,0:km),u2(0:im,0:jm,0:km),u3(0:im,0:jm,0:km) )
    allocate( du1(1:3,0:im,0:jm,0:km),du2(1:3,0:im,0:jm,0:km),         &
              du3(1:3,0:im,0:jm,0:km) )
    !
    m_iim_jj  = 0.d0
    m_ijm_ji  = 0.d0
    m_ijm_ij  = 0.d0
    m_iim_jjm_kk = 0.d0
    m_ijm_jim_kk = 0.d0
    m_ijm_ijm_kk = 0.d0
    m_ijm_jkm_ki = 0.d0
    m_ijm_kjm_ki = 0.d0
    !
    interval=255
    !
    print*,'iref=',iref,'interval=',interval
    !
    do n=nsta,nend
      !
      write(fname,'(i4.4)')n
      !
      filename=outfolder//'flowfield'//fname//'.h5'
      print*,' >> data:',trim(filename)
      !
      ! call H5ReadArray(ro,im,jm,km,'ro',filename)
      call H5ReadArray(u1,im,jm,km,'u1',filename)
      call H5ReadArray(u2,im,jm,km,'u2',filename)
      call H5ReadArray(u3,im,jm,km,'u3',filename)
      !
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
      !$OMP DO
      do j=0,jm
      do i=0,im
        u1(i,j,:)=u1(i,j,:)-u1m(i,j)
        u2(i,j,:)=u2(i,j,:)-u2m(i,j)
      enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !
      if(n==nsta) then
        du1=grad_3d(u1,x,y,z)
      else
        du1=grad_3d(u1)
      endif
      du2=grad_3d(u2)
      du3=grad_3d(u3)
      !
      do j=0,jm
      do k=1,km
        !
        do i=iref-interval,iref+interval
        ! ros=ro(iref,j,k)
        ros=1.d0
        du11=du1(1,i,j,k); du21=du2(1,i,j,k); du31=du3(1,i,j,k)
        du12=du1(2,i,j,k); du22=du2(2,i,j,k); du32=du3(2,i,j,k)
        du13=du1(3,i,j,k); du23=du2(3,i,j,k); du33=du3(3,i,j,k)
        !
        div=du11+du22+du33
        !
        var1  = div*div
        var2  = du11*du11 + du12*du21 + du13*du31 +     &
                du21*du12 + du22*du22 + du23*du32 +     &
                du31*du13 + du32*du23 + du33*du33
        var3  = du11*du11 + du12*du12 + du13*du13 +     &
                du21*du21 + du22*du22 + du23*du23 +     &
                du31*du31 + du32*du32 + du33*du33
        !
        m_iim_jj(j)  = m_iim_jj(j) + var1 
        m_ijm_ji(j)  = m_ijm_ji(j) + var2     
        m_ijm_ij(j)  = m_ijm_ij(j) + var3
        !
        m_iim_jjm_kk(j) = m_iim_jjm_kk(j)+var1*div
        m_ijm_jim_kk(j) = m_ijm_jim_kk(j)+var2*div
        m_ijm_ijm_kk(j) = m_ijm_ijm_kk(j)+var3*div
        m_ijm_jkm_ki(j) = m_ijm_jkm_ki(j) +                            &
                 du11*du11*du11+  du11*du12*du21 + du11*du13*du31 +    &
                 du12*du21*du11+  du12*du22*du21 + du12*du23*du31 +    &
                 du13*du31*du11+  du13*du32*du21 + du13*du33*du31 +    &
                 du21*du11*du12+  du21*du12*du22 + du21*du13*du32 +    &
                 du22*du21*du12+  du22*du22*du22 + du22*du23*du32 +    &
                 du23*du31*du12+  du23*du32*du22 + du23*du33*du32 +    &
                 du31*du11*du13+  du31*du12*du23 + du31*du13*du33 +    &
                 du32*du21*du13+  du32*du22*du23 + du32*du23*du33 +    &
                 du33*du31*du13+  du33*du32*du23 + du33*du33*du33 
        m_ijm_kjm_ki(j) = m_ijm_kjm_ki(j) +                            &
                 du11*du11*du11+  du11*du21*du21 + du11*du31*du31 +    &
                 du12*du12*du11+  du12*du22*du21 + du12*du32*du31 +    &
                 du13*du13*du11+  du13*du23*du21 + du13*du33*du31 +    &
                 du21*du11*du12+  du21*du21*du22 + du21*du31*du32 +    &
                 du22*du12*du12+  du22*du22*du22 + du22*du32*du32 +    &
                 du23*du13*du12+  du23*du23*du22 + du23*du33*du32 +    &
                 du31*du11*du13+  du31*du21*du23 + du31*du31*du33 +    &
                 du32*du12*du13+  du32*du22*du23 + du32*du32*du33 +    &
                 du33*du13*du13+  du33*du23*du23 + du33*du33*du33 
        enddo
      end do
      end do
      !
    enddo
    !
    m_iim_jj     = m_iim_jj    /dble((2*interval+1)*km*(nend-nsta+1)) !/rom(iref,:)
    m_ijm_ji     = m_ijm_ji    /dble((2*interval+1)*km*(nend-nsta+1)) !/rom(iref,:)
    m_ijm_ij     = m_ijm_ij    /dble((2*interval+1)*km*(nend-nsta+1)) !/rom(iref,:)
    m_iim_jjm_kk = m_iim_jjm_kk/dble((2*interval+1)*km*(nend-nsta+1)) !/rom(iref,:)
    m_ijm_jim_kk = m_ijm_jim_kk/dble((2*interval+1)*km*(nend-nsta+1)) !/rom(iref,:)
    m_ijm_ijm_kk = m_ijm_ijm_kk/dble((2*interval+1)*km*(nend-nsta+1)) !/rom(iref,:)
    m_ijm_jkm_ki = m_ijm_jkm_ki/dble((2*interval+1)*km*(nend-nsta+1)) !/rom(iref,:)
    m_ijm_kjm_ki = m_ijm_kjm_ki/dble((2*interval+1)*km*(nend-nsta+1)) !/rom(iref,:)
    !
    allocate(u11(0:im,0:jm),u22(0:im,0:jm),u33(0:im,0:jm))
    call H5ReadArray(u11,im,jm,'u11','Results/2order.fav.zm.h5')
    call H5ReadArray(u22,im,jm,'u22','Results/2order.fav.zm.h5')
    call H5ReadArray(u33,im,jm,'u33','Results/2order.fav.zm.h5')
    !
    do j=0,jm
      var1=sqrt(abs(u11(iref,j)+u22(iref,j)+u33(iref,j)))
      Marms(j)=var1/sqrt(tm(iref,j))*Mach
    enddo
    !
    open(18,file='m_Betchov.dat')
    write(18,"(10(1X,A15))")'y','m_iim_jj','m_ijm_ji','m_ijm_ij',      &
                        'm_iim_jjm_kk','m_ijm_jim_kk','m_ijm_ijm_kk',  &
                        'm_ijm_jkm_ki','m_ijm_kjm_ki','Marms'
    write(18,"(10(1X,E15.7E3))")(y(iref,j,0),m_iim_jj(j),m_ijm_ji(j),  &
                         m_ijm_ij(j),m_iim_jjm_kk(j),m_ijm_jim_kk(j),  &
                     m_ijm_ijm_kk(j),m_ijm_jkm_ki(j),m_ijm_kjm_ki(j),  &
                     Marms(j),j=0,jm)
    close(18)
    print*,' << m_Betchov.dat ... done !'

    !
  end subroutine momtana
  !
  subroutine sonicline
    !
    use commvardefine,only: im,jm,km,Reynolds,Mach,pinf
    use interpolation
    !
    integer,parameter :: ims=2180
    integer :: i,j,k,i1
    !
    real(8),allocatable :: z1d(:) 
    real(8),allocatable,dimension(:,:,:) :: x,y,z,u,w,p
    real(8),allocatable,dimension(:,:) :: xf,yf,zf,uf,wf,pf,uzm,pzm
    !
    real(8) :: xsin(0:ims),ysin(0:ims),rosin(0:ims),usin(0:ims),       &
               vsin(0:ims),psin(0:ims),tsin(0:ims)
    real(8) :: xs(0:im),ys(0:im),ros(0:im),us(0:im),               &
               vs(0:im),ps(0:im),ts(0:im)
    !
    open(12,file='Results/sonicline.dat')
    read(12,'(////////////)')
    do i=0,ims
      read(12,*)xsin(i),ysin(i),rosin(i),usin(i),vsin(i),psin(i),tsin(i)
    enddo
    close(12)
    print*,' >> Results/sonicline.dat'
    !
    ! sort from small to large
    do i=0,ims
      !
      do i1=i,ims
        !
        if(xsin(i)>xsin(i1)) then
          call swap(xsin(i),xsin(i1))
          call swap(ysin(i),ysin(i1))
          call swap(rosin(i),rosin(i1))
          call swap(usin(i),usin(i1))
          call swap(vsin(i),vsin(i1))
          call swap(psin(i),psin(i1))
          call swap(tsin(i),tsin(i1))
        endif
        !
      enddo
      !
    enddo
    print*,' ** array sorted from small to large'
    !
    allocate(z1d(0:km))
    call H5ReadSubset(z1d,im,jm,km,'z',trim(gridfile),islice=0,jslice=0)
    !
    allocate(x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km))
    call H5ReadArray(x,im,jm,km,'x',trim(gridfile))
    call H5ReadArray(y,im,jm,km,'y',trim(gridfile))
    call H5ReadArray(z,im,jm,km,'z',trim(gridfile))
    !
    do i=0,im
      xs(i)=x(i,0,0)
    enddo
    !
    call regularlinearinterp(xsin,ysin,rosin,usin,vsin,psin,tsin,ims,  &
                                 xs,ys,ros,us,vs,ps,ts,im)
    !
    open(18,file='Results/sonicline_ordered.dat')
    write(18,"(7(1X,E15.7E3))")(xs(i),ys(i),ros(i),us(i),vs(i),ps(i),ts(i),i=0,im)
    close(18)
    print*,' << Results/sonicline_ordered.dat'
    !
    allocate(u(0:im,0:jm,0:km),w(0:im,0:jm,0:km),p(0:im,0:jm,0:km))
    call H5ReadArray(u,im,jm,km,'u1','outdat/flowfield0219.h5')
    call H5ReadArray(w,im,jm,km,'u3','outdat/flowfield0219.h5')
    call H5ReadArray(p,im,jm,km, 'p','outdat/flowfield0219.h5')
    !
    allocate(uzm(0:im,0:jm),pzm(0:im,0:jm))
    call H5ReadArray(pzm, im,jm, 'p','Results/mean.fav.zm.h5')
    call H5ReadArray(uzm, im,jm, 'u1','Results/mean.fav.zm.h5')
    !
    do k=0,km
      u(:,:,k)=u(:,:,k)-uzm
      p(:,:,k)=p(:,:,k)-pzm
    enddo
    !
    call writetecbin('Results/tecj6.plt', x(:,6,:),'x',y(:,6,:),'y',&
                           z(:,6,:),'z',u(:,6,:),'u"',w(:,6,:),'w"',&
                           p(:,6,:),"p'",im,km)
    !
    stop
    !
    allocate(xf(0:im,0:km),yf(0:im,0:km),zf(0:im,0:km),             &
             uf(0:im,0:km),wf(0:im,0:km),pf(0:im,0:km))
    !
    do k=0,km
      do i=0,im
        xf(i,k)=xs(i)
        yf(i,k)=ys(i)
        zf(i,k)=z1d(k)
      enddo
    enddo 
    !
    call binterp(x(:,:,0),y(:,:,0),u(:,:,:),w(:,:,:),p(:,:,:),       &
                  xf(:,0), yf(:,0), uf(:,:), wf(:,:), pf(:,:),km )
    !
    do k=0,km
      do i=0,im
        uf(i,k)=uf(i,k)-us(i)
        pf(i,k)=pf(i,k)-ps(i)
      enddo
    enddo 
    !
    call writetecbin('Results/tecsonic.plt', xf,'x',yf,'y',zf,'z',     &
                                       uf,'u"',wf,'w"',pf,"p'",im,km)
    !

    !
  end subroutine sonicline
  !
  subroutine shearcore
    !
    use commvardefine,only: im,jm,km,Reynolds,Mach,pinf
    use interpolation
    !
    integer,parameter :: ims=1749
    integer :: i,j,k,i1
    !
    real(8),allocatable :: z1d(:) 
    real(8),allocatable,dimension(:,:,:) :: x,y,z,u,w,p
    real(8),allocatable,dimension(:,:) :: xf,yf,zf,uf,wf,pf,uzm,pzm
    !
    real(8) :: xsin(0:ims),ysin(0:ims),rosin(0:ims),usin(0:ims),       &
               vsin(0:ims),psin(0:ims),tsin(0:ims)
    real(8) :: xs(0:im),ys(0:im),ros(0:im),us(0:im),               &
               vs(0:im),ps(0:im),ts(0:im)
    !
    open(12,file='Results/shear-core-smooth.dat')
    do i=0,ims
      read(12,*)xsin(i),ysin(i),rosin(i),usin(i),vsin(i),psin(i),tsin(i)
    enddo
    close(12)
    print*,' >> Results/shear-core-smooth.dat'
    !
    allocate(z1d(0:km))
    call H5ReadSubset(z1d,im,jm,km,'z',trim(gridfile),islice=0,jslice=0)
    !
    allocate(x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km))
    call H5ReadArray(x,im,jm,km,'x',trim(gridfile))
    call H5ReadArray(y,im,jm,km,'y',trim(gridfile))
    call H5ReadArray(z,im,jm,km,'z',trim(gridfile))
    !
    allocate(u(0:im,0:jm,0:km),w(0:im,0:jm,0:km),p(0:im,0:jm,0:km))
    call H5ReadArray(u,im,jm,km,'u1','outdat/flowfield0219.h5')
    call H5ReadArray(w,im,jm,km,'u3','outdat/flowfield0219.h5')
    call H5ReadArray(p,im,jm,km, 'p','outdat/flowfield0219.h5')
    !
    allocate(uzm(0:im,0:jm),pzm(0:im,0:jm))
    call H5ReadArray(pzm, im,jm, 'p','Results/mean.fav.zm.h5')
    call H5ReadArray(uzm, im,jm, 'u1','Results/mean.fav.zm.h5')
    !
    do k=0,km
      u(:,:,k)=u(:,:,k)-uzm
      p(:,:,k)=p(:,:,k)-pzm
    enddo
    !
    !
    call writetecbin('Results/tecfluc_0219..plt', x(:,0:178,:),'x',     &
                                                 y(:,0:178,:),'y',     &
                                                 z(:,0:178,:),'z',     &
                              u(:,0:178,:),'u"',w(:,0:178,:),'w"',     &
                              p(:,0:178,:),"p'",im,178,km)
    !
    ! allocate(xf(0:ims,0:km),yf(0:ims,0:km),zf(0:ims,0:km),             &
    !          uf(0:ims,0:km),wf(0:ims,0:km),pf(0:ims,0:km))
    ! !
    ! do k=0,km
    !   do i=0,ims
    !     xf(i,k)=xsin(i)
    !     yf(i,k)=ysin(i)
    !     zf(i,k)=z1d(k)
    !   enddo
    ! enddo 
    ! !
    ! call binterp(x(:,:,0),y(:,:,0),u(:,:,:),w(:,:,:),p(:,:,:),       &
    !               xf(:,0), yf(:,0), uf(:,:), wf(:,:), pf(:,:),km )
    ! !
    ! do k=0,km
    !   do i=0,ims
    !     uf(i,k)=uf(i,k)-usin(i)
    !     pf(i,k)=pf(i,k)-psin(i)
    !   enddo
    ! enddo 
    ! !
    ! call writetecbin('Results/tecshearcore.plt', xf,'x',yf,'y',zf,'z',     &
    !                                    uf,'u"',wf,'w"',pf,"p'",ims,km)
    !

    !
  end subroutine shearcore
  !
  subroutine swap(var1,var2)
    !
    real(8) :: var1,var2
    real(8) :: var3
    !
    var3=var1
    var1=var2
    var2=var3
    !
    return
    !
  end subroutine swap
  !
  subroutine wall3d
    !
    use commvardefine,only: im,jm,km,Reynolds,Mach,pinf
    !
    real(8),allocatable,dimension(:,:,:) :: x,y,z,dudy,dwdy,p
    real(8),allocatable,dimension(:,:) :: pzm
    integer :: k
    !
    allocate(x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km))
    allocate(dudy(0:im,0:jm,0:km),dwdy(0:im,0:jm,0:km),p(0:im,0:jm,0:km))
    allocate(pzm(0:im,0:jm))
    !
    call H5ReadArray(x,im,jm,km,'x',trim(gridfile))
    call H5ReadArray(y,im,jm,km,'y',trim(gridfile))
    call H5ReadArray(z,im,jm,km,'z',trim(gridfile))
    !
    call H5ReadArray(pzm, im,jm, 'p','Results/mean.fav.zm.h5')
    call H5ReadArray(p,   im,jm,km,    'p','Results/mean.fav.h5')
    call H5ReadArray(dudy,im,jm,km,'du1dy','Results/mean.fav.h5')
    call H5ReadArray(dwdy,im,jm,km,'du3dy','Results/mean.fav.h5')
    !
    do k=0,km
      p(:,:,k)=p(:,:,k)-pzm(:,:)
    enddo
    !
    call writetecbin('Results/tecwall.plt', x(:,0,:),'x',       &
                                            y(:,0,:),'y',       &
                                            z(:,0,:),'z',       &
                                            p(:,0,:),'p',       &
                                            dudy(:,0,:),'cfx',  &
                                            dwdy(:,0,:),'cfz',im,km)
    !
  end subroutine wall3d
  !
  subroutine swtbli(iref)
    !
    use commvardefine,only: im,jm,km,Reynolds,Mach,pinf
    use basicfunction
    use interpolation 
    use writetec
    use boundarylayer
    use gradsolver
    !
    ! argument
    integer,intent(in) :: iref
    !
    ! local data
    integer :: i,j,jeir,jrh,ns,n,js,jsub
    integer,allocatable :: is(:),jedge(:)
    ! jeir: edge of the rotational part of the boundary layer (δe) is 
    ! defined as the point where the mean spanwise vorticity becomes 
    ! less than a suitable threshold value (here set to 0.005 u∞/δin).
    !
    real(8),allocatable,dimension(:,:) :: x,y,ro,u1,u2,p,t,dudy,dvdx,omegaz,ma,dtdy
    real(8),allocatable,dimension(:,:) :: u11,u22,u33,u12,tke,tt,pp
    real(8),allocatable,dimension(:) :: z,yplus,uplus,bl_delta99,      &
                                        bl_dstar,bl_theta,cf,shapefac, &
                                        cp,pwall,qw,st
    real(8),allocatable,dimension(:,:) :: convect,pre_ace,pre_dil,     &
                                          pre_tra,product,tur_tra,     &
                                          vis_ace,vis_dif,dissipa,balance, &
                                          productxx,productxy,productyy
    !
    real(8),allocatable,dimension(:) :: roin,uin,xfs,yfs,yfsup,yfsub,omzax,yuu,yvv,yww,yuv,yk,uup,vvp,wwp,uvp,tkep
    real(8),allocatable,dimension(:) :: cfv,cft,cfg,cfg1,cfg2,cfg3,yedge,cfg1u,cfg1v,cfg2res,cfg2vis
    real(8),allocatable,dimension(:,:) :: tauxx,tauxy,dudx,dvdy,ww,dpz2,c_dpdz_w
    real(8),allocatable,dimension(:) :: urev
    real(8),allocatable,dimension(:) :: xsonic,ysonic,v1son,v2son,v3son, &
                                        v4son,v5son,v6son
    real(8),allocatable :: temp1(:,:),dtemp1(:,:,:),temp2(:,:),dtemp2(:,:,:),dp(:,:,:)
    real(8) :: utaw,lvis,omegaz_ref,miu,var1,var2,lref,x0,xs,xr,       &
               lsep,betai,ximp,dxmin,ttotal,dy,uedge
    character(len=5) :: irefname,iname
    !
    print*,' ** get statistics at the location i=',iref
    if(iref<=0) stop
    write(irefname,'(1hi,I4.4)')iref
    !
    allocate(x(0:im,0:jm),y(0:im,0:jm),z(0:km))
    call H5ReadSubset(x,im,jm,km,'x',trim(gridfile),kslice=0)
    call H5ReadSubset(y,im,jm,km,'y',trim(gridfile),kslice=0)
    call H5ReadSubset(z,im,jm,km,'z',trim(gridfile),islice=iref,jslice=0)
    !
    ns=9
    allocate(is(ns))
    is(1)=500
    is(2)=708
    is(3)=769
    is(4)=818
    is(5)=863
    is(6)=926
    is(7)=991
    is(8)=1126
    is(9)=1456
    !
    write(*,*)'------- i sample -------'
    do n=1,ns
      write(*,'(2(A,I0),A,F12.6)')'  is',n,'= ',is(n),', x=',x(is(n),0)
    enddo
    write(*,*)'------------------------'
    !
    allocate(roin(0:jm),uin(0:jm))
    open(13,file='./datin/inlet.prof')
    read(13,*)
    read(13,*)
    read(13,*)
    read(13,*)
    do j=0,jm
      read(13,*)roin(j),uin(j)
    enddo
    close(13)
    print*,' >> ./datin/inlet.prof'
    !
    allocate(xsonic(0:im),ysonic(0:im))
    open(12,file='Results/sonicline_ordered.dat')
    do i=0,im
      read(12,*)xsonic(i),ysonic(i)
    enddo
    close(12)
    print*,' >> Results/sonicline_ordered.dat'
    !
    do j=1,jm
      if(uin(j)<0.95d0*uin(j-1)) then
        jrh=j
        exit
      endif
    enddo
    !
    print*,' ** RH jump is at j=',jrh,'y=',y(0,jrh)
    !
    betai=22.76d0
    ximp=y(0,jrh+1)/tan(betai/180.d0*pi)
    !
    ! x0=ximp
    x0=0.d0
    !
    print*,' ** inviscid impinging point is at x=',ximp
    !
    allocate(ro(0:im,0:jm),u1(0:im,0:jm),p(0:im,0:jm),                 &
             t(0:im,0:jm),dudy(0:im,0:jm),dvdx(0:im,0:jm),             &
             omegaz(0:im,0:jm),u2(0:im,0:jm),ma(0:im,0:jm),            &
             dtdy(0:im,0:jm),pp(0:im,0:jm))
    !
    call H5ReadArray(ro,im,jm,'ro','Results/mean.fav.zm.h5')
    call H5ReadArray(u1,im,jm,'u1','Results/mean.fav.zm.h5')
    call H5ReadArray(u2,im,jm,'u2','Results/mean.fav.zm.h5')
    call H5ReadArray(p, im,jm, 'p','Results/mean.fav.zm.h5')
    call H5ReadArray(t, im,jm, 't','Results/mean.fav.zm.h5')
    call H5ReadArray(dudy,im,jm,'du1dy','Results/mean.fav.zm.h5')
    call H5ReadArray(dvdx,im,jm,'du2dx','Results/mean.fav.zm.h5')
    call H5ReadArray(dtdy,im,jm, 'dtdy','Results/mean.fav.zm.h5')
    call H5ReadArray( pp,im,jm, 'pp','Results/2order.fav.zm.h5')
    !
    allocate(cf(0:im),cp(0:im),pwall(0:im),qw(0:im),st(0:im))
    ttotal=1.d0+0.5d0*(gamma-1.d0)*Mach**2
    do i=0,im
      miu=miucal(t(i,0))/Reynolds
      cf(i)=miu*dudy(i,0)*2.d0
      cp(i)=p(i,0)*2.d0
      pwall(i)=p(i,0)/pinf
      !
      qw(i)=miu*dtdy(i,0)
      st(i)=qw(i)/((ttotal-t(i,0)))/Prandtl
    enddo
    !
    do i=1,im
      if(cf(i-1)>=0.d0 .and. cf(i)<=0.d0) then
        xs=linear1d(cf(i-1),cf(i),x(i-1,j),x(i,j),0.d0)
        exit
      endif
    enddo
    !
    do i=im,1,-1
      if(cf(i)>=0.d0 .and. cf(i-1)<=0.d0) then
        xr=linear1d(cf(i-1),cf(i),x(i-1,j),x(i,j),0.d0)
        exit
      endif
    enddo
    lsep=xr-xs
    !
    print*,' ** separation point:',xs
    print*,' ** reattchment point:',xr
    print*,' ** separation size:',lsep
    !
    ! lref=lsep
    lref=1.d0
    !
    allocate(bl_delta99(0:im),bl_dstar(0:im),bl_theta(0:im), &
             shapefac(0:im),jedge(0:im),yedge(0:im))
    !
    do i=0,im
      bl_delta99(i)=0.d0
      do j=1,jm
        if(u1(i,j-1)<=0.99d0 .and. u1(i,j)>=0.99d0) then
          bl_delta99(i)=linear1d(u1(i,j-1),u1(i,j),y(i,j-1),y(i,j),0.99d0)
          exit
        endif
      enddo
    enddo
    omegaz_ref=0.005d0/bl_delta99(iref)
    !
    omegaz=dudy-dvdx
    ma=sqrt(u1**2+u2**2)/sqrt(t)*mach
    !
    call writetecbin('Results/tecmeanflow.plt',x,'x',y,'y',            &
                              ro,'ro',u1,'u',u2,'v', p,'p',            &
                              t,'T',ma,'Ma',-omegaz,'omegaz',im,jm) 
    do i=0,im
      !
      do j=1,jm
        if( abs(omegaz(i,j))<=0.01d0 .and. y(i,j)>0.5d0) then
          jedge(i)=j
          exit
        endif
      enddo
      !
      yedge(i)=linear1d( abs(omegaz(i,jedge(i))),abs(omegaz(i,jedge(i)-1)), &
                         y(i,jedge(i)), y(i,jedge(i)-1),0.01d0)
      !
      bl_dstar(i)=0.d0
      bl_theta(i)=0.d0
      !
      do j=1,jedge(i)
        var1=0.5d0*(u1(i,j)*ro(i,j)+u1(i,j-1)*ro(i,j-1))
        bl_dstar(i)=bl_dstar(i)+(1.d0-var1)*(y(i,j)-y(i,j-1))
        !
        var1=0.5d0*(u1(i,j)*ro(i,j)+u1(i,j-1)*ro(i,j-1))
        var2=1.d0-0.5d0*(u1(i,j)+u1(i,j-1))
        bl_theta(i)=bl_theta(i)+var1*var2*(y(i,j)-y(i,j-1))
      enddo
      !
      shapefac(i)=bl_dstar(i)/bl_theta(i)
      !
    enddo
    !
    write(*,"(A,F7.3,A)")'  ---------- BL thickness at x=',x(iref,0),' ----------'
    write(*,"(A)")'            δ          δ*           θ           H'
    write(*,"(1X,4(F12.7))")bl_delta99(iref),bl_dstar(iref),           &
                                           bl_theta(iref),shapefac(iref)
    write(*,"(A)")'  -----------------------------------------------'
    !
    open(18,file='Results/bledge.dat')
    ! write(18,"(3(1X,A15))")'x','yedge','jedge'
    write(18,"(A)")'TITLE     = "BL edge"'
    write(18,"(A)")'VARIABLES = "x" "y" "j"'
    write(18,"(A)")'ZONE T="ZONE 001"'
    write(18,"(A,I0,A)")'I=',im-99,', J=1, K=1, ZONETYPE=Ordered DATAPACKING=POINT'
    write(18,"(2(1X,E15.7E3),1X,I15)")(x(i,0),yedge(i),jedge(i),i=0,im-100)
    close(18)
    print*,' << Results/bledge.dat ... done. '
    !
    allocate(yplus(0:jm),uplus(0:jm))
    !
    dxmin=1.d10
    do i=1,im
      dxmin=min(dxmin,x(i,0)-x(i-1,0))
    enddo
    !
    i=iref
    !
    open(18,file='Results/profile.'//irefname//'.dat')
    write(18,"(6(1X,A15))")'y/delta','u','v','ro','p','t'
    write(18,"(6(1X,E15.7E3))")(y(i,j)/bl_delta99(i),u1(i,j),u2(i,j),  &
                                           ro(i,j),p(i,j),t(i,j),j=0,jm)
    close(18)
    print*,' << Results/profile.',irefname,'.delta ... done. '
    !
    call upluscal(uplus=uplus,yplus=yplus,u=u1(i,:),y=y(i,:),          &
                                         ro=ro(i,:),tw=t(i,0),utaw=utaw)
    !
    print*,' ** utaw=',utaw
    !
    open(18,file='Results/uplus.'//irefname//'.dat')
    write(18,"(2(1X,A15))")'yplus','uplus'
    write(18,"(2(1X,E15.7E3))")(yplus(j),uplus(j),j=0,jm)
    close(18)
    print*,' << uplus.',irefname,'.dat ... done. '
    !
    lvis=ro(iref,0)*utaw/miucal(T(iref,0))*Reynolds
    open(18,file='Results/mesh_report.'//irefname//'.dat')
    write(18,"(A,F12.7,A,I4)")'mesh resolution at x=',x(i,0),', i=',i
    write(18,"(A)")'-----------------------------------------------------------------'
    write(18,"(A)")'          Δx+       Δxmin+         Δye+          Δz+          Lz+'
    write(18,"(A)")'-----------------------------------------------------------------'
    write(18,"(6(1X,F12.7))")0.5d0*(x(i+1,0)-x(i-1,0))*lvis,           &
                                                 dxmin*lvis,           &
                                       (y(i,1)-y(i,0))*lvis,           &
                             (y(i,jedge(i))-y(i,jedge(i)-1))*lvis,     &
                                           (z(1)-z(0))*lvis,           &
                                          (z(km)-z(0))*lvis
    write(18,"(A)")'-----------------------------------------------------------------'
    close(18)
    print*,' << mesh_report.',irefname,'.dat ... done !'
    !
    open(18,file='Results/uprofile.is')
    do j=0,jm
    write(18,"(11(1X,E15.7E3))")y(iref,j)/bl_delta99(iref),yplus(j),(u1(is(n),j),n=1,9)
    enddo
    close(18)
    print*,' << Results/uprofile.is '
    !
    open(18,file='Results/omegazprofile.is')
    do j=0,jm
    write(18,"(11(1X,E15.7E3))")y(iref,j)/bl_delta99(iref),yplus(j),(omegaz(is(n),j),n=1,9)
    enddo
    close(18)
    print*,' << Results/omegazprofile.is '
    !
    ! free shear-layer core
    allocate(xfs(0:im),yfs(0:im),yfsup(0:im),yfsub(0:im),omzax(0:im))
    xfs(:)=x(:,0)
    yfs(:)=0.d0
    omzax=0.d0
    do i=1,im-100
      !
      var1=0.d0
      js=0
      !
      if(i<=708) then
        do j=1,178
          !
          if(omegaz(i,j)>=0.001d0 .and. &
             omegaz(i,j)>omegaz(i,j+1) .and. &
             omegaz(i,j)>omegaz(i,j-1) ) then
            !
            if(omegaz(i,j)>var1) then
              var1=omegaz(i,j)
              yfs(i)=y(i,j)
              js=j
            endif
            !
          endif
          !
        enddo
      else
        do j=15,178
          !
          if(omegaz(i,j)>=0.001d0 .and. &
             omegaz(i,j)>omegaz(i,j+1) .and. &
             omegaz(i,j)>omegaz(i,j-1) ) then
            !
            if(omegaz(i,j)>var1) then
              var1=omegaz(i,j)
              yfs(i)=y(i,j)
              js=j
            endif
            !
          endif
          !
        enddo
      endif
      !
      ! print*,i,yfs(i),omegaz(i,js)
      call quadfit(y(i,js-1),omegaz(i,js-1), &
                   y(i,js),  omegaz(i,js),   &
                   y(i,js+1),omegaz(i,js+1), yfs(i),omzax(i))
      !
      do j=js,jedge(i)+20
        !
        if(omegaz(i,j)<0.1d0*omzax(i)) then
          yfsup(i)=linear1d( omegaz(i,j-1),omegaz(i,j),  &
                                  y(i,j-1),     y(i,j),0.1d0*omzax(i))
          exit
        endif
        !
      enddo
      !
      var1=1.d10
      do j=1,js
        !
        if(omegaz(i,j)<var1) then
          var1=omegaz(i,j)
          jsub=j
        endif
        !
      enddo
      call quadfit(y(i,jsub-1),omegaz(i,jsub-1), &
                   y(i,jsub),  omegaz(i,jsub),   &
                   y(i,jsub+1),omegaz(i,jsub+1), yfsub(i),var1)
      !
      do j=js,1,-1
        !
        if(omegaz(i,j)<0.1d0*omzax(i)) then
          yfsub(i)=linear1d( omegaz(i,j),omegaz(i,j+1),  &
                                  y(i,j),     y(i,j+1),0.1d0*omzax(i))
          exit
        endif
        !
      enddo
      !
      ! print*,i,yfsub(i),yfs(i),yfsup(i)
      ! print*,'-----------------------------------'
      !
    enddo
    !
    open(18,file='Results/sharelayer_core.dat')
    write(18,"(A)")'VARIABLES = "x" "y" "omegaz"'
    write(18,"(A)")'ZONE T="ZONE 001"'
    write(18,"(A,I0,A)")' I=',im-99,',ZONETYPE=Ordered DATAPACKING=POINT'
    write(18,"(3(1X,E15.7E3))")(xfs(i),yfs(i),omzax(i),i=0,im-100)
    close(18)
    print*,' << Results/sharelayer_core.dat'
    !
    open(18,file='Results/shear-layer.dat')
    write(18,"(5(1X,A15))")'x','ybe','yc','yue','-omegaz'
    write(18,"(5(1X,E15.7E3))")(xfs(i),yfsub(i),yfs(i),yfsup(i),omzax(i),i=0,im-100)
    close(18)
    print*,' << Results/shear-layer.dat'
    !
    stop
    !
    ! do n=1,ns
    !   !
    !   i=is(n)
    !   write(iname,'(1hi,I4.4)')i
    !   !
    !   open(18,file='Results/profile.'//iname//'.dat')
    !   write(18,"(7(1X,A15))")'y/delta','u','v','ro','p','t','omegaz'
    !   write(18,"(7(1X,E15.7E3))")(y(i,j)/bl_delta99(is(1)),u1(i,j),u2(i,j),  &
    !                                          ro(i,j),p(i,j),t(i,j),omegaz(i,j),j=0,jm)
    !   close(18)
    !   print*,' << Results/profile.',iname,'.delta ... done. '
    !   !
    !   call upluscal(uplus=uplus,yplus=yplus,u=u1(i,:),y=y(i,:),              &
    !                                                      ro=ro(i,:),tw=t(i,0))
    !   !
    !   open(18,file='Results/uplus.'//iname//'.dat')
    !   write(18,"(2(1X,A15))")'yplus','uplus'
    !   write(18,"(2(1X,E15.7E3))")(yplus(j),uplus(j),j=0,jm)
    !   close(18)
    !   print*,' << Results/uplus.',iname,'.delta ... done. '
    !   !
    ! enddo
    !
    ! lref=bl_delta99(iref)
    !
    ! x=(x-x0)/lref
    ! y=y/lref
    ! z=z/lref
    !
    allocate(urev(0:im))
    do i=0,im
      !
      urev(i)=0.d0
      do j=0,jm
        !
        if(u1(i,j)<0.d0) then
          urev(i)=min(urev(i),u1(i,j))
        endif
        !
      enddo
      !
    enddo
    !
    open(18,file='Results/BLParameter.dat')
    write(18,"(9(1X,A15))")'x','Cf','pw','delta','delta*','theta','H','St','pwrms'
    write(18,"(9(1X,E15.7E3))")((x(i,0)-ximp)/lsep,cf(i),pwall(i),bl_delta99(i),     &
                                 bl_dstar(i),bl_theta(i),shapefac(i),st(i),   &
                                 sqrt(pp(i,0)/pp(is(1),0)),i=0,im-100)
    close(18)
    print*,' << BLParameter.dat ... done !'
    !
    open(18,file='Results/pwrms.dat')
    write(18,"(3(1X,A15))")'x','pwrms','ur'
    write(18,"(3(1X,E15.7E3))")(x(i,0),sqrt(pp(i,0)/pp(is(1),0)),urev(i),i=0,im-100)
    close(18)
    print*,' << pwrms.dat ... done !'
    !
    allocate(u11(0:im,0:jm),u22(0:im,0:jm),u33(0:im,0:jm),             &
             u12(0:im,0:jm),tke(0:im,0:jm),tt(0:im,0:jm))
    !
    call H5ReadArray(u11,im,jm,'u11','Results/2order.fav.zm.h5')
    call H5ReadArray(u22,im,jm,'u22','Results/2order.fav.zm.h5')
    call H5ReadArray(u33,im,jm,'u33','Results/2order.fav.zm.h5')
    call H5ReadArray(u12,im,jm,'u12','Results/2order.fav.zm.h5')
    call H5ReadArray( tt,im,jm, 'tt','Results/2order.fav.zm.h5')
    !
    tke=0.5d0*(u11+u22+u33)
    !
    call writetecbin('Results/tecrestress.plt',x,'x',y,'y',             &
                                 u11/utaw/utaw,'uu',u22/utaw/utaw,'vv', &
                                 u33/utaw/utaw,'ww',u12/utaw/utaw,'uv', &
                                 tke/utaw/utaw,'TKE',2.d0*sqrt(pp),'prms',im,jm)
    !
    i=iref
    open(18,file='Results/restress.'//irefname//'.dat')
    write(18,"(7(1X,A15))")'y/delta','yplus','uu','vv','ww','uv','k'
    do j=0,jm
      var1=ro(i,j)/ro(i,0)/utaw/utaw
      write(18,"(7(1X,E15.7E3))")y(i,j)/bl_delta99(i),yplus(j),u11(i,j)*var1,u22(i,j)*var1,u33(i,j)*var1,u12(i,j)*var1,tke(i,j)*var1
    enddo
    close(18)
    print*,' << Results/restress.',irefname,'.delta ... done. '
    !
    allocate(yuu(0:im),yvv(0:im),yww(0:im),yuv(0:im),yk(0:im),uup(0:im),vvp(0:im),wwp(0:im),uvp(0:im),tkep(0:im))
    !
    xfs(:)=x(:,0)
    !
    yuu=0.d0
    yvv=0.d0
    yww=0.d0
    yuv=0.d0
    !
    uup=0.d0
    vvp=0.d0
    wwp=0.d0
    uvp=0.d0
    !
    jeir=180
    do i=0,im-100
      !
      var1=0.d0
      js=0
      do j=1,jeir
        !
        if( u11(i,j)>var1 ) then
          !
          var1=u11(i,j)
          yuu(i)=y(i,j)
          js=j
          !
        endif
        !
      enddo
      !
      call quadfit(y(i,js-1),u11(i,js-1), &
                   y(i,js),  u11(i,js),   &
                   y(i,js+1),u11(i,js+1), yuu(i),uup(i))
      !
      var1=0.d0
      js=0
      do j=1,jeir
        !
        if(u22(i,j)>var1 ) then
          !
          var1=u22(i,j)
          yvv(i)=y(i,j)
          js=j
          !
        endif
        !
      enddo
      !
      call quadfit(y(i,js-1),u22(i,js-1), &
                   y(i,js),  u22(i,js),   &
                   y(i,js+1),u22(i,js+1), yvv(i),vvp(i))
      !
      var1=0.d0
      js=0
      do j=1,jeir
        !
        if(u33(i,j)>var1 ) then
          !
          var1=u33(i,j)
          yww(i)=y(i,j)
          js=j
          !
        endif
        !
      enddo
      !
      call quadfit(y(i,js-1),u33(i,js-1), &
                   y(i,js),  u33(i,js),   &
                   y(i,js+1),u33(i,js+1), yww(i),wwp(i))
      !
      var1=0.d0
      js=0
      do j=1,jeir
        !
        if(abs(u12(i,j))>var1 ) then
          !
          var1=abs(u12(i,j))
          yuv(i)=y(i,j)
          js=j
          !
        endif
        !
      enddo
      !
      call quadfit(y(i,js-1),u12(i,js-1), &
                   y(i,js),  u12(i,js),   &
                   y(i,js+1),u12(i,js+1), yuv(i),uvp(i))
      !
      var1=0.d0
      js=0
      do j=1,jeir
        !
        if(tke(i,j)>var1 ) then
          !
          var1=tke(i,j)
          yk(i)=y(i,j)
          js=j
          !
        endif
        !
      enddo
      !
      call quadfit(y(i,js-1),tke(i,js-1), &
                   y(i,js),  tke(i,js),   &
                   y(i,js+1),tke(i,js+1), yk(i),tkep(i))
      !
    enddo
    !
    var1=bl_delta99(iref)
    open(18,file='Results/restress_peak.dat')
    ! write(18,"(A)")'VARIABLES = "x" "yuu" "yvv" "yww" "yuv" uup" "vvp" "wwp" "uvp" '
    ! write(18,"(A)")'ZONE T="ZONE 001"'
    ! write(18,"(A,I0,A)")' I=',im-99,',ZONETYPE=Ordered DATAPACKING=POINT'

    write(18,"(11(1X,A))")'x','yuu','yvv','yww','yuv','ytke','uv','uu','vv','ww','uv','tke'
    write(18,"(11(1X,E15.7E3))")(xfs(i),yuu(i)/var1,yvv(i)/var1,yww(i)/var1,yuv(i)/var1,yk(i)/var1, &
                          uup(i)/uup(iref),vvp(i)/vvp(iref),wwp(i)/wwp(iref),uvp(i)/uvp(iref),tkep(i)/tkep(iref),i=0,im-100)
    close(18)
    print*,' << Results/restress_peak.dat'
    !
    open(18,file='Results/uuprofile.is')
    do j=0,jm
    write(18,"(11(1X,E15.7E3))")y(iref,j)/bl_delta99(iref),yplus(j),(u11(is(n),j)/utaw/utaw,n=1,9)
    enddo
    close(18)
    print*,' << Results/uuprofile.is '
    !
    open(18,file='Results/vvprofile.is')
    do j=0,jm
    write(18,"(11(1X,E15.7E3))")y(iref,j)/bl_delta99(iref),yplus(j),(u22(is(n),j)/utaw/utaw,n=1,9)
    enddo
    close(18)
    print*,' << Results/vvprofile.is '
    !
    open(18,file='Results/wwprofile.is')
    do j=0,jm
    write(18,"(11(1X,E15.7E3))")y(iref,j)/bl_delta99(iref),yplus(j),(u33(is(n),j)/utaw/utaw,n=1,9)
    enddo
    close(18)
    print*,' << Results/wwprofile.is '
    !
    open(18,file='Results/uvprofile.is')
    do j=0,jm
    write(18,"(11(1X,E15.7E3))")y(iref,j)/bl_delta99(iref),yplus(j),(u12(is(n),j)/utaw/utaw,n=1,9)
    enddo
    close(18)
    print*,' << Results/uvprofile.is '
    !
    ! reset jedge
    open(12,file='Results/yedge.dat')
    do i=0,im-100
      read(12,*)var1,yedge(i)
    enddo
    close(12)
    print*,' >> Results/yedge.dat'
    !
    do i=0,im-100
      do j=1,jm
        if( y(i,j)>=yedge(i)) then
          jedge(i)=j
          exit
        endif
      enddo
    enddo
    !
    ! drag decomposition
    allocate(cfv(0:im),cft(0:im),cfg(0:im),cfg1(0:im),cfg2(0:im),cfg3(0:im))
    allocate(cfg1u(0:im),cfg1v(0:im),cfg2res(0:im),cfg2vis(0:im))
    allocate(tauxy(0:im,0:jm),tauxx(0:im,0:jm),dudx(0:im,0:jm),dvdy(0:im,0:jm))
    !
    call H5ReadArray(tauxx,im,jm,'sgmam11','Results/budget.zm.h5')
    call H5ReadArray(tauxy,im,jm,'sgmam12','Results/budget.zm.h5')
    call H5ReadArray(dudx,im,jm,'du1dx','Results/mean.fav.zm.h5')
    call H5ReadArray(dvdy,im,jm,'du2dy','Results/mean.fav.zm.h5')
    !
    allocate(temp1(0:im,0:jm),dtemp1(1:2,0:im,0:jm), &
             temp2(0:im,0:jm),dtemp2(1:2,0:im,0:jm),dp(1:2,0:im,0:jm))
    !
    temp1=ro*u11
    dtemp1=grad(temp1,x,y)
    temp2=-tauxx
    dtemp2=grad(temp2,x,y)
    dp=grad(p,x,y)
    !
    cfv=0.d0
    cft=0.d0
    cfg1=0.d0
    cfg1u=0.d0
    cfg1v=0.d0
    cfg2=0.d0
    cfg3=0.d0
    cfg2res=0.d0
    cfg2vis=0.d0
    !
    do i=0,im-100
      !
      uedge=1.d0
      !
      do j=1,jedge(i)
        !
        dy=y(i,j)-y(i,j-1)
        !
        var1=tauxy(i,j)*dudy(i,j)
        var2=tauxy(i,j-1)*dudy(i,j-1)
        cfv(i)=cfv(i)+0.5d0*(var1+var2)*dy
        !
        var1=-ro(i,j)  *u12(i,j)  *dudy(i,j)
        var2=-ro(i,j-1)*u12(i,j-1)*dudy(i,j-1)
        cft(i)=cft(i)+0.5d0*(var1+var2)*dy
        !
        var1=(u1(i,j)-uedge)*ro(i,j)*(u1(i,j)*dudx(i,j)+u2(i,j)*dudy(i,j))
        var2=(u1(i,j-1)-uedge)*ro(i,j-1)*(u1(i,j-1)*dudx(i,j-1)+u2(i,j-1)*dudy(i,j-1))
        cfg1(i)=cfg1(i)+0.5d0*(var1+var2)*dy
        !
        var1=(u1(i,j)-uedge)*ro(i,j)*(u1(i,j)*dudx(i,j))
        var2=(u1(i,j-1)-uedge)*ro(i,j-1)*(u1(i,j-1)*dudx(i,j-1))
        cfg1u(i)=cfg1u(i)+0.5d0*(var1+var2)*dy
        !
        var1=(u1(i,j)-uedge)*ro(i,j)*(u2(i,j)*dudy(i,j))
        var2=(u1(i,j-1)-uedge)*ro(i,j-1)*(u2(i,j-1)*dudy(i,j-1))
        cfg1v(i)=cfg1v(i)+0.5d0*(var1+var2)*dy
        !
        var1=(u1(i,j)-uedge)*(dtemp1(1,i,j)+dtemp2(1,i,j))
        var2=(u1(i,j-1)-uedge)*(dtemp1(1,i,j-1)+dtemp2(1,i,j-1))
        cfg2(i)=cfg2(i)+0.5d0*(var1+var2)*dy
        !
        var1=(u1(i,j)-uedge)*(dtemp1(1,i,j))
        var2=(u1(i,j-1)-uedge)*(dtemp1(1,i,j-1))
        cfg2res(i)=cfg2res(i)+0.5d0*(var1+var2)*dy
        !
        var1=(u1(i,j)-uedge)*(dtemp2(1,i,j))
        var2=(u1(i,j-1)-uedge)*(dtemp2(1,i,j-1))
        cfg2vis(i)=cfg2vis(i)+0.5d0*(var1+var2)*dy
        !
        var1=(u1(i,j)-uedge)*dp(1,i,j)
        var2=(u1(i,j-1)-uedge)*dp(1,i,j-1)
        !
        cfg3(i)=cfg3(i)+0.5d0*(var1+var2)*dy
        !
      enddo
      !
      cfg(i)=cfg1(i)+cfg2(i)+cfg3(i)
      !
    enddo
    !
    cfv=2.d0*cfv
    cft=2.d0*cft
    cfg=2.d0*cfg
    cfg1=2.d0*cfg1
    cfg2=2.d0*cfg2
    cfg3=2.d0*cfg3
    !
    cfg1u=2.d0*cfg1u
    cfg1v=2.d0*cfg1v
    !
    cfg2res=2.d0*cfg2res
    cfg2vis=2.d0*cfg2vis
    !
    open(18,file='Results/cf_decomp.dat')
    write(18,"(13(1X,A15))")'x','Cf','Cfv','Cft','Cfg','Cfg1','Cfg2','Cfg3',    &
                           'cfsum','cfg1u','cfg1v','cfg2res','cfg2vis'
    write(18,"(13(1X,E15.7E3))")((x(i,0)-ximp)/lsep,cf(i),cfv(i),cft(i),cfg(i), &
                         cfg1(i),cfg2(i),cfg3(i),cfv(i)+cft(i)+cfg(i),cfg1u(i), &
                         cfg1v(i),cfg2res(i),cfg2vis(i),i=0,im-100)
    close(18)
    print*,' << Results/cf_decomp.dat '
    !
    ! budget terms
    allocate(convect(0:im,0:jm),pre_ace(0:im,0:jm),pre_dil(0:im,0:jm), &
             pre_tra(0:im,0:jm),product(0:im,0:jm),tur_tra(0:im,0:jm), &
             vis_ace(0:im,0:jm),vis_dif(0:im,0:jm),dissipa(0:im,0:jm), &
             balance(0:im,0:jm) )
    allocate(productxx(0:im,0:jm),productxy(0:im,0:jm), &
             productyy(0:im,0:jm))
    !
    call H5ReadArray(convect,im,jm,'convection','Results/tke_budget.zm.h5')
    call H5ReadArray(dissipa,im,jm,'dissipation','Results/tke_budget.zm.h5')
    call H5ReadArray(pre_ace,im,jm,'pressure_accelaration','Results/tke_budget.zm.h5')
    call H5ReadArray(pre_dil,im,jm,'pressure_dilatation','Results/tke_budget.zm.h5')
    call H5ReadArray(pre_tra,im,jm,'pressure_transport','Results/tke_budget.zm.h5')
    call H5ReadArray(product,im,jm,'production','Results/tke_budget.zm.h5')
    call H5ReadArray(tur_tra,im,jm,'turbulent_transport','Results/tke_budget.zm.h5')
    call H5ReadArray(vis_ace,im,jm,'viscous_accelaration','Results/tke_budget.zm.h5')
    call H5ReadArray(vis_dif,im,jm,'viscous_diffusion','Results/tke_budget.zm.h5')
    !
    balance=convect+pre_ace+pre_dil+pre_tra+product+tur_tra+vis_ace+vis_dif+dissipa
    !
    dissipa=dissipa-0.7d0*balance
    balance=convect+pre_ace+pre_dil+pre_tra+product+tur_tra+vis_ace+vis_dif+dissipa
    !
    miu=miucal(t(iref,0))/Reynolds
    var1=ro(iref,0)**2*utaw**4/miu
    !
    call writetecbin('Results/tecbudget.plt',x,'x',y,'y',        &
                                       convect/var1,'Ck',      &
                             (tur_tra+pre_tra)/var1,'Tk',      &
                                       product/var1,'Pk',      &
                                       vis_dif/var1,'Vk',      &
                     (pre_ace+pre_dil+vis_ace)/var1,'Πk',      &
                                       dissipa/var1,'ρε',      &
                                       balance/var1,'balance',im,jm)
    !
    do n=1,ns
      !
      i=is(n)
      write(iname,'(1hi,I4.4)')i
      !
      open(18,file='Results/budget_profile.'//iname//'.dat')
      write(18,"(9(1X,A15))")'y/delta','yplus','Ck','Tk','Pk',          &
                                                'Vk','Πk','ρε','balance'
      do j=0,jm
        write(18,"(9(1X,E15.7E3))")y(iref,j)/bl_delta99(iref),yplus(j),  &
                                                      convect(i,j)/var1, &
                                       (tur_tra(i,j)+pre_tra(i,j))/var1, &
                                                      product(i,j)/var1, &
                                                      vis_dif(i,j)/var1, &
                          (pre_ace(i,j)+pre_dil(i,j)+vis_ace(i,j))/var1, &
                                                      dissipa(i,j)/var1, &
                                                      balance(i,j)/var1
      end do
      close(18)
      print*,' << budget_profile.',iname,'.dat ... done. '
      !
    enddo
    !
    allocate(v1son(0:im),v2son(0:im),v3son(0:im),v4son(0:im),          &
             v5son(0:im),v6son(0:im))
    call binterp(x,y,convect,tur_tra+pre_tra,product,       &
                 xsonic,ysonic, v1son,v2son,v3son )
    call binterp(x,y,vis_dif,pre_ace+pre_dil+vis_ace,dissipa,       &
                 xsonic,ysonic, v4son,v5son,v6son )
    !
    open(18,file='Results/budget_sonic.dat')
    write(18,"(8(1X,A15))")'x','y','Ck','Tk','Pk','Vk','Πk','ρε'
    do i=0,im
      write(18,"(8(1X,E15.7E3))")xsonic(i),ysonic(i),  &
                                                    v1son(i)/var1, &
                                                    v2son(i)/var1, &
                                                    v3son(i)/var1, &
                                                    v4son(i)/var1, &
                                                    v5son(i)/var1, &
                                                    v6son(i)/var1 
    end do
    close(18)
    print*,' << Results/budget_sonic.dat ... done. '
    !
    do j=0,jm
    do i=0,im
      productxx(i,j)=-ro(i,j)*u11(i,j)*dudx(i,j)
      productyy(i,j)=-ro(i,j)*u22(i,j)*dvdy(i,j)
      productxy(i,j)=-ro(i,j)*u12(i,j)*(dudy(i,j)+dvdx(i,j))
    enddo
    enddo
    !
    call writetecbin('Results/tecproduction.plt',x,'x',y,'y',     &
                                       product/var1,'Pk',         &
                                       productxx/var1,'Pkx',      &
                                       productyy/var1,'Pky',      &
                                       productxy/var1,'Pks', im,jm)
    do n=1,ns
      !
      i=is(n)
      write(iname,'(1hi,I4.4)')i
      !
      open(18,file='Results/production_profile.'//iname//'.dat')
      write(18,"(6(1X,A15))")'y/delta','yplus','Pk','Pkx','Pky','Pks'
      do j=0,jm
        write(18,"(6(1X,E15.7E3))")y(iref,j)/bl_delta99(iref),yplus(j),  &
                                                      product(i,j)/var1, &
                                                    productxx(i,j)/var1, &
                                                    productyy(i,j)/var1, &
                                                    productxy(i,j)/var1
      end do
      close(18)
      print*,' << production_profile.',iname,'.dat ... done. '
      !
    enddo
    !
    allocate(c_dpdz_w(0:im,0:jm),ww(0:im,0:jm),dpz2(0:im,0:jm))
    !
    call H5ReadArray(c_dpdz_w,im,jm,'c_dpdz_w','Results/c_dpdz_w.h5')
    call H5ReadArray(      ww,im,jm,      'ww','Results/c_dpdz_w.h5')
    call H5ReadArray(    dpz2,im,jm,    'dpz2','Results/c_dpdz_w.h5')
    !
    do j=0,jm
    do i=0,im
      c_dpdz_w(i,j)=c_dpdz_w(i,j)/sqrt(ww(i,j)*dpz2(i,j))
    enddo
    enddo
    !
    open(18,file='Results/c_dpdz_w_profile.is')
    do j=0,jm
    write(18,"(11(1X,E15.7E3))")y(iref,j)/bl_delta99(iref),yplus(j),   &
                                               (c_dpdz_w(is(n),j),n=1,9)
    enddo
    close(18)
    print*,' << Results/c_dpdz_w_profile.is '
    !
    stop
    !
  end subroutine swtbli
  !
  subroutine tke_budget_swtbli(iref)
    !
    use commvardefine, only: x,y,im,jm,u1,ro,t
    !
    integer,intent(in) :: iref
    !
    ! real(8),allocatable,dimension(:,:) :: convect,pre_ace,pre_dil, &
    !                                       pre_tra,product,tur_tra, &
    !                                       vis_ace,vis_dif,dissipa, &
    !                                       balance
    ! real(8) :: yplus(0:jm),uplus(0:jm)
    ! real(8) :: budnorm,utaw
    ! !
    ! ! budget terms
    ! allocate(convect(0:im,0:jm),pre_ace(0:im,0:jm),pre_dil(0:im,0:jm), &
    !          pre_tra(0:im,0:jm),product(0:im,0:jm),tur_tra(0:im,0:jm), &
    !          vis_ace(0:im,0:jm),vis_dif(0:im,0:jm),dissipa(0:im,0:jm), &
    !          balance(0:im,0:jm) )
    ! !
    ! call H5ReadArray(convect,im,jm,'convection','Results/tke_budget.zm.h5')
    ! call H5ReadArray(dissipa,im,jm,'dissipation','Results/tke_budget.zm.h5')
    ! call H5ReadArray(pre_ace,im,jm,'pressure_accelaration','Results/tke_budget.zm.h5')
    ! call H5ReadArray(pre_dil,im,jm,'pressure_dilatation','Results/tke_budget.zm.h5')
    ! call H5ReadArray(pre_tra,im,jm,'pressure_transport','Results/tke_budget.zm.h5')
    ! call H5ReadArray(product,im,jm,'production','Results/tke_budget.zm.h5')
    ! call H5ReadArray(tur_tra,im,jm,'turbulent_transport','Results/tke_budget.zm.h5')
    ! call H5ReadArray(vis_ace,im,jm,'viscous_accelaration','Results/tke_budget.zm.h5')
    ! call H5ReadArray(vis_dif,im,jm,'viscous_diffusion','Results/tke_budget.zm.h5')
    ! !
    ! balance=convect+pre_ace+pre_dil+pre_tra+product+tur_tra+vis_ace+vis_dif+dissipa
    ! !
    ! dissipa=dissipa-0.5d0*balance
    ! balance=convect+pre_ace+pre_dil+pre_tra+product+tur_tra+vis_ace+vis_dif+dissipa
    ! !
    ! call upluscal(uplus=uplus,yplus=yplus,u=u1(i,:),y=y(i,:),          &
    !                                       ro=ro(i,:),tw=t(i,0),utaw=utaw)
    ! !
    ! print*,' ** utaw=',utaw
    ! !
    ! miu=miucal(t(iref,0))/Reynolds
    ! budnorm=ro(iref,0)**2*utaw**4/miu
    ! !
    ! call writetecbin('Results/tecbudget.plt',x,'x',y,'y',        &
    !                                    convect/budnorm,'Ck',      &
    !                          (tur_tra+pre_tra)/budnorm,'Tk',      &
    !                                    product/budnorm,'Pk',      &
    !                                    vis_dif/budnorm,'Vk',      &
    !                  (pre_ace+pre_dil+vis_ace)/budnorm,'Πk',      &
    !                                            dissipa,'ρε',      &
    !                                            balance,'balance',im,jm)
    ! !
    ! stop
    !
    ! i=iref
    ! open(18,file='Results/budget_profile.'//irefname//'.dat')
    ! write(18,"(12(1X,A15))")'y/delta','yplus','convection',            &
    !         'pressure_accelaration','pressure_dilatation',             &
    !         'pressure_transport','production','turbulent_transport',   &
    !         'viscous_accelaration','viscous_diffusion','dissipation',  &
    !         'balance'
    ! print*,' ** BL thickness=',bl_delta99(iref)
    ! do j=0,jm
    !   miu=miucal(t(iref,0))/Reynolds
    !   var1=ro(iref,0)**2*utaw**4/miu
    !   write(18,"(12(1X,E15.7E3))")y(iref,j)/bl_delta99(iref),yplus(j), &
    !             convect(i,j)/var1,pre_ace(i,j)/var1,pre_dil(i,j)/var1, &
    !             pre_tra(i,j)/var1,product(i,j)/var1,tur_tra(i,j)/var1, &
    !             vis_ace(i,j)/var1,vis_dif(i,j)/var1,dissipa(i,j)/var1, &
    !             balance(i,j)/var1
    ! end do
    ! close(18)
    ! print*,' << budget_profile.',irefname,'.dat ... done. '
    ! !
    ! i=739
    ! write(iname,'(1hi,I4.4)')i
    ! open(18,file='Results/budget_profile.'//iname//'.dat')
    ! write(18,"(12(1X,A15))")'y/delta','yplus','convection',            &
    !         'pressure_accelaration','pressure_dilatation',             &
    !         'pressure_transport','production','turbulent_transport',   &
    !         'viscous_accelaration','viscous_diffusion','dissipation',  &
    !         'balance'
    ! print*,' ** BL thickness=',bl_delta99(iref)
    ! do j=0,jm
    !   miu=miucal(t(iref,0))/Reynolds
    !   var1=ro(iref,0)**2*utaw**4/miu
    !   write(18,"(12(1X,E15.7E3))")y(iref,j)/bl_delta99(iref),yplus(j), &
    !             convect(i,j)/var1,pre_ace(i,j)/var1,pre_dil(i,j)/var1, &
    !             pre_tra(i,j)/var1,product(i,j)/var1,tur_tra(i,j)/var1, &
    !             vis_ace(i,j)/var1,vis_dif(i,j)/var1,dissipa(i,j)/var1, &
    !             balance(i,j)/var1
    ! end do
    ! close(18)
    ! print*,' << budget_profile.',iname,'.dat ... done. '
    ! !
    ! i=992
    ! write(iname,'(1hi,I4.4)')i
    ! open(18,file='Results/budget_profile.'//iname//'.dat')
    ! write(18,"(12(1X,A15))")'y/delta','yplus','convection',            &
    !         'pressure_accelaration','pressure_dilatation',             &
    !         'pressure_transport','production','turbulent_transport',   &
    !         'viscous_accelaration','viscous_diffusion','dissipation',  &
    !         'balance'
    ! print*,' ** BL thickness=',bl_delta99(iref)
    ! do j=0,jm
    !   miu=miucal(t(iref,0))/Reynolds
    !   var1=ro(iref,0)**2*utaw**4/miu
    !   write(18,"(12(1X,E15.7E3))")y(iref,j)/bl_delta99(iref),yplus(j), &
    !             convect(i,j)/var1,pre_ace(i,j)/var1,pre_dil(i,j)/var1, &
    !             pre_tra(i,j)/var1,product(i,j)/var1,tur_tra(i,j)/var1, &
    !             vis_ace(i,j)/var1,vis_dif(i,j)/var1,dissipa(i,j)/var1, &
    !             balance(i,j)/var1
    ! end do
    ! close(18)
    ! print*,' << budget_profile.',iname,'.dat ... done. '
    ! !
    ! i=1338
    ! write(iname,'(1hi,I4.4)')i
    ! open(18,file='Results/budget_profile.'//iname//'.dat')
    ! write(18,"(12(1X,A15))")'y/delta','yplus','convection',            &
    !         'pressure_accelaration','pressure_dilatation',             &
    !         'pressure_transport','production','turbulent_transport',   &
    !         'viscous_accelaration','viscous_diffusion','dissipation',  &
    !         'balance'
    ! print*,' ** BL thickness=',bl_delta99(iref)
    ! do j=0,jm
    !   miu=miucal(t(iref,0))/Reynolds
    !   var1=ro(iref,0)**2*utaw**4/miu
    !   write(18,"(12(1X,E15.7E3))")y(iref,j)/bl_delta99(iref),yplus(j), &
    !             convect(i,j)/var1,pre_ace(i,j)/var1,pre_dil(i,j)/var1, &
    !             pre_tra(i,j)/var1,product(i,j)/var1,tur_tra(i,j)/var1, &
    !             vis_ace(i,j)/var1,vis_dif(i,j)/var1,dissipa(i,j)/var1, &
    !             balance(i,j)/var1
    ! end do
    ! close(18)
    ! print*,' << budget_profile.',iname,'.dat ... done. '
    ! !
    ! i=iref
    ! open(18,file='Results/mesh_report.'//irefname//'.dat',position='append')
    ! !
    ! write(18,"(A)")' ratio of mesh to kolmogorov scale at wall          '
    ! write(18,"(A)")'----------------------------------------------------'
    ! write(18,"(A)")'         Δx/η         Δy/η         Δz/η         Δe/η'
    ! write(18,"(A)")'----------------------------------------------------'
    ! !
    ! miu=miucal(t(i,0))/Reynolds
    ! var1=sqrt(sqrt(miu**3/ro(i,0)**2/(-dissipa(i,0))))
    ! write(18,"(4(1X,F12.7))")0.5d0*(x(i+1,0)-x(i-1,0))/var1,           &
    !                                    (y(i,1)-y(i,0))/var1,           &
    !                                        (z(1)-z(0))/var1,           &
    !                      cube_root(0.5d0*(x(i+1,0)-x(i-1,0))*          &
    !                             (y(i,1)-y(i,0))*(z(1)-z(0)))/var1
    ! write(18,"(A)")'----------------------------------------------------'
    ! !
    ! write(18,"(A)")' ratio of mesh to kolmogorov scale at BL edge       '
    ! write(18,"(A)")'----------------------------------------------------'
    ! write(18,"(A)")'         Δx/η         Δy/η         Δz/η         Δe/η'
    ! write(18,"(A)")'----------------------------------------------------'
    ! !
    ! miu=miucal(t(i,jedge(i)))/Reynolds
    ! var1=sqrt(sqrt(miu**3/ro(i,jedge(i))**2/(-dissipa(i,jedge(i)))))
    ! write(18,"(4(1X,F12.7))")0.5d0*(x(i+1,jedge(i))-x(i-1,jedge(i)))/var1,   &
    !                          0.5d0*(y(i,jedge(i)+1)-y(i,jedge(i)-1))/var1,   &
    !                                              (z(1)-z(0))/var1,     &
    !                    cube_root(0.25d0*(x(i+1,jedge(i))-x(i-1,jedge(i)))*   &
    !                        (y(i,jedge(i)+1)-y(i,jedge(i)-1))*(z(1)-z(0)))/var1
    ! write(18,"(A)")'----------------------------------------------------'
    ! close(18)
    ! print*,' << mesh_report.',irefname,'.dat ... done !'
    !

    !
  end subroutine tke_budget_swtbli
  !
  subroutine comexpsl
    !
    use commvardefine,only: im,jm,km,Mach,Reynolds,pinf,pi,isp
    use interpolation 
    !
    integer :: is1,is2,i,j,k
    real(8) :: hstep,delta,xc
    real(8),allocatable,dimension(:,:) :: xsl1,ysl1,zsl1,xsl2,ysl2,zsl2,x2d,y2d
    real(8),allocatable,dimension(:,:) :: u1_zm,u2_zm,p_zm
    real(8),allocatable,dimension(:,:,:) :: u1_m,u2_m,p_m
    real(8),allocatable,dimension(:,:,:) :: u,v,w,p,x,y,z
    real(8),allocatable,dimension(:,:) :: u_sl1,v_sl1,w_sl1,p_sl1,   &
                                            u_sl2,v_sl2,w_sl2,p_sl2
    !
    hstep=2.0435385976288067
    xc=17.997071067811866
    delta=0.8488107d0
    !
    is1=1680
    is2=1310
    allocate( xsl1(0:is1,0:km),ysl1(0:is1,0:km),zsl1(0:is1,0:km),      &
              xsl2(0:is2,0:km),ysl2(0:is2,0:km),zsl2(0:is2,0:km) )
    !
    open(12,file='Results/stream_line1.dat')
    do i=1,17
      read(12,*)
    enddo
    do i=0,is1
      read(12,*)xsl1(i,0),ysl1(i,0)
      do j=1,9
        read(12,*)
      enddo
    enddo
    close(12)
    print*,' >> Results/stream_line1.dat'
    open(12,file='Results/stream_line2.dat')
    do i=1,16
      read(12,*)
    enddo
    do i=0,is2
      read(12,*)xsl2(i,0),ysl2(i,0)
      do j=1,9
        read(12,*)
      enddo
    enddo
    close(12)
    print*,' >> Results/stream_line2.dat'
    !
    ! allocate(x2d(0:im,0:jm),y2d(0:im,0:jm),z(0:km))
    ! call H5ReadSubset(x2d,im,jm,km,'x',gridfile,kslice=0)
    ! call H5ReadSubset(y2d,im,jm,km,'y',gridfile,kslice=0)
    ! call H5ReadSubset(z,im,jm,km,'z','datin/grid.h5',islice=0,jslice=0)
    !
    allocate(x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km))
    call H5ReadArray(x,im,jm,km,'x','datin/grid.h5')
    call H5ReadArray(y,im,jm,km,'y','datin/grid.h5')
    call H5ReadArray(z,im,jm,km,'z','datin/grid.h5')
    !
    x=(x-xc)/hstep
    y=y/hstep
    z=z/hstep
    !
    do i=0,is1
      xsl1(i,1:km)=xsl1(i,0)
      ysl1(i,1:km)=ysl1(i,0)
      !
      zsl1(i,:)=z(0,0,:)
    enddo
    !
    do i=0,is2
      xsl2(i,1:km)=xsl2(i,0)
      ysl2(i,1:km)=ysl2(i,0)
      !
      zsl2(i,:)=z(0,0,:)
    enddo
    !
    allocate( u(0:im,0:jm,0:km),v(0:im,0:jm,0:km),               &
              w(0:im,0:jm,0:km),p(0:im,0:jm,0:km),               &
              u1_zm(0:im,0:jm),u2_zm(0:im,0:jm),p_zm(0:im,0:jm) )
    allocate( u_sl1(0:is1,0:km),v_sl1(0:is1,0:km),w_sl1(0:is1,0:km),   &
              p_sl1(0:is1,0:km),u_sl2(0:is2,0:km),v_sl2(0:is2,0:km),   &
              w_sl2(0:is2,0:km),p_sl2(0:is2,0:km) )
    !
    call H5ReadArray(u1_zm,im,jm,'u1','Results/mean.fav.zm.h5')
    call H5ReadArray(u2_zm,im,jm,'u2','Results/mean.fav.zm.h5')
    call H5ReadArray(p_zm, im,jm, 'p','Results/mean.fav.zm.h5')
    
    ! call H5ReadArray(u,im,jm,km,'u1','Outdat/flowfield0222.h5')
    ! call H5ReadArray(v,im,jm,km,'u2','Outdat/flowfield0222.h5')
    ! call H5ReadArray(w,im,jm,km,'u3','Outdat/flowfield0222.h5')
    ! call H5ReadArray(p,im,jm,km, 'p','Outdat/flowfield0222.h5')
    call H5ReadArray(u,im,jm,km,'u1','Results/mean.fav.h5')
    call H5ReadArray(v,im,jm,km,'u2','Results/mean.fav.h5')
    call H5ReadArray(w,im,jm,km,'u3','Results/mean.fav.h5')
    call H5ReadArray(p,im,jm,km, 'p','Results/mean.fav.h5')
    !
    do k=0,km
      u(:,:,k)=u(:,:,k)-u1_zm(:,:)
      v(:,:,k)=v(:,:,k)-u2_zm(:,:)
      p(:,:,k)=p(:,:,k)- p_zm(:,:)
    enddo
    p=p/pinf
    !
    ! call writetecbin('Results/tecfluc.plt',x(996-100:996+100,:,:),'x', &
    !                                        y(996-100:996+100,:,:),'y', &
    !                                        z(996-100:996+100,:,:),'z', &
    !                                        u(996-100:996+100,:,:),'u"',&
    !                                        v(996-100:996+100,:,:),'v"',&
    !                                        w(996-100:996+100,:,:),'w"',&
    !                                        p(996-100:996+100,:,:),"p'",&
    !                                        200,jm,km)
    !
    call InvdisInterp2D( x(:,:,0),y(:,:,0),u,v,w,p,                    &
                         xsl1(:,0),ysl1(:,0),u_sl1,v_sl1,w_sl1,p_sl1 )
    !
    call writetecbin('Results/tec_dis_SL1.plt',xsl1,'x',ysl1,'y',zsl1,'z',  &
                     u_sl1,'u"',v_sl1,'v"',w_sl1,'w"',p_sl1,"p'",is1,km)
    !
    call InvdisInterp2D( x(:,:,0),y(:,:,0),u,v,w,p,                    &
                         xsl2(:,0),ysl2(:,0),u_sl2,v_sl2,w_sl2,p_sl2 )
    !
    call writetecbin('Results/tec_dis_SL2.plt',xsl2,'x',ysl2,'y',zsl2,'z',  &
                     u_sl2,'u"',v_sl2,'v"',w_sl2,'w"',p_sl2,"p'",is2,km)
    !
  end subroutine comexpsl
  !
  subroutine comexp3d
    !
    use commvardefine,only: im,jm,km,Mach,Reynolds,pinf,pi,isp
    use interpolation 
    use gradsolver
    !
    real(8),allocatable,dimension(:,:,:) :: x,y,z,u,v,w,p,t,dwdy,dvdz,   &
                                            um,vm,wm,pm,tm,ma,divp
    real(8),allocatable,dimension(:,:,:) :: x_s,y_s,z_s,u_s,v_s,w_s,   &
                                            p_s,dwdy_s,dvdz_s,omegax,divp_s
    real(8),allocatable,dimension(:,:,:,:) :: du1,du2,du3,gradp
    real(8),allocatable,dimension(:,:) :: u1_zm,u2_zm,p_zm,t_zm
    real(8),allocatable,dimension(:) :: dy
    real(8) :: hstep,delta,xc,ymin,ymax,var1
    integer :: nis,i1,j,k,js
    integer,allocatable :: is(:)
    ! 
    hstep=2.0435385976288067
    xc=17.997071067811866
    delta=0.8488107d0
    !
    allocate( dwdy(0:im,0:jm,0:km),dvdz(0:im,0:jm,0:km),               &
              u(0:im,0:jm,0:km),v(0:im,0:jm,0:km),                     &
              w(0:im,0:jm,0:km),p(0:im,0:jm,0:km),                     &
              um(0:im,0:jm,0:km),vm(0:im,0:jm,0:km),                   &
              wm(0:im,0:jm,0:km),pm(0:im,0:jm,0:km),                   &
              t(0:im,0:jm,0:km),tm(0:im,0:jm,0:km),ma(0:im,0:jm,0:km), &
              x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km),   &
              u1_zm(0:im,0:jm),u2_zm(0:im,0:jm),p_zm(0:im,0:jm),t_zm(0:im,0:jm) )
    !
    call H5ReadArray(x,im,jm,km,'x','datin/grid.h5')
    call H5ReadArray(y,im,jm,km,'y','datin/grid.h5')
    call H5ReadArray(z,im,jm,km,'z','datin/grid.h5')
    !
    x=(x-xc)/hstep
    y=y/hstep
    z=z/hstep
    !
    call H5ReadArray(u1_zm,im,jm,'u1','Results/mean.fav.zm.h5')
    call H5ReadArray(u2_zm,im,jm,'u2','Results/mean.fav.zm.h5')
    call H5ReadArray(p_zm, im,jm, 'p','Results/mean.fav.zm.h5')
    call H5ReadArray(t_zm, im,jm, 't','Results/mean.fav.zm.h5')
    !
    ! call H5ReadArray(dwdy,im,jm,km,'du3dy','Results/mean.fav.h5')
    ! call H5ReadArray(dvdz,im,jm,km,'du2dz','Results/mean.fav.h5')
    ! call H5ReadArray(um,im,jm,km,'u1','Results/mean.fav.h5')
    ! call H5ReadArray(vm,im,jm,km,'u2','Results/mean.fav.h5')
    ! call H5ReadArray(wm,im,jm,km,'u3','Results/mean.fav.h5')
    ! call H5ReadArray(pm,im,jm,km, 'p','Results/mean.fav.h5')
    ! call H5ReadArray(tm,im,jm,km, 't','Results/mean.fav.h5')
    !
    !
    call H5ReadArray(u,im,jm,km,'u1','Outdat/flowfield0222.h5')
    call H5ReadArray(v,im,jm,km,'u2','Outdat/flowfield0222.h5')
    call H5ReadArray(w,im,jm,km,'u3','Outdat/flowfield0222.h5')
    call H5ReadArray(p,im,jm,km, 'p','Outdat/flowfield0222.h5')
    !
    allocate(gradp(3,0:im,0:jm,0:km),divp(0:im,0:jm,0:km))
    gradp=grad(p,x,y,z)
    divp=sqrt(gradp(1,:,:,:)**2+gradp(1,:,:,:)**2+gradp(3,:,:,:)**2)
    ! p=p/pinf
    ! allocate(du2(1:3,0:im,0:jm,0:km),du3(1:3,0:im,0:jm,0:km))

    ! du2=grad_3d(v,x,y,z)
    ! du3=grad_3d(w)
    ! !
    ! dwdy=du3(2,:,:,:)
    ! dvdz=du2(3,:,:,:)
    !
    ma=sqrt((um**2+vm**2)/tm)*mach
    do k=0,km
      u(:,:,k)=u(:,:,k)-u1_zm
      v(:,:,k)=v(:,:,k)-u2_zm
      p(:,:,k)=(p(:,:,k)-p_zm)/pinf
    enddo
    !
    call writetecbin('Results/tec_disp_k160.plt',              &
                                  x(:,:,160),'x',              &
                                  y(:,:,160),'y',              &
                                  z(:,:,160),'z',              &
                                  um(:,:,160),'<u>',           &
                                  vm(:,:,160),'<v>',           &
                                  wm(:,:,160),'<w>',           &
                                  pm(:,:,160),'p',             &
                                  ma(:,:,160),'Ma',            &
                                  u(:,:,160),'u<',             &
                                  v(:,:,160),'v<',             &
                                  w(:,:,160),'w<',             &
                                  p(:,:,160),"p<",im,jm)
    call writetecbin('Results/tec_disp_k420.plt',              &
                                  x(:,:,420),'x',              &
                                  y(:,:,420),'y',              &
                                  z(:,:,420),'z',              &
                                  um(:,:,420),'<u>',           &
                                  vm(:,:,420),'<v>',           &
                                  wm(:,:,420),'<w>',           &
                                  pm(:,:,420),'p',             &
                                  ma(:,:,420),'Ma',            &
                                  u(:,:,420),'u<',             &
                                  v(:,:,420),'v<',             &
                                  w(:,:,420),'w<',             &
                                  p(:,:,420),"p<",im,jm)
    ! u=u-um
    ! v=v-vm
    ! w=w-wm
    ! p=(p-pm)/pinf
    !
    ! call writetecbin('Results/tecinst_dplus10.plt',                    &
    !                               x(700:im-100,10,:),'x',              &
    !                               y(700:im-100,10,:),'y',              &
    !                               z(700:im-100,10,:),'z',              &
    !                               u(700:im-100,10,:),'u"',             &
    !                               v(700:im-100,10,:),'v"',             &
    !                               w(700:im-100,10,:),'w"',             &
    !                               p(700:im-100,10,:),"p'",im-800,km)
    !
    nis=10
    js=240
    !
    allocate(is(0:nis))
    allocate(dy(0:js))
    allocate( x_s(0:nis,0:js,0:km),y_s(0:nis,0:js,0:km),               &
              z_s(0:nis,0:js,0:km),u_s(0:nis,0:js,0:km),               &
              v_s(0:nis,0:js,0:km),w_s(0:nis,0:js,0:km),               &
              p_s(0:nis,0:js,0:km),dwdy_s(0:nis,0:js,0:km),            &
              dvdz_s(0:nis,0:js,0:km),divp_s(0:nis,0:js,0:km),omegax(0:nis,0:js,0:km))
    !
    is(0)=996
    is(1)=1124
    is(2)=1377
    is(3)=1528
    is(4)=1635
    is(5)=1800
    is(6)=1895
    ! is(7)=1988
    is(7)=2020
    is(8)=2089
    is(9)=2182
    is(10)=2608
    !
    x_s(0,0:js,:)=x(is(0),0:js,:)
    y_s(0,0:js,:)=y(is(0),0:js,:)
    z_s(0,0:js,:)=z(is(0),0:js,:)
    !
    dy=y_s(0,:,0)
    !
    do i1=1,nis
      x_s(i1,0,:)=x(is(i1),0,:)
      y_s(i1,0,:)=y(is(i1),0,:)
      z_s(i1,0:js,:)=z(is(i1),0:js,:)
    enddo
    !
    do i1=1,nis
      !
      ! if( is(i1)>=1800 .and. is(i1)<=2089 ) then
      !   ! on the step
      !   print*,' ** stations on step',i1,x_s(i1,0,0),y_s(i1,0,0)
      !   !
      !   do j=1,js
      !     x_s(i1,j,:)=x(is(i1),0,:)-dy(j)*cos(0.25d0*pi)
      !     y_s(i1,j,:)=y(is(i1),0,:)+dy(j)*cos(0.25d0*pi)
      !   enddo
      ! else
        do j=1,js
          x_s(i1,j,:)=x(is(i1),0,:)
          y_s(i1,j,:)=y(is(i1),0,:)+dy(j)
        enddo
      ! endif
      !
    enddo
    !
    call InvdisInterp2D5_km(x(:,:,0),y(:,:,0),u,v,w,p,divp,         &
                            x_s(:,:,0),y_s(:,:,0),u_s,v_s,w_s,p_s,divp_s)
    !
    ! omegax=(dwdy_s-dvdz_s)*delta
    !
    call writetecbin('Results/tec3d_samp_fluc2.plt',x_s,'x',y_s,'y',z_s,'z',u_s,'u',   &
                                            v_s,'v',w_s,'w',p_s,'p',divp_s,'dp',nis,js,km)
    !

    !
  end subroutine comexp3d
  !
  subroutine cavity(iref)
    !
    use commvardefine,only: im,jm,km,Mach,Reynolds,pinf,pi,isp,gamma,gridfile,Mach,Prandtl
    use h5readwrite
    use basicfunction
    use interpolation 
    use boundarylayer

    use gradsolver
    !
    integer,intent(in) :: iref
    !
    real(8),allocatable,dimension(:,:) :: x,y,ro,u1,u2,u3,p,t,rns,     &
                                          dudx,dudy,dvdx,dvdy,omegaz,  &
                                          ma,divu,divp,u11,u22,u33,u12,&
                                          u13,u23,tt,pp,lumley_xi,     &
                                          lumley_eta,dtdx,dtdy
    real(8),allocatable,dimension(:,:) :: convect,pre_ace,pre_dil,     &
                                          pre_tra,product,tur_tra,     &
                                          vis_ace,vis_dif,dissipa,     &
                                          balance,wallnorm,pwxy_upper, &
                                          pwxy_lower,xlower,ylower,    &
                                          zlower,xupper,yupper,zupper, &
                                          prods,prodd
    real(8),allocatable,dimension(:) :: z,yplus,uplus,yinl,pw_upper,   &
                                        pw_lower,cf_upper,cf_lower,    &
                                        dp_lower,dp_upper,prms_upper,  &
                                        prms_lower
    real(8),allocatable :: gradp(:,:,:)
    real(8),allocatable,dimension(:,:,:) :: p3d,dudx_3d,dudy_3d,dudz_3d, &
                                           dvdx_3d,dvdy_3d,dvdz_3d,    &
                                           dwdx_3d,dwdy_3d,dwdz_3d,    &
                                           cf3d_upper,cf3d_lower
    real(8),allocatable :: xs(:,:),ys(:,:),yref(:),us(:,:),vs(:,:),ts(:,:),ptot(:,:)
    real(8),allocatable,dimension(:,:) :: u11_samp,u22_samp,u33_samp,prodc_samp,prods_samp,prodd_samp
    real(8),allocatable :: xint(:,:),yint(:,:),ptint(:,:),ptavg(:),area(:),ptcen(:)
    !
    integer :: jcav,icav,i,j,k,ii,jedge,i_cav_front,i_cav_rear,ixi,imxi
    real(8) :: utaw,lvis,bl_delta99,bl_dstar,bl_theta,var1,miu,     &
               shapefac,Retau,div,s11,s12,s13,s22,s23,s33,f11,f12,f13, &
               f22,f23,f33,tawx,tawy,ptref,ymin,ymax,a11,a12,a13,a22,a23,a33
    real(8) :: x_cav_front,x_cav_rear,Ttotal,wcav,dcav
    real(8),allocatable :: cf(:),xsi(:),chsi(:),twsi(:),st(:),qst(:),s(:)
    character(len=5) :: irefname,iname
    logical :: lfilalive
    !
    print*,' ** iref=',iref
    !
    ! x_cav_front=0.13812d0
    x_cav_front=0.14094d0
    x_cav_rear=0.1486d0
    !
    allocate(x(0:im,0:jm),y(0:im,0:jm),z(0:km),wallnorm(2,0:im))
    call H5ReadSubset(x,im,jm,km,'x',gridfile,kslice=0)
    call H5ReadSubset(y,im,jm,km,'y',gridfile,kslice=0)
    call H5ReadSubset(z,im,jm,km,'z',gridfile,islice=iref,jslice=0)
    !
    call gridjacobian_xy(x,y,wallnormal=wallnorm)
    !
    allocate(ro(0:im,0:jm),u1(0:im,0:jm),u2(0:im,0:jm),p(0:im,0:jm),   &
             t(0:im,0:jm),dudx(0:im,0:jm),dudy(0:im,0:jm),             &
             dvdx(0:im,0:jm),dvdy(0:im,0:jm),omegaz(0:im,0:jm),        &
             divu(0:im,0:jm),divp(0:im,0:jm),u3(0:im,0:jm),            &
             ptot(0:im,0:jm),dtdx(0:im,0:jm),dtdy(0:im,0:jm))
    !
    call H5ReadArray(ro,im,jm,'ro','Results/mean.fav.zm.h5')
    call H5ReadArray(u1,im,jm,'u1','Results/mean.fav.zm.h5')
    call H5ReadArray(u2,im,jm,'u2','Results/mean.fav.zm.h5')
    call H5ReadArray(u3,im,jm,'u3','Results/mean.fav.zm.h5')
    call H5ReadArray(p, im,jm, 'p','Results/mean.fav.zm.h5')
    call H5ReadArray(t, im,jm, 't','Results/mean.fav.zm.h5')
    call H5ReadArray(dudx,im,jm,'du1dx','Results/mean.fav.zm.h5')
    call H5ReadArray(dudy,im,jm,'du1dy','Results/mean.fav.zm.h5')
    call H5ReadArray(dvdx,im,jm,'du2dx','Results/mean.fav.zm.h5')
    call H5ReadArray(dvdy,im,jm,'du2dy','Results/mean.fav.zm.h5')
    call H5ReadArray(dtdx,im,jm,'dtdx','Results/mean.fav.zm.h5')
    call H5ReadArray(dtdy,im,jm,'dtdy','Results/mean.fav.zm.h5')
    !
    allocate(u11(0:im,0:jm),u22(0:im,0:jm),u33(0:im,0:jm),u12(0:im,0:jm))
    allocate(u13(0:im,0:jm),u23(0:im,0:jm),lumley_xi(0:im,0:jm),lumley_eta(0:im,0:jm))
    allocate(tt(0:im,0:jm),pp(0:im,0:jm))
    !
    ! call H5ReadArray(u11,im,jm,'u11','Results/2order.fav.zm.h5')
    ! call H5ReadArray(u22,im,jm,'u22','Results/2order.fav.zm.h5')
    ! call H5ReadArray(u33,im,jm,'u33','Results/2order.fav.zm.h5')
    ! call H5ReadArray(u12,im,jm,'u12','Results/2order.fav.zm.h5')
    ! call H5ReadArray(u13,im,jm,'u13','Results/2order.fav.zm.h5')
    ! call H5ReadArray(u23,im,jm,'u23','Results/2order.fav.zm.h5')
    ! !
    ! call H5ReadArray(tt,im,jm,'tt','Results/2order.fav.zm.h5')
    ! call H5ReadArray(pp,im,jm,'pp','Results/2order.fav.zm.h5')
    !
    allocate(gradp(2,0:im,0:jm),ma(0:im,0:jm))
    gradp=grad_xy(p,x,y)
    divp=sqrt(gradp(1,:,:)**2+gradp(2,:,:)**2)
    !
    var1=maxval(divp)
    divp=divp/var1
    !
    ma=sqrt((u1**2+u2**2)/t)*mach
    !
    p=p/pinf
    !
    ptot=p*(1.d0+0.5d0*(gamma-1.d0)*ma**2)**(gamma/(gamma-1.d0))
    ptref=(1.d0+0.5d0*(gamma-1.d0)*mach**2)**(gamma/(gamma-1.d0))
    !
    ptot=ptot/ptref
    !
    omegaz=dvdx(:,:)-dudy(:,:)
    !
    isp=im-20
    !
    call writetecbin('Results/tecmean.plt',x(0:isp,:),'x', y(0:isp,:),'y',  &
                                     u1(0:isp,:),'u',u2(0:isp,:),'v',  &
                                    ro(0:isp,:),'ro', p(0:isp,:),'p',  &
                                     t(0:isp,:),'T', ma(0:isp,:),'Ma', &
                  divp(0:isp,:),'divp',omegaz(0:isp,:),'omegaz',ptot(0:isp,:),'Pt',isp,jm)
    !
    do j=0,jm
    do i=0,im

      var1=u11(i,j)+u22(i,j)+u33(i,j)
      a11=u11(i,j)/var1-1.d0/3.d0
      a12=u12(i,j)/var1
      a13=u13(i,j)/var1
      a22=u22(i,j)/var1-1.d0/3.d0
      a23=u23(i,j)/var1
      a33=u33(i,j)/var1-1.d0/3.d0
      !
      lumley_xi(i,j)= a11*a11*a11+a12*a11*a12+a13*a11*a13+   &
                      a12*a12*a11+a22*a12*a12+a23*a12*a13+   &
                      a13*a13*a11+a23*a13*a12+a33*a13*a13+   &
                      a11*a12*a12+a12*a12*a22+a13*a12*a23+   &
                      a12*a22*a12+a22*a22*a22+a23*a22*a23+   &
                      a13*a23*a12+a23*a23*a22+a33*a23*a23+   &
                      a11*a13*a13+a12*a13*a23+a13*a13*a33+   &
                      a12*a23*a13+a22*a23*a23+a23*a23*a33+   &
                      a13*a33*a13+a23*a33*a23+a33*a33*a33
      lumley_xi(i,j)= cube_root(lumley_xi(i,j)/6.d0)
      !
      lumley_eta(i,j)= a11*a11+a12*a12+a13*a13+a12*a12+a22*a22+a23*a23+a13*a13+a23*a23+a33*a33
      !
      lumley_eta(i,j)= sqrt(lumley_eta(i,j)/6.d0)
      !
    enddo
    enddo
    !
    call writetecbin('Results/lumley_tri.plt',x,'x', y,'y',lumley_xi,'ξ',lumley_eta,'η',im,jm)
    !
    ! do i=0,im
    !   if((7.32d0-x(i,0))<1.466d0) then
    !     iref=i
    !     exit
    !   endif
    ! enddo
    print*,' ** the reference location is at i=',iref,'x=',x(iref,0)
    !
    i=0
    !
    do j=0,jm
      !
      if(y(i,j)>=0.d0) then
        jcav=j
        exit
      endif
      !
    enddo
    !
    do i=1,im
      !
      if(x(i-1,jcav)<=x_cav_front .and. x(i,jcav)>x_cav_front) then
        i_cav_front=i-1
      endif
      !
      if(x(i-1,jcav)<x_cav_rear .and. x(i,jcav)>=x_cav_rear) then
        i_cav_rear=i
      endif
      !
    enddo
    !
    print*,' ** the front wall of the cavity is at    i=',i_cav_front,'x=',x(i_cav_front,jcav)
    print*,' ** the rear wall of the cavity is at     i=',i_cav_rear, 'x=',x(i_cav_rear,jcav)
    print*,' ** the upper surface of the cavity is at j=',jcav,       'y=',y(i_cav_front,jcav)
    print*,' ** the bottom wall of the cavity is at   j=',0,          'y=',y(i_cav_front,0)
    !
    wcav=x(i_cav_rear,jcav)-x(i_cav_front,jcav)
    dcav=y(i_cav_front,jcav)-y(i_cav_front,0)
    !
    print*,' ** w=',wcav,'d=',dcav,'w/d=',wcav/dcav
    !
    imxi=-1
    do j=jcav,0,-1
      imxi=imxi+1
    enddo
    !
    do i=i_cav_front+1,i_cav_rear
      imxi=imxi+1
    enddo
    !
    do j=1,jcav
      imxi=imxi+1
    enddo
    !
    print*,' ** imxi=',imxi
    !
    allocate(xsi(0:imxi),chsi(0:imxi),twsi(0:imxi),st(0:imxi),cf(0:imxi),qst(0:imxi),s(0:imxi))
    !
    ttotal=1.d0+0.5d0*(gamma-1.d0)*Mach**2
    !
    ixi=-1
    xsi=0.d0
    ! front wall
    i=i_cav_front
    do j=jcav,0,-1
      ixi=ixi+1
      !
      miu=miucal(t(i,j))/Reynolds
      !
      if(j==jcav) then 
        xsi(ixi)=0.d0
      else
         xsi(ixi)=xsi(ixi-1)+y(i,j+1)-y(i,j)
      endif
      twsi(ixi)=t(i,j)
      !
      chsi(ixi)=2.d0*miu*dtdx(i,j)/(Mach**2*Prandtl*(gamma-1))
        st(ixi)=miu*dtdx(i,j)/(ttotal-twsi(ixi))/Prandtl

        cf(ixi)=2.d0*miu*dvdx(i,j)
      !
    enddo
    !
    ! bottom wall
    j=0
    do i=i_cav_front+1,i_cav_rear
      ixi=ixi+1
      !
      miu=miucal(t(i,j))/Reynolds
      !
      twsi(ixi)=t(i,j)
      !
       xsi(ixi)=xsi(ixi-1)+x(i,j)-x(i-1,j)
      chsi(ixi)=2.d0*miu*dtdy(i,j)/(Mach**2*Prandtl*(gamma-1))
        st(ixi)=miu*dtdy(i,j)/(ttotal-twsi(ixi))/Prandtl

        cf(ixi)=2.d0*miu*dudy(i,j)
      !
    enddo
    !
    ! rear wall
    i=i_cav_rear
    do j=1,jcav
      ixi=ixi+1
      !
      miu=miucal(t(i,j))/Reynolds
      !
      twsi(ixi)=t(i,j)
       xsi(ixi)=xsi(ixi-1)+y(i,j)-y(i,j-1)
      chsi(ixi)=-2.d0*miu*dtdx(i,j)/(Mach**2*Prandtl*(gamma-1))
        st(ixi)= -miu*dtdx(i,j)/(ttotal-twsi(ixi))/Prandtl
      !
        cf(ixi)=-2.d0*miu*dvdx(i,j)
      !
    enddo
    !
    do i=1,imxi
      s(i)=xsi(imxi)-xsi(i)
      qst(i)=qstheory(s=s(i),w=wcav,d=dcav)
      ! print*,i,s(i),qst(i)
    enddo
    !
    open(18,file='Results/qwall.dat')
    write(18,"(6(1X,A15))")'Xs','tw','ch','st','cf','q_theory'
    write(18,"(6(1X,E15.7E3))")(xsi(i)/xsi(imxi),twsi(i),chsi(i),st(i),cf(i),qst(i),i=0,imxi)
    close(18)
    print*,' << Results/qwall.dat'
    !
    deallocate(xsi,chsi,twsi,st,cf)
    !
    allocate(xsi(0:im),chsi(0:im),twsi(0:im),st(0:im),cf(0:im))
    !
    ! j=jcav
    j=0
    do i=0,im
      !
      miu=miucal(t(i,j))/Reynolds
      twsi(i)=t(i,j)
      chsi(i)=2.d0*miu*dtdy(i,j)/(Mach**2*Prandtl*(gamma-1))
        st(i)=     miu*dtdy(i,j)/(ttotal-twsi(i))/Prandtl
        cf(i)=2.d0*miu*dudy(i,j)
      !
    enddo
    !
    open(18,file='Results/BLParameter.dat')
    write(18,"(5(1X,A15))")'s','tw','ch','st','cf'
    write(18,"(5(1X,E15.7E3))")(x(i,j),twsi(i),chsi(i),st(i),cf(i),i=0,im)
    close(18)
    print*,' << Results/BLParameter.dat'
    !
    open(18,file='Results/profile.dat')
    write(18,"(6(1X,A15))")'y','u','v','T','ro','p'
    do j=0,jm
      write(18,"(6(1X,E15.7E3))")y(iref,j)-y(iref,jcav),u1(iref,j),u2(iref,j),t(iref,j),ro(iref,j),p(iref,j)
    end do
    close(18)
    print*,' << Results/profile.dat'
    !

    stop
    !
    allocate( xlower(0:im+jcav+1,0:km),ylower(0:im+jcav+1,0:km),       &
              zlower(0:im+jcav+1,0:km),pwxy_lower(0:im+jcav+1,0:km),   &
              cf3d_lower(3,0:im+jcav+1,0:km))
    !
    allocate( p3d(0:im,0:jm,0:km),dudx_3d(0:im,0:jm,0:km),             &
              dudy_3d(0:im,0:jm,0:km),dudz_3d(0:im,0:jm,0:km),         &
              dvdx_3d(0:im,0:jm,0:km),dvdy_3d(0:im,0:jm,0:km),         &
              dvdz_3d(0:im,0:jm,0:km),dwdx_3d(0:im,0:jm,0:km),         &
              dwdy_3d(0:im,0:jm,0:km),dwdz_3d(0:im,0:jm,0:km))
    !
    call H5ReadArray(p3d,    im,jm,km,'p','Results/mean.fav.h5')
    call H5ReadArray(dudx_3d,im,jm,km,'du1dx','Results/mean.fav.h5')
    call H5ReadArray(dudy_3d,im,jm,km,'du1dy','Results/mean.fav.h5')
    call H5ReadArray(dudz_3d,im,jm,km,'du1dz','Results/mean.fav.h5')
    call H5ReadArray(dvdx_3d,im,jm,km,'du2dx','Results/mean.fav.h5')
    call H5ReadArray(dvdy_3d,im,jm,km,'du2dy','Results/mean.fav.h5')
    call H5ReadArray(dvdz_3d,im,jm,km,'du2dz','Results/mean.fav.h5')
    call H5ReadArray(dwdx_3d,im,jm,km,'du3dx','Results/mean.fav.h5')
    call H5ReadArray(dwdy_3d,im,jm,km,'du3dy','Results/mean.fav.h5')
    call H5ReadArray(dwdz_3d,im,jm,km,'du3dz','Results/mean.fav.h5')
    !
    p3d=p3d/pinf
    !
    xlower(0:icav,0)=x(0:icav,jcav)
    ylower(0:icav,0)=y(0:icav,jcav)
    !
    pwxy_lower(0:icav,:)=p3d(0:icav,jcav,:)
    !
    print*,' ** part 0'
    !
    do i=0,icav
      !
      miu=miucal(t(i,jcav))/Reynolds
      !
      do k=0,km
        !
        div=(dudx_3d(i,jcav,k)+dvdy_3d(i,jcav,k)+dwdz_3d(i,jcav,k))/3.d0
        !
        s11=dudx_3d(i,jcav,k)-div
        s12=0.5d0*(dudy_3d(i,jcav,k)+dvdx_3d(i,jcav,k))
        s13=0.5d0*(dudz_3d(i,jcav,k)+dwdx_3d(i,jcav,k))
        !
        s22=dvdy_3d(i,jcav,k)-div
        s23=0.5d0*(dvdz_3d(i,jcav,k)+dwdy_3d(i,jcav,k))
        !
        s33=dwdz_3d(i,jcav,k)-div
        !
        f11=2.d0*miu*s11; f12=2.d0*miu*s12; f13=2.d0*miu*s13
                          f22=2.d0*miu*s22; f23=2.d0*miu*s23
                                            f33=2.d0*miu*s33
        !
        cf3d_lower(1,i,k)=(f11*wallnorm(1,i)+f12*wallnorm(2,i))*2.d0
        cf3d_lower(2,i,k)=(f12*wallnorm(1,i)+f22*wallnorm(2,i))*2.d0
        cf3d_lower(3,i,k)=(f13*wallnorm(1,i)+f23*wallnorm(2,i))*2.d0
        !
      enddo
      !
    enddo
    !
    print*,' ** part 1'
    !
    ii=0
    do i=icav+1,icav+jcav+1
      ii=ii+1
      !
      xlower(i,0)=x(icav,jcav)
      ylower(i,0)=y(icav,jcav-ii)
      !
      pwxy_lower(i,:)=p3d(icav,jcav-ii,:)
      !
      miu=miucal(t(i,jcav-ii))/Reynolds
      !
      do k=0,km
      !   !
      !   div=(dudx_3d(icav,jcav-ii,k)+dvdy_3d(icav,jcav-ii,k)+dwdz_3d(icav,jcav-ii,k))/3.d0
      !   !
      !   s11=dudx_3d(icav,jcav-ii,k)-div
      !   s12=0.5d0*(dudy_3d(icav,jcav-ii,k)+dvdx_3d(icav,jcav-ii,k))
      !   s13=0.5d0*(dudz_3d(icav,jcav-ii,k)+dwdx_3d(icav,jcav-ii,k))
      !   !
      !   s22=dvdy_3d(icav,jcav-ii,k)-div
      !   s23=0.5d0*(dvdz_3d(icav,jcav-ii,k)+dwdy_3d(icav,jcav-ii,k))
      !   !
      !   s33=dwdz_3d(icav,jcav-ii,k)-div
      !   !
      !   f11=2.d0*miu*s11; f12=2.d0*miu*s12; f13=2.d0*miu*s13
      !                     f22=2.d0*miu*s22; f23=2.d0*miu*s23
      !                                       f33=2.d0*miu*s33
      !   !
        cf3d_lower(1,i,k)=0.d0 !f11*2.d0
        cf3d_lower(2,i,k)=0.d0 !f12*2.d0
        cf3d_lower(3,i,k)=0.d0 !f13*2.d0
      !   !
      enddo
      !
    enddo
    !
    print*,' ** part 2'
    !
    do i=icav+1,im
      !
      xlower(i+jcav+1,0)=x(i,0)
      ylower(i+jcav+1,0)=y(i,0)
      !
      pwxy_lower(i+jcav+1,:)=p3d(i,0,:)
      !
      miu=miucal(t(i,0))/Reynolds
      !
      do k=0,km
        !
        div=(dudx_3d(i,0,k)+dvdy_3d(i,0,k)+dwdz_3d(i,0,k))/3.d0
        !
        s11=dudx_3d(i,0,k)-div
        s12=0.5d0*(dudy_3d(i,0,k)+dvdx_3d(i,0,k))
        s13=0.5d0*(dudz_3d(i,0,k)+dwdx_3d(i,0,k))
        !
        s22=dvdy_3d(i,0,k)-div
        s23=0.5d0*(dvdz_3d(i,0,k)+dwdy_3d(i,0,k))
        !
        s33=dwdz_3d(i,0,k)-div
        !
        f11=2.d0*miu*s11; f12=2.d0*miu*s12; f13=2.d0*miu*s13
                          f22=2.d0*miu*s22; f23=2.d0*miu*s23
                                            f33=2.d0*miu*s33
        !
        cf3d_lower(1,i+jcav+1,k)=(f11*wallnorm(1,i)+f12*wallnorm(2,i))*2.d0
        cf3d_lower(2,i+jcav+1,k)=(f12*wallnorm(1,i)+f22*wallnorm(2,i))*2.d0
        cf3d_lower(3,i+jcav+1,k)=(f13*wallnorm(1,i)+f23*wallnorm(2,i))*2.d0
        !
      enddo
      !
    enddo
    !
    print*,' ** part 3'
    !
    do k=0,km
      !
      xlower(:,k)=xlower(:,0)
      ylower(:,k)=ylower(:,0)
      zlower(:,k)=z(k)
      !
    enddo
    !
    allocate( xupper(0:im,0:km),yupper(0:im,0:km),zupper(0:im,0:km),   &
              pwxy_upper(0:im,0:km),cf3d_upper(1:3,0:im,0:km))
    !
    do i=0,im
      xupper(i,:)=x(i,jm)
      yupper(i,:)=y(i,jm)
      zupper(i,:)=z(:)
      !
      pwxy_upper(i,:)=p3d(i,jm,:)
      !
      do k=0,km
        !
        div=(dudx_3d(i,jm,k)+dvdy_3d(i,jm,k)+dwdz_3d(i,jm,k))/3.d0
        !
        s11=dudx_3d(i,jm,k)-div
        s12=0.5d0*(dudy_3d(i,jm,k)+dvdx_3d(i,jm,k))
        s13=0.5d0*(dudz_3d(i,jm,k)+dwdx_3d(i,jm,k))
        !
        s22=dvdy_3d(i,jm,k)-div
        s23=0.5d0*(dvdz_3d(i,jm,k)+dwdy_3d(i,jm,k))
        !
        s33=dwdz_3d(i,jm,k)-div
        !
        f11=2.d0*miu*s11; f12=2.d0*miu*s12; f13=2.d0*miu*s13
                          f22=2.d0*miu*s22; f23=2.d0*miu*s23
                                            f33=2.d0*miu*s33
        !
        cf3d_upper(1,i,k)=(-f12*wallnorm(2,i))*2.d0
        cf3d_upper(2,i,k)=(-f22*wallnorm(2,i))*2.d0
        cf3d_upper(3,i,k)=(-f23*wallnorm(2,i))*2.d0
        !
      enddo
      !
    enddo
    !
    call writetecbin('Results/tec_lowerwall.plt',xlower,'x',           &
                                     ylower,'y',zlower,'z',            &
                                     pwxy_lower,'p',                   &
                                     cf3d_lower(1,:,:),'cfx',          &
                                     cf3d_lower(2,:,:),'cfy',          &
                                     cf3d_lower(3,:,:),'cfz',im+jcav+1,km)
    call writetecbin('Results/tec_upperwall.plt',xupper,'x',           &
                                     yupper,'y',zupper,'z',            &
                                     pwxy_upper,'p',                   &
                                     cf3d_upper(1,:,:),'cfx',          &
                                     cf3d_upper(2,:,:),'cfy',          &
                                     cf3d_upper(3,:,:),'cfz',im,km)
    !
    allocate(pw_upper(0:im),pw_lower(0:im),           &
             cf_upper(0:im),cf_lower(0:im),dp_lower(0:im),             &
             dp_upper(0:im),prms_upper(0:im),prms_lower(0:im))
    !
    do i=0,im
      !
      if(i<=icav) then
        !
        j=jcav
        !
        miu=miucal(t(i,j))/Reynolds
        cf_lower(i)=2.d0*miu*dudy(i,j)
        dp_lower(i)=gradp(1,i,j)
        !
        prms_lower(i)=2.d0*sqrt(pp(i,j))
      else
        !
        j=0
        !
        miu=miucal(t(i,j))/Reynolds
        !
        div=0.333333333333333d0*(dudx(i,j)+dvdy(i,j))
        s11=dudx(i,j)-div
        s12=0.5d0*(dudy(i,j)+dvdx(i,j))
        s22=dvdy(i,j)-div
        !
        f11=2.d0*miu*s11
        f12=2.d0*miu*s12
        f22=2.d0*miu*s22
        !
        tawx=f11*wallnorm(1,i)+f12*wallnorm(2,i)
        tawy=f12*wallnorm(1,i)+f22*wallnorm(2,i)
        !
        ! taw(i)=sqrt(tawx**2+tawy**2)
        cf_lower(i)=2.d0*(tawx*wallnorm(2,i)-tawy*wallnorm(1,i))
        !
        dp_lower(i)=gradp(1,i,j)*wallnorm(2,i)-gradp(2,i,j)*wallnorm(1,i)
        !
      endif
      !
      pw_lower(i)=p(i,j)
      prms_lower(i)=2.d0*sqrt(pp(i,j))
      !
      miu=miucal(t(i,jm))/Reynolds
      cf_upper(i)=2.d0*miu*dudy(i,jm)
      pw_upper(i)=p(i,jm)
      !
      dp_upper(i)=gradp(1,i,jm)
      !
      prms_upper(i)=2.d0*sqrt(pp(i,jm))
      !
    enddo
    !
    open(18,file='Results/wall_variable.dat')
    write(18,"(9(1X,A15))")'x','pw_lower','pw_upper','cf_lower',       &
                           'cf_upper','dp/dx_lower','dp/dx_upper',     &
                           'prms_lower','prms_upper'
    write(18,"(9(1X,E15.7E3))")(x(i,0),pw_lower(i),pw_upper(i),        &
                                         cf_lower(i),-cf_upper(i),     &
                                         dp_lower(i),dp_upper(i),      &
                                         prms_lower(i),prms_upper(i),i=0,im)
    close(18)
    print*,' << Results/wall_variable.dat'
    !
    !
    i=100
    write(irefname,'(1hi,I4.4)')i
    !
    bl_delta99=blthickness(y(i,:),u1(i,:),jedge=jedge,btype='nominal')
    bl_delta99=bl_delta99-y(i,jcav)
    print*,' ** the boundary layer thickness is : ',bl_delta99
    !
    bl_dstar=blthickness(y(i,jcav:jcav+100),u1(i,jcav:jcav+100),ro=ro(i,jcav:jcav+100),btype='displacement')
    bl_theta=blthickness(y(i,jcav:jcav+100),u1(i,jcav:jcav+100),ro=ro(i,jcav:jcav+100),btype='momentum')
    !
    print*,' ** the edge of the boundary layer is at j=',jedge,'y=',y(i,jedge)
    !
    allocate(uplus(0:100),yplus(0:100),yinl(0:100))
    !
    yinl(0:100)=y(i,jcav:jcav+100)-y(i,jcav)
    !
    call upluscal(uplus=uplus,yplus=yplus,u=u1(i,jcav:jcav+100), &
                                             y=yinl,                   &
                                            ro=ro(i,jcav:jcav+100), &
                                         tw=t(i,jcav),utaw=utaw)
    !
    miu=miucal(t(i,jcav))/Reynolds
    cf=miu*dudy(i,jcav)*2.d0
    shapefac=bl_dstar/bl_theta
    Retau=ro(i,jcav)*utaw*bl_delta99/miu
    !
    open(18,file='Results/BLParameter'//irefname//'.dat')
    write(18,"(A,F12.7,A,I4)")'boundary layer paramters at x=',x(i,0),', i=',i
    write(18,"(10(1X,A15))")'Cf','delta','delta*','theta','H','Re_delta','Re_delta*','Re_theta','Re_tau','utau'
    write(18,"(10(1X,E15.7E3))")cf,bl_delta99,bl_dstar,bl_theta,shapefac,     &
                               Reynolds*bl_delta99,Reynolds*bl_dstar,Reynolds*bl_theta,Retau,utaw
    close(18)
    print*,' << BLParameter.dat ... done !'
    !
    !
    open(18,file='Results/uplus.dat')
    write(18,"(2(1X,A15))")'yplus','uplus'
    write(18,"(2(1X,E15.7E3))")(yplus(j),uplus(j),j=0,100)
    close(18)
    print*,' << uplus.dat ... done. '
    !
    open(18,file='Results/uprofile.dat')
    write(18,"(2(1X,E15.7E3))")(yinl(j)/bl_delta99,u1(iref,j+jcav),j=0,100)
    close(18)
    print*,' << Results/uprofile.dat'
    !
    lvis=ro(i,jcav)*utaw/miucal(t(i,jcav))*Reynolds
    !
    open(18,file='Results/mesh_report.'//irefname//'.dat')
    write(18,"(A,F12.7,A,I4)")'mesh resolution at x=',x(i,0),', i=',i
    write(18,"(A)")'----------------------------------------------------'
    write(18,"(A)")'          Δx+         Δy1+         Δye+          Δz+'
    write(18,"(A)")'----------------------------------------------------'
    write(18,"(4(1X,F12.7))")0.5d0*(x(i+1,0)-x(i,0))*lvis,           &
                                     (y(i,jcav+1)-y(i,jcav))*lvis,           &
                               (y(i,jedge)-y(i,jedge-1))*lvis,         &
                                           (z(1)-z(0))*lvis
    write(18,"(A)")'----------------------------------------------------'
    close(18)
    print*,' << mesh_report.',irefname,'.dat ... done !'
    !
    inquire(file='Results/tke_budget.zm.h5',exist=lfilalive)
    !
    if(lfilalive) then
      !
      allocate(dissipa(0:im,0:jm))
      !
      call H5ReadArray(dissipa,im,jm,'dissipation','Results/tke_budget.zm.h5')
      !
      open(18,file='Results/mesh_report.'//irefname//'.dat',position='append')
      !
      write(18,"(A)")' ratio of mesh to kolmogorov scale at wall          '
      write(18,"(A)")'----------------------------------------------------'
      write(18,"(A)")'         Δx/η         Δy/η         Δz/η         Δe/η'
      write(18,"(A)")'----------------------------------------------------'
      !
      miu=miucal(t(i,0))/Reynolds
      var1=sqrt(sqrt(miu**3/ro(i,jcav)**2/(-dissipa(i,jcav))))
      write(18,"(4(1X,F12.7))")0.5d0*(x(i+1,jcav)-x(i-1,jcav))/var1,           &
                                       (y(i,jcav+1)-y(i,jcav))/var1,           &
                                                   (z(1)-z(0))/var1,           &
                           cube_root(0.5d0*(x(i+1,jcav)-x(i-1,jcav))*          &
                                           (y(i,jcav+1)-y(i,jcav))*(z(1)-z(0)))/var1
      write(18,"(A)")'----------------------------------------------------'
      !
      write(18,"(A)")' ratio of mesh to kolmogorov scale at BL edge       '
      write(18,"(A)")'----------------------------------------------------'
      write(18,"(A)")'         Δx/η         Δy/η         Δz/η         Δe/η'
      write(18,"(A)")'----------------------------------------------------'
      !
      miu=miucal(t(i,jedge))/Reynolds
      var1=sqrt(sqrt(miu**3/ro(i,jedge)**2/(-dissipa(i,jedge))))
      write(18,"(4(1X,F12.7))")0.5d0*(x(i+1,jedge)-x(i-1,jedge))/var1,     &
                               0.5d0*(y(i,jedge+1)-y(i,jedge-1))/var1,     &
                                                     (z(1)-z(0))/var1,     &
                         cube_root(0.25d0*(x(i+1,jedge)-x(i-1,jedge))*     &
                             (y(i,jedge+1)-y(i,jedge-1))*(z(1)-z(0)))/var1
      write(18,"(A)")'----------------------------------------------------'
      close(18)
      print*,' << mesh_report.',irefname,'.dat ... done !'
      !
    endif
    !
    do j=0,jm
    do i=0,im
      var1=max(abs(pp(i,j)),1.d-10)
      pp(i,j)=2.d0*sqrt(var1)
    enddo
    enddo
    !
    call writetecbin('Results/tecres.plt',x(0:isp,:),'x', y(0:isp,:),'y',  &
                                     u11(0:isp,:)/utaw/utaw,'uu',           &
                                     u22(0:isp,:)/utaw/utaw,'vv',           &
                                     u33(0:isp,:)/utaw/utaw,'ww',           &
                                     u12(0:isp,:)/utaw/utaw,'uv',           &
                     0.5d0*(u11(0:isp,:)+u22(0:isp,:)+u33(0:isp,:))/utaw/utaw,'K',    &
                     tt(0:isp,:),'TT',pp(0:isp,:),'prms',isp,jm)
    !
    open(18,file='Results/restress.dat')
    write(18,"(2(1X,A15))")'y/delta','yplus','urms','vrms','wrms','uv'
    write(18,"(6(1X,E15.7E3))")(yinl(j)/bl_delta99,yplus(j),sqrt(u11(iref,j+jcav))/utaw, &
                                                            sqrt(u22(iref,j+jcav))/utaw, &
                                                            sqrt(u33(iref,j+jcav))/utaw, &
                                                            u12(iref,j+jcav)/utaw/utaw,j=0,100)
    close(18)
    print*,' << Results/restress.dat'
    !
    open(18,file='Results/lumley_profile.dat')
    write(18,"(3(1X,A15))")'y/delta','xi','eta'
    write(18,"(3(1X,E15.7E3))")(yinl(j)/bl_delta99,lumley_xi(iref,j+jcav),lumley_eta(iref,j+jcav),j=0,100)
    close(18)
    print*,' << Results/lumley_profile.dat'
    !
    ! budget terms
    allocate(convect(0:im,0:jm),pre_ace(0:im,0:jm),pre_dil(0:im,0:jm), &
             pre_tra(0:im,0:jm),product(0:im,0:jm),tur_tra(0:im,0:jm), &
             vis_ace(0:im,0:jm),vis_dif(0:im,0:jm),balance(0:im,0:jm) )
    !
    call H5ReadArray(convect,im,jm,'convection','Results/tke_budget.zm.h5')
    call H5ReadArray(dissipa,im,jm,'dissipation','Results/tke_budget.zm.h5')
    call H5ReadArray(pre_ace,im,jm,'pressure_accelaration','Results/tke_budget.zm.h5')
    call H5ReadArray(pre_dil,im,jm,'pressure_dilatation','Results/tke_budget.zm.h5')
    call H5ReadArray(pre_tra,im,jm,'pressure_transport','Results/tke_budget.zm.h5')
    call H5ReadArray(product,im,jm,'production','Results/tke_budget.zm.h5')
    call H5ReadArray(tur_tra,im,jm,'turbulent_transport','Results/tke_budget.zm.h5')
    call H5ReadArray(vis_ace,im,jm,'viscous_accelaration','Results/tke_budget.zm.h5')
    call H5ReadArray(vis_dif,im,jm,'viscous_diffusion','Results/tke_budget.zm.h5')
    !
    balance=convect+pre_ace+pre_dil+pre_tra+product+tur_tra+vis_ace+vis_dif+dissipa
    !
    open(18,file='Results/budget_profile.'//irefname//'.dat')
    write(18,"(12(1X,A15))")'y/delta','yplus','convection',            &
            'pressure_accelaration','pressure_dilatation',             &
            'pressure_transport','production','turbulent_transport',   &
            'viscous_accelaration','viscous_diffusion','dissipation',  &
            'balance'
    do j=0,100
      miu=miucal(t(i,jcav))/Reynolds
      var1=ro(i,jcav)**2*utaw**4/miu
      write(18,"(12(1X,E15.7E3))")yinl(j)/bl_delta99,yplus(j),         &
                convect(i,j+jcav)/var1,pre_ace(i,j+jcav)/var1,pre_dil(i,j+jcav)/var1, &
                pre_tra(i,j+jcav)/var1,product(i,j+jcav)/var1,tur_tra(i,j+jcav)/var1, &
                vis_ace(i,j+jcav)/var1,vis_dif(i,j+jcav)/var1,dissipa(i,j+jcav)/var1, &
                balance(i,j+jcav)/var1
    end do
    close(18)
    print*,' << budget_profile.',irefname,'.dat ... done.'
    !
    allocate(prods(0:im,0:jm),prodd(0:im,0:jm))
    miu=miucal(t(iref,jcav))/Reynolds
    var1=ro(iref,jcav)**2*utaw**4/miu
    do j=0,jm
    do i=0,im
      product(i,j)=product(i,j)/var1
      prods(i,j)=-ro(i,j)*u12(i,j)*(dudy(i,j)+dvdx(i,j))/var1
      prodd(i,j)=-ro(i,j)*(u11(i,j)*dudx(i,j)+u22(i,j)*dvdy(i,j))/var1
    enddo
    enddo
    !
    call writetecbin('Results/tec_production.plt',x,'x',y,'y',     &
                          product,'PK',prods,'PS',prodd,'PD',im,jm)
    !
    allocate(xs(8,0:jm),ys(8,0:jm),yref(0:jm))
    !
    do j=0,jm
      yref(j)=(y(im,j)-y(im,0))/(y(im,jm)-y(im,0))
    enddo
    !
    xs(1,:)=5.89644d0
    ys(1,0)=3.05
    xs(2,:)=14.d0
    ys(2,0)=0.d0
    xs(3,:)=22.8258d0
    ys(3,0)=2.17871d0
    xs(4,:)=33.0248d0
    ys(4,0)=3.05d0
    xs(5,:)=40.3764d0
    ys(5,0)=3.05d0
    xs(6,:)=48.401d0
    ys(6,0)=3.05d0
    xs(7,:)=56.6327d0
    ys(7,0)=3.05d0
    xs(8,:)=64.2949d0
    ys(8,0)=3.05d0
    !
    do i=1,9
     do j=1,jm
       ys(i,j)=(17.71d0-ys(i,0))*yref(j)+ys(i,0)
     enddo
    enddo
    !
    allocate(us(8,0:jm),vs(8,0:jm),ts(8,0:jm))
    allocate(u11_samp(8,0:jm),u22_samp(8,0:jm),u33_samp(8,0:jm))
    allocate(prodc_samp(8,0:jm),prods_samp(8,0:jm),prodd_samp(8,0:jm))
    !
    call InvdisInterp2D(x,  y,   u1,u2, t,     u11,     u22,    u33,im,jm, &
                        xs,ys,   us,vs,ts,u11_samp,u22_samp,u33_samp,7,jm)
    !
    call InvdisInterp2D(x,  y,   product,prods,prodd,im,jm, &
                        xs,ys,   prodc_samp,prods_samp,prodd_samp,7,jm)
    !
    open(18,file='Results/uprofile.dat')
    do j=0,jm
      write(18,"(16(1X,E15.7E3))")(ys(i,j),us(i,j),i=1,8)
    enddo
    close(18)
    print*,' << Results/uprofile.dat'
    !
    open(18,file='Results/Tprofile.dat')
    do j=0,jm
      write(18,"(16(1X,E15.7E3))")(ys(i,j),ts(i,j),i=1,8)
    enddo
    close(18)
    print*,' << Results/Tprofile.dat'
    !
    open(18,file='Results/u11_profile.dat')
    do j=0,jm
      write(18,"(16(1X,E15.7E3))")(ys(i,j),u11_samp(i,j)/utaw/utaw,i=1,8)
    enddo
    close(18)
    print*,' << Results/u11_profile.dat'
    !
    open(18,file='Results/u22_profile.dat')
    do j=0,jm
      write(18,"(16(1X,E15.7E3))")(ys(i,j),u22_samp(i,j)/utaw/utaw,i=1,8)
    enddo
    close(18)
    print*,' << Results/u22_profile.dat'
    !
    open(18,file='Results/u33_profile.dat')
    do j=0,jm
      write(18,"(16(1X,E15.7E3))")(ys(i,j),u33_samp(i,j)/utaw/utaw,i=1,8)
    enddo
    close(18)
    print*,' << Results/u33_profile.dat'
    !
    open(18,file='Results/production_profile.dat')
    write(18,"(32(1X,A15))")('y','PK','PS','PD',i=1,8)
    do j=0,jm
      write(18,"(32(1X,E15.7E3))")(ys(i,j),prodc_samp(i,j),prods_samp(i,j),prodd_samp(i,j),i=1,8)
    enddo
    close(18)
    print*,' << Results/production_profile.dat'
    !
    stop
    !
    allocate(xint(0:im,0:jm),yint(0:im,0:jm),ptint(0:im,0:jm),ptavg(0:im),ptcen(0:im),area(0:im))
    do i=0,im
      !
      if(i<=icav) then
        ymin=3.05d0
      else
        ymin=y(i,0)
      endif
      ymax=y(i,jm)
      !
      xint(i,:)=x(i,0)
      do j=0,jm
        yint(i,j)=(ymax-ymin)*yref(j)+ymin
      enddo
      !
    enddo
    !
    call InvdisInterp2D(x,      y, ptot,  im,jm,        &
                        xint,yint, ptint, im,jm)
    !
    call writetecbin('Results/tecint.plt',xint,'x',yint,'y',ptint,'Pt',im,jm)
    !
    ptavg=0.d0
    area=0.d0
    do i=0,im
      !
      do j=1,jm
        var1=abs(yint(i,j)-yint(i,j-1))
        area(i)=area(i)+var1
        ptavg(i)=ptavg(i)+0.5d0*(ptint(i,j)+ptint(i,j-1))*var1
      enddo
      !
    enddo
    !
    open(18,file='Results/pt_avg.dat')
    write(18,"(4(1X,E15.7E3))")(xint(i,0),ptavg(i),area(i),ptavg(i)/area(i),i=0,im)
    close(18)
    print*,' << Results/pt_avg.hstep'
    !
    contains
    !
    function qstheory(s,w,d) result(qs)
      !
      real(8) :: qs,s,w,d
      !
      qs=0.21d0/sqrt(1.d0+d/w)*(xi(0.5d0*s/(w+d))-xi(0.5d0*(s+w)/(w+d)))
      !
    end function qstheory
    !
    function xi(v)
      !
      real(8) :: v,xi
      !
      xi=1.d0/sqrt(v)+ 0.803323d0 -3.89728d0*(v+1.d0) +  2.55002d0*(v+1.d0)**2 + &
         -1.19121*(v+1.d0)**3+0.308284d0*(v+1.d0)**4 -0.0335024*(v+1.d0)**5
      !
    end function xi
    !
  end subroutine cavity
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used as the post-process for                   |
  !|  ecompression-expansion corner case.                              |
  !+-------------------------------------------------------------------+
  subroutine comexp(iref)
    !
    use commvardefine,only: im,jm,km,Mach,Reynolds,pinf,pi,isp
    use h5readwrite
    use basicfunction
    use interpolation 
    use boundarylayer

    use gradsolver
    !
    ! argument
    integer,intent(in) :: iref
    !
    real(8),parameter :: epslion=1.d-9
    ! local data
    !
    real(8),allocatable,dimension(:,:) :: x,y,ro,u1,u2,u3,p,t,rns,     &
                                          dudx,dudy,dvdx,dvdy,omegaz,ma,divu,divp
    real(8),allocatable,dimension(:) :: z,xs,yplus,uplus,taw,pwall,cf,ppw,dy,dpw
    real(8),allocatable,dimension(:,:) :: wallnorm
    real(8),allocatable,dimension(:) :: yref,uref,xexp,xexp2
    real(8),allocatable,dimension(:,:) :: usamp,vsamp,rosamp,xsamp,    &
                                          ysamp,usamp_max,usamp_min,   &
                                          masamp
    real(8),allocatable,dimension(:,:,:) :: u1_tm,usamp_tm,gradp
    integer,allocatable,dimension(:) :: is
    real(8),allocatable,dimension(:,:) :: xwall,ywall,zwall,twall,     &
                                    dudx_wall,dudy_wall,dudz_wall,     &
                                    dvdx_wall,dvdy_wall,dvdz_wall,     &
                                    dwdx_wall,dwdy_wall,dwdz_wall,     &
                                    taw_tan,pres_wall,prms_wall
    real(8),allocatable,dimension(:,:,:) :: taw_wall,dp
    real(8),allocatable,dimension(:,:) :: convect,pre_ace,pre_dil,     &
                                          pre_tra,product,tur_tra,     &
                                          vis_ace,vis_dif,dissipa,     &
                                          balance,prodxx,prodyy,prodxy,&
                                          prods,prodd,prodc,curv
    real(8),allocatable,dimension(:,:) :: u11,u22,u33,u12,tke,us,un,   &
                                          uss,unn,usn,dusds,dusdn,pp,  &
                                          u11_dis,u22_dis,u33_dis,     &
                                          u11_hom,u22_hom,u33_hom,     &
                                          u11_samp,u22_samp,u33_samp,  &
                                          u12_samp,uss_samp,unn_samp,  &
                                          usn_samp,pu1,pu2,pu3,        &
                                          pu1_samp,pu2_samp,pu3_samp
    real(8),allocatable,dimension(:,:,:) :: u11_3d,u22_3d,u33_3d
    real(8),allocatable,dimension(:,:,:) :: dus
    integer :: i,j,k,ic,ie,jedge,i1,i2,i3,i4,i5,nis,n,js
    real(8) :: xc,yc,xe,ye,alfa,hstep,lstep,ymin,ymax,xsep,xatt,lsep,  &
               xxsep,xxatt,lxsep,ppwmax
    real(8) :: bl_delta99,bl_dstar,bl_theta,omegaz_ref,utaw,lvis,      &
               pwref,cfref,shapefac,rops1,ups1,vps1,pps1,Tps1
    real(8) :: div,s11,s12,s13,s22,s23,s33,f11,f12,f13,f22,f23,f33,    &
               var1,var2,miu,tawx,tawy,delta,delta1,delta_exp,p1,p2
    real(8) :: delta99(0:im),pinv(0:im)
    real(8) :: prodmax(0:im),xprodmax(0:im),yprodmax(0:im),            &
               prodsmax(0:im),prodcmax(0:im),proddmax(0:im),           &
               x_tkemax(0:im),y_tkemax(0:im),tkemax(0:im)
    real(8) :: dgrid(0:jm)
    real(8) :: rotmat(2,2),bmat(2,2),budref
    integer :: iu,ju,iv,jv,iw,jw,ik,jk
    real(8) :: uumax,vvmax,wwmax,kmax,uumax_ref,vvmax_ref,wwmax_ref,   &
               kmax_ref
    !
    character(len=5) :: irefname,iname
    !
    !+-------+
    !| begin |
    !+-------+
    !
    print*,' ** processing data of a compression-expansion corner flow.'
    !
    allocate(x(0:im,0:jm),y(0:im,0:jm),z(0:km))
    call H5ReadSubset(x,im,jm,km,'x','datin/grid.h5',kslice=0)
    call H5ReadSubset(y,im,jm,km,'y','datin/grid.h5',kslice=0)
    call H5ReadSubset(z,im,jm,km,'z','datin/grid.h5',islice=iref,jslice=0)
    !
    ! to find the compression corner
    ! alfa=45.d0/180.d0*pi
    !
    do i=1,im
      if(y(i,0)>=0.d0 .and. abs(y(i-1,0))<=epslion) then
        xc=x(i-1,0)
        yc=y(i-1,0)
        ic=i-1
      endif
      if(y(i,0)>y(i-1,0)) then
        xe=x(i,0)
        ye=y(i,0)
        ie=i
      endif
    enddo
    !
    hstep=ye-yc
    lstep=xe-xc
    !
    alfa=atan(hstep/lstep)
    !
    print*,' ** compression corner at: (',xc,yc,') at i= ',ic
    print*,' ** expansion  corner  at: (',xe,ye,') at i= ',ie
    print*,' ** height of the step: ',hstep
    print*,' ** x-length of the step: ',lstep
    print*,' ** length of the step: ',lstep/cos(alfa),hstep/sin(alfa)
    print*,' ** deflection angle  : ',alfa/pi*180.d0
    !
    allocate(xs(0:im))
    xs(0)=0.d0
    do i=1,im
      var1=sqrt((x(i,0)-x(i-1,0))**2+(y(i,0)-y(i-1,0))**2)
      xs(i)=xs(i-1)+var1
    enddo
    xs=xs-xc
    !
    allocate(ro(0:im,0:jm),u1(0:im,0:jm),u2(0:im,0:jm),p(0:im,0:jm),   &
             t(0:im,0:jm),dudx(0:im,0:jm),dudy(0:im,0:jm),             &
             dvdx(0:im,0:jm),dvdy(0:im,0:jm),omegaz(0:im,0:jm),        &
             divu(0:im,0:jm),divp(0:im,0:jm),u3(0:im,0:jm))
    !
    call H5ReadArray(ro,im,jm,'ro','Results/mean.fav.zm.h5')
    call H5ReadArray(u1,im,jm,'u1','Results/mean.fav.zm.h5')
    call H5ReadArray(u2,im,jm,'u2','Results/mean.fav.zm.h5')
    call H5ReadArray(u3,im,jm,'u3','Results/mean.fav.zm.h5')
    call H5ReadArray(p, im,jm, 'p','Results/mean.fav.zm.h5')
    call H5ReadArray(t, im,jm, 't','Results/mean.fav.zm.h5')
    call H5ReadArray(dudx,im,jm,'du1dx','Results/mean.fav.zm.h5')
    call H5ReadArray(dudy,im,jm,'du1dy','Results/mean.fav.zm.h5')
    call H5ReadArray(dvdx,im,jm,'du2dx','Results/mean.fav.zm.h5')
    call H5ReadArray(dvdy,im,jm,'du2dy','Results/mean.fav.zm.h5')
    !
    allocate( u11(0:im,0:jm),u22(0:im,0:jm),u33(0:im,0:jm),            &
              u12(0:im,0:jm),tke(0:im,0:jm),pp(0:im,0:jm),             &
              uss(0:im,0:jm),unn(0:im,0:jm) )
    call H5ReadArray(u11,im,jm,'u11','Results/2order.fav.zm.h5')
    call H5ReadArray(u12,im,jm,'u12','Results/2order.fav.zm.h5')
    call H5ReadArray(u22,im,jm,'u22','Results/2order.fav.zm.h5')
    call H5ReadArray(u33,im,jm,'u33','Results/2order.fav.zm.h5')
    !
    allocate(pu1(0:im,0:jm),pu2(0:im,0:jm),pu3(0:im,0:jm))
    call H5ReadArray(u11,im,jm,'u11','Results/2order.rey.zm.h5')
    call H5ReadArray(u12,im,jm,'u12','Results/2order.rey.zm.h5')
    call H5ReadArray(u22,im,jm,'u22','Results/2order.rey.zm.h5')
    call H5ReadArray( pp,im,jm, 'pp','Results/2order.rey.zm.h5')
    call H5ReadArray(pu1,im,jm,'pu1','Results/2order.rey.zm.h5')
    call H5ReadArray(pu2,im,jm,'pu2','Results/2order.rey.zm.h5')
    call H5ReadArray(pu3,im,jm,'pu3','Results/2order.rey.zm.h5')
    !
    do j=1,jm
    do i=0,im
      pu1(i,j)=pu1(i,j)/sqrt(pp(i,j)*u11(i,j))
      pu2(i,j)=pu2(i,j)/sqrt(pp(i,j)*u22(i,j))
      pu3(i,j)=pu3(i,j)/sqrt(pp(i,j)*u33(i,j))
    enddo
    enddo
    !
    call writetecbin('Results/tecpu.plt',(x-xc)/hstep,'x',y/hstep,'y', &
                                       pu1,'pu',pu2,'pv',pu3,'pw',im,jm)
    ! allocate(u11_3d(0:im,0:jm,0:km),u22_3d(0:im,0:jm,0:km), &
    !          u33_3d(0:im,0:jm,0:km),u11_dis(0:im,0:jm),     &
    !          u22_dis(0:im,0:jm),u33_dis(0:im,0:jm),u11_hom(0:im,0:jm), &
    !          u22_hom(0:im,0:jm),u33_hom(0:im,0:jm))
    ! !
    ! u11_hom=0.d0
    ! u22_hom=0.d0
    ! u33_hom=0.d0
    ! call H5ReadArray(u11_3d,im,jm,km,'u11','Results/2order.fav.h5')
    ! do k=1,km
    !   u11_hom=u11_hom+u11_3d(:,:,k)
    ! enddo
    ! call H5ReadArray(u11_3d,im,jm,km,'u22','Results/2order.fav.h5')
    ! do k=1,km
    !   u22_hom=u22_hom+u11_3d(:,:,k)
    ! enddo
    ! call H5ReadArray(u11_3d,im,jm,km,'u33','Results/2order.fav.h5')
    ! do k=1,km
    !   u33_hom=u33_hom+u11_3d(:,:,k)
    ! enddo
    ! u11_hom=u11_hom/dble(km)
    ! u22_hom=u22_hom/dble(km)
    ! u33_hom=u33_hom/dble(km)
    ! !
    ! u11_dis=u11-u11_hom
    ! u22_dis=u22-u22_hom
    ! u33_dis=u33-u33_hom
    !
    omegaz=dvdx-dudy
    divu=(dudx+dvdy)
    !
    allocate(dp(2,0:im,0:jm))
    dp=grad_xy(p,x,y)
    divp=sqrt(dp(1,:,:)**2+dp(2,:,:)**2)
    !
    allocate(rns(0:im,0:jm))
    rns=rnsxy(ro)
    !
    allocate(ma(0:im,0:jm))
    ma=sqrt((u1**2+u2**2)/t)*mach
    !
    i=iref
    !
    write(irefname,'(1hi,I4.4)')iref
    !
    bl_delta99=blthickness(y(i,:),u1(i,:),jedge=jedge,btype='nominal')
    !
    bl_dstar=blthickness(y(i,:),u1(i,:),ro=ro(i,:),omegaz=omegaz(i,:), &
                                                   btype='displacement')
    bl_theta=blthickness(y(i,:),u1(i,:),ro=ro(i,:),omegaz=omegaz(i,:), &
                                                       btype='momentum')
    !
    shapefac=bl_dstar/bl_theta
    write(*,"(A,F7.3,A)")'  ---------- BL thickness at x=',x(iref,0),' -----------'
    write(*,"(A)")'            δ          δ*           θ           H'
    write(*,"(1X,4(F12.7))")bl_delta99,bl_dstar,bl_theta,shapefac
    write(*,"(A)")'  -----------------------------------------------'
    write(*,"(A)")'         δ*/δ        θ/δ         h/δ'
    write(*,"(1X,3(F12.7))")bl_dstar/bl_delta99,bl_theta/bl_delta99,hstep/bl_delta99
    write(*,"(A)")'  -----------------------------------------------'
    !
    print*,' ** reference location at i: ',i
    print*,' **          distance to cc: ', (x(iref,0)-xc)/hstep,'/h'
    !
    print*,' ** length of the domain: ',x(isp,0)/bl_delta99,x(im,0)/bl_delta99
    print*,' ** height of the domain: ',y(0,jm)/bl_delta99
    print*,' ** height of the step: ',hstep/bl_delta99
    print*,' ** length of the step: ',lstep/bl_delta99
    print*,' ** compression corner at: (',xc/bl_delta99,yc/bl_delta99,')'
    print*,' ** expansion  corner  at: (',xe/bl_delta99,ye/bl_delta99,')'
    !
    call writetecbin('Results/tecmean_iref.plt',(x-xc)/bl_delta99,'x',y/bl_delta99,'y', &
                                       ro,'ro',u1,'u',u2,'v',p,'p',t,'T',rns,'rns',im,jm)
    !
    divu=divu*bl_delta99
    divp=divp/0.5d0*bl_delta99
    !
    omegaz=omegaz*bl_delta99
    !
    call writetecbin('Results/tecmean.plt',(x-xc)/hstep,'x',y/hstep,'y',      &
                                       ro,'ro',u1,'u',u2,'v',p/pinf,'p',      &
                                       t,'T',ma,'Ma',rns,'rns',divu,'div',    &
                                       divp,'dP',omegaz,'omegaz',im,jm)
    !
    allocate(yplus(0:jm),uplus(0:jm))
    !
    call upluscal(uplus=uplus,yplus=yplus,u=u1(i,:),y=y(i,:),          &
                                         ro=ro(i,:),tw=t(i,0),utaw=utaw)
    !
    print*,' ** utaw=',utaw
    print*,' ** miu=',miucal(t(i,0))/Reynolds
    print*,' ** row=',ro(i,0)
    !
    write(*,"(A,F7.3,A)")'  ---------- Reynolds number at x=',x(iref,0),' -----------'
    write(*,"(A)")'                 Reδ               Reδ*                Reθ                Reτ'
    write(*,"(1X,4(1X,F18.7))")Reynolds*bl_delta99,Reynolds*bl_dstar,Reynolds*bl_theta, &
                               Reynolds*ro(iref,0)*bl_delta99*utaw/miucal(t(iref,0))
    write(*,"(A)")'  -----------------------------------------------'
    !
    open(18,file='Results/uplus.'//irefname//'.dat')
    write(18,"(2(1X,A15))")'yplus','uplus'
    write(18,"(2(1X,E15.7E3))")(yplus(j),uplus(j),j=0,jm)
    close(18)
    print*,' << uplus.',irefname,'.dat ... done. '
    !
    open(18,file='Results/profile.'//irefname//'.dat')
    write(18,"(6(1X,A15))")'y/delta','u','v','ro','p','t'
    write(18,"(6(1X,E15.7E3))")(y(i,j)/bl_delta99,u1(i,j),u2(i,j),ro(i,j),p(i,j),t(i,j),j=0,jm)
    close(18)
    print*,' << Results/profile.',irefname,'.delta ... done. '
    !
    call writeupluslay('uplus.'//irefname//'.dat',                     &
                                              'Results/uplus-yplus.lay')
    !
    lvis=ro(i,0)*utaw/miucal(t(i,0))*Reynolds
    open(18,file='Results/mesh_report.'//irefname//'.dat')
    write(18,"(A,F12.7,A,I4)")'mesh resolution at x=',x(i,0),', i=',i
    write(18,"(A)")'----------------------------------------------------'
    write(18,"(A)")'          Δx+         Δy1+         Δye+          Δz+'
    write(18,"(A)")'----------------------------------------------------'
    write(18,"(4(1X,F12.7))")0.5d0*(x(i+1,0)-x(i-1,0))*lvis,           &
                                       (y(i,1)-y(i,0))*lvis,           &
                               (y(i,jedge)-y(i,jedge-1))*lvis,         &
                                           (z(1)-z(0))*lvis
    write(18,"(A)")'----------------------------------------------------'
    close(18)
    print*,' << mesh_report.',irefname,'.dat ... done !'
    !
    tke=0.5d0*(u11+u22+u33)
    uss=0.d0
    unn=0.d0
    uumax=0.d0
    vvmax=0.d0
    wwmax=0.d0
    kmax=0.d0
    do j=1,jm
    do i=0,im
      var1=u1(i,j)/sqrt(u1(i,j)**2+u2(i,j)**2)
      var2=u2(i,j)/sqrt(u1(i,j)**2+u2(i,j)**2)
      !
      uss(i,j)=u11(i,j)*var1**2+u22(i,j)*var2**2+2.d0*u12(i,j)*var1*var2
      unn(i,j)=u11(i,j)*var2**2+u22(i,j)*var1**2-2.d0*u12(i,j)*var1*var2
      !
      if(u11(i,j)>uumax) then
        uumax=u11(i,j)
        iu=i
        ju=j
      endif
      if(u22(i,j)>vvmax) then
        vvmax=u22(i,j)
        iv=i
        jv=j
      endif
      if(u33(i,j)>wwmax) then
        wwmax=u33(i,j)
        iw=i
        jw=j
      endif
      if(tke(i,j)>kmax) then
        kmax=tke(i,j)
        ik=i
        jk=j
      endif
      !
    enddo
    enddo
    !
    i=iref
    uumax_ref=0.d0
    vvmax_ref=0.d0
    wwmax_ref=0.d0
    kmax_ref=0.d0
    do j=1,jm
      uumax_ref=max(uumax_ref,u11(i,j))
      vvmax_ref=max(vvmax_ref,u22(i,j))
      wwmax_ref=max(wwmax_ref,u33(i,j))
      kmax_ref=max(kmax_ref,tke(i,j))
    enddo
    !
    write(*,"(A)")'-----------  Restress amplification  ---------------'
    write(*,"(3(A,E15.7E3))")' uu amplified:',uumax/uumax_ref,        &
                     ' max uu at:(',(x(iu,ju)-xc)/hstep,',',y(iu,ju)/hstep
    write(*,"(3(A,E15.7E3))")' vv amplified:',vvmax/vvmax_ref,        &
                     ' max vv at:(',(x(iv,jv)-xc)/hstep,',',y(iv,jv)/hstep
    write(*,"(3(A,E15.7E3))")' ww amplified:',wwmax/wwmax_ref,        &
                     ' max ww at:(',(x(iw,jw)-xc)/hstep,',',y(iw,jw)/hstep
    write(*,"(3(A,E15.7E3))")' TKE amplified:',kmax/kmax_ref,     &
                     ' max tke at:(',(x(ik,jk)-xc)/hstep,',',y(ik,jk)/hstep
    write(*,"(A)")'----------------------------------------------------'
    !

    !
    call writetecbin('Results/tec2order.plt',(x-xc)/hstep,'x',         &
                                    y/hstep,'y',u1,'u',u2,'v',         &
                                     u11/utaw/utaw,'uu',               &
                                     u22/utaw/utaw,'vv',               &
                                     u33/utaw/utaw,'ww',               &
                                     tke/utaw/utaw,'TKE',              &
                                     u12/utaw/utaw,'uv',               &
                                     uss/utaw/utaw,'usus',             &
                                     unn/utaw/utaw,'unun',             &
                                     2.d0*sqrt(pp),'prms',im,jm)
    ! call writetecbin('Results/tecredis.plt',(x-xc)/hstep,'x',          &
    !                                 y/hstep,'y',u1,'u',u2,'v',         &
    !                                  u11/utaw/utaw,'uu',             &
    !                                  u22/utaw/utaw,'vv',             &
    !                                  u33/utaw/utaw,'ww',             &
    !                                  u11_hom/utaw/utaw,'uu_>',       &
    !                                  u22_hom/utaw/utaw,'vv_>',       &
    !                                  u33_hom/utaw/utaw,'ww_>',       &
    !                                  u11_dis/utaw/utaw,'uu_<',       &
    !                                  u22_dis/utaw/utaw,'vv_<',       &
    !                                  u33_dis/utaw/utaw,'ww_<',im,jm)
    !
    do i=0,im
      !
      tkemax(i)=0.d0
      do j=0,200
        if(tke(i,j)>=tkemax(i)) then
          tkemax(i)=tke(i,j)
          x_tkemax(i)=(x(i,j)-xc)/hstep
          y_tkemax(i)=y(i,j)/hstep
        endif
        !
      enddo
      !
      !
    enddo
    !
    open(18,file='Results/tkemax.dat')
    write(18,"(4(1X,A15))")'x','y','tke','amp'
    write(18,"(4(1X,E15.7E3))")(x_tkemax(i),y_tkemax(i),tkemax(i),     &
                                          tkemax(i)/tkemax(iref),i=0,im)
    close(18)
    print*,' << Results/tkemax.dat'
    !
    i=iref
    open(18,file='Results/restress.'//irefname//'.dat')
    write(18,"(7(1X,A15))")'y/delta','yplus','uu','vv','ww','uv','k'
    do j=0,jm
      var1=ro(i,j)/ro(i,0)/utaw/utaw
      write(18,"(7(1X,E15.7E3))")y(i,j)/bl_delta99,yplus(j),u11(i,j)*var1,u22(i,j)*var1,u33(i,j)*var1,u12(i,j)*var1,tke(i,j)*var1
    enddo
    close(18)
    print*,' << Results/restress.',irefname,'.delta ... done. '
    !
    allocate( wallnorm(2,0:im),taw(0:im),pwall(0:im),cf(0:im),         &
              ppw(0:im),dpw(0:im),gradp(1:2,0:im,0:jm) )
    call gridjacobian_xy(x,y,wallnormal=wallnorm)
    !
    pwall(:)=p(:,0)
    !
    gradp=grad_xy(p,x,y)
    do i=0,im
      dpw(i)=gradp(1,i,0)*wallnorm(1,i)+gradp(2,i,0)*wallnorm(2,i)
    enddo
    dpw=dpw*2.d0*bl_delta99
    !
    ppw(:)  =sqrt(pp(:,0))
    !
    call PostShockCal(1.d0,1.d0,0.d0,pinf,1.d0,                     &
                                rops1,ups1,vps1,pps1,Tps1,32.d0,Mach)
    !
    write(*,'(2X,A)')'--------------------  post shock --------------------'
    write(*,'(2X,A)')'       ro        u        v        p        t'
    write(*,'(2X,5(F9.3))')rops1,ups1,vps1,pps1,Tps1
    write(*,'(2X,A)')'-----------------------------------------------------'
    !
    pwref=pwall(iref)
    pwall=pwall/pwref
    ppw=ppw/(0.5*roinf*uinf*uinf)
    !
    pwref=ppw(iref)
    ppwmax=0.d0
    do i=0,im-60
      ppwmax=max(ppwmax,ppw(i))
    enddo
    !
    print*,' ** max r.m.s. pressure fluctuation: ',ppwmax
    print*,' **              pp amplified ratio: ',ppwmax/pwref
    print*,' **              pp amplified ratio: ',20.d0*log((ppwmax/pwref)),'dB'
    !
    do i=0,im
      !
      miu=miucal(t(i,0))/Reynolds
      !
      div=0.333333333333333d0*(dudx(i,0)+dvdy(i,0))
      s11=dudx(i,0)-div
      s12=0.5d0*(dudy(i,0)+dvdx(i,0))
      s22=dvdy(i,0)-div
      !
      f11=2.d0*miu*s11
      f12=2.d0*miu*s12
      f22=2.d0*miu*s22
      !
      tawx=f11*wallnorm(1,i)+f12*wallnorm(2,i)
      tawy=f12*wallnorm(1,i)+f22*wallnorm(2,i)
      !
      ! taw(i)=sqrt(tawx**2+tawy**2)
      taw(i)=tawx*wallnorm(2,i)-tawy*wallnorm(1,i)
      !
      cf(i)=taw(i)*2.d0
      !
      delta99(i)=blthickness(y(i,:),u1(i,:),jedge=j,btype='nominal')
      !
    enddo
    !
    cfref=cf(iref)
    !
    cf=cf/cfref
    !
    do i=1,im-1
      !
      if(cf(i)>=0.d0 .and. cf(i+1)<=0.d0) then
        xsep=linear1d(cf(i),cf(i+1),xs(i),xs(i+1),0.d0)
        exit
      endif
    enddo
    do i=im-1,1,-1
      if(cf(i)<=0.d0 .and. cf(i+1)>=0.d0) then
        xatt=linear1d(cf(i),cf(i+1),xs(i),xs(i+1),0.d0)
        exit
      endif
    enddo
    lsep=xatt-xsep
    print*,' ** separation xs: ',xsep,'reattchment xs: ',xatt
    print*,' ** separation xs: ',xsep/hstep,'reattchment xs: ',xatt/hstep,'/h'
    print*,' ** separation xs: ',xsep/bl_delta99,'reattchment xs: ',xatt/bl_delta99,'/δ'
    print*,' ** separation size: ',lsep,lsep/hstep,'/h',lsep/bl_delta99,'/δ'
    !
    open(18,file='Results/cf-pw.dat')
    write(18,"(5(1X,A15))")'xs/h','pwall','cf','pp','dp'
    write(18,"(5(1X,E15.7E3))")(xs(i)/hstep,pwall(i),cf(i),ppw(i),dpw(i),i=0,im-60)
    close(18)
    print*,' << Results/cf-pw.dat'
    !
    allocate(xexp(0:4),xexp2(0:4))
    !
    ! input experimental data
    allocate(yref(0:18),uref(0:18))
    !
    open(12,file='Results/exp_profile.dat')
    read(12,*)
    read(12,*)(xexp(i1),i1=0,4)
    read(12,*)
    read(12,*)xxsep,xxatt
    read(12,*)
    do j=0,18
      read(12,*)yref(j),uref(j)
    enddo
    close(12)
    print*,' >> Results/exp_profile.dat'
    ! !
    lxsep=xxatt-xxsep
    print*,' ** Lsep in experiment: ',lxsep
    delta_exp=blthickness(yref,uref,btype='nominal')
    print*,' ** δ of the first station in experiment: ',delta_exp,'(mm)'
    print*,' ** h/δ: ',15.d0/delta_exp
    print*,' ** δ of the first station in experiment: ',delta_exp/15.d0,'/h'
    !
    nis=10
    allocate(is(0:nis))
    is(0)=996
    is(1)=1124
    is(2)=1377
    is(3)=1528
    is(4)=1635
    is(5)=1800
    is(6)=1895
    ! is(7)=1988
    is(7)=2020
    is(8)=2089
    is(9)=2182
    is(10)=2608
    ! scale with height of step
    ! print*,' +---------------+'
    ! print*,' | scaled with h |'
    ! print*,' +---------------+'
    ! !
    ! do i=1,im
    !   !
    !   if(xs(i)/hstep>=xexp(0) .and. xs(i-1)/hstep<xexp(0)) then
    !     is(0)=i
    !   endif
    !   if(xs(i)/hstep>=xexp(1) .and. xs(i-1)/hstep<xexp(1)) then
    !     is(1)=i
    !   endif
    !   if(xs(i)/hstep>=xexp(2) .and. xs(i-1)/hstep<xexp(2)) then
    !     is(2)=i
    !   endif
    !   if(xs(i)/hstep>=xexp(3) .and. xs(i-1)/hstep<xexp(3)) then
    !     is(3)=i
    !   endif
    !   if(xs(i)/hstep>=xexp(4) .and. xs(i-1)/hstep<xexp(4)) then
    !     is(4)=i
    !   endif
    !   !
    ! enddo
    ! ! !
    ! ! is(4)=min(2950,is(4))
    ! ! !
    do i1=0,nis
      write(*,"(A,I0,A,I0,(A,F10.5))")'  ** sample(',i1,               &
                       '): i= ',is(i1),'; x/h=',xs(is(i1))/hstep
    enddo
    !
    js=240
    allocate(usamp(0:nis,0:js),vsamp(0:nis,0:js),rosamp(0:nis,0:js),   &
             xsamp(0:nis,0:js),ysamp(0:nis,0:js),dy(0:js),             &
             us(0:nis,0:js),un(0:nis,0:js),masamp(0:nis,0:js))
    allocate(u11_samp(0:nis,0:js),u22_samp(0:nis,0:js),                &
             u33_samp(0:nis,0:js),u12_samp(0:nis,0:js),                &
             uss_samp(0:nis,0:js),unn_samp(0:nis,0:js),                &
             usn_samp(0:nis,0:js))
    allocate(pu1_samp(0:nis,0:js),pu2_samp(0:nis,0:js),pu3_samp(0:nis,0:js))
    !
    xsamp(0,:)=x(is(0),0)
    ysamp(0,:)=y(is(0),0:js)
    dy=ysamp(0,:)
    !
    do i1=1,nis
      xsamp(i1,0)=x(is(i1),0)
      ysamp(i1,0)=y(is(i1),0)
    enddo
    !
    do i1=1,nis
      !
      ! if( is(i1)>=1800 .and. is(i1)<=2089 ) then
      !   ! on the step
      !   print*,' ** stations on step',i1,xsamp(i1,0),ysamp(i1,0)
      !   do j=1,js
      !     xsamp(i1,j)=x(is(i1),0)-dy(j)*cos(0.25d0*pi)
      !     ysamp(i1,j)=y(is(i1),0)+dy(j)*cos(0.25d0*pi)
      !   enddo
      ! else
        do j=1,js
          xsamp(i1,j)=x(is(i1),0)
          ysamp(i1,j)=y(is(i1),0)+dy(j)
        enddo
      ! endif
      !
    enddo
    ! !
    call InvdisInterp2D(x,y,u1,u2,ma,im,jm,        &
                        xsamp,ysamp,usamp,vsamp,masamp,nis,js)
    call InvdisInterp2D(x,y,u11,u22,u33,u12,im,jm, &
                        xsamp,ysamp,u11_samp,u22_samp,u33_samp,u12_samp,nis,js)
    call InvdisInterp2D(x,y,pu1,pu2,pu3,im,jm,     &
                        xsamp,ysamp,pu1_samp,pu2_samp,pu3_samp,nis,js)
    !
    do i1=0,nis
      !
      if( is(i1)>1800 .and. is(i1)<=2089 ) then
        ! on the step
        us(i1,:)=(usamp(i1,:)+vsamp(i1,:))*cos(0.25d0*pi)
        un(i1,:)=(vsamp(i1,:)-usamp(i1,:))*cos(0.25d0*pi)
        uss_samp(i1,:)=u11_samp(i1,:)*cos(0.25d0*pi)**2+  &
                       u22_samp(i1,:)*sin(0.25d0*pi)**2+  &
                       2.d0*u12_samp(i1,:)*cos(0.25d0*pi)*sin(0.25d0*pi)
        unn_samp(i1,:)=u11_samp(i1,:)*sin(0.25d0*pi)**2+  &
                       u22_samp(i1,:)*cos(0.25d0*pi)**2-  &
                       2.d0*u12_samp(i1,:)*cos(0.25d0*pi)*sin(0.25d0*pi)
        usn_samp(i1,:)=-u11_samp(i1,:)*cos(0.25d0*pi)*sin(0.25d0*pi) + &
                        u22_samp(i1,:)*sin(0.25d0*pi)*cos(0.25d0*pi)
      else
        us(i1,:)=usamp(i1,:)
        un(i1,:)=vsamp(i1,:)
        !
        uss_samp(i1,:)=u11_samp(i1,:)
        unn_samp(i1,:)=u22_samp(i1,:)
        usn_samp(i1,:)=u12_samp(i1,:)
      endif
      !
    enddo
    !
    u11_samp=u11_samp/utaw/utaw
    u22_samp=u22_samp/utaw/utaw
    u33_samp=u33_samp/utaw/utaw
    u12_samp=u12_samp/utaw/utaw
    uss_samp=uss_samp/utaw/utaw
    unn_samp=unn_samp/utaw/utaw
    usn_samp=usn_samp/utaw/utaw
    !
    call writetecbin('Results/tecsamp.plt',(xsamp-xc)/hstep,'x',       &
                                       ysamp/hstep,'y',usamp,'u',nis,js)
    !
    open(18,file='Results/uprofile.hstep')
    do j=0,js
      write(18,"(12(1X,E15.7E3))")dy(j)/hstep,(usamp(i1,j),i1=0,nis)
    enddo
    close(18)
    print*,' << Results/uprofile.hstep'
    !
    open(18,file='Results/vprofile.hstep')
    do j=0,js
      write(18,"(12(1X,E15.7E3))")dy(j)/hstep,(vsamp(i1,j),i1=0,nis)
    enddo
    close(18)
    print*,' << Results/vprofile.hstep'
    !
    open(18,file='Results/usprofile.hstep')
    do j=0,js
      write(18,"(12(1X,E15.7E3))")dy(j)/hstep,(us(i1,j),i1=0,nis)
    enddo
    close(18)
    print*,' << Results/usprofile.hstep'
    !
    open(18,file='Results/unprofile.hstep')
    do j=0,js
      write(18,"(12(1X,E15.7E3))")dy(j)/hstep,(un(i1,j),i1=0,nis)
    enddo
    close(18)
    print*,' << Results/unprofile.hstep'
    !
    open(18,file='Results/maprofile.hstep')
    do j=0,js
      write(18,"(12(1X,E15.7E3))")dy(j)/hstep,(masamp(i1,j),i1=0,nis)
    enddo
    close(18)
    print*,' << Results/maprofile.hstep'
    !
    open(18,file='Results/ususprofile.hstep')
    do j=0,js
      write(18,"(12(1X,E15.7E3))")dy(j)/hstep,(uss_samp(i1,j),i1=0,nis)
    enddo
    close(18)
    print*,' << Results/ususprofile.hstep'
    !
    open(18,file='Results/ununprofile.hstep')
    do j=0,js
      write(18,"(12(1X,E15.7E3))")dy(j)/hstep,(unn_samp(i1,j),i1=0,nis)
    enddo
    close(18)
    print*,' << Results/ununprofile.hstep'
    !
    open(18,file='Results/usunprofile.hstep')
    do j=0,js
      write(18,"(12(1X,E15.7E3))")dy(j)/hstep,(usn_samp(i1,j),i1=0,nis)
    enddo
    close(18)
    print*,' << Results/usunprofile.hstep'
    !
    open(18,file='Results/u11profile.hstep')
    do j=0,js
      write(18,"(12(1X,E15.7E3))")dy(j)/hstep,(u11_samp(i1,j),i1=0,nis)
    enddo
    close(18)
    print*,' << Results/u11profile.hstep'
    !
    open(18,file='Results/u22profile.hstep')
    do j=0,js
      write(18,"(12(1X,E15.7E3))")dy(j)/hstep,(u22_samp(i1,j),i1=0,nis)
    enddo
    close(18)
    print*,' << Results/u22profile.hstep'
    !
    open(18,file='Results/wwprofile.hstep')
    do j=0,js
      write(18,"(12(1X,E15.7E3))")dy(j)/hstep,(u33_samp(i1,j),i1=0,nis)
    enddo
    close(18)
    print*,' << Results/wwprofile.hstep'
    !
    open(18,file='Results/puprofile.hstep')
    do j=0,js
      write(18,"(12(1X,E15.7E3))")dy(j)/hstep,(pu1_samp(i1,j),i1=0,nis)
    enddo
    close(18)
    print*,' << Results/puprofile.hstep'
    !
    open(18,file='Results/pvprofile.hstep')
    do j=0,js
      write(18,"(12(1X,E15.7E3))")dy(j)/hstep,(pu2_samp(i1,j),i1=0,nis)
    enddo
    close(18)
    print*,' << Results/pvprofile.hstep'
    !
    open(18,file='Results/pwprofile.hstep')
    do j=0,js
      write(18,"(12(1X,E15.7E3))")dy(j)/hstep,(pu3_samp(i1,j),i1=0,nis)
    enddo
    close(18)
    print*,' << Results/pwprofile.hstep'
    !
    stop
    !
    ! delta=blthickness(ysamp(0,:),usamp(0,:),btype='nominal')
    ! print*,' ** δ at the first station: ',delta
    ! ! 
    ! do i1=0,nis
    !   write(irefname,'(1hi,I4.4)')is(i1)
    !   open(18,file='Results/uprofile.'//irefname//'.hstep')
    !   write(18,"(6(1X,A15))")'y/delta','u','y/h','u','y/lsep','u'
    !   write(18,"(6(1X,E15.7E3))")((ysamp(i1,j)-ysamp(i1,0))/delta,usamp(i1,j),       &
    !                               (ysamp(i1,j)-ysamp(i1,0))/hstep,usamp(i1,j),       &
    !                               (ysamp(i1,j)-ysamp(i1,0))/lsep, usamp(i1,j),j=0,jm)
    !   close(18)
    !   print*,' << Results/uprofile.',irefname,'.hstep ... done. '
    ! enddo
    !
    ! allocate(u1_tm(0:im,0:jm,0:km),usamp_tm(0:nis,0:jm,0:km),usamp_max(0:nis,0:jm),usamp_min(0:nis,0:jm))
    ! call H5ReadArray(u1_tm,im,jm,km,'u1','Results/mean.fav.h5')
    ! !
    ! do k=0,km
    !   call InvdisInterp2D(x,y,u1_tm(:,:,k),im,jm, xsamp,ysamp,usamp_tm(:,:,k),nis,jm)
    ! enddo
    ! !
    ! do i1=0,nis
    ! do j=0,jm
    !   !
    !   usamp_max(i1,j)=-1.d10
    !   usamp_min(i1,j)= 1.d10
    !   do k=1,km
    !     usamp_max(i1,j)=max(usamp_max(i1,j),usamp_tm(i1,j,k))
    !     usamp_min(i1,j)=min(usamp_min(i1,j),usamp_tm(i1,j,k))
    !   enddo
    !   !
    ! enddo
    ! enddo
    ! !
    ! do i1=0,nis
    !   write(irefname,'(1hi,I4.4)')is(i1)
    !   open(18,file='Results/uprofile_min-max.'//irefname//'.hstep')
    !   write(18,"(4(1X,A15))")'y/h','u','umin','umax'
    !   write(18,"(4(1X,E15.7E3))")((ysamp(i1,j)-ysamp(i1,0))/hstep,usamp(i1,j),usamp_min(i1,j),usamp_max(i1,j),j=0,jm)
    !   close(18)
    !   print*,' << Results/uprofile_min-max.',irefname,'.hstep ... done. '
    ! enddo
    !
    ! dgrid(:)=y(0,:) !/y(0,jm)*2.d0
    ! !
    ! do i1=0,nis
    !   !
    !   print*,i1,x(is(i1),0),xc,xe
    !   !
    !   ymin=y(is(i1),0)
    !   ymax=y(is(i1),jm)
    !   !
    !   if(x(is(i1),0)>=xc .and. x(is(i1),0)<=xe) then
    !     ! on the step
    !     !
    !     xsamp(i1,:)=x(is(i1),0)-dgrid(:)*sin(alfa)
    !     ysamp(i1,:)=y(is(i1),0)+dgrid(:)*cos(alfa)
    !     !
    !   else
    !     ! not on the step
    !     xsamp(i1,:)=x(is(i1),0)
    !     ysamp(i1,:)=y(is(i1),0)+dgrid(:)
    !     !
    !   endif
    !   !
    !   !
    ! enddo
    ! !
    ! call InvdisInterp2D(x,y,u1,u2,ro,im,jm,xsamp,ysamp,usamp,vsamp,rosamp,nis,jm)
    ! !
    ! call writetecbin('Results/tecsamp2.plt',(xsamp-xc)/hstep,'x',      &
    !                                    ysamp/hstep,'y',usamp,'u',      &
    !                                    vsamp,'v',rosamp,'ro',nis,jm)
    ! !
    ! do i1=0,nis
    !   !
    !   write(irefname,'(1hi,I4.4)')is(i1)
    !   !
    !   if(x(is(i1),0)>=xc .and. x(is(i1),0)<=xe) then
    !     ! on the step
    !     !
    !     usamp(i1,:)=usamp(i1,:)*cos(alfa)+vsamp(i1,:)*sin(alfa)
    !     !
    !   else
    !     ! not on the step
    !     !
    !   endif
    !   !
    !   call upluscal(uplus=uplus,yplus=yplus,u=usamp(i1,:),y=dgrid,     &
    !                   ro=rosamp(i1,:),tw=t(is(i1),0),taw=taw(is(i1)),utaw=utaw)
    !   !
    !   open(18,file='Results/uplus.'//irefname//'.dat')
    !   write(18,"(2(1X,A15))")'yplus','uplus'
    !   write(18,"(2(1X,E15.7E3))")(yplus(j),uplus(j),j=0,jm)
    !   close(18)
    !   print*,' << uplus.',irefname,'.dat ... done. '
    !   !
    ! enddo
    !
    ! ! scale with bl thickness
    ! print*,' +---------------+'
    ! print*,' | scaled with δ |'
    ! print*,' +---------------+'
    ! n=0
    ! do while(n<1)
    !   !
    !   n=n+1
    !   !
    !   xexp2(0)=xexp(0)/delta_exp
    !   xexp2(1)=xexp(1)/delta_exp
    !   xexp2(2)=xexp(2)/delta_exp
    !   xexp2(3)=xexp(3)/delta_exp
    !   xexp2(4)=xexp(4)/delta_exp
    !   !
    !   delta=bl_delta99
    !   !
    !   do i=1,im
    !     !
    !     if(xs(i)/delta>=xexp2(0) .and. xs(i-1)/delta<xexp2(0)) then
    !       is(0)=i
    !     endif
    !     if(xs(i)/delta>=xexp2(1) .and. xs(i-1)/delta<xexp2(1)) then
    !       is(1)=i
    !     endif
    !     if(xs(i)/delta>=xexp2(2) .and. xs(i-1)/delta<xexp2(2)) then
    !       is(2)=i
    !     endif
    !     if(xs(i)/delta>=xexp2(3) .and. xs(i-1)/delta<xexp2(3)) then
    !       is(3)=i
    !     endif
    !     if(xs(i)/delta>=xexp2(4) .and. xs(i-1)/delta<xexp2(4)) then
    !       is(4)=i
    !     endif
    !     !
    !   enddo
    !   !
    !   is(0)=min(590,is(0))
    !   is(4)=min(2950,is(4))
    !   !
    !   do i1=0,nis
    !     write(*,"(A,I1,A,I4,2(A,F10.5))")'  ** sample(',i1,              &
    !                      '): i= ',is(i1),'; x/δ',xs(is(i1))/delta,       &
    !                      '; x/δ (exp)=',xexp2(i1)
    !   enddo
    !   !
    !   xsamp(0,:)=x(is(0),0)
    !   ysamp(0,:)=y(is(0),:)
    !   !
    !   do i1=1,nis
    !     !
    !     ymin=y(is(i1),0)
    !     ymax=y(is(i1),jm)
    !     xsamp(i1,:)=x(is(i1),0)
    !     ysamp(i1,:)=y(is(0),:)/y(is(0),jm)*(ymax-ymin)+ymin
    !     !
    !   enddo
    !   !
    !   call InvdisInterp2D(x,y,u1,im,jm, xsamp,ysamp,usamp,nis,jm)
    !   !
    !   ! call writetecbin('Results/tecsamp.plt',(xsamp-xc)/hstep,'x',       &
    !   !                                    ysamp/hstep,'y',usamp,'u',nis,jm)
    !   !
    !   delta=blthickness(ysamp(0,:),usamp(0,:),btype='nominal')
    !   delta1=delta
    !   print*,' ** δ at the first station: ',delta
    !   ! 
    !   do i1=0,nis
    !     write(irefname,'(1hi,I4.4)')is(i1)
    !     open(18,file='Results/uprofile.'//irefname//'.delta')
    !     write(18,"(6(1X,A15))")'y/delta','u','y/h','u','y/lsep','u'
    !     write(18,"(6(1X,E15.7E3))")((ysamp(i1,j)-ysamp(i1,0))/delta,usamp(i1,j),       &
    !                                 (ysamp(i1,j)-ysamp(i1,0))/hstep,usamp(i1,j),       &
    !                                 (ysamp(i1,j)-ysamp(i1,0))/lsep, usamp(i1,j),j=0,jm)
    !     close(18)
    !     print*,' << Results/uprofile.',irefname,'.delta ... done. '
    !   enddo
    !   !
    ! enddo
    ! !
    ! ! scale with size of separation zone
    ! print*,' +------------------+'
    ! print*,' | scaled with Lspe |'
    ! print*,' +------------------+'
    ! !
    ! xexp2(0)=xexp(0)/lxsep
    ! xexp2(1)=xexp(1)/lxsep
    ! xexp2(2)=xexp(2)/lxsep
    ! xexp2(3)=xexp(3)/lxsep
    ! xexp2(4)=xexp(4)/lxsep
    ! !
    ! do i=1,im
    !   !
    !   if(xs(i)/lsep>=xexp2(0) .and. xs(i-1)/lsep<xexp2(0)) then
    !     is(0)=i
    !   endif
    !   if(xs(i)/lsep>=xexp2(1) .and. xs(i-1)/lsep<xexp2(1)) then
    !     is(1)=i
    !   endif
    !   if(xs(i)/lsep>=xexp2(2) .and. xs(i-1)/lsep<xexp2(2)) then
    !     is(2)=i
    !   endif
    !   if(xs(i)/lsep>=xexp2(3) .and. xs(i-1)/lsep<xexp2(3)) then
    !     is(3)=i
    !   endif
    !   if(xs(i)/lsep>=xexp2(4) .and. xs(i-1)/lsep<xexp2(4)) then
    !     is(4)=i
    !   endif
    !   !
    ! enddo
    ! !
    ! is(4)=min(2950,is(4))
    ! !
    ! print*,' ** Lsep/h(exp) ',lxsep,'Lsep/δ(exp) ',lxsep/delta_exp
    ! print*,' ** Lsep/h(DNS) ',lsep/hstep,'Lsep/δ(DNS) ',lsep/delta

    ! do i1=0,nis
    !   write(*,"(A,I1,A,I4,2(A,F10.5))")'  ** sample(',i1,              &
    !                    '): i= ',is(i1),'; x/Lsep',xs(is(i1))/lsep,     &
    !                    '; x/lsep (exp)=',xexp2(i1)
    ! enddo
    ! !
    ! xsamp(0,:)=x(is(0),0)
    ! ysamp(0,:)=y(is(0),:)
    ! !
    ! do i1=1,nis
    !   !
    !   ymin=y(is(i1),0)
    !   ymax=y(is(i1),jm)
    !   xsamp(i1,:)=x(is(i1),0)
    !   ysamp(i1,:)=y(is(0),:)/y(is(0),jm)*(ymax-ymin)+ymin
    !   !
    ! enddo
    ! !
    ! call InvdisInterp2D(x,y,u1,im,jm, xsamp,ysamp,usamp,nis,jm)
    ! !
    ! delta=blthickness(ysamp(0,:),usamp(0,:),btype='nominal')
    ! print*,' ** δ at the first station: ',delta
    ! ! 
    ! do i1=0,nis
    !   write(irefname,'(1hi,I4.4)')is(i1)
    !   open(18,file='Results/uprofile.'//irefname//'.lsep')
    !   write(18,"(6(1X,A15))")'y/delta','u','y/h','u','y/lsep','u'
    !   write(18,"(6(1X,E15.7E3))")((ysamp(i1,j)-ysamp(i1,0))/delta,usamp(i1,j),       &
    !                               (ysamp(i1,j)-ysamp(i1,0))/hstep,usamp(i1,j),       &
    !                               (ysamp(i1,j)-ysamp(i1,0))/lsep, usamp(i1,j),j=0,jm)
    !   close(18)
    !   print*,' << Results/uprofile.',irefname,'.lsep ... done. '
    ! enddo
    ! !
    ! ! call shockexpcom(beta=alfa/pi*180.d0,pcc=p1,pec=p2)
    ! !
    ! do i=0,im
    !   if(i<ic) then
    !     pinv(i)=1.d0
    !   elseif(i>=ic .and. i<=ie) then
    !     pinv(i)=p1/pinf
    !   elseif(i>ie) then
    !     pinv(i)=p2/pinf
    !   else
    !     stop ' error 1 at comexp'
    !   endif
    ! enddo
    ! !
    !
    ! 3-D wall tecwall-skin-friction line
    allocate(xwall(0:im,0:km),ywall(0:im,0:km),zwall(0:im,0:km))
    !
    do k=0,km
    do i=0,im
      xwall(i,k)=x(i,0)
      ywall(i,k)=y(i,0)
      zwall(i,k)=z(k)
    enddo
    enddo
    xwall=(xwall-xc)/hstep
    ywall=ywall/hstep
    zwall=zwall/hstep
    !
    allocate(dudx_wall(0:im,0:km),dudy_wall(0:im,0:km),                &
             dudz_wall(0:im,0:km),dvdx_wall(0:im,0:km),                &
             dvdy_wall(0:im,0:km),dvdz_wall(0:im,0:km),                &
             dwdx_wall(0:im,0:km),dwdy_wall(0:im,0:km),                &
             dwdz_wall(0:im,0:km),pres_wall(0:im,0:km),                &
             twall(0:im,0:km),taw_wall(1:3,0:im,0:km),                 &
             taw_tan(0:im,0:km),prms_wall(0:im,0:km))
    !
    call H5ReadSubset(dudx_wall,im,jm,km,'du1dx','Results/mean.fav.h5',jslice=0)
    call H5ReadSubset(dudy_wall,im,jm,km,'du1dy','Results/mean.fav.h5',jslice=0)
    call H5ReadSubset(dudz_wall,im,jm,km,'du1dz','Results/mean.fav.h5',jslice=0)
    call H5ReadSubset(dvdx_wall,im,jm,km,'du2dx','Results/mean.fav.h5',jslice=0)
    call H5ReadSubset(dvdy_wall,im,jm,km,'du2dy','Results/mean.fav.h5',jslice=0)
    call H5ReadSubset(dvdz_wall,im,jm,km,'du2dz','Results/mean.fav.h5',jslice=0)
    call H5ReadSubset(dwdx_wall,im,jm,km,'du3dx','Results/mean.fav.h5',jslice=0)
    call H5ReadSubset(dwdy_wall,im,jm,km,'du3dy','Results/mean.fav.h5',jslice=0)
    call H5ReadSubset(dwdz_wall,im,jm,km,'du3dz','Results/mean.fav.h5',jslice=0)
    call H5ReadSubset(pres_wall,im,jm,km,'p','Results/mean.fav.h5',jslice=0)
    call H5ReadSubset(prms_wall,im,jm,km,'pp','Results/2order.fav.h5',jslice=0)
    call H5ReadSubset(twall,im,jm,km,'t','Results/mean.fav.h5',jslice=0)
    !
    do k=0,km
    do i=0,im
      !
      miu=miucal(twall(i,k))/Reynolds
      !
      div=(dudx_wall(i,k)+dvdy_wall(i,k)+dwdz_wall(i,k))/3.d0
      !
      s11=dudx_wall(i,k)-div
      s12=0.5d0*(dudy_wall(i,k)+dvdx_wall(i,k))
      s13=0.5d0*(dudz_wall(i,k)+dwdx_wall(i,k))
      !
      s22=dvdy_wall(i,k)-div
      s23=0.5d0*(dvdz_wall(i,k)+dwdy_wall(i,k))
      !
      s33=dwdz_wall(i,k)-div
      !
      f11=2.d0*miu*s11; f12=2.d0*miu*s12; f13=2.d0*miu*s13
                        f22=2.d0*miu*s22; f23=2.d0*miu*s23
                                          f33=2.d0*miu*s33
      !
      taw_wall(1,i,k)=(f11*wallnorm(1,i)+f12*wallnorm(2,i))*2.d0/cfref
      taw_wall(2,i,k)=(f12*wallnorm(1,i)+f22*wallnorm(2,i))*2.d0/cfref
      taw_wall(3,i,k)=(f13*wallnorm(1,i)+f23*wallnorm(2,i))*2.d0/cfref
      !
      taw_tan(i,k)=taw_wall(1,i,k)*wallnorm(2,i)-taw_wall(2,i,k)*wallnorm(1,i)
      !
      pres_wall(i,k)=pres_wall(i,k)/pinf
      prms_wall(i,k)=sqrt(prms_wall(i,k)/pp(iref,0))
    enddo
    enddo
    !
    call writetecbin('Results/tecwall-skin-friction.plt',              &
                           xwall,'x',ywall,'y',zwall,'z',              &
                     taw_wall(1,:,:),'tawx',taw_wall(2,:,:),'tawy',    &
                     taw_wall(3,:,:),'tawz',taw_tan,'taw',pres_wall,'p',prms_wall,'prms',im,km)
    !
    ! budget terms
    ! allocate(product(0:im,0:jm),prodxx(0:im,0:jm),prodyy(0:im,0:jm),   &
    !          prodxy(0:im,0:jm),prods(0:im,0:jm),prodc(0:im,0:jm),    &
    !          prodd(0:im,0:jm))
    ! allocate(us(0:im,0:jm),un(0:im,0:jm),uss(0:im,0:jm),unn(0:im,0:jm),&
    !          usn(0:im,0:jm),curv(0:im,0:jm),dus(1:2,0:im,0:jm),        &
    !          dusds(0:im,0:jm),dusdn(0:im,0:jm))
    ! !
    ! call H5ReadArray(product,im,jm,'production','Results/tke_budget.zm.h5')
    ! !
    ! miu=miucal(t(iref,0))/Reynolds
    ! budref=ro(iref,0)**2*utaw**4/miu
    ! !
    ! do j=0,jm
    ! do i=0,im
    !   curv(i,j)=-curvature(u1(i,j),u2(i,j),dudx(i,j),dudy(i,j),dvdx(i,j),dvdy(i,j))
    !   us(i,j)=sqrt(u1(i,j)**2+u2(i,j)**2)
    ! enddo
    ! enddo
    ! !
    ! dus=grad(us,x,y)
    ! !
    ! !
    ! do j=0,jm
    ! do i=0,im
    !   !
    !   if(j==0) then
    !     rotmat(1,1)=1.d0
    !     rotmat(1,2)=0.d0
    !     rotmat(2,1)=0.d0
    !     rotmat(2,2)=1.d0
    !   else
    !     var1=sqrt(u1(i,j)**2+u2(i,j)**2)
    !     rotmat(1,1)= u1(i,j)/var1
    !     rotmat(1,2)= u2(i,j)/var1
    !     rotmat(2,1)=-u2(i,j)/var1
    !     rotmat(2,2)= u1(i,j)/var1
    !   endif
    !   uss(i,j)=rotmat(1,1)**2*u11(i,j)+rotmat(1,2)**2*u22(i,j) +       &
    !            2.d0*rotmat(1,1)*rotmat(1,2)*u12(i,j)
    !   unn(i,j)=rotmat(2,1)**2*u11(i,j)+rotmat(2,2)**2*u22(i,j) +       &
    !            2.d0*rotmat(2,1)*rotmat(1,2)*u12(i,j)
    !   usn(i,j)=rotmat(1,1)*rotmat(2,1)*u11(i,j)+                       &
    !            rotmat(1,2)*rotmat(2,2)*u22(i,j)+                       &
    !            (rotmat(1,1)*rotmat(2,2)+rotmat(1,2)*rotmat(2,1))*u12(i,j)
    !   !
    !   dusds(i,j)=rotmat(1,1)*dus(1,i,j)+rotmat(1,2)*dus(2,i,j)
    !   dusdn(i,j)=rotmat(2,1)*dus(1,i,j)+rotmat(2,2)*dus(2,i,j)
    !   !
    !   prodd(i,j)=-ro(i,j)*uss(i,j)*dusds(i,j)
    !   prods(i,j)=-ro(i,j)*usn(i,j)*dusdn(i,j)
    !   prodc(i,j)= ro(i,j)*usn(i,j)*us(i,j)*curv(i,j)
    !   !
    ! enddo
    ! enddo
    ! !
    ! prodxx=-ro*u11*dudx
    ! prodyy=-ro*u22*dvdy
    ! prodxy=-ro*u12*(dudy+dvdx)
    ! !
    ! product=product/budref
    ! prodxx =prodxx /budref
    ! prodyy =prodyy /budref
    ! prodxy =prodxy /budref
    ! !
    ! prodd  =prodd  /budref
    ! prods  =prods  /budref
    ! prodc  =prodc  /budref
    ! !
    ! do j=0,jm
    ! do i=0,im
    !   !
    !   if(j==0) then
    !     rotmat(1,1)=1.d0
    !     rotmat(1,2)=0.d0
    !     rotmat(2,1)=0.d0
    !     rotmat(2,2)=1.d0
    !   else
    !     var1=sqrt(u1(i,j)**2+u2(i,j)**2)
    !     rotmat(1,1)= u1(i,j)/var1
    !     rotmat(1,2)= u2(i,j)/var1
    !     rotmat(2,1)=-u2(i,j)/var1
    !     rotmat(2,2)= u1(i,j)/var1
    !   endif
    !   !
    !   us(i,j)=rotmat(1,1)*u1(i,j)+rotmat(1,2)*u2(i,j)
    !   un(i,j)=rotmat(2,1)*u1(i,j)+rotmat(2,2)*u2(i,j)
    !   !
    !   bmat(1,1)=rotmat(1,1)*dudx(i,j)+rotmat(1,2)*dvdx(i,j)  !dus/dx
    !   bmat(1,2)=rotmat(2,1)*dudx(i,j)+rotmat(2,2)*dvdx(i,j)  !dun/dx
    !   bmat(2,1)=rotmat(1,1)*dudy(i,j)+rotmat(1,2)*dvdy(i,j)  !dus/dy
    !   bmat(2,2)=rotmat(2,1)*dudy(i,j)+rotmat(2,2)*dvdy(i,j)  !dun/dy
    !   !
    !   dus(1,i,j)=rotmat(1,1)*bmat(1,1)+rotmat(1,2)*bmat(2,1) ! dus/ds
    !   dus(2,i,j)=rotmat(2,1)*bmat(1,1)+rotmat(2,2)*bmat(2,1) ! dus/dn
    !   dun(1,i,j)=rotmat(1,1)*bmat(1,2)+rotmat(1,2)*bmat(2,2) ! dun/ds
    !   dun(2,i,j)=rotmat(2,1)*bmat(1,2)+rotmat(2,2)*bmat(2,2) ! dun/dn
    !   !
    !   uss(i,j)=rotmat(1,1)**2*u11(i,j)+rotmat(1,2)**2*u22(i,j) +       &
    !            2.d0*rotmat(1,1)*rotmat(1,2)*u12(i,j)
    !   unn(i,j)=rotmat(2,1)**2*u11(i,j)+rotmat(2,2)**2*u22(i,j) +       &
    !            2.d0*rotmat(2,1)*rotmat(1,2)*u12(i,j)
    !   usn(i,j)=rotmat(1,1)*rotmat(2,1)*u11(i,j)+                       &
    !            rotmat(1,2)*rotmat(2,2)*u22(i,j)+                       &
    !            (rotmat(1,1)*rotmat(2,2)+rotmat(1,2)*rotmat(2,1))*u12(i,j)
    ! enddo
    ! enddo
    ! !
    ! i=iref
    ! miu=miucal(t(i,0))/Reynolds
    ! var1=ro(i,0)**2*utaw**4/miu
    ! !
    ! prodss = -ro*uss*dus(1,:,:) /var1
    ! prodnn = -ro*unn*dun(2,:,:) /var1
    ! prodsn = -ro*usn*(dus(2,:,:)+dun(1,:,:)) /var1
    !
    ! call writetecbin('Results/tecprod.plt',(x-xc)/hstep,'x',           &
    !               y/hstep,'y',u1,'u',u2,'v',product,'P',prodxx,'Pxx',  &
    !               prodyy,'Pyy',prodxy,'Pxy',prods,'Ps',prodd,'Pd',     &
    !               prodc,'Pc',curv,'curvature',im,jm)
    ! do i=0,im
    !   !
    !   prodmax(i)=0.d0
    !   do j=0,140
    !     !
    !     if(product(i,j)>=prodmax(i)) then
    !       prodmax(i)=product(i,j)
    !       prodsmax(i)=prods(i,j)
    !       prodcmax(i)=prodc(i,j)
    !       proddmax(i)=prodd(i,j)
    !       !
    !       xprodmax(i)=(x(i,j)-xc)/hstep
    !       yprodmax(i)=y(i,j)/hstep
    !     endif
    !     !
    !   enddo
    !   !
    ! enddo
    ! !
    ! open(18,file='Results/maxproduction.dat')
    ! write(18,"(6(1X,A15))")'x','y','P','Ps','Pd','Pc'
    ! write(18,"(6(1X,E15.7E3))")(xprodmax(i),yprodmax(i),prodmax(i),    &
    !                          prodsmax(i),proddmax(i),prodcmax(i),i=0,im)
    ! close(18)
    ! print*,' << Results/maxproduction.dat'
    !
    allocate(convect(0:im,0:jm),pre_ace(0:im,0:jm),pre_dil(0:im,0:jm), &
             pre_tra(0:im,0:jm),product(0:im,0:jm),tur_tra(0:im,0:jm), &
             vis_ace(0:im,0:jm),vis_dif(0:im,0:jm),dissipa(0:im,0:jm), &
             balance(0:im,0:jm) )
    ! !
    call H5ReadArray(convect,im,jm,'convection','Results/tke_budget.zm.h5')
    call H5ReadArray(dissipa,im,jm,'dissipation','Results/tke_budget.zm.h5')
    call H5ReadArray(pre_ace,im,jm,'pressure_accelaration','Results/tke_budget.zm.h5')
    call H5ReadArray(pre_dil,im,jm,'pressure_dilatation','Results/tke_budget.zm.h5')
    call H5ReadArray(pre_tra,im,jm,'pressure_transport','Results/tke_budget.zm.h5')
    call H5ReadArray(product,im,jm,'production','Results/tke_budget.zm.h5')
    call H5ReadArray(tur_tra,im,jm,'turbulent_transport','Results/tke_budget.zm.h5')
    call H5ReadArray(vis_ace,im,jm,'viscous_accelaration','Results/tke_budget.zm.h5')
    call H5ReadArray(vis_dif,im,jm,'viscous_diffusion','Results/tke_budget.zm.h5')
    !
    balance=convect+pre_ace+pre_dil+pre_tra+product+tur_tra+vis_ace+vis_dif+dissipa
    !
    do j=0,jm
    do i=0,im
      !
      if(balance(i,j)>0.d0) then
        dissipa(i,j)=dissipa(i,j)-0.7d0*balance(i,j)
      endif
      !
    enddo
    enddo
    balance=convect+pre_ace+pre_dil+pre_tra+product+tur_tra+vis_ace+vis_dif+dissipa
    !
    i=iref
    !
    write(irefname,'(1hi,I4.4)')i
    open(18,file='Results/mesh_report.'//irefname//'.dat',position='append')
    !
    write(18,"(A)")' ratio of mesh to kolmogorov scale at wall          '
    write(18,"(A)")'----------------------------------------------------'
    write(18,"(A)")'         Δx/η         Δy/η         Δz/η         Δe/η'
    write(18,"(A)")'----------------------------------------------------'
    !
    miu=miucal(t(i,0))/Reynolds
    var1=sqrt(sqrt(miu**3/ro(i,0)**2/(-dissipa(i,0))))
    write(18,"(4(1X,F12.7))")0.5d0*(x(i+1,0)-x(i-1,0))/var1,           &
                                       (y(i,1)-y(i,0))/var1,           &
                                           (z(1)-z(0))/var1,           &
                         cube_root(0.5d0*(x(i+1,0)-x(i-1,0))*          &
                                (y(i,1)-y(i,0))*(z(1)-z(0)))/var1
    write(18,"(A)")'----------------------------------------------------'
    !
    write(18,"(A)")' ratio of mesh to kolmogorov scale at BL edge       '
    write(18,"(A)")'----------------------------------------------------'
    write(18,"(A)")'         Δx/η         Δy/η         Δz/η         Δe/η'
    write(18,"(A)")'----------------------------------------------------'
    !
    miu=miucal(t(i,jedge))/Reynolds
    var1=sqrt(sqrt(miu**3/ro(i,jedge)**2/(-dissipa(i,jedge))))
    write(18,"(4(1X,F12.7))")0.5d0*(x(i+1,jedge)-x(i-1,jedge))/var1,   &
                             0.5d0*(y(i,jedge+1)-y(i,jedge-1))/var1,   &
                                                 (z(1)-z(0))/var1,     &
                       cube_root(0.25d0*(x(i+1,jedge)-x(i-1,jedge))*   &
                           (y(i,jedge+1)-y(i,jedge-1))*(z(1)-z(0)))/var1
    write(18,"(A)")'----------------------------------------------------'
    close(18)
    print*,' << mesh_report.',irefname,'.dat ... done !'
    ! 
    ! is(1)=1
    ! is(2)=2
    ! is(3)=3
    ! is(4)=4
    !
    do i1=0,4
      !
      i=is(i1)
      !
      write(irefname,'(1hi,I4.4)')i
      !
      call upluscal(uplus=uplus,yplus=yplus,u=u1(i,:),y=y(i,:),        &
                              ro=ro(i,:),tw=t(i,0),taw=taw(i),utaw=utaw)
      !
      open(18,file='Results/budget_profile.'//irefname//'.dat')
      write(18,"(A,E15.7E3)")'# x=',x(i,0)/bl_delta99
      write(18,"(12(1X,A15))")'y/delta','yplus','convection',          &
              'pressure_accelaration','pressure_dilatation',           &
              'pressure_transport','production','turbulent_transport', &
              'viscous_accelaration','viscous_diffusion','dissipation',&
              'balance'
      do j=0,jm
        miu=miucal(t(i,0))/Reynolds
        var1=ro(i,0)**2*utaw**4/miu
        write(18,"(12(1X,E15.7E3))")y(i,j)/delta99(i),yplus(j),          &
                  convect(i,j)/var1,pre_ace(i,j)/var1,pre_dil(i,j)/var1, &
                  pre_tra(i,j)/var1,product(i,j)/var1,tur_tra(i,j)/var1, &
                  vis_ace(i,j)/var1,vis_dif(i,j)/var1,dissipa(i,j)/var1, &
                  balance(i,j)/var1
      end do
      close(18)
      print*,' << budget_profile.',irefname,'.dat ... done. '
      !
    enddo
    !
    call writebudgetlay('budget_profile.'//irefname//'.dat',           &
                        'Results/plot_budget.lay')
    !
    ! delta=blthickness(y(i1,:),u1(i1,:),jedge=jedge,btype='nominal')
    ! !
    ! write(irefname,'(1hi,I4.4)')i1
    ! open(18,file='Results/profile.'//irefname//'.dat')
    ! write(18,"(5(1X,A15))")'y/delta','ro','u','T','du/dy'
    ! write(18,"(5(1X,E15.7E3))")(y(i1,j)/hstep,ro(i1,j),                &
    !                                  u1(i1,j),t(i1,j),dudy(i1,j),j=0,jm)
    ! close(18)
    ! print*,' << profile.',irefname,'.dat ... done. '
    ! !
    ! write(irefname,'(1hi,I4.4)')i2
    ! open(18,file='Results/profile.'//irefname//'.dat')
    ! write(18,"(5(1X,A15))")'y/delta','ro','u','T','du/dy'
    ! write(18,"(5(1X,E15.7E3))")(y(i2,j)/hstep,ro(i2,j),                &
    !                                  u1(i2,j),t(i2,j),dudy(i2,j),j=0,jm)
    ! close(18)
    ! print*,' << profile.',irefname,'.dat ... done. '
    ! !
    ! write(irefname,'(1hi,I4.4)')i3
    ! open(18,file='Results/profile.'//irefname//'.dat')
    ! write(18,"(5(1X,A15))")'y/delta','ro','u','T','du/dy'
    ! write(18,"(5(1X,E15.7E3))")(y(i3,j)/hstep,ro(i3,j),                &
    !                                  u1(i3,j),t(i3,j),dudy(i3,j),j=0,jm)
    ! close(18)
    ! print*,' << profile.',irefname,'.dat ... done. '
    ! !
    ! write(irefname,'(1hi,I4.4)')i4
    ! open(18,file='Results/profile.'//irefname//'.dat')
    ! write(18,"(5(1X,A15))")'y/delta','ro','u','T','du/dy'
    ! write(18,"(5(1X,E15.7E3))")(y(i4,j)/hstep,ro(i4,j),                &
    !                                  u1(i4,j),t(i4,j),dudy(i4,j),j=0,jm)
    ! close(18)
    ! print*,' << profile.',irefname,'.dat ... done. '
    ! !
    ! write(irefname,'(1hi,I4.4)')i5
    ! open(18,file='Results/profile.'//irefname//'.dat')
    ! write(18,"(5(1X,A15))")'y/delta','ro','u','T','du/dy'
    ! write(18,"(5(1X,E15.7E3))")(y(i5,j)/hstep,ro(i5,j),                &
    !                                  u1(i5,j),t(i5,j),dudy(i5,j),j=0,jm)
    ! close(18)
    ! print*,' << profile.',irefname,'.dat ... done. '
    ! !
    !
    ! write(18,*) 'variables=x,y,u,v'
    ! write(18,*) 'zone f=point,i=',im+1,'j=',1
    ! write(18,"(4(1X,E15.7E3))")(x(i,0),y(i,0),wallnorm(1,i),           &
    !                                           wallnorm(2,i),i=0,im)
    ! close(18)
    ! print*,' << tecwall.dat'
    ! !
    ! do i=0,im
    !   print*,wallnorm(1:2,i)
    ! enddo
    !
  end subroutine comexp
  !+-------------------------------------------------------------------+
  !| The end of the subroutine comexp.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| The subroutine is used to analyse flow past a cylinder.           |
  !+-------------------------------------------------------------------+
  !| Writen by Jian Fang, 2020-06-19.                                  |
  !+-------------------------------------------------------------------+
  subroutine cylinder
    !
    use commvardefine,only: im,jm,km,Reynolds,pinf,pi
    use h5readwrite
    use basicfunction
    use interpolation 
    use writetec
    use boundarylayer
    use gradsolver
    !
    real(8),allocatable,dimension(:,:) :: x,y,ro,u1,u2,p,t,            &
                                          dudx,dudy,dvdx,dvdy,omegaz
    real(8),allocatable,dimension(:) :: z,taw,pwall,cf,utau,lplus,alfa,&
                                        ds,dr,dz
    real(8),allocatable,dimension(:,:) :: wallnorm
    real(8) :: div,s11,s12,s13,s22,s23,s33,f11,f12,f13,f22,f23,f33,     &
               var1,miu,tawx,tawy,delta,pwref
    integer :: i,j,k
    !
    !+-------+
    !| begin |
    !+-------+
    !
    print*,' ** processing data of a flow past a cylinder'
    !
    allocate(x(0:im,0:jm),y(0:im,0:jm),z(0:km))
    call H5ReadSubset(x,im,jm,km,'x','datin/grid.h5',kslice=0)
    call H5ReadSubset(y,im,jm,km,'y','datin/grid.h5',kslice=0)
    call H5ReadSubset(z,im,jm,km,'z','datin/grid.h5',islice=0,jslice=0)
    !
    allocate(ro(0:im,0:jm),u1(0:im,0:jm),u2(0:im,0:jm),p(0:im,0:jm),   &
             t(0:im,0:jm),dudx(0:im,0:jm),dudy(0:im,0:jm),             &
             dvdx(0:im,0:jm),dvdy(0:im,0:jm),omegaz(0:im,0:jm))
    !
    call H5ReadArray(ro,im,jm,'ro','Results/mean.fav.zm.h5')
    call H5ReadArray(u1,im,jm,'u1','Results/mean.fav.zm.h5')
    call H5ReadArray(u2,im,jm,'u2','Results/mean.fav.zm.h5')
    call H5ReadArray(p, im,jm, 'p','Results/mean.fav.zm.h5')
    call H5ReadArray(t, im,jm, 't','Results/mean.fav.zm.h5')
    call H5ReadArray(dudx,im,jm,'du1dx','Results/mean.fav.zm.h5')
    call H5ReadArray(dudy,im,jm,'du1dy','Results/mean.fav.zm.h5')
    call H5ReadArray(dvdx,im,jm,'du2dx','Results/mean.fav.zm.h5')
    call H5ReadArray(dvdy,im,jm,'du2dy','Results/mean.fav.zm.h5')
    !
    print*,' ** calcualte wall variables'
    !
    allocate(wallnorm(2,0:im),taw(0:im),pwall(0:im),cf(0:im))
    allocate(utau(0:im),lplus(0:im),alfa(0:im))
    !
    call gridjacobian_xy(x,y,wallnormal=wallnorm)
    !
    pwall(:)=p(:,0)
    !
    pwref=pwall(0)
    pwall=pwall/pwref
    !
    do i=0,im
      !
      alfa(i)=asin(y(i,0)/sqrt(x(i,0)**2+y(i,0)**2))
      !
      miu=miucal(t(i,0))/Reynolds
      !
      div=0.333333333333333d0*(dudx(i,0)+dvdy(i,0))
      s11=dudx(i,0)-div
      s12=0.5d0*(dudy(i,0)+dvdx(i,0))
      s22=dvdy(i,0)-div
      !
      f11=2.d0*miu*s11
      f12=2.d0*miu*s12
      f22=2.d0*miu*s22
      !
      tawx=f11*wallnorm(1,i)+f12*wallnorm(2,i)
      tawy=f12*wallnorm(1,i)+f22*wallnorm(2,i)
      !
      ! taw(i)=sqrt(tawx**2+tawy**2)
      taw(i)=tawx*wallnorm(2,i)-tawy*wallnorm(1,i)
      !
      utau(i)=sqrt(taw(i)/ro(i,0))
      lplus(i)=ro(i,0)*utau(i)/miu
      !
      cf(i)=taw(i)*2.d0
      !
    enddo
    !
    print*,' ** calcualte mesh resolution'
    allocate(ds(1:im-1),dr(1:im-1),dz(1:im-1))
    !
    do i=1,im-1
      !
      ds(i)=0.5d0*sqrt((x(i+1,0)-x(i-1,0))**2+(y(i+1,0)-y(i-1,0))**2)
      dr(i)=sqrt((x(i,1)-x(i,0))**2+(y(i,1)-y(i,0))**2)
      dz(i)=z(1)-z(0)
      !
      ds(i)=ds(i)*lplus(i)
      dr(i)=dr(i)*lplus(i)
      dz(i)=dz(i)*lplus(i)
      !
    enddo
    !
    open(18,file='Results/mesh_resolution.dat')
    write(18,"(6(1X,A15))")'alfa','ds+','dr+','dz+','l+','rhow'
    write(18,"(6(1X,E15.7E3))")(alfa(i),ds(i),dr(i),dz(i),1.d0/lplus(i),ro(i,0),i=1,im-1)
    close(18)
    print*,' << Results/mesh_resolution.dat'
    !
  end subroutine cylinder
  !+-------------------------------------------------------------------+
  !| The end of the subroutine comexp.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to analyse islice.                        |
  !+-------------------------------------------------------------------+
  subroutine isliceanalys(nsta,nend,mode)
    !
    use commvardefine,only: im,jm,km,Reynolds,Mach,pinf,pi,const2,gridfile
    use h5readwrite
    use writetec
    use interpolation
    use basicfunction
    !
    integer,intent(in) :: nsta,nend
    character(len=*),intent(in) :: mode
    !
    real(8),allocatable,dimension(:) :: rom,u1m,u2m,u3m,pm,tm,uu,vv,ww,&
                                        uv,uw,vw,tke,yplus,uplus,mat,marms
    real(8),allocatable,dimension(:,:) :: rotm,u1tm,u2tm,u3tm,ptm,ttm
    real(8),allocatable,dimension(:,:) :: ro,u1,u2,u3,p,t,x,y,z,       &
                                          zcoru,zcorv,zcorw,           &
                                          tcoru,tcorv,tcorw
    real(8),allocatable,dimension(:,:,:) :: ro_i,u1_i,u2_i,u3_i,p_i,t_i
    real(8),allocatable :: rsamp(:)
    integer,allocatable :: nsamp(:)
    character(len=5) :: fname
    character(len=3) :: jname
    character(len=255) :: filename
    integer :: n,i,j,k,k1,k2,js(4),nlen,ncoun,ntc,nic,n1,n2,kslec
    real(8) :: time,time1,time2,deltat,var1,var2,bl_delta99,bl_dstar,  &
               bl_theta,utaw,lvis,tvis,nsamples,shapefac
    logical :: lfilalive
    !
    if(mode=='check') then
      !
      do n=nsta,nend
        !
        write(fname,'(I5.5)')n
        !
        filename='islice/islice'//fname//'.h5'
        !
        inquire(file=trim(filename),exist=lfilalive)
        !
        if(lfilalive) then
          !
          write(*,'(1A1,3(A),I4,A2,$)')char(13),'  ** checking file: ', &
                       trim(filename),' : ',100*(n-nsta+1)/(nend-nsta+1),' %'
          !
        else
          print*,' ** file: ',trim(filename),' is missed'
        endif
        !
      enddo
      !
      write(*,*)
      if(n==nend+1) then
        print*,' ** all the slice files are checked'
      endif
      !
      return
      !
    else
      !
      allocate(x(0:jm,0:km),y(0:jm,0:km),z(0:jm,0:km))
      call h5_read2dfrom3d(x,im,jm,km,'x',trim(gridfile),islice=0)
      call h5_read2dfrom3d(y,im,jm,km,'y',trim(gridfile),islice=0)
      call h5_read2dfrom3d(z,im,jm,km,'z',trim(gridfile),islice=0)
      ! call H5ReadArray(x,jm,km,'x',trim(gridfile))
      ! call H5ReadArray(y,jm,km,'y',trim(gridfile))
      ! call H5ReadArray(z,jm,km,'z',trim(gridfile))
      !
      allocate(rom(0:jm),u1m(0:jm),u2m(0:jm),u3m(0:jm),pm(0:jm),tm(0:jm))
      allocate(rotm(0:jm,0:km),u1tm(0:jm,0:km),u2tm(0:jm,0:km),          &
               u3tm(0:jm,0:km),ptm(0:jm,0:km),ttm(0:jm,0:km))
      !
    endif
    !
    if(mode=='bl') then
      !
      call H5ReadArray(rom,jm,'ro','Results/meanprofile.h5')
      call H5ReadArray(u1m,jm,'u1','Results/meanprofile.h5')
      call H5ReadArray(u2m,jm,'u2','Results/meanprofile.h5')
      call H5ReadArray(u3m,jm,'u3','Results/meanprofile.h5')
      call H5ReadArray( tm,jm, 't','Results/meanprofile.h5')
      call H5ReadArray( pm,jm, 'p','Results/meanprofile.h5')
      !
      do j=1,jm
        if(u1m(j-1)<=0.99d0 .and. u1m(j)>=0.99d0) then
          bl_delta99=linear1d(u1m(j-1),u1m(j),y(j-1,0),y(j,0),0.99d0)
          exit
        endif
      enddo
      !
      print*,' **    δ=',bl_delta99
      print*,' ** Re_δ=',bl_delta99*Reynolds
      !
      bl_dstar=0.d0
      bl_theta=0.d0
      do j=1,jm
        var1=0.5d0*(u1m(j)*rom(j)+u1m(j-1)*rom(j-1))
        bl_dstar=bl_dstar+(1.d0-var1)*(y(j,0)-y(j-1,0))
        !
        var1=0.5d0*(u1m(j)*rom(j)+u1m(j-1)*rom(j-1))
        var2=1.d0-0.5d0*(u1m(j)+u1m(j-1))
        bl_theta=bl_theta+var1*var2*(y(j,0)-y(j-1,0))
      enddo
      !
      shapefac=bl_dstar/bl_theta
      !
      open(18,file='Results/profile.dat')
      write(18,"(6(1X,A15))")'y/δ','ro','u','v','w','T'
      write(18,"(6(1X,E15.7E3))")(y(j,0)/bl_delta99,rom(j),u1m(j),     &
                                            u2m(j),u3m(j),tm(j),j=0,jm)
      close(18)
      print*,' << profile.dat ... done. '
      !
      !
      allocate(yplus(0:jm),uplus(0:jm))
      call upluscal(uplus=uplus,yplus=yplus,u=u1m,y=y(:,0),ro=rom,     &
                                                     tw=tm(0),utaw=utaw)
      !
      ! call H5WriteArray(uplus,jm,'uplus','Results/meanprofile.h5')
      ! call H5WriteArray(yplus,jm,'yplus','Results/meanprofile.h5')
      ! call H5WriteArray(utaw,'utaw','Results/meanprofile.h5')
      !
      write(*,"(A,F7.3,A)")'  ---------------- BL thickness  ----------------'
      write(*,"(A)")'            δ          δ*           θ           H'
      write(*,"(1X,4(F12.7))")bl_delta99,bl_dstar,bl_theta,shapefac
      write(*,"(A)")'         Re_δ       Re_δ*        Re_θ        Re_taw'
      write(*,"(1X,4(E12.6E1))")bl_delta99*Reynolds,bl_dstar*Reynolds,bl_theta*Reynolds, &
                              rom(0)*utaw*bl_delta99/miucal(tm(0))*Reynolds
      write(*,"(A)")'  -----------------------------------------------'
      !
      open(16,file='Results/meanprofile.dat')    
      write(16,"(4(1X,A14))")'delta','delta*','theter','utaw'
      write(16,"(4(1X,E14.7E2))")bl_delta99,bl_dstar,bl_theta,utaw
      write(16,"(6(1X,A14))")'y','ro','u1','u2','t','p'
      write(16,"(6(1X,E14.7E2))")(y(j,0),rom(j),u1m(j),u2m(j),tm(j),pm(j),j=0,jm)
      close(16)
      print*,' <<< meanprofile.dat ... done!'
      !
      open(18,file='Results/uplus.dat')
      write(18,"(2(1X,A15))")'yplus','uplus'
      write(18,"(2(1X,E15.7E3))")(yplus(j),uplus(j),j=0,jm)
      close(18)
      print*,' << uplus.dat ... done. '
      !
      allocate(uu(0:jm),vv(0:jm),ww(0:jm),uv(0:jm),uw(0:jm),vw(0:jm),  &
               tke(0:jm),mat(0:jm),marms(0:jm))
      !
      call H5ReadArray(uu,jm,'uu','Results/restress.h5')
      call H5ReadArray(vv,jm,'vv','Results/restress.h5')
      call H5ReadArray(ww,jm,'ww','Results/restress.h5')
      call H5ReadArray(uv,jm,'uv','Results/restress.h5')
      call H5ReadArray(uw,jm,'uw','Results/restress.h5')
      call H5ReadArray(vw,jm,'vw','Results/restress.h5')
      !
      tke=0.5d0*(uu+vv+ww)
      !
      mat=sqrt(uu+vv+ww)/sqrt(tm)*Mach
      !
      open(18,file='Results/Restress_profile.dat')
      write(18,"(8(1X,A15))")'y/delta','yplus','uu','vv','ww','tke','uv','Mt'
      do j=0,jm
        ! var1=rom(j)/rom(0)/utaw/utaw
        var1=1.d0/utaw/utaw
        write(18,"(8(1X,E15.7E3))")y(j,0)/bl_delta99,yplus(j),         &
                 uu(j)*var1,vv(j)*var1,ww(j)*var1,tke(j)*var1,         &
                 uv(j)*var1,mat(j)
      end do
      close(18)
      print*,' << Restress_profile.dat ... done. '
      !
      allocate(zcoru(0:jm,0:km/2),zcorv(0:jm,0:km/2),zcorw(0:jm,0:km/2))
      call H5ReadArray(zcoru,jm,km/2,'coru','Results/correlation.h5')
      call H5ReadArray(zcorv,jm,km/2,'corv','Results/correlation.h5')
      call H5ReadArray(zcorw,jm,km/2,'corw','Results/correlation.h5')
      !
      js(1)=7
      js(2)=27
      js(3)=62
      js(4)=84
      do j=1,4
        write(jname,'(I3.3)')js(j)
        !
        lvis=rom(0)*utaw/miucal(tm(0))*Reynolds
        !
        open(18,file='Results/correlation'//jname//'.dat')
        write(18,'(2(A,F8.4))')'# yplus=',yplus(js(j)),' y/δ=',        &
                                           y(js(j),0)/bl_delta99
        write(18,'(4(1X,A15))')'zplus','corr_u','corr_v','corr_w'
        write(18,'(4(1X,E15.7E3))')(z(js(j),k)*lvis,zcoru(js(j),k),    &
                                  zcorv(js(j),k),zcorw(js(j),k),k=0,km/2)
        close(18)
        print*,' << Results/correlation',jname,'.dat ... done. '
      enddo
      !
      nlen=200
      allocate(tcoru(0:jm,0:nlen),tcorv(0:jm,0:nlen),tcorw(0:jm,0:nlen))
      call H5ReadArray(tcoru,jm,nlen,'tcoru','Results/tcorrelation.h5')
      call H5ReadArray(tcorv,jm,nlen,'tcorv','Results/tcorrelation.h5')
      call H5ReadArray(tcorw,jm,nlen,'tcorw','Results/tcorrelation.h5')
      !
      write(fname,'(I5.5)')nsta
      filename='islice/islice'//fname//'.h5'
      call H5ReadArray(time1,'time',trim(filename))
      write(fname,'(I5.5)')nsta+1
      filename='islice/islice'//fname//'.h5'
      call H5ReadArray(time2,'time',trim(filename))
      !
      deltat=time2-time1
      !
      do j=1,4
        write(jname,'(I3.3)')js(j)
        !
        lvis=rom(0)*utaw/miucal(tm(0))*Reynolds
        tvis=lvis/utaw
        !
        open(18,file='Results/tcorrelation'//jname//'.dat')
        write(18,'(2(A,F8.4))')'# yplus=',yplus(js(j)),' y/δ=',        &
                                           y(js(j),0)/bl_delta99
        write(18,'(6(1X,A15))')'t','t_delta','tplus','corr_u',         &
                               'corr_v','corr_w'
        write(18,'(6(1X,E15.7E3))')(deltat*n,deltat*n/bl_delta99,      &
                     deltat*n/tvis,tcoru(js(j),n),tcorv(js(j),n),      &
                                   tcorw(js(j),n),n=0,nlen)
        close(18)
        print*,' << Results/tcorrelation',jname,'.dat ... done. '
      enddo
      !
      print*,' ** isliceanalys job done !' 
      !
      stop 
      !
    endif
    !
    if(mode=='mean') then
      !
      rom=0.d0
      u1m=0.d0
      u2m=0.d0
      u3m=0.d0
       pm=0.d0
       tm=0.d0
      !
      rotm=0.d0
      u1tm=0.d0
      u2tm=0.d0
      u3tm=0.d0
       ptm=0.d0
       ttm=0.d0
      !
    endif
    !
    if(mode=='horder') then
      !
      call H5ReadArray(rom,jm,'ro','Results/meanprofile.h5')
      call H5ReadArray(u1m,jm,'u1','Results/meanprofile.h5')
      call H5ReadArray(u2m,jm,'u2','Results/meanprofile.h5')
      call H5ReadArray(u3m,jm,'u3','Results/meanprofile.h5')
      call H5ReadArray( tm,jm, 't','Results/meanprofile.h5')
      call H5ReadArray( pm,jm, 'p','Results/meanprofile.h5')
      !
      allocate(uu(0:jm),vv(0:jm),ww(0:jm),uv(0:jm),uw(0:jm),vw(0:jm))
      !
      uu=0.d0; vv=0.d0; ww=0.d0
      uv=0.d0; uw=0.d0; vw=0.d0
      !
      allocate(zcoru(0:jm,0:km/2),zcorv(0:jm,0:km/2),zcorw(0:jm,0:km/2))
      !
      zcoru=0.d0
      zcorv=0.d0
      zcorw=0.d0
      !
    endif
    !
    if(mode=='tcorr') then
      !
      call H5ReadArray(rom,jm,'ro','Results/meanprofile.h5')
      call H5ReadArray(u1m,jm,'u1','Results/meanprofile.h5')
      call H5ReadArray(u2m,jm,'u2','Results/meanprofile.h5')
      call H5ReadArray(u3m,jm,'u3','Results/meanprofile.h5')
      call H5ReadArray( tm,jm, 't','Results/meanprofile.h5')
      call H5ReadArray( pm,jm, 'p','Results/meanprofile.h5')
      !
      nlen=200
      !
      allocate(ro_i(0:jm,0:km,1:nlen),u1_i(0:jm,0:km,1:nlen),          &
               u2_i(0:jm,0:km,1:nlen),u3_i(0:jm,0:km,1:nlen),          &
                p_i(0:jm,0:km,1:nlen), t_i(0:jm,0:km,1:nlen))
      !
      allocate(tcoru(0:jm,0:nlen),tcorv(0:jm,0:nlen),tcorw(0:jm,0:nlen))
      allocate(rsamp(0:nlen),nsamp(1:nlen))
      !
      kslec=64
      print*,' ** selected k position=',kslec
      !
      tcoru=0.d0
      tcorv=0.d0
      tcorw=0.d0
      !
      rsamp=0.d0
      ncoun=0
      !
    endif
    !
    allocate(ro(0:jm,0:km),u1(0:jm,0:km),u2(0:jm,0:km),u3(0:jm,0:km),  &
              p(0:jm,0:km), t(0:jm,0:km))
    do n=nsta,nend
      !
      write(fname,'(I5.5)')n
      !
      ! filename='inflow/islice'//fname//'.h5'
      filename='islice/islice'//fname//'.h5'
      !
      inquire(file=trim(filename),exist=lfilalive)
      !
      if(lfilalive) then
        !
        call H5ReadArray(time,'time',trim(filename),explicit=.false.)
        call H5ReadArray(p,jm,km, 'p',trim(filename),explicit=.false.)
        call H5ReadArray(u1,jm,km,'u1',trim(filename),explicit=.false.)
        call H5ReadArray(u2,jm,km,'u2',trim(filename),explicit=.false.)
        call H5ReadArray(u3,jm,km,'u3',trim(filename),explicit=.false.)
        call H5ReadArray( t,jm,km, 't',trim(filename),explicit=.false.)
        !
        ! p=ro*t/const2
        ro=p/t*const2
        !
        nsamples=nsamples+1.d0
        !
        write(*,'(1A1,3(A),I0,A,$)')char(13),'  >> ',trim(filename),   &
                                   ' ... ',(n-nsta)*100/(nend-nsta),' %'
      else
        print*,' ** ',trim(filename),' is missing.'
        cycle
      endif
      !
      if(mode=='mean') then
        do k=1,km
          rom(:)=rom(:)+ro(:,k)
          u1m(:)=u1m(:)+ro(:,k)*u1(:,k)
          u2m(:)=u2m(:)+ro(:,k)*u2(:,k)
          u3m(:)=u3m(:)+ro(:,k)*u3(:,k)
           tm(:)= tm(:)+ro(:,k)* t(:,k)
           pm(:)= pm(:)+ p(:,k)
        enddo
        !
        rotm=rotm+ro
        u1tm=u1tm+u1
        u2tm=u2tm+u2
        u3tm=u3tm+u3
         ptm= ptm+ t
         ttm= ttm+ p
        !
      elseif(mode=='horder') then
        !
        do k=0,km
          u1(:,k)=u1(:,k)-u1m(:)
          u2(:,k)=u2(:,k)-u2m(:)
          u3(:,k)=u3(:,k)-u3m(:)
        enddo
        !
        do k=1,km
          uu=uu+ro(:,k)*u1(:,k)*u1(:,k)
          vv=vv+ro(:,k)*u2(:,k)*u2(:,k)
          ww=ww+ro(:,k)*u3(:,k)*u3(:,k)
          uv=uv+ro(:,k)*u1(:,k)*u2(:,k)
          uw=uw+ro(:,k)*u1(:,k)*u3(:,k)
          vw=vw+ro(:,k)*u2(:,k)*u3(:,k)
        enddo
        !
        do k=0,km/2
          do k1=0,km
            k2=k1+k
            if(k2>km) k2=k2-km
            zcoru(:,k)=zcoru(:,k)+u1(:,k1)*u1(:,k2)
            zcorv(:,k)=zcorv(:,k)+u2(:,k1)*u2(:,k2)
            zcorw(:,k)=zcorw(:,k)+u3(:,k1)*u3(:,k2)
          enddo
        enddo
        !
      elseif(mode=='tcorr') then
        !
        ! u1=u1-u1tm
        ! u2=u2-u2tm
        ! u3=u3-u3tm
        !
        do k=0,km
          u1(:,k)=u1(:,k)-u1m(:)
          u2(:,k)=u2(:,k)-u2m(:)
          u3(:,k)=u3(:,k)-u3m(:)
        enddo
        !
        if(ncoun<nlen) then
          !
          ncoun=ncoun+1
          !
          nsamp(ncoun)=n
          !
          u1_i(:,:,ncoun)=u1(:,:)
          u2_i(:,:,ncoun)=u2(:,:)
          u3_i(:,:,ncoun)=u3(:,:)
          !
        elseif(ncoun==nlen) then
          !
          ! cal. coorelation
          ncoun=ncoun+1
          !
          do n1=1,nlen
            do n2=n1,nlen
              !
              ntc=nsamp(n2)-nsamp(n1)
              !
              ! do k=kslec,kslec
              do k=1,km
              do j=0,jm
                tcoru(j,ntc)=tcoru(j,ntc)+u1_i(j,k,n1)*u1_i(j,k,n2)
                tcorv(j,ntc)=tcorv(j,ntc)+u2_i(j,k,n1)*u2_i(j,k,n2)
                tcorw(j,ntc)=tcorw(j,ntc)+u3_i(j,k,n1)*u3_i(j,k,n2)
              enddo
              enddo
              !
              rsamp(ntc)=rsamp(ntc)+1.d0
              !
            enddo
          enddo
          !
          ! for the newly loaded sample
          do n1=1,nlen
            !
            ntc=n-nsamp(n1)
            !
            ! do k=kslec,kslec
            do k=1,km
            do j=0,jm
              tcoru(j,ntc)=tcoru(j,ntc)+u1_i(j,k,n1)*u1(j,k)
              tcorv(j,ntc)=tcorv(j,ntc)+u2_i(j,k,n1)*u2(j,k)
              tcorw(j,ntc)=tcorw(j,ntc)+u3_i(j,k,n1)*u3(j,k)
            enddo
            enddo
            !
            rsamp(ntc)=rsamp(ntc)+1.d0
            !
          enddo
          !
          ! replace the data
          !
          do n1=1,nlen-1
            u1_i(:,:,n1)=u1_i(:,:,n1+1)
            u2_i(:,:,n1)=u2_i(:,:,n1+1)
            u3_i(:,:,n1)=u3_i(:,:,n1+1)
            !
            nsamp(n1)=nsamp(n1+1)
          enddo
          !
          u1_i(:,:,nlen)=u1(:,:)
          u2_i(:,:,nlen)=u2(:,:)
          u3_i(:,:,nlen)=u3(:,:)
          !
          nsamp(nlen)=n
          !
        elseif(ncoun>nlen) then
          !
          ncoun=ncoun+1
          !
          ! for the newly loaded sample
          do n1=1,nlen
            !
            ntc=n-nsamp(n1)
            !
            ! do k=kslec,kslec
            do k=1,km
            do j=0,jm
              tcoru(j,ntc)=tcoru(j,ntc)+u1_i(j,k,n1)*u1(j,k)
              tcorv(j,ntc)=tcorv(j,ntc)+u2_i(j,k,n1)*u2(j,k)
              tcorw(j,ntc)=tcorw(j,ntc)+u3_i(j,k,n1)*u3(j,k)
            enddo
            enddo
            !
            rsamp(ntc)=rsamp(ntc)+1.d0
            !
          enddo
          !
          ! replace the data
          !
          do n1=1,nlen-1
            u1_i(:,:,n1)=u1_i(:,:,n1+1)
            u2_i(:,:,n1)=u2_i(:,:,n1+1)
            u3_i(:,:,n1)=u3_i(:,:,n1+1)
            !
            nsamp(n1)=nsamp(n1+1)
          enddo
          !
          u1_i(:,:,nlen)=u1(:,:)
          u2_i(:,:,nlen)=u2(:,:)
          u3_i(:,:,nlen)=u3(:,:)
          !
          nsamp(nlen)=n
          !
          !
        endif
        !
      else
        stop ' !! mode not defined !!'
      endif
      !
      ! print*,' ** time=',time
      !
    enddo
    !
    call writetecbin('Results/tecyz'//fname//'.plt',x,'x',y,'y',z,'z', &
                         ro,'ro',u1,'u',u2,'v',u3,'w',p,'p',t,'t',jm,km)
    !
    if(mode=='mean') then
      u1m=u1m/rom
      u2m=u2m/rom
      u3m=u3m/rom
       tm= tm/rom
      !
      var1=nsamples*dble(km)
      rom=rom/var1
       pm= pm/var1
      !
       rotm=rotm/nsamples
       u1tm=u1tm/nsamples
       u2tm=u2tm/nsamples
       u3tm=u3tm/nsamples
        ptm= ptm/nsamples
        ttm= ttm/nsamples
       !
      inquire(file='Results/meanprofile.h5',exist=lfilalive)
      if(lfilalive) call system('mv  -v Results/meanprofile.h5 Results/meanprofile.bak')
      !
      call H5WriteArray(rom,jm,'ro','Results/meanprofile.h5')
      call H5WriteArray(u1m,jm,'u1','Results/meanprofile.h5')
      call H5WriteArray(u2m,jm,'u2','Results/meanprofile.h5')
      call H5WriteArray(u3m,jm,'u3','Results/meanprofile.h5')
      call H5WriteArray( tm,jm, 't','Results/meanprofile.h5')
      call H5WriteArray( pm,jm, 'p','Results/meanprofile.h5')
      !
      inquire(file='Results/meanflow.h5',exist=lfilalive)
      if(lfilalive) call system('mv  -v Results/meanflow.h5 Results/meanflow.bak')
      !
      call H5WriteArray(rotm,jm,km,'ro','Results/meanflow.h5')
      call H5WriteArray(u1tm,jm,km,'u1','Results/meanflow.h5')
      call H5WriteArray(u2tm,jm,km,'u2','Results/meanflow.h5')
      call H5WriteArray(u3tm,jm,km,'u3','Results/meanflow.h5')
      call H5WriteArray( ttm,jm,km, 't','Results/meanflow.h5')
      call H5WriteArray( ptm,jm,km, 'p','Results/meanflow.h5')
      !
      ! do j=0,jm
      !   write(*,*)y(j,0),rom(j),u1m(j),tm(j)
      ! enddo
      !
    elseif(mode=='horder') then
      !
      var1=nsamples*dble(km)
      !
      uu=uu/rom/var1
      vv=vv/rom/var1
      ww=ww/rom/var1
      uv=uv/rom/var1
      uw=uw/rom/var1
      vw=vw/rom/var1
      !
      inquire(file='Results/restress.h5',exist=lfilalive)
      if(lfilalive) call system('mv  -v Results/restress.h5 Results/restress.bak')
      !
      call H5WriteArray(uu,jm,'uu','Results/restress.h5')
      call H5WriteArray(vv,jm,'vv','Results/restress.h5')
      call H5WriteArray(ww,jm,'ww','Results/restress.h5')
      call H5WriteArray(uv,jm,'uv','Results/restress.h5')
      call H5WriteArray(uw,jm,'uw','Results/restress.h5')
      call H5WriteArray(vw,jm,'vw','Results/restress.h5')
      !
      do j=0,jm
        var1=zcoru(j,0)
        zcoru(j,:)=zcoru(j,:)/var1
        var1=zcorv(j,0)
        zcorv(j,:)=zcorv(j,:)/var1
        var1=zcorw(j,0)
        zcorw(j,:)=zcorw(j,:)/var1
      enddo
      !
      inquire(file='Results/correlation.h5',exist=lfilalive)
      if(lfilalive) call system('mv  -v Results/correlation.h5 Results/correlation.bak')
      !
      call H5WriteArray(zcoru,jm,km/2,'coru','Results/correlation.h5')
      call H5WriteArray(zcorv,jm,km/2,'corv','Results/correlation.h5')
      call H5WriteArray(zcorw,jm,km/2,'corw','Results/correlation.h5')
      !
    elseif(mode=='tcorr') then
      !
      do ntc=0,nlen
        tcoru(:,ntc)=tcoru(:,ntc)/rsamp(ntc)/dble(km)
        tcorv(:,ntc)=tcorv(:,ntc)/rsamp(ntc)/dble(km)
        tcorw(:,ntc)=tcorw(:,ntc)/rsamp(ntc)/dble(km)
      enddo
      !
      do j=0,jm
        var1=tcoru(j,0)
        tcoru(j,:)=tcoru(j,:)/var1
        var1=tcorv(j,0)
        tcorv(j,:)=tcorv(j,:)/var1
        var1=tcorw(j,0)
        tcorw(j,:)=tcorw(j,:)/var1
      enddo
      !
      inquire(file='Results/tcorrelation.h5',exist=lfilalive)
      if(lfilalive) call system('mv  -v Results/tcorrelation.h5 Results/tcorrelation.bak')
      !
      call H5WriteArray(tcoru,jm,nlen,'tcoru','Results/tcorrelation.h5')
      call H5WriteArray(tcorv,jm,nlen,'tcorv','Results/tcorrelation.h5')
      call H5WriteArray(tcorw,jm,nlen,'tcorw','Results/tcorrelation.h5')
      !
    endif
    !
  end subroutine isliceanalys
  !+-------------------------------------------------------------------+
  !| The end of the subroutine isliceanalys.                           |
  !+-------------------------------------------------------------------+
  !
  subroutine shockjump
    !
    use commvardefine,only: Reynolds,mach,roinf,uinf,vinf,roinf,tinf,  &
                            pinf,pi,gamma,const2
    use basicfunction,only: fxsolver
    !
    real(8) :: ro0,u0,v0,p0,t0,ro1,u1,v1,p1,t1,ro2,u2,v2,p2,t2,ro3,u3,v3,p3,t3 
    real(8) :: delta1,beta0,beta1,delta2,beta2,beta1_target,delta1_target
    real(8) :: var1,var2
    !
    !
    ! u0 =1.d0 !1.0018d0
    ! v0 =0.d0
    ! ro0=1.d0 !1.017d0
    ! t0 =1.d0/ro0
    ! p0 =ro0*t0/const2 
    ! print*,' ** mach0=',u0/sqrt(t0)*mach
    ! !
    ! var1=5.d0
    ! var2=25.d0
    ! beta0=fxsolver(xmin=var1,xmax=var2,ftarget=0.01030631d0,fx=ledv)
    ! !
    ! print*,' ** leading shock angle=',beta0
    ! print*,' **       deflect angle=',atan(v1/u1)/pi*180.d0
    ! print*,' ** ro0=',ro0,'ro1',ro1
    ! print*,' **  u0=', u0,' u1', u1
    ! print*,' **  v0=', v0,' v1', v1
    ! print*,' **  T0=', T0,' T1', T1
    ! print*,' **  p0=', p0,' p1', p1,'p1/p0=',p1/p0
    !
    u1 =0.9995405d00 
    v1 =0.1030631d-1 
    t1 =0.1003848d1  
    ro1=0.9964236d0  
    ! u1 =1.d0
    ! v1 =0.d0
    ! t1 =1.d0
    ! ro1=1.d0
    p1 =ro1*t1/const2 
    print*,' ** mach1=',u1/sqrt(t1)*mach
    !
    !
    var1=-20.d0
    var2=-30.d0
    !0.12548336709962035
    !4.3919178484867123
    beta1=fxsolver(xmin=var1,xmax=var2,ftarget=-13.44d0,fx=impang)
    ! beta1=fxsolver(xmin=var1,xmax=var2,ftarget=4.3919178484867123*p1,fx=postprs)
    !
    print*,' ** impinging angle=',beta1
    print*,' **   deflect angle=',atan(v2/u2)/pi*180.d0
    print*,' ** ro1=',ro1,'ro2',ro2
    print*,' **  u1=', u1,' u2', u2
    print*,' **  v1=', v1,' v2', v2
    print*,' **  T1=', T1,' T2', T2
    print*,' **  p1=', p1,' p2', p2,'p2/p1=',p2/p1
    !
    var1=10.d0
    var2=30.d0
    !
    beta2=fxsolver(xmin=var1,xmax=var2,ftarget=0.d0,fx=refang)
    !
    print*,' ** reflectio angle=',beta2
    print*,' **   deflect angle=',atan(v3/u3)/pi*180.d0
    print*,' ** ro3',ro3
    print*,' **  u3', u3
    print*,' **  v3', v3
    print*,' **  T3', T3
    print*,' **  p3', p3,'p3/p1=',p3/p1,'p3/p2=',p3/p2
    !
    contains
    !
    function ledv(beta)
      real(8) :: ledv,beta
      !
      call PostShockCal(ro0,u0,v0,p0,t0,                            &
                        ro1,u1,v1,p1,t1,beta,mach)
      ledv=v1
    end function ledv
    !
    function impang(beta)
      real(8) :: impang,beta
      !
      call PostShockCal(ro1,u1,v1,p1,t1,                            &
                        ro2,u2,v2,p2,t2,beta,mach)
      impang=atan(v2/u2)/pi*180.d0
    end function impang
    !
    function refang(beta)
      real(8) :: refang,beta
      !
      call PostShockCal(ro2,u2,v2,p2,t2,                            &
                        ro3,u3,v3,p3,t3,beta,mach)
      refang=atan(v3/u3)/pi*180.d0
    end function refang
    !
    function postprs(beta)
      real(8) :: postprs,beta
      !
      call PostShockCal(ro1,u1,v1,p1,t1,                            &
                        ro2,u2,v2,p2,t2,beta,mach)
      postprs=p2
    end function postprs
    !
  end subroutine shockjump
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to to an inviscid analysis of flow before |
  !| after shock in the compression-expansion corner flow.             |
  !+-------------------------------------------------------------------+
  subroutine shockexpcom(rocc,ucc,vcc,pcc,tcc,roec,uec,vec,pec,tec)
    !
    use commvardefine,only: Reynolds,mach,roinf,uinf,vinf,roinf,tinf,  &
                            pinf,pi,gamma,const2
    !
    real(8),optional,intent(out) :: rocc,ucc,vcc,pcc,tcc,roec,uec,vec,pec,tec      
    !
    real(8) :: beta
    integer :: j,ncoun
    real(8) :: beter,beter1,beter2,error,beta_tar
    real(8) :: rops0,ups0,vps0,pps0,Tps0,rops1,ups1,vps1,pps1,Tps1,    &
               roexp,uexp,vexp,pexp,texp,rops2,ups2,vps2,pps2,Tps2
    real(8) :: dbeter,err,vsave,var1,var2,theter,machps0,machps2,machpe,dmach
    real(8) :: xiw1,xiw2,ly,dvp,dvpsave
    !
    print*,' >> input target angle 1 of flow: '
    read(*,*)beta
    !
    ! beta=-23.33d0
    ! beta=-22.76d0
    !
    beter=beta
    !      
    ! rops0=0.9964236E+00      
    ! ups0=0.9995405E+00
    ! vps0=0.1030631E-01
    ! Tps0=0.1003848E+01
    !
    rops0=1.d0
    ups0=1.d0
    vps0=0.d0
    Tps0=1.d0
    !
    pps0=rops0*Tps0/const2 
    !
    call PostShockCal(rops0,ups0,vps0,pps0,Tps0,                       &
                                   rops1,ups1,vps1,pps1,Tps1,beter,mach)
    !
    print*,' ** rops1=',rops1
    print*,' **  ups1=', ups1
    print*,' **  vps1=', vps1
    print*,' **  Tps1=', Tps1
    print*,' **  pps1=', pps1
    !
    ! beter=beta_tar+3.d0
    ! error=1.d10
    ! dbeter=0.01d0
    ! ncoun=0
    ! var2=0.d0
    ! do while(abs(error)>1.d-3 )
    !   !
    !   ncoun=ncoun+1
    !   !
    !   call PostShockCal(rops0,ups0,vps0,pps0,Tps0,                     &
    !                                     rops1,ups1,vps1,pps1,Tps1,beter)
    !   !
    !   beter1=atan(vps0/ups0)/pi*180.d0
    !   beter2=atan(vps1/ups1)/pi*180.d0
    !   !
    !   var1=beter2-beta_tar
    !   if(error*var1<0) dbeter=0.5d0*dbeter
    !   error=var1
    !   !
    !   if(error>0.d0) then
    !     beter=beter-dbeter
    !   elseif(error<0.d0) then
    !     beter=beter+dbeter
    !   else
    !     stop 'error 1 @ shockangle'
    !   endif
    !   !
    !   print*,' ** ', ncoun,beter,beter2,error
    !   !
    !   var2=beter2
    !   !
    ! enddo
    !
    ! call PostShockCal(rops0,ups0,vps0,pps0,Tps0,                     &
    !                                     rops1,ups1,vps1,pps1,Tps1,beter,mach)
    beter1=atan(vps0/ups0)/pi*180.d0
    beter2=atan(vps1/ups1)/pi*180.d0
    !
    print*,' **     shock angle:',beter,' deg'
    print*,' ** deflection angle:',beter2,' deg'
    write(*,'(2X,A)')'-------------------- before shock --------------------'
    write(*,'(2X,A)')'       ro        u        v        p        t        θ'
    write(*,'(2X,6(F9.3))')rops0,ups0,vps0,pps0,Tps0,beter1
    write(*,'(2X,A)')'--------------------  after shock --------------------'
    write(*,'(2X,A)')'       ro        u        v        p        t        θ'
    write(*,'(2X,6(F9.3))')rops1,ups1,vps1,pps1,Tps1,beter2
    write(*,'(2X,A)')'------------------------------------------------------'
    !
    machps0=sqrt((ups0**2+vps0**2)/Tps0)*Mach
    machps2=sqrt((ups1**2+vps1**2)/Tps1)*Mach
    print*,' ** mach number before shock',machps0
    print*,' ** mach number after shock',machps2
    print*,' ** pressure ratio before/after shock',pps1/pps0
    !
    print*,' >> input target angle 2 of flow: '
    read(*,*)beter
    !
    call PostShockCal(rops1,ups1,vps1,pps1,Tps1,                     &
                                        rops2,ups2,vps2,pps2,Tps2,beter,mach)
    beter2=atan(vps2/ups2)/pi*180.d0
    !
    print*,' **     shock angle:',beter,' deg'
    print*,' ** deflection angle:',beter2,' deg'
    write(*,'(2X,A)')'-------------------- before shock --------------------'
    write(*,'(2X,A)')'       ro        u        v        p        t        θ'
    write(*,'(2X,6(F9.3))')rops1,ups1,vps1,pps1,Tps1,beter1
    write(*,'(2X,A)')'--------------------  after shock --------------------'
    write(*,'(2X,A)')'       ro        u        v        p        t        θ'
    write(*,'(2X,6(F9.3))')rops2,ups2,vps2,pps2,Tps2
    write(*,'(2X,A)')'------------------------------------------------------'
    !
    machps2=sqrt((ups2**2+vps2**2)/Tps2)*Mach
    print*,' ** mach number after shock',machps2
    !
    print*,' ** pressure ratio before/after shock',pps2/pps1
    print*,' ** totao pressure ratio ',pps2/pps0
    !
    stop
    !
    ! var1=(1.d0+0.5d0*(gamma-1.d0)*machps2**2)**(gamma/(gamma-1.d0))
    ! print*,' ** totao pressure ratio ',pps1/pps0*var1
    !
    machpe=2.9d0
    error=1.d10
    dmach=0.1d0
    ncoun=0
    do while(abs(error)>1.d-3 .and. ncoun<15)
      !
      ncoun=ncoun+1
      theter=PrandtlMeyer(machpe)-PrandtlMeyer(machps2)
      theter=theter/pi*180.d0
      !
      var1=theter-beta_tar
      if(error*var1<0) dmach=0.5d0*dmach
      error=var1
      !
      if(abs(error)<1.d-3) exit
      !
      if(error>0.d0) then
        machpe=machpe-dmach
      elseif(error<0.d0) then
        machpe=machpe+dmach
      else
        stop 'error 2 @ shockangle'
      endif
      !
      ! print*,' ** ', ncoun,machpe,dmach,theter
      !
    enddo
    print*,' ** turnning angle number after expansion wave',theter,', Mach=',machpe
    !
    var1=1.d0+0.5d0*(gamma-1.d0)*machps2**2
    var2=1.d0+0.5d0*(gamma-1.d0)*machpe**2
    pexp=(var1/var2)**(gamma/(gamma-1.d0))*pps1
    roexp=(var1/var2)**(1.d0/(gamma-1.d0))*rops1
    texp=(var1/var2)*Tps1
    !
    uexp=sqrt(texp)/mach*machpe*cos(theter/180.d0*pi)
    vexp=sqrt(texp)/mach*machpe*sin(theter/180.d0*pi)
    !
    write(*,'(2X,A)')'----------------  after expansion fan ----------------'
    write(*,'(2X,A)')'       ro        u        v        p        t        θ'
    write(*,'(2X,6(F9.3))')roexp,uexp,vexp,pexp,texp,theter
    write(*,'(2X,A)')'------------------------------------------------------'
    !
    print*,' ** pressure ratio after Prandtl-Meyer expansion',pexp/pps0
    !
    if(present(pcc)) pcc=pps1
    if(present(pec)) pec=pexp
    !
  end subroutine shockexpcom
  !+-------------------------------------------------------------------+
  !| The end of the subroutine shockexpcom.                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to post-process flow pass a blunt body.   |
  !+-------------------------------------------------------------------+
  subroutine nbody(ipost)
    !
    use commvardefine,only: im,jm,km,Reynolds,pinf,pi,mach,Prandtl,gamma
    use h5readwrite
    use basicfunction
    use interpolation 
    use writetec
    use boundarylayer
    use gradsolver
    !
    integer,intent(in),optional :: ipost
    !
    real(8),allocatable,dimension(:,:) :: x,y,ro,u1,u2,p,t,            &
                                          dudx,dudy,dvdx,dvdy,omegaz,  &
                                          dtdx,dtdy,ma
    real(8),allocatable,dimension(:) :: z,taw,pwall,cf,utau,lplus,alfa,&
                                        qw,ch,cp,st,tw,delta_alfa,ttau
    real(8),allocatable,dimension(:,:) :: wallnorm
    real(8) :: div,s11,s12,s13,s22,s23,s33,f11,f12,f13,f22,f23,f33,    &
               var1,miu,tawx,tawy,delta,pwref,ttotal,alfamin,ds,dr,dz
    real(8),allocatable :: total_enthalpy(:),total_energy(:)
    real(8),allocatable,dimension(:,:) :: convect,pre_ace,pre_dil,     &
                                          pre_tra,product,tur_tra,     &
                                          vis_ace,vis_dif,dissipa,balance
    real(8),allocatable :: du1(:,:,:),du2(:,:,:),dt(:,:,:)
    real(8),allocatable :: theter(:,:),radius(:,:),us(:),dis2wall(:),  &
                           yplus(:),uplus(:)
    real(8),allocatable,dimension(:,:) :: u11,u22,u33,u12,tke,tt,pp,tu1,tu2
    real(8),allocatable,dimension(:) :: uss,unn,usn,tus,tun
    real(8) :: vec_s(2),vec_n(2),dis,lkhomo
    integer :: i,j,k,i0,iref,n,jedge
    character(len=4) :: fname

    !
    !+-------+
    !| begin |
    !+-------+
    !
    print*,' ** processing data of a flow past a cylinder'
    !
    allocate(x(0:im,0:jm),y(0:im,0:jm),z(0:km))
    call H5ReadSubset(x,im,jm,km,'x','datin/grid.h5',kslice=0)
    call H5ReadSubset(y,im,jm,km,'y','datin/grid.h5',kslice=0)
    call H5ReadSubset(z,im,jm,km,'z','datin/grid.h5',islice=0,jslice=0)
    !
    allocate(theter(0:im,0:jm),radius(0:im,0:jm))
    do j=0,jm
    do i=0,im
      radius(i,j)=sqrt((x(i,j)-1.d0)**2+y(i,j)**2)
      theter(i,j)=asin(y(i,j)/radius(i,j))/pi*180.d0
    enddo
    enddo
    !
    allocate(ro(0:im,0:jm),u1(0:im,0:jm),u2(0:im,0:jm),p(0:im,0:jm),   &
             t(0:im,0:jm),dudx(0:im,0:jm),dudy(0:im,0:jm),             &
             dvdx(0:im,0:jm),dvdy(0:im,0:jm),omegaz(0:im,0:jm),        &
             dtdx(0:im,0:jm),dtdy(0:im,0:jm),ma(0:im,0:jm),            &
             qw(0:im),ch(0:im),cp(0:im),st(0:im),tw(0:im))
    !
    call H5ReadArray(ro,im,jm,'ro','Results/mean.fav.zm.h5')
    call H5ReadArray(u1,im,jm,'u1','Results/mean.fav.zm.h5')
    call H5ReadArray(u2,im,jm,'u2','Results/mean.fav.zm.h5')
    call H5ReadArray(p, im,jm, 'p','Results/mean.fav.zm.h5')
    call H5ReadArray(t, im,jm, 't','Results/mean.fav.zm.h5')
    call H5ReadArray(dudx,im,jm,'du1dx','Results/mean.fav.zm.h5')
    call H5ReadArray(dudy,im,jm,'du1dy','Results/mean.fav.zm.h5')
    call H5ReadArray(dvdx,im,jm,'du2dx','Results/mean.fav.zm.h5')
    call H5ReadArray(dvdy,im,jm,'du2dy','Results/mean.fav.zm.h5')
    call H5ReadArray(dtdx,im,jm,'dtdx','Results/mean.fav.zm.h5')
    call H5ReadArray(dtdy,im,jm,'dtdy','Results/mean.fav.zm.h5')
    !
    do j=0,jm
      do i=0,im
        ma(i,j)=sqrt(u1(i,j)**2+u2(i,j)**2)/sqrt(t(i,j))*mach
        omegaz(i,j)=dudy(i,j)-dvdx(i,j)
      enddo
    enddo
    !
    call writetecbin('Results/tecmeanflow.plt',x,'x',y,'y',            &
                              ro,'ro',u1,'u',u2,'v', p,'p',t,'T',      &
                              ma,'Mach',im,jm)
    ! fname='mean'
    !
    !
    allocate(du1(2,0:im,0:jm),du2(2,0:im,0:jm),dt(2,0:im,0:jm))
    allocate(wallnorm(2,0:im),taw(0:im),pwall(0:im),cf(0:im))
    allocate(utau(0:im),ttau(0:im),lplus(0:im),alfa(0:im))
    !
    ! do n=64,94,10
      
    !   write(fname,'(I4.4)')n
    !   !
    !   call H5ReadSubset(ro,im,jm,km,'ro',outfolder//'flowfield'//fname//'.h5',kslice=256)
    !   call H5ReadSubset(u1,im,jm,km,'u1',outfolder//'flowfield'//fname//'.h5',kslice=256)
    !   call H5ReadSubset(u2,im,jm,km,'u2',outfolder//'flowfield'//fname//'.h5',kslice=256)
    !   call H5ReadSubset (p,im,jm,km, 'p',outfolder//'flowfield'//fname//'.h5',kslice=256)
    !   call H5ReadSubset( t,im,jm,km, 't',outfolder//'flowfield'//fname//'.h5',kslice=256)
    !   !
    !   du1=grad_xy(u1,x,y)
    !   du2=grad_xy(u2)
    !   dt =grad_xy(t)
    !   dudx=du1(1,:,:)
    !   dudy=du1(2,:,:)
    !   dvdx=du2(1,:,:)
    !   dvdy=du2(2,:,:)
    !   dtdx=dt(1,:,:)
    !   dtdy=dt(2,:,:)
      !
      ! print*,' ** calcualte wall variables'
      !
      call gridjacobian_xy(x,y,wallnormal=wallnorm)
      !
      pwall(:)=p(:,0)
      !
      pwref=pwall(0)
      pwall=pwall/pwref
      !
      alfamin=1.d10
      do i=0,im
        !
        alfa(i)=asin(y(i,0)/sqrt((x(i,0)-1.d0)**2+y(i,0)**2))/pi*180.d0
        !
        tw(i)=t(i,0)
        !
        miu=miucal(t(i,0))/Reynolds
        !
        div=0.333333333333333d0*(dudx(i,0)+dvdy(i,0))
        s11=dudx(i,0)-div
        s12=0.5d0*(dudy(i,0)+dvdx(i,0))
        s22=dvdy(i,0)-div
        !
        f11=2.d0*miu*s11
        f12=2.d0*miu*s12
        f22=2.d0*miu*s22
        !
        tawx=f11*wallnorm(1,i)+f12*wallnorm(2,i)
        tawy=f12*wallnorm(1,i)+f22*wallnorm(2,i)
        !
        ! taw(i)=sqrt(tawx**2+tawy**2)
        taw(i)=tawx*wallnorm(2,i)-tawy*wallnorm(1,i)
        !
        utau(i)=sqrt(taw(i)/ro(i,0))
        lplus(i)=ro(i,0)*utau(i)/miu
        !
        cf(i)=taw(i)*2.d0
        cp(i)=p(i,0)*2.d0
        !
        qw(i)=miu*(dtdx(i,0)*wallnorm(1,i)+dtdy(i,0)*wallnorm(2,i))
        ch(i)=2.d0*qw(i)/(Mach**2*Prandtl*(gamma-1))
        ttotal=1.d0+0.5d0*(gamma-1.d0)*Mach**2
        st(i)=qw(i)/((ttotal-t(i,0)))/Prandtl
        ttau(i)=qw(i)/utau(i)/ro(i,0)/Prandtl
        !
        if(abs(alfa(i))<alfamin) then
          alfamin=abs(alfa(i))
          i0=i
        endif
        !
      enddo
      !
      print*,' ** the stagnation line at i=',i0,'alfa=',alfa(i0)
      !
      ! if(present(ipost)) then
      !   iref=ipost
      ! else
        do i=1,im
          !
          if(alfa(i)>=60.d0 .and. alfa(i-1)<60.d0) then
            iref=i
          endif
          !
        enddo
        !
      ! endif
      !
      write(fname,'(I4.4)')iref
      !
      print*,' ** the reference location at i=',iref,'alfa=',alfa(iref)
      print*,' ** taw(iref)=',taw(iref),dudx(iref,0),dudy(iref,0)
      print*,' ** utau(iref)=',utau(iref)
      !
      do j=0,jm
        dis=sqrt((x(i,j)-x(i,0))**2+(y(i,j)-y(i,0))**2)
        !
        if(dis>0.02d0) then
          jedge=j
          exit
        endif
      enddo
      !
      print*,' ** the edge of BL at j=',jedge,'d=',dis
      !
      open(18,file='Results/surface_coefficient.dat')
      write(18,"(8(1X,A15))")'x','alfa','cf','cp','ch','st','cf/ch','twall'
      write(18,"(8(1X,E15.7E3))")(x(i,0),alfa(i),cf(i),cp(i),ch(i),st(i),&
                                                 cf(i)/ch(i),tw(i),i=0,im)
      close(18)
      print*,' << surface_coefficient.dat ... done !'
      !
      allocate(total_enthalpy(0:jm),total_energy(0:jm))
      i=i0
      do j=0,jm
        total_energy(j)=p(i,j)/ro(i,j)/(gamma-1.d0)+                      &
                                             0.5d0*(u1(i,j)**2+u2(i,j)**2)
        total_enthalpy(j)=total_energy(j)+p(i,j)/ro(i,j)
      enddo
      !
      open(18,file='Results/stagnation_line.dat')
      write(18,"(5(1X,A15))")'x','alfa','T','total_energy','total_enthalpy'
      write(18,"(5(1X,E15.7E3))")(x(i,j),alfa(i),t(i,j),total_energy(j),&
                                                 total_enthalpy(j),j=0,jm)
      close(18)
      print*,' << stagnation_line.dat ... done !'
      !
      deallocate(total_enthalpy,total_energy)
      !
    ! enddo
    !
    open(18,file='Results/mesh_resolution.dat')
    write(18,"(5(1X,A15))")'alfa','ds+','dr+','dz+','delta_alfa'
    do i=1,im
      ds=sqrt((x(i,0)-x(i-1,0))**2+(y(i,0)-y(i-1,0))**2)*lplus(i)
      dr=sqrt((x(i,1)-x(i,0))**2+(y(i,1)-y(i,0))**2)*lplus(i)
      dz=(z(1)-z(0))*lplus(i)
      ! write(*,*)alfa(i),ds,dr,dz,alfa(i)-alfa(i-1)
      write(18,*)alfa(i),ds,dr,dz,alfa(i)-alfa(i-1)
    enddo
    close(18)
    print*,' << mesh_resolution.dat'
    !
    allocate(u11(0:im,0:jm),u22(0:im,0:jm),u33(0:im,0:jm),             &
             u12(0:im,0:jm),tke(0:im,0:jm),tt(0:im,0:jm),pp(0:im,0:jm),&
             tu1(0:im,0:jm),tu2(0:im,0:jm))
    !
    call H5ReadArray(u11,im,jm,'u11','Results/2order.fav.zm.h5')
    call H5ReadArray(u22,im,jm,'u22','Results/2order.fav.zm.h5')
    call H5ReadArray(u33,im,jm,'u33','Results/2order.fav.zm.h5')
    call H5ReadArray(u12,im,jm,'u12','Results/2order.fav.zm.h5')
    call H5ReadArray( tt,im,jm, 'tt','Results/2order.fav.zm.h5')
    call H5ReadArray(tu1,im,jm, 'tu1','Results/2order.fav.zm.h5')
    call H5ReadArray(tu2,im,jm, 'tu2','Results/2order.fav.zm.h5')
    call H5ReadArray( pp,im,jm, 'pp','Results/2order.fav.zm.h5')
    !
    tke=0.5d0*(u11+u22+u33)
    !
    do j=0,jm
    do i=0,im
      var1=ro(i,j)/ro(iref,0)/utau(iref)/utau(iref)
      u11(i,j)=u11(i,j)*var1
      u22(i,j)=u22(i,j)*var1
      u33(i,j)=u33(i,j)*var1
      u12(i,j)=u12(i,j)*var1
      tke(i,j)=tke(i,j)*var1
      tt(i,j)=tt(i,j)/ttau(iref)/ttau(iref)
      pp(i,j)=pp(i,j)/pinf/pinf
    end do
    end do
    !
    call writetecbin('Results/tec2order.plt',x,'x',y,'y',              &
                         theter,'theter',radius,'radius',              &
                         u11,'uu',u22,'vv',u33,'ww',u12,'uv',tke,'TKE',&
                         tt,'tt',pp,'pp',im,jm)
    !
    i=iref
    print*,' ** writing profile at i=',i
    !
    allocate(us(0:jm),dis2wall(0:jm))
    allocate(yplus(0:jm),uplus(0:jm))
    !
    do j=0,jm
      dis2wall(j)=sqrt((x(i,j)-x(i,0))**2+(y(i,j)-y(i,0))**2)
      !
      us(j)=sin(alfa(i)/180.d0*pi)*u1(i,j)+cos(alfa(i)/180.d0*pi)*u2(i,j)
      !
    enddo
    call upluscal(uplus=uplus,yplus=yplus,u=us,y=dis2wall,             &
                                        ro=ro(i,:),tw=t(i,0))
    !
    open(18,file='Results/uplus.alfa60.dat')
    write(18,"(2(1X,A15))")'yplus','uplus'
    write(18,"(2(1X,E15.7E3))")(yplus(j),uplus(j),j=0,jm)
    close(18)
    print*,' << uplus.alfa60.dat ... done. '
    !
    !
    vec_s(1)=sin(alfa(i)/180.d0*pi); vec_n(1)=-cos(alfa(i)/180.d0*pi)
    vec_s(2)=cos(alfa(i)/180.d0*pi); vec_n(2)= sin(alfa(i)/180.d0*pi)
    allocate(uss(0:jm),unn(0:jm),usn(0:jm),tus(0:jm),tun(0:jm))
    do j=0,jm
      uss(j)=vec_s(1)**2*u11(i,j)+vec_s(2)**2*u22(i,j) +               &
                                         2.d0*vec_s(1)*vec_s(2)*u12(i,j)
      unn(j)=vec_n(1)**2*u11(i,j)+vec_n(2)**2*u22(i,j) +               &
                                         2.d0*vec_n(1)*vec_n(2)*u12(i,j)
      usn(j)=vec_s(1)*vec_n(1)*u11(i,j)+                       &
             vec_s(2)*vec_n(2)*u22(i,j)+                       &
             (vec_s(1)*vec_n(2)+vec_s(2)*vec_n(1))*u12(i,j)
      tus(j)=vec_s(1)*tu1(i,j)+vec_s(2)*tu2(i,j)
      tun(j)=vec_n(1)*tu1(i,j)+vec_n(2)*tu2(i,j)
    enddo
    !
    open(18,file='Results/profile.alfa60.dat')
    write(18,"(14(1X,A15))")'d','yplus','us','t','p','uss','unn','ww', &
                                      'tke','usun','tt','tus','tun','pp'
    do j=0,jm
      write(18,"(14(1X,E15.7E3))")dis2wall(j),yplus(j),us(j),t(i,j),   &
                              p(i,j),uss(j),unn(j),u33(i,j),           &
                              tke(i,j),usn(j),tt(i,j),tus(j),        &
                              tun(j),pp(i,j)
    end do
    close(18)
    print*,' << profile.alfa60.dat ... done. '
    !
    ! budget terms
    allocate(convect(0:im,0:jm),pre_ace(0:im,0:jm),pre_dil(0:im,0:jm), &
             pre_tra(0:im,0:jm),product(0:im,0:jm),tur_tra(0:im,0:jm), &
             vis_ace(0:im,0:jm),vis_dif(0:im,0:jm),dissipa(0:im,0:jm), &
             balance(0:im,0:jm) )
    !
    call H5ReadArray(convect,im,jm,'convection','Results/tke_budget.zm.h5')
    call H5ReadArray(dissipa,im,jm,'dissipation','Results/tke_budget.zm.h5')
    call H5ReadArray(pre_ace,im,jm,'pressure_accelaration','Results/tke_budget.zm.h5')
    call H5ReadArray(pre_dil,im,jm,'pressure_dilatation','Results/tke_budget.zm.h5')
    call H5ReadArray(pre_tra,im,jm,'pressure_transport','Results/tke_budget.zm.h5')
    call H5ReadArray(product,im,jm,'production','Results/tke_budget.zm.h5')
    call H5ReadArray(tur_tra,im,jm,'turbulent_transport','Results/tke_budget.zm.h5')
    call H5ReadArray(vis_ace,im,jm,'viscous_accelaration','Results/tke_budget.zm.h5')
    call H5ReadArray(vis_dif,im,jm,'viscous_diffusion','Results/tke_budget.zm.h5')
    !
    balance=convect+pre_ace+pre_dil+pre_tra+product+tur_tra+vis_ace+vis_dif+dissipa
    !
    dissipa=dissipa+0.6d0*balance
    balance=balance*0.4d0
    !
    i=iref
    print*,' ** writing budget at i=',i
    open(18,file='Results/budget_profile.alfa60.dat')
    write(18,"(10(1X,A15))")'d','y+','Ck','Tk','Pk','Vk','Pik','Ak','epslion','balance'
    do j=0,jm
      
      miu=miucal(t(i,0))/Reynolds
      var1=ro(i,0)**2*utau(i)**4/miu
      write(18,"(10(1X,E15.7E3))")dis2wall(j),yplus(j),                &
                                 convect(i,j)/var1,                    &
                                 (tur_tra(i,j)+pre_tra(i,j))/var1,     &
                                 product(i,j)/var1,                    &
                                 vis_dif(i,j)/var1,                    &
                                 pre_dil(i,j)/var1,                    &
                                 (vis_ace(i,j)+pre_ace(i,j))/var1,     &
                                 dissipa(i,j)/var1,balance(i,j)/var1
    end do
    close(18)
    print*,' << budget_profile.alfa60.dat ... done. '
    !
    open(18,file='Results/mesh_report.'//fname//'.dat')
    write(18,"(A,F12.7,A,I4,A)")'mesh resolution at x=',x(iref,0),' i=',iref,' j=0'

    i=iref
    ds=sqrt((x(i,0)-x(i-1,0))**2+(y(i,0)-y(i-1,0))**2)
    dr=sqrt((x(i,1)-x(i,0))**2+(y(i,1)-y(i,0))**2)
    dz=(z(1)-z(0))
    !
    miu=miucal(t(i,0))/Reynolds
    lkhomo=sqrt(sqrt(miu**3/ro(i,0)**2/(-dissipa(i,0))))
    !
    write(18,"(A)")'-------------------------------------------------'
    write(18,"(A)")'          Δs+         Δr+         Δz+         Δ/η'
    write(18,"(A)")'-------------------------------------------------'
    !
    write(18,*)ds*lplus(i),dr*lplus(i),dz*lplus(i),cube_root(ds*dr*dz)/lkhomo 
    write(18,"(A)")'-------------------------------------------------'
    write(18,"(A,F12.7,A,I4,A,I0)")'mesh resolution at x=',x(iref,0),' i=',iref,' j=',jedge

    i=iref
    j=jedge
    ds=sqrt((x(i,j)-x(i-1,j))**2+(y(i,j)-y(i-1,j))**2)
    dr=sqrt((x(i,j)-x(i,j-1))**2+(y(i,j)-y(i,j-1))**2)
    dz=(z(1)-z(0))
    !
    miu=miucal(t(i,j))/Reynolds
    lkhomo=sqrt(sqrt(miu**3/ro(i,j)**2/(-dissipa(i,j))))
    !
    write(18,"(A)")'-------------------------------------------------'
    write(18,"(A)")'          Δs+         Δr+         Δz+         Δ/η'
    write(18,"(A)")'-------------------------------------------------'
    !
    write(18,*)ds*lplus(i),dr*lplus(i),dz*lplus(i),cube_root(ds*dr*dz)/lkhomo 
    write(18,"(A)")'-------------------------------------------------'
    close(18)
    print*,' << mesh_report.',fname,'.dat ... done !'
    !
    !
    ! open(18,file='tecwallvec.dat')
    ! write(18,*) 'variables=x,y,nx,ny'
    ! write(18,*) 'zone f=point,i=',im+1,'j=1'
    ! write(18,"(9(1X,E15.7E3))")(x(i,0),y(i,0),wallnorm(1,i),wallnorm(2,i),i=0,im)
    ! close(18)
    ! print*,' << tecwallvec.dat ... done'
    ! 
    !
  end subroutine nbody
  !+-------------------------------------------------------------------+
  !| The end of the subroutine nbody.                                  |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to analyse mixing layer flow.                  |
  !+-------------------------------------------------------------------+
  subroutine mixinglayer(iref)
    !
    use commvardefine,only: im,jm,km,Reynolds,pinf,mach
    use h5readwrite
    use basicfunction
    use interpolation 
    use writetec
    !
    integer,intent(in) :: iref
    !
    ! local data
    real(8),allocatable,dimension(:) :: z,theter,yc,delta_omega
    real(8),allocatable,dimension(:,:) :: x,y,ro,u1,p,t,dvdx,dudy,marms,malo
    real(8),allocatable,dimension(:,:) :: u11,u22,u33,u12,tke
    real(8),allocatable,dimension(:,:) :: convect,pre_ace,pre_dil,     &
                                          pre_tra,product,tur_tra, &
                                          vis_ace,vis_dif,dissipa, &
                                          balance
    real(8),allocatable,dimension(:) :: convect_pro,pre_ace_pro,       &
                                        pre_dil_pro,pre_tra_pro,       &
                                        product_pro,tur_tra_pro,       &
                                        vis_ace_pro,vis_dif_pro,       &
                                        dissipa_pro,balance_pro
    real(8),allocatable,dimension(:,:) :: m_iim_jj,m_ijm_ji,m_ijm_ij,  &
                              m_iim_jjm_kk,m_ijm_jim_kk,m_ijm_ijm_kk,  &
                              m_ijm_jkm_ki,m_ijm_kjm_ki
    real(8),allocatable,dimension(:) :: m_iim_jj_prof,m_ijm_ji_prof,   &
                                  m_ijm_ij_prof,m_iim_jjm_kk_prof,     &
                                  m_ijm_jim_kk_prof,m_ijm_ijm_kk_prof, &
                                  m_ijm_jkm_ki_prof,m_ijm_kjm_ki_prof
    real(8),allocatable,dimension(:,:) :: eta,unorm,vtmp
    integer,allocatable,dimension(:) :: isamp
    integer :: i,j,n,ii
    real(8) :: miu,u_1,u_2,d_u,u_c,dy,var1,var2,var3,growthrate,lkhomo
    character(len=5) :: irefname
    real(8),allocatable,dimension(:,:) :: a1111,a2222,a3333, &
                                          a111111,a222222,a333333, &
                                          a121211,a212122,a131311, &
                                          a313133,a323233,a232322, &
                                          a111122,a111133,a222211, &
                                          a222233,a333311,a333322, &
                                          a113232,a112323,a223131, &
                                          a221313,a331212,a332121
    real(8),allocatable,dimension(:) :: a1111_prof,a2222_prof,a3333_prof, &
                                        a111111_prof,a222222_prof,a333333_prof, &
                                        a121211_prof,a212122_prof,a131311_prof, &
                                        a313133_prof,a323233_prof,a232322_prof, &
                                        a111122_prof,a111133_prof,a222211_prof, &
                                        a222233_prof,a333311_prof,a333322_prof, &
                                        a113232_prof,a112323_prof,a223131_prof, &
                                        a221313_prof,a331212_prof,a332121_prof
    !
    allocate( x(0:im,0:jm),y(0:im,0:jm),ro(0:im,0:jm),u1(0:im,0:jm),   &
              p(0:im,0:jm),t(0:im,0:jm),dvdx(0:im,0:jm),               &
              dudy(0:im,0:jm),z(0:km) )
    !
    call H5ReadSubset(x,im,jm,km,'x','datin/grid.h5',kslice=0)
    call H5ReadSubset(y,im,jm,km,'y','datin/grid.h5',kslice=0)
    call H5ReadSubset(z,im,jm,km,'z','datin/grid.h5',islice=iref,jslice=0)
    !
    call H5ReadArray(ro,im,jm,'ro','Results/mean.fav.zm.h5')
    call H5ReadArray(u1,im,jm,'u1','Results/mean.fav.zm.h5')
    call H5ReadArray(p, im,jm, 'p','Results/mean.fav.zm.h5')
    call H5ReadArray(t, im,jm, 't','Results/mean.fav.zm.h5')
    call H5ReadArray(dudy,im,jm,'du1dy','Results/mean.fav.zm.h5')
    call H5ReadArray(dvdx,im,jm,'du2dx','Results/mean.fav.zm.h5')
    !
    allocate(u11(0:im,0:jm),u22(0:im,0:jm),u33(0:im,0:jm),             &
             u12(0:im,0:jm),tke(0:im,0:jm))
    !
    call H5ReadArray(u11,im,jm,'u11','Results/2order.fav.zm.h5')
    call H5ReadArray(u22,im,jm,'u22','Results/2order.fav.zm.h5')
    call H5ReadArray(u33,im,jm,'u33','Results/2order.fav.zm.h5')
    call H5ReadArray(u12,im,jm,'u12','Results/2order.fav.zm.h5')
    !
    allocate(marms(0:im,0:jm),malo(0:im,0:jm))
    do j=0,jm
    do i=0,im
      var1=sqrt(abs(u11(i,j)+u22(i,j)+u33(i,j)))
      marms(i,j)=var1/sqrt(t(i,j))*Mach
      malo(i,j)=u1(i,j)/sqrt(t(i,j))*Mach
    enddo
    enddo
    !
    u_1=0.5d0
    u_2=2.5d0
    d_u=u_2-u_1
    u_c=0.5d0*(u_1+u_2)
    !
    do j=1,jm
      if(y(0,j-1)<=0.d0 .and. y(0,j)>=0.d0) then
        print*,' ** Re_midd=',ro(0,j)/miucal(t(0,j)),'at j=',j,'y=',y(0,j)
      endif
    enddo
    !
    allocate(theter(0:im),yc(0:im),delta_omega(0:im))
    !
    do i=0,im
      !
      do j=1,jm
        !
        if(u1(i,j)>=u_c .and. u1(i,j-1)<=u_c) then
          yc(i)=linear1d(u1(i,j-1),u1(i,j),y(i,j-1),y(i,j),u_c)
          !
          exit
        endif
        !
      enddo
      !
      var1=0.d0
      do j=1,jm
        !
        dy=y(i,j)-y(i,j-1)
        var2=0.5d0*(u1(i,j-1)+u1(i,j))
        var3=0.5d0*(ro(i,j-1)+ro(i,j))
        var1=var1+var3*(u_2-var2)*(var2-u_1)*dy
        !
      enddo
      theter(i)=var1/(d_u**2)
      !
      var1=0.d0
      do j=1,jm
        var1=max(var1,dudy(i,j))
      enddo
      delta_omega(i)=d_u/var1
      !
    enddo
    !
    do i=100,im-100
      !
      var1=0.d0
      do n=-50,50
        var1=var1+delta_omega(i+n)
      enddo
      delta_omega(i)=var1/101.d0
      !
    enddo
    !
    open(18,file='Results/BLParameter.dat')
    write(18,"(4(1X,A15))")'x','yc','theter','delta_omega'
    write(18,"(4(1X,E15.7E3))")(x(i,0),yc(i),theter(i),delta_omega(i),i=0,im)
    close(18)
    print*,' << Results/BLParameter.dat'
    !
    growthrate=0.d0
    n=0
    do i=0,im
      if(x(i,0)>=350.d0 .and. x(i,0)<=450.d0) then
        var1=(theter(i+1)-theter(i-1))/(x(i+1,0)-x(i-1,0))
        growthrate=growthrate+var1
        n=n+1
      endif
    enddo
    growthrate=growthrate/dble(n)
    print*,' ** growth rate averaged between 350-450 is: ',growthrate
    print*,' ** normlized growth rate is:',growthrate/0.032d0/((u_2-u_1)/(u_2+u_1))
    !
    allocate(eta(0:im,0:jm),unorm(5,0:jm),isamp(5))
    !
    do i=0,im
      do j=0,jm
        eta(i,j)=(y(i,j)-yc(i))/theter(i)
      enddo
    enddo
    !
    do n=1,4
      !
      do i=1,im
        if(x(i,0)>=50.d0*n .and. x(i-1,0)<=50.d0*n) then
          do j=0,jm
            unorm(n,j)=(u1(i,j)-u_c)/(u_2-u_1)
          enddo
          !
          isamp(n)=i
          print*,' ** x=',x(i,0),'i=',i
          exit
          !
        endif
      enddo
      !
    enddo
    !
    n=5
    do i=1,im
      if(x(i,0)>=185.d0 .and. x(i-1,0)<=185.d0) then
      ! if(x(i,0)>=430.d0 .and. x(i-1,0)<=430.d0) then
        do j=0,jm
          eta(i,j)=(y(i,j)-yc(i))/theter(i)
          unorm(n,j)=(u1(i,j)-u_c)/(u_2-u_1)
        enddo
        !
        isamp(n)=i
        print*,' ** x=',x(i,0)
        exit
        !
      endif
    enddo
    !
    open(18,file='Results/uprofile.dat')
    write(18,"(10(1X,A15))")'y','u_100','y','u_200','y','u_300',       &
                            'y','u_400','y','u_430'
    write(18,"(10(1X,E15.7E3))")((y(isamp(n),j),u1(isamp(n),j),n=1,5),j=0,jm)
    close(18)
    print*,' << Results/uprofile.dat'
    !
    open(18,file='Results/uprofile_norm.dat')
    write(18,"(10(1X,A15))")'eta','u_100','eta','u_200','eta','u_300', &
                            'eta','u_400','eta','u_430'
    write(18,"(10(1X,E15.7E3))")((eta(isamp(n),j),unorm(n,j),n=1,5),j=0,jm)
    close(18)
    print*,' << Results/uprofile_norm.dat'
    !
    open(18,file='Results/restress_profile_norm.dat')
    write(18,"(25(1X,A15))")'eta_100','uu','vv','ww','uv',             &
                            'eta_200','uu','vv','ww','uv',             &
                            'eta_300','uu','vv','ww','uv',             &
                            'eta_400','uu','vv','ww','uv',             &
                            'eta_430','uu','vv','ww','uv'
    write(18,"(25(1X,E15.7E3))")((eta(isamp(n),j),u11(isamp(n),j),     &
         u22(isamp(n),j),u33(isamp(n),j),u12(isamp(n),j),n=1,5),j=0,jm)
    close(18)
    print*,' << Results/restress_profile_norm.dat'
    !
    !
    i=iref
    print*,' ** reference location at:',x(i,0),'i=',i
    write(irefname,'(1hi,I4.4)')i
    !
    open(18,file='Results/profile'//irefname//'.dat')
    write(18,"(8(1X,A15))")'y','u','T','p','rho','Mt','M','dudy'
    do j=0,jm
      write(18,"(8(1X,E15.7E3))")y(i,j),u1(i,j),t(i,j),p(i,j),ro(i,j), &
                                          marms(i,j),malo(i,j),dudy(i,j)
    end do
    close(18)
    print*,' << Results/profile',irefname,'.dat'
    !
    allocate(convect(0:im,0:jm),pre_ace(0:im,0:jm),pre_dil(0:im,0:jm), &
             pre_tra(0:im,0:jm),product(0:im,0:jm),tur_tra(0:im,0:jm), &
             vis_ace(0:im,0:jm),vis_dif(0:im,0:jm),dissipa(0:im,0:jm), &
             balance(0:im,0:jm) )
    allocate(convect_pro(0:jm),pre_ace_pro(0:jm),pre_dil_pro(0:jm),    &
             pre_tra_pro(0:jm),product_pro(0:jm),tur_tra_pro(0:jm),    &
             vis_ace_pro(0:jm),vis_dif_pro(0:jm),dissipa_pro(0:jm),    &
             balance_pro(0:jm))
    !
    call H5ReadArray(convect,im,jm,'convection','Results/tke_budget.zm.h5')
    call H5ReadArray(dissipa,im,jm,'dissipation','Results/tke_budget.zm.h5')
    call H5ReadArray(pre_ace,im,jm,'pressure_accelaration','Results/tke_budget.zm.h5')
    call H5ReadArray(pre_dil,im,jm,'pressure_dilatation','Results/tke_budget.zm.h5')
    call H5ReadArray(pre_tra,im,jm,'pressure_transport','Results/tke_budget.zm.h5')
    call H5ReadArray(product,im,jm,'production','Results/tke_budget.zm.h5')
    call H5ReadArray(tur_tra,im,jm,'turbulent_transport','Results/tke_budget.zm.h5')
    call H5ReadArray(vis_ace,im,jm,'viscous_accelaration','Results/tke_budget.zm.h5')
    call H5ReadArray(vis_dif,im,jm,'viscous_diffusion','Results/tke_budget.zm.h5')
    !
    balance=convect+pre_ace+pre_dil+pre_tra+product+tur_tra+vis_ace+vis_dif+dissipa
    !
    convect_pro=0.d0
    pre_ace_pro=0.d0
    pre_dil_pro=0.d0
    pre_tra_pro=0.d0
    product_pro=0.d0
    tur_tra_pro=0.d0
    vis_ace_pro=0.d0
    vis_dif_pro=0.d0
    dissipa_pro=0.d0
    balance_pro=0.d0
    !
    ! allocate(vtmp(9,0:jm))
    do ii=i-20,i+20
    !   !
    !   call regularlinearinterp(eta(ii,:),convect(ii,:),pre_ace(ii,:),     &
    !                            pre_dil(ii,:),pre_tra(ii,:),product(ii,:), &
    !                            tur_tra(ii,:),vis_ace(ii,:),vis_dif(ii,:), &
    !                            dissipa(ii,:),jm,                          &
    !                            eta(i,:),vtmp(1,:),vtmp(2,:),vtmp(3,:),    &
    !                            vtmp(4,:),vtmp(5,:),vtmp(6,:),vtmp(7,:),   &
    !                            vtmp(8,:),vtmp(9,:),jm)
    !   !
      do j=0,jm
        convect_pro(j)=convect_pro(j)+convect(ii,j)
        pre_ace_pro(j)=pre_ace_pro(j)+pre_ace(ii,j)
        pre_dil_pro(j)=pre_dil_pro(j)+pre_dil(ii,j)
        pre_tra_pro(j)=pre_tra_pro(j)+pre_tra(ii,j)
        product_pro(j)=product_pro(j)+product(ii,j)
        tur_tra_pro(j)=tur_tra_pro(j)+tur_tra(ii,j)
        vis_ace_pro(j)=vis_ace_pro(j)+vis_ace(ii,j)
        vis_dif_pro(j)=vis_dif_pro(j)+vis_dif(ii,j)
        dissipa_pro(j)=dissipa_pro(j)+dissipa(ii,j)
      enddo
    !   !
    enddo
    ! deallocate(vtmp)
    !
    convect_pro=convect_pro/21.d0
    pre_ace_pro=pre_ace_pro/21.d0
    pre_dil_pro=pre_dil_pro/21.d0
    pre_tra_pro=pre_tra_pro/21.d0
    product_pro=product_pro/21.d0
    tur_tra_pro=tur_tra_pro/21.d0
    vis_ace_pro=vis_ace_pro/21.d0
    vis_dif_pro=vis_dif_pro/21.d0
    dissipa_pro=dissipa_pro/21.d0
    !
    ! convect_pro(:)=convect(i,:)
    ! pre_ace_pro(:)=pre_ace(i,:)
    ! pre_dil_pro(:)=pre_dil(i,:)
    ! pre_tra_pro(:)=pre_tra(i,:)
    ! product_pro(:)=product(i,:)
    ! tur_tra_pro(:)=tur_tra(i,:)
    ! vis_ace_pro(:)=vis_ace(i,:)
    ! vis_dif_pro(:)=vis_dif(i,:)
    ! dissipa_pro(:)=dissipa(i,:)
    !
    balance_pro=convect_pro+pre_ace_pro+pre_dil_pro+pre_tra_pro+       &
                product_pro+tur_tra_pro+vis_ace_pro+vis_dif_pro+dissipa_pro
    !
    ! convect_pro=convect_pro-balance_pro*0.3d0
    ! balance_pro=balance_pro*0.7d0
    !
    open(18,file='Results/budget_profile.'//irefname//'.dat')
    write(18,"(A,E15.7E3)")'# x=',x(i,0)
    write(18,"(8(1X,A15))")'eta','C','T','P','V','K','epslion','balance'
    do j=0,jm
      write(18,"(8(1X,E15.7E3))")eta(i,j),                           &
                            convect_pro(j),                          &
                            tur_tra_pro(j)+pre_tra_pro(j),           &
                            product_pro(j),                          &
                            vis_dif_pro(j),                          &
                            pre_dil_pro(j)+pre_ace_pro(j)+vis_ace_pro(j), &
                            dissipa_pro(j),                          &
                            balance_pro(j)
    end do
    close(18)
    print*,' << budget_profile.',irefname,'.dat ... done. '
    ! !
    open(18,file='Results/mesh_report.'//irefname//'.dat')
    write(18,"(A,F12.7)")'# ratio of mesh to kolmogorov scale at x=',x(i,0)
    write(18,"(5(1X,A12))")'eta','delta_x/elta','delta_y/elta','delta_z/elta','delta_z/elta'
    !
    lkhomo=1000.d0
    do j=1,jm-1
      miu=miucal(1.d0)/Reynolds
      var1=sqrt(sqrt(miu**3/1.d0**2/(-dissipa(i,j))))
      var2=cube_root(0.25d0*(x(i+1,j)-x(i-1,j))*(y(i,j+1)-y(i,j-1))*   &
                                                (z(1)-z(0)))
      write(18,"(6(1X,F12.7))")0.5d0*(eta(i,j+1)+eta(i,j-1)),          &
                               0.5d0*(x(i+1,j)-x(i-1,j))/var1,         &
                               0.5d0*(y(i,j+1)-y(i,j-1))/var1,         &
                                             (z(1)-z(0))/var1,         &
                                             var2/var1,                &
                                             var1
      lkhomo=min(lkhomo,var1)
    enddo
    close(18)
    !
    print*,' ** min lkhomo=',lkhomo
    !
    print*,' << Results/mesh_report.',irefname,'.dat'
    !
    allocate(m_iim_jj(0:im,0:jm),m_ijm_ji(0:im,0:jm),            &
             m_ijm_ij(0:im,0:jm),m_iim_jjm_kk(0:im,0:jm),        &
             m_ijm_jim_kk(0:im,0:jm),m_ijm_ijm_kk(0:im,0:jm),    &
             m_ijm_jkm_ki(0:im,0:jm),m_ijm_kjm_ki(0:im,0:jm))
    allocate(a1111(0:im,0:jm),a2222(0:im,0:jm),a3333(0:im,0:jm), &
             a111111(0:im,0:jm),a222222(0:im,0:jm),a333333(0:im,0:jm), &
             a121211(0:im,0:jm),a212122(0:im,0:jm),a131311(0:im,0:jm), &
             a313133(0:im,0:jm),a323233(0:im,0:jm),a232322(0:im,0:jm), &
             a111122(0:im,0:jm),a111133(0:im,0:jm),a222211(0:im,0:jm), &
             a222233(0:im,0:jm),a333311(0:im,0:jm),a333322(0:im,0:jm), &
             a113232(0:im,0:jm),a112323(0:im,0:jm),a223131(0:im,0:jm), &
             a221313(0:im,0:jm),a331212(0:im,0:jm),a332121(0:im,0:jm) )
    !
    call H5ReadArray(    m_iim_jj,im,jm,    'm_iim_jj','Results/betchov.zm.h5')     
    call H5ReadArray(    m_ijm_ji,im,jm,    'm_ijm_ji','Results/betchov.zm.h5')      
    call H5ReadArray(    m_ijm_ij,im,jm,    'm_ijm_ij','Results/betchov.zm.h5')     
    call H5ReadArray(m_iim_jjm_kk,im,jm,'m_iim_jjm_kk','Results/betchov.zm.h5')     
    call H5ReadArray(m_ijm_jim_kk,im,jm,'m_ijm_jim_kk','Results/betchov.zm.h5')    
    call H5ReadArray(m_ijm_ijm_kk,im,jm,'m_ijm_ijm_kk','Results/betchov.zm.h5')    
    call H5ReadArray(m_ijm_jkm_ki,im,jm,'m_ijm_jkm_ki','Results/betchov.zm.h5')    
    call H5ReadArray(m_ijm_kjm_ki,im,jm,'m_ijm_kjm_ki','Results/betchov.zm.h5')
    !
    call H5ReadArray(    m_iim_jj,im,jm,    'm_iim_jj','Results/betchov.zm.h5')     
    call H5ReadArray(    m_ijm_ji,im,jm,    'm_ijm_ji','Results/betchov.zm.h5')      
    call H5ReadArray(    m_ijm_ij,im,jm,    'm_ijm_ij','Results/betchov.zm.h5')     
    call H5ReadArray(m_iim_jjm_kk,im,jm,'m_iim_jjm_kk','Results/betchov.zm.h5')     
    call H5ReadArray(m_ijm_jim_kk,im,jm,'m_ijm_jim_kk','Results/betchov.zm.h5')    
    call H5ReadArray(m_ijm_ijm_kk,im,jm,'m_ijm_ijm_kk','Results/betchov.zm.h5')    
    call H5ReadArray(m_ijm_jkm_ki,im,jm,'m_ijm_jkm_ki','Results/betchov.zm.h5')    
    call H5ReadArray(m_ijm_kjm_ki,im,jm,'m_ijm_kjm_ki','Results/betchov.zm.h5')
    !
    call H5ReadArray(a1111  ,im,jm,'A1111',  'Results/betchov_a.zm.h5')
    call H5ReadArray(a2222  ,im,jm,'A2222',  'Results/betchov_a.zm.h5')
    call H5ReadArray(a3333  ,im,jm,'A3333',  'Results/betchov_a.zm.h5')
    call H5ReadArray(a111111,im,jm,'A111111','Results/betchov_a.zm.h5')
    call H5ReadArray(a222222,im,jm,'A222222','Results/betchov_a.zm.h5')
    call H5ReadArray(a333333,im,jm,'A333333','Results/betchov_a.zm.h5')
    call H5ReadArray(a121211,im,jm,'A121211','Results/betchov_a.zm.h5')
    call H5ReadArray(a212122,im,jm,'A212122','Results/betchov_a.zm.h5')
    call H5ReadArray(a131311,im,jm,'A131311','Results/betchov_a.zm.h5')
    call H5ReadArray(a313133,im,jm,'A313133','Results/betchov_a.zm.h5')
    call H5ReadArray(a323233,im,jm,'A323233','Results/betchov_a.zm.h5')
    call H5ReadArray(a232322,im,jm,'A232322','Results/betchov_a.zm.h5')
    call H5ReadArray(a111122,im,jm,'A111122','Results/betchov_a.zm.h5')
    call H5ReadArray(a111133,im,jm,'A111133','Results/betchov_a.zm.h5')
    call H5ReadArray(a222211,im,jm,'A222211','Results/betchov_a.zm.h5')
    call H5ReadArray(a222233,im,jm,'A222233','Results/betchov_a.zm.h5')
    call H5ReadArray(a333311,im,jm,'A333311','Results/betchov_a.zm.h5')
    call H5ReadArray(a333322,im,jm,'A333322','Results/betchov_a.zm.h5')
    call H5ReadArray(a113232,im,jm,'A113232','Results/betchov_a.zm.h5')
    call H5ReadArray(a112323,im,jm,'A112323','Results/betchov_a.zm.h5')
    call H5ReadArray(a223131,im,jm,'A223131','Results/betchov_a.zm.h5')
    call H5ReadArray(a221313,im,jm,'A221313','Results/betchov_a.zm.h5')
    call H5ReadArray(a331212,im,jm,'A331212','Results/betchov_a.zm.h5')
    call H5ReadArray(a332121,im,jm,'A332121','Results/betchov_a.zm.h5')
    !
    allocate(m_iim_jj_prof(0:jm),    m_ijm_ji_prof(0:jm),         &
             m_ijm_ij_prof(0:jm),    m_iim_jjm_kk_prof(0:jm),     &
             m_ijm_jim_kk_prof(0:jm),m_ijm_ijm_kk_prof(0:jm),     &
             m_ijm_jkm_ki_prof(0:jm),m_ijm_kjm_ki_prof(0:jm) )
    allocate(a1111_prof(0:jm),a2222_prof(0:jm),a3333_prof(0:jm),       &
             a111111_prof(0:jm),a222222_prof(0:jm),a333333_prof(0:jm), &
             a121211_prof(0:jm),a212122_prof(0:jm),a131311_prof(0:jm), &
             a313133_prof(0:jm),a323233_prof(0:jm),a232322_prof(0:jm), &
             a111122_prof(0:jm),a111133_prof(0:jm),a222211_prof(0:jm), &
             a222233_prof(0:jm),a333311_prof(0:jm),a333322_prof(0:jm), &
             a113232_prof(0:jm),a112323_prof(0:jm),a223131_prof(0:jm), &
             a221313_prof(0:jm),a331212_prof(0:jm),a332121_prof(0:jm) )
    !
        m_iim_jj_prof=0.d0
        m_ijm_ji_prof=0.d0
        m_ijm_ij_prof=0.d0
    m_iim_jjm_kk_prof=0.d0
    m_ijm_jim_kk_prof=0.d0
    m_ijm_ijm_kk_prof=0.d0
    m_ijm_jkm_ki_prof=0.d0
    m_ijm_kjm_ki_prof=0.d0
    !
      a1111_prof=0.d0
      a2222_prof=0.d0
      a3333_prof=0.d0
    a111111_prof=0.d0
    a222222_prof=0.d0
    a333333_prof=0.d0
    a121211_prof=0.d0
    a212122_prof=0.d0
    a131311_prof=0.d0
    a313133_prof=0.d0
    a323233_prof=0.d0
    a232322_prof=0.d0
    a111122_prof=0.d0
    a111133_prof=0.d0
    a222211_prof=0.d0
    a222233_prof=0.d0
    a333311_prof=0.d0
    a333322_prof=0.d0
    a113232_prof=0.d0
    a112323_prof=0.d0
    a223131_prof=0.d0
    a221313_prof=0.d0
    a331212_prof=0.d0
    a332121_prof=0.d0
    !
    allocate(vtmp(8,0:jm))
    do ii=iref-100,iref+100
      !
      ! call regularlinearinterp(eta(ii,:),m_iim_jj(ii,:),m_ijm_ji(ii,:),  &
      !             m_ijm_ij(ii,:),m_iim_jjm_kk(ii,:),m_ijm_jim_kk(ii,:),  &
      !             m_ijm_ijm_kk(ii,:),m_ijm_jkm_ki(ii,:),                 &
      !             m_ijm_kjm_ki(ii,:),jm,                                 &
      !             eta(i,:),vtmp(1,:),vtmp(2,:),vtmp(3,:),              &
      !             vtmp(4,:),vtmp(5,:),vtmp(6,:),vtmp(7,:),             &
      !             vtmp(8,:),jm)
      !
      do j=0,jm
            m_iim_jj_prof(j)=    m_iim_jj_prof(j)+    m_iim_jj(ii,j)
            m_ijm_ji_prof(j)=    m_ijm_ji_prof(j)+    m_ijm_ji(ii,j)
            m_ijm_ij_prof(j)=    m_ijm_ij_prof(j)+    m_ijm_ij(ii,j)
        m_iim_jjm_kk_prof(j)=m_iim_jjm_kk_prof(j)+m_iim_jjm_kk(ii,j)
        m_ijm_jim_kk_prof(j)=m_ijm_jim_kk_prof(j)+m_ijm_jim_kk(ii,j)
        m_ijm_ijm_kk_prof(j)=m_ijm_ijm_kk_prof(j)+m_ijm_ijm_kk(ii,j)
        m_ijm_jkm_ki_prof(j)=m_ijm_jkm_ki_prof(j)+m_ijm_jkm_ki(ii,j)
        m_ijm_kjm_ki_prof(j)=m_ijm_kjm_ki_prof(j)+m_ijm_kjm_ki(ii,j)
      enddo
      !
      ! call regularlinearinterp(eta(ii,:),a1111(ii,:),a2222(ii,:),  &
      !             a3333(ii,:),a111111(ii,:),a222222(ii,:),         &
      !             a333333(ii,:),a121211(ii,:),a212122(ii,:),jm,    &
      !             eta(i,:),vtmp(1,:),vtmp(2,:),vtmp(3,:),              &
      !             vtmp(4,:),vtmp(5,:),vtmp(6,:),vtmp(7,:),             &
      !             vtmp(8,:),jm)
      !
      do j=0,jm
          a1111_prof(j)=  a1111_prof(j)+  a1111(ii,j)
          a2222_prof(j)=  a2222_prof(j)+  a2222(ii,j)
          a3333_prof(j)=  a3333_prof(j)+  a3333(ii,j)
        a111111_prof(j)=a111111_prof(j)+a111111(ii,j)
        a222222_prof(j)=a222222_prof(j)+a222222(ii,j)
        a333333_prof(j)=a333333_prof(j)+a333333(ii,j)
        a121211_prof(j)=a121211_prof(j)+a121211(ii,j)
        a212122_prof(j)=a212122_prof(j)+a212122(ii,j)
      enddo
      !
      ! call regularlinearinterp(eta(ii,:),a131311(ii,:),a313133(ii,:),  &
      !             a323233(ii,:),a232322(ii,:),a111122(ii,:),           &
      !             a111133(ii,:),a222211(ii,:),a222233(ii,:),jm,        &
      !             eta(i,:),vtmp(1,:),vtmp(2,:),vtmp(3,:),              &
      !             vtmp(4,:),vtmp(5,:),vtmp(6,:),vtmp(7,:),             &
      !             vtmp(8,:),jm)
      !
      do j=0,jm
        a131311_prof(j)=a131311_prof(j)+a131311(ii,j)
        a313133_prof(j)=a313133_prof(j)+a313133(ii,j)
        a323233_prof(j)=a323233_prof(j)+a323233(ii,j)
        a232322_prof(j)=a232322_prof(j)+a232322(ii,j)
        a111122_prof(j)=a111122_prof(j)+a111122(ii,j)
        a111133_prof(j)=a111133_prof(j)+a111133(ii,j)
        a222211_prof(j)=a222211_prof(j)+a222211(ii,j)
        a222233_prof(j)=a222233_prof(j)+a222233(ii,j)
      enddo
      !
      ! call regularlinearinterp(eta(ii,:),a333311(ii,:),a333322(ii,:),  &
      !             a113232(ii,:),a112323(ii,:),a223131(ii,:),           &
      !             a221313(ii,:),a331212(ii,:),a332121(ii,:),jm,        &
      !             eta(i,:),vtmp(1,:),vtmp(2,:),vtmp(3,:),              &
      !             vtmp(4,:),vtmp(5,:),vtmp(6,:),vtmp(7,:),             &
      !             vtmp(8,:),jm)
      !
      do j=0,jm
        a333311_prof(j)=a333311_prof(j)+a333311(ii,j)
        a333322_prof(j)=a333322_prof(j)+a333322(ii,j)
        a113232_prof(j)=a113232_prof(j)+a113232(ii,j)
        a112323_prof(j)=a112323_prof(j)+a112323(ii,j)
        a223131_prof(j)=a223131_prof(j)+a223131(ii,j)
        a221313_prof(j)=a221313_prof(j)+a221313(ii,j)
        a331212_prof(j)=a331212_prof(j)+a331212(ii,j)
        a332121_prof(j)=a332121_prof(j)+a332121(ii,j)
      enddo
      !
    enddo
    deallocate(vtmp)
    !
        m_iim_jj_prof=    m_iim_jj_prof/201.d0
        m_ijm_ji_prof=    m_ijm_ji_prof/201.d0
        m_ijm_ij_prof=    m_ijm_ij_prof/201.d0
    m_iim_jjm_kk_prof=m_iim_jjm_kk_prof/201.d0
    m_ijm_jim_kk_prof=m_ijm_jim_kk_prof/201.d0
    m_ijm_ijm_kk_prof=m_ijm_ijm_kk_prof/201.d0
    m_ijm_jkm_ki_prof=m_ijm_jkm_ki_prof/201.d0
    m_ijm_kjm_ki_prof=m_ijm_kjm_ki_prof/201.d0
    !
      a1111_prof=  a1111_prof/201.d0
      a2222_prof=  a2222_prof/201.d0
      a3333_prof=  a3333_prof/201.d0
    a111111_prof=a111111_prof/201.d0
    a222222_prof=a222222_prof/201.d0
    a333333_prof=a333333_prof/201.d0
    a121211_prof=a121211_prof/201.d0
    a212122_prof=a212122_prof/201.d0
    !
    a131311_prof=a131311_prof/201.d0
    a313133_prof=a313133_prof/201.d0
    a323233_prof=a323233_prof/201.d0
    a232322_prof=a232322_prof/201.d0
    a111122_prof=a111122_prof/201.d0
    a111133_prof=a111133_prof/201.d0
    a222211_prof=a222211_prof/201.d0
    a222233_prof=a222233_prof/201.d0
    !
    a333311_prof=a333311_prof/201.d0
    a333322_prof=a333322_prof/201.d0
    a113232_prof=a113232_prof/201.d0
    a112323_prof=a112323_prof/201.d0
    a223131_prof=a223131_prof/201.d0
    a221313_prof=a221313_prof/201.d0
    a331212_prof=a331212_prof/201.d0
    a332121_prof=a332121_prof/201.d0
    !
    open(18,file='m_Betchov'//irefname//'.dat')
    write(18,"(9(1X,A15))")'y','m_iim_jj','m_ijm_ji','m_ijm_ij',       &
                        'm_iim_jjm_kk','m_ijm_jim_kk','m_ijm_ijm_kk',  &
                        'm_ijm_jkm_ki','m_ijm_kjm_ki'
    write(18,"(9(1X,E15.7E3))")(eta(i,j),m_iim_jj_prof(j),m_ijm_ji_prof(j), &
                           m_ijm_ij_prof(j),m_iim_jjm_kk_prof(j),           &
                           m_ijm_jim_kk_prof(j),m_ijm_ijm_kk_prof(j),       &
                           m_ijm_jkm_ki_prof(j),m_ijm_kjm_ki_prof(j),j=0,jm)
    close(18)
    print*,' << m_Betchov',irefname,'.dat ... done !'
    !
    open(18,file='A_Betchov'//irefname//'.dat')
    write(18,"(25(1X,A15))")'y','A1111','A2222','A3333','A111111',      &
              'A222222','A333333','A121211','A212122','A131311',       &
              'A313133','A323233','A232322','A111122','A111133',       &
              'A222211','A222233','A333311','A333322','A113232',       &
              'A112323','A223131','A221313','A331212','A332121'
    write(18,"(25(1X,E15.7E3))")(eta(i,j),a1111_prof(j),a2222_prof(j),  &
       a3333_prof(j),a111111_prof(j),a222222_prof(j),a333333_prof(j),  &
       a121211_prof(j),a212122_prof(j),a131311_prof(j),a313133_prof(j), &
       a323233_prof(j),a232322_prof(j),a111122_prof(j),a111133_prof(j), &
       a222211_prof(j),a222233_prof(j),a333311_prof(j),a333322_prof(j), &
       a113232_prof(j),a112323_prof(j),a223131_prof(j),a221313_prof(j), &
       a331212_prof(j),a332121_prof(j),j=0,jm)
    close(18)
    print*,' << m_Betchov',irefname,'.dat ... done !'
    !
  end subroutine mixinglayer
  !+-------------------------------------------------------------------+
  !| The end of the subroutine mixinglayer.                            |
  !+-------------------------------------------------------------------+
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to cauculate the post-shock flow paramters,
  ! by using pre-shock flow paramters and shock angle and invisid R-H
  ! Relation.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2010-10-23.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine PostShockCal(ropre,upre,vpre,ppre,Tpre,ropos,upos,vpos,   &
                                              ppos,Tpos,ang,mach,lshock)
    !
    use commvardefine, only: const2,gamma,pi
    !
    real(8) :: ropre,upre,vpre,ppre,Tpre,ropos,upos,vpos,ppos,Tpos,ang,mach
    logical,intent(in),optional :: lshock
    real(8) :: CssPre,uNpre,uTPre,MachNPre,uNpos,uTPos,SRatio,ang2
    logical :: lc
    !
    if(present(lshock)) then
        lc=lshock
    else
        lc=.true.
    endif
    !
    if(lc) then
      ang2=ang/180.d0*pi
      !
      CssPre=dsqrt(Tpre)/Mach
      !
      uNpre=upre*dsin(ang2)-vpre*dcos(ang2)
      uTPre=upre*dcos(ang2)+vpre*dsin(ang2)
      MachNPre=uNpre/CssPre
      ! print*,'MachNPre=',MachNPre
      ! print*,'uNpre=',uNpre
      ! print*,'CssPre=',CssPre
      SRatio=2.d0*gamma/(gamma+1.d0)*MachNPre**2-(gamma-1.d0)/(gamma+1.d0)
      !
      ppos=ppre*SRatio
      ropos=(SRatio*(gamma+1.d0)+gamma-1.d0)/                            &
                                    (gamma+1.d0+SRatio*(gamma-1.d0))*ropre
      Tpos=ppos/ropos*const2
      !
      uNpos=ropre*uNpre/ropos
      uTpos=uTpre
      upos=uNpos*dsin(ang2)+uTpos*dcos(ang2)
      vpos=-uNpos*dcos(ang2)+uTpos*dsin(ang2)
      !
    else
      stop
    endif
    !
  end subroutine PostShockCal
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine PostShockCal.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !+-------------------------------------------------------------------+
  !| This Prandtl–Meyer function describes the angle through which a   |
  !| flow turns isentropically from sonic velocity (M=1) to a Mach (M) |
  !| number greater than 1. The maximum angle through which a sonic    |
  !| (M = 1) flow .                                                    |
  !+-------------------------------------------------------------------+
  real(8) function PrandtlMeyer(ma)
    !
    use commvardefine, only: gamma
    !
    real(8),intent(in) :: ma
    real(8) :: var1,var2,var3
    !
    var1=sqrt((gamma+1.d0)/(gamma-1.d0))
    var2=sqrt((gamma-1.d0)/(gamma+1.d0)*(Ma**2-1.d0))
    var3=sqrt(Ma**2-1.d0)
    PrandtlMeyer=var1*atan(var2)-atan(var3)
    !
  end function PrandtlMeyer
  !+-------------------------------------------------------------------+
  !| The end of the subroutine PrandtlMeyer.                           |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function returns the local streamline curvature according to |
  !! local velocity and velocity gradient.                             |
  !+-------------------------------------------------------------------+
  real(8) function curvature(u,v,dudx,dudy,dvdx,dvdy)
    !
    real(8),intent(in) :: u,v,dudx,dudy,dvdx,dvdy
    !
    real(8),parameter :: eps=1.d-15
    real(8) :: slop,dslop(2),var1,var2,dx,dy
    !
    slop=v/(u+eps)
    !
    ! the derivative of slop
    dslop(1)=-v/(u**2+eps)*dudx+1.d0/(u+eps)*dvdx
    dslop(2)=-v/(u**2+eps)*dudy+1.d0/(u+eps)*dvdy
    !
    var1=sqrt(u**2+v**2)
    dx=u/var1
    dy=v/var1
    !
    ! project to the s direction
    var2=dslop(1)*dx+dslop(2)*dy
    var2=var2/dx ! derivative in the x direction
    !
    curvature=var2/sqrt((1.d0+slop**2)**3)
    !
    !
    if(abs(curvature)>1.d9 .or. isnan(curvature)) then
      curvature=0.d0
    endif
    
    !
    ! if(abs(var2)>1.d-15 .and. abs(var2)<1.d15) then
    !   curvature=var2/sqrt((1.d0+slop**2)**3)
    ! else
    !   print*,slop,dslop(1),dslop(2)
    ! endif
    !
    return
    !
  end function curvature
  !+-------------------------------------------------------------------+
  !| The end of the subroutine curvature.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| this subroutine is used to integrate velocity field to get        |
  !| streamline.                                                       |
  !+-------------------------------------------------------------------+
  subroutine streamlinecal(iref,jref)
    !
    use interpolation
    !
    integer,intent(in) :: iref,jref
    real(8),allocatable,dimension(:,:) :: x,y,u1,u2,dudx,dvdx,dudy,    &
                                          dvdy,cur,rad,udv
    real(8),allocatable,dimension(:,:,:) :: du1,du2,uddv
    real(8),allocatable,dimension(:) :: xs,ys,us,vs
    integer :: is,js,m,i,j,n
    real(8) :: dx,dy,ds,xmax,ymax,xmin,ymin,var1,var2,eps
    !
    allocate(x(0:im,0:jm),y(0:im,0:jm))
    call H5ReadSubset(x,im,jm,km,'x','datin/grid.h5',kslice=0)
    call H5ReadSubset(y,im,jm,km,'y','datin/grid.h5',kslice=0)
    ! do j=0,jm
    ! do i=0,im
    !   x(i,j)=10.d0/im*i
    !   y(i,j)=10.d0/jm*j
    ! enddo
    ! enddo
    !
    xmax=maxval(x); ymax=maxval(y)
    xmin=minval(x); ymin=minval(y)
    !
    allocate(u1(0:im,0:jm),u2(0:im,0:jm))
    call H5ReadArray(u1,im,jm,'u1','Results/mean.fav.zm.h5')
    call H5ReadArray(u2,im,jm,'u2','Results/mean.fav.zm.h5')
    ! do j=0,jm
    ! do i=0,im
    !   var1=(x(i,j)-5.d0)**2+(y(i,j)-5.d0)**2
    !   u1(i,j)=-1.d0*(y(i,j)-5.d0)/1.d0*exp(-0.5d0*var1)
    !   u2(i,j)= 1.d0*(x(i,j)-5.d0)/1.d0*exp(-0.5d0*var1)
    ! enddo
    ! enddo
    !
    allocate(du1(1:2,0:im,0:jm),du2(1:2,0:im,0:jm))
    allocate(cur(0:im,0:jm),rad(0:im,0:jm))
    !
    call H5ReadArray(du1(1,:,:),im,jm,'du1dx','Results/mean.fav.zm.h5')
    call H5ReadArray(du2(1,:,:),im,jm,'du2dx','Results/mean.fav.zm.h5')
    call H5ReadArray(du1(2,:,:),im,jm,'du1dy','Results/mean.fav.zm.h5')
    call H5ReadArray(du2(2,:,:),im,jm,'du2dy','Results/mean.fav.zm.h5')
    ! du1=grad_xy(u1,x,y)
    ! du2=grad_xy(u2)
    !
    do j=0,jm
    do i=0,im
      cur(i,j)=curvature(u1(i,j),u2(i,j),du1(1,i,j),du1(2,i,j),du2(1,i,j),du2(2,i,j))
    enddo
    enddo
    !
    call writetecbin('Results/tecuv.plt',x,'x',y,'y',u1,'u',u2,'v',    &
                                                  cur,'curvature',im,jm) 
    !
    ! m=im*5
    ! allocate(xs(0:m),ys(0:m),us(0:m),vs(0:m))
    ! !
    ! is=iref
    ! js=jref
    ! !
    ! xs=0.d0
    ! ys=0.d0
    ! !
    ! xs(0)=x(is,js)
    ! ys(0)=y(is,js)
    ! us(0)=u1(is,js)
    ! vs(0)=u2(is,js)
    ! !
    ! ds=sqrt((x(is+1,js+1)-x(is,js))**2+(y(is+1,js+1)-y(is,js))**2)
    ! ds=0.2d0*ds
    ! !
    ! do n=1,m
    !   !
    !   dx=ds*us(n-1)/sqrt(us(n-1)**2+vs(n-1)**2)
    !   dy=dx*vs(n-1)/sqrt(us(n-1)**2+vs(n-1)**2)
    !   !
    !   xs(n)=xs(n-1)+dx
    !   ys(n)=ys(n-1)+dy
    !   !
    !   if(xs(n)>=xmax .or. xs(n)<=xmin .or.                             &
    !      ys(n)>=ymax .or. ys(n)<=ymin ) exit
    !   !
    !   call pointinterp(x,y,u1,u2,im,jm,  xs(n),ys(n),us(n),vs(n))
    !   !
    !   if(mod(n,50)==0)                                                 &
    !   write(*,'(1A1,A,I0,A,$)')char(13),'  ** calcualte the ',n,'th point'
    !   !
    ! enddo
    ! !
    ! open(18,file='Results/tecstreamline.dat')
    ! write(18,*) 'variables=x,y,u,v'
    ! write(18,*) 'zone f=point,i=',m+1
    ! write(18,"(4(1X,E15.7E3))")(xs(i),ys(i),us(i),vs(i),i=0,m)
    ! close(18)
    ! print*,' << Results/tecstreamline.dat ... done !'
    !
    ! eps=1.d-15
    ! ! eps=0.d0
    ! !
    ! allocate( dudx(0:im,0:jm),dvdx(0:im,0:jm),                         &
    !           dudy(0:im,0:jm),dvdy(0:im,0:jm),udv(0:im,0:jm),          &
    !           uddv(1:2,0:im,0:jm))
    ! allocate(du1(1:2,0:im,0:jm),du2(1:2,0:im,0:jm))
    ! call H5ReadArray(du1(1,:,:),im,jm,'du1dx','Results/mean.fav.zm.h5')
    ! call H5ReadArray(du2(1,:,:),im,jm,'du2dx','Results/mean.fav.zm.h5')
    ! call H5ReadArray(du1(2,:,:),im,jm,'du1dy','Results/mean.fav.zm.h5')
    ! call H5ReadArray(du2(2,:,:),im,jm,'du2dy','Results/mean.fav.zm.h5')
    ! !
    ! udv=u2/(u1+eps)
    ! !
    ! ! du1=grad_xy(u1,x,y)
    ! ! du2=grad_xy(u2)
    ! !
    ! uddv(1,:,:)=-u2/(u1**2+eps)*du1(1,:,:)+1.d0/(u1+eps)*du2(1,:,:)
    ! uddv(2,:,:)=-u2/(u1**2+eps)*du1(2,:,:)+1.d0/(u1+eps)*du2(2,:,:)

    ! ! du2=grad_xy(u2)
    ! ! dudx=du1(1,:,:)
    ! ! dvdx=du2(1,:,:)
    ! !
    ! !
    ! !
    ! allocate(cur(0:im,0:jm),rad(0:im,0:jm))
    ! cur=0.d0
    ! rad=0.d0
    ! do j=1,jm
    !   do i=0,im
    !     var1=sqrt(u1(i,j)**2+u2(i,j)**2)
    !     dx=u1(i,j)/var1
    !     dy=u2(i,j)/var1
    !     !
    !     var1=dy/(dx+eps)
    !     !
    !     var2=uddv(1,i,j)+uddv(2,i,j)*var1
    !     ! var2=-u2(i,j)/(u1(i,j)**2+eps)*dudx(i,j)+1.d0/(u1(i,j)+eps)*dvdx(i,j)
    !     !
    !     if(abs(var2)>1.d-15 .and. abs(var2)<1.d15) then
    !       cur(i,j)=var2/sqrt((1.d0+var1**2)**3)
    !       rad(i,j)=1.d0/cur(i,j) !sqrt((x(i,j)-5.d0)**2+(y(i,j)-5.d0)**2) !sqrt((1.d0+var1**2)**3)/var2
    !     else
    !       print*,var1,var2,uddv(1,i,j),uddv(2,i,j)
    !     endif
    !     !
    !   enddo
    ! enddo
    !
  end subroutine streamlinecal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine streamlinecal.                          |
  !+-------------------------------------------------------------------+
  !
  subroutine MVG3D
    !
    use commvardefine,only: im,jm,km,x,y,z,Reynolds,Mach,lihomo,lkhomo
    !
    real(8),allocatable,dimension(:,:,:) :: u1,u2,u3,u1f
    integer :: i,j,k
    !
    allocate(x(0:im,0:jm,0:km),  y(0:im,0:jm,0:km), z(0:im,0:jm,0:km), &
             u1(0:im,0:jm,0:km),u2(0:im,0:jm,0:km),u3(0:im,0:jm,0:km), &
             u1f(0:im,0:jm,0:km))
    call H5ReadArray(x,im,jm,km,'x','datin/grid.h5')
    call H5ReadArray(y,im,jm,km,'y','datin/grid.h5')
    call H5ReadArray(z,im,jm,km,'z','datin/grid.h5')
    !
    call H5ReadArray(u1,im,jm,km,'u1','Results/mean.fav.h5')
    call H5ReadArray(u2,im,jm,km,'u2','Results/mean.fav.h5')
    call H5ReadArray(u3,im,jm,km,'u3','Results/mean.fav.h5')
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      u1f(i,j,k)=u1(i,j,k)-u1(0,j,k)
    enddo
    enddo
    enddo
    !
    call writetecbin('Results/tecmeanflow.plt',x,'x',y,'y',z,'z',u1,'u',u2,'v',u3,'w',u1f,'u"',im,jm,km) 
    !
  end subroutine MVG3D
  !
  subroutine fxtest
    !
    use basicfunction, only: fxsolver
    !
    real(8) :: x,f
    !
    x=fxsolver(xmin=0.d0,xmax=3.d0,ftarget=0d0,fx=xmcosx)
    !
  end subroutine fxtest
  !
  function square_root2(x)
    real(8) :: x,square_root2
    square_root2=sqrt(x)
  end function square_root2
  function xmcosx(x)
    real(8) :: x,xmcosx
    xmcosx=x-cos(x)
  end function xmcosx
  !
end module flowanalyse