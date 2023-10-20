module chem
  !
#ifdef COMB
  use cantera
#endif 
  !
  implicit none
  !
  integer :: num_species
  !
#ifdef COMB
  !
  integer :: ncstep,nbody=0,ngibb=0,nlind=0,ntroe=0,nsrif=0
  real(8),parameter :: rguniv=8.3142d3
  real(8),parameter :: alamdc=2.58d-5,rlamda=7.0d-1,tlamda=2.98d2
  real(8) :: prefgb,alamda
  !
  character(len=10),allocatable :: spcsym(:),bdysym(:)
  character(len=256) :: chemxmlfile
  integer,allocatable :: &
    ntint(:),ncofcp(:,:),nsslen(:),nsspec(:,:),nrslen(:),nrspec(:,:),npslen(:) &
    ,npspec(:,:),nrclen(:),nrcpec(:,:),npclen(:),npcpec(:,:),mblist(:)         &
    ,mglist(:),mllist(:),mtlist(:),mslist(:),ncpoly(:,:),ncpom1(:,:)           &
    ,ncenth(:,:),ncenpy(:,:)
  real(8),allocatable :: &
    wmolar(:),clewis(:),tintlo(:,:),tinthi(:,:),amolcp(:,:,:),Arrhenius(:,:)   &
    ,crspec(:,:),cpspec(:,:),diffmu(:,:),effy3b(:,:),rclind(:,:),rctroe(:,:)   &
    ,rcsrif(:,:),diffmw(:,:),ovwmol(:),rgspec(:),amascp(:,:,:),amascv(:,:,:)   &
    ,amasct(:,:,:),amasch(:,:,:),amasce(:,:,:),amascs(:,:,:),amolgb(:,:,:)     &
    ,olewis(:),wirate(:)
    !
  type(phase_t) :: mixture
  !
#endif
  !
  character(len=5) :: tranmod='mixav'
  logical :: ctrflag=.true.
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine reads the chemistry data from cantera format file |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 23-Feb-2021  | Created by Z.X. Chen @ Cambridge                   |
  !+-------------------------------------------------------------------+
  subroutine chemread(cheminfile)
    !
    ! arguments
    character(len=*),intent(in) :: cheminfile
    !
#ifdef COMB
    ! local data
    integer :: &
      js,is,it,icp,jr,ir,j3b,i3b,ilind,jlind,itroe,jtroe,isrif,jsrif,  &
      ireac,iprod,icol,ispa,jreac,jprod
    real(8), allocatable :: &
      rcoeff(:,:),pcoeff(:,:),Arrhenius0(:,:),troeparams(:,:), &
      sriparams(:,:)
    real(8) :: dum
    character(len=10) :: spcnm,eunit,reactype,reverse,fallofftype,phase_id
    character(len=100) :: stringline
    !
    phase_id='gas'
    !---CANTERA---
    mixture=importPhase(cheminfile,trim(phase_id))
    if(speciesIndex(mixture,'O')<0 .and. &
    speciesIndex(mixture,'o')<0) &
      stop '!! Species "O" not exist, check phase_id !!'
    ncstep=nReactions(mixture)
    num_species=nSpecies(mixture)
    prefgb=refPressure(mixture)
    call setPressure(mixture,prefgb)
    chemxmlfile=cheminfile
    !
    print*,' ** num_species:',num_species
    !
    !===============ALLOCATE=============================
    allocate(spcsym(num_species),wmolar(num_species),                  &
             clewis(num_species),ntint(num_species),                   &
             ovwmol(num_species),rgspec(num_species))
    !
    !---CANTERA---
    call getMolecularWeights(mixture,wmolar)
    !
    write(*,'(2X,A12)')'------------'
    write(*,'(2X,A12)')'   chemistry'
    write(*,'(2X,A12)')'----+-------'
    do js=1,num_species
      call getSpeciesName(mixture,js,spcsym(js))
      
      write(*,'(2X,I3,A3,A5)')js,' | ',spcsym(js)
    enddo
    write(*,'(2X,A12)')'----+-------'
    !
    clewis(:)=1.d0
    !
    !RECIPROCAL OF MOLAR MASS
    ovwmol(:)=1.0d0/wmolar(:)
    !SPECIFIC mixture CONSTANT
    rgspec(:)=rguniv*ovwmol(:)
    !
#endif
    !
  end subroutine chemread
  !+-------------------------------------------------------------------+
  !| The end of the subroutine chemread_ctr                            |
  !+-------------------------------------------------------------------+
  !+-------------------------------------------------------------------+
  !| This subroutine prints the chemistry data for display.            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-Aug-2020  | Created by Z.X. Chen @ Cambridge                   |
  !+-------------------------------------------------------------------+
  subroutine chemrep(filename)
    ! arguments
    character(len=*),intent(in) :: filename
    !
#ifdef COMB
    ! local data
    integer :: js,icp,jcp,kcp,jr,istr1,istr2,ks
    ! real(8) :: ttemp(5),fornow,ttold(5)
    ! character*132 char132
    ! character*10 char10
    ! character*5 char5
    ! character*4 char4
    ! character*1 char1
    ! !
    open(16,file=filename)
    !
    ! SPECIES LIST, 43 char length per line
    write(16,'(A)')' +------------ Chemical Data -------------+'
    write(16,'(A,I5,A18)')' Number of species:',num_species,''
    write(16,'(A7,3X,A7,3X,A8,4X,A9,A2)') &
    ' Index','Species','Mol.Mass','Lewis No.',''
    do js=1,num_species
      write(16,'(A2,I5,3X,A7,3X,1PE9.3,3X,1PE9.3,A2)') &
      ' ',js,spcsym(js),wmolar(js),clewis(js),''
    enddo
    !
    close(16)
    !
    print*,' << ',filename
    !
#endif
    !
  end subroutine chemrep
  !+-------------------------------------------------------------------+
  !| The end of the function chemrep.                                  |
  !+-------------------------------------------------------------------+
  !!+-------------------------------------------------------------------+
  !| This subroutine is used to compute thermodynamic quantities       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16-Aug-2020: Created by Z.X. Chen @ Cambridge                     |
  !+-------------------------------------------------------------------+
  subroutine thermdyn
    !
#ifdef COMB
    !
    ! local data
    integer :: is,it,itnm1,jc,icp
    !
    real(8) :: &
      tbreak,twidth,ovtwid,ovtwrp,ttemp1,cpmol1,ttemp2,cpmol2,deltcp   &
      ,giblet
    !
    !===============ALLOCATE=============================
    allocate(  olewis(num_species),wirate(num_species)   )
    !
    !RECIPROCAL OF LEWIS NUMBER
    do is=1,num_species
      olewis(is)=1.0d0/clewis(is)
    enddo
    !
    !CONDUCTIVITY COEFFICIENT
    alamda=alamdc*exp(-rlamda*log(tlamda))
    !
#endif
    !
  end subroutine thermdyn
  !+-------------------------------------------------------------------+
  !| The end of the subroutine cheminit.                               |
  !+-------------------------------------------------------------------+
  !+-------------------------------------------------------------------+
  !| This subroutine computes heat release rate.                       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 19-Nov-2020  | Created by Z.X. Chen @ Cambridge                   |
  !+-------------------------------------------------------------------+
  real(8) function heatrate(den,tmp,spc)
    !
    ! arguments
    real(8),intent(in) :: den,tmp,spc(:)
    !
#ifdef COMB
    ! local data
    real(8) :: hi(num_species)
    !
    call enthpy(tmp,hi)
    call chemrate(den,tmp,spc(:))
    !
    heatrate=-1.d0*sum(wirate(:)*hi(:))
    !
#endif
    !
  end function heatrate
  !+-------------------------------------------------------------------+
  !| The end of the subroutine aceval.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine computes local species enthalpy.                  |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-Oct-2020  | Created by Z.X. Chen @ Cambridge                   |
  !+-------------------------------------------------------------------+
  subroutine enthpy(tmp,hi)
    !
    ! arguments
    real(8),intent(in) :: tmp
    real(8),intent(out) :: hi(:)
    !
#ifdef COMB
    ! local data
    integer :: is,it,jt,icp
    real(8) :: fornow
    !
    if(.true. .and. ctrflag) then
      !
      call setTemperature(mixture,tmp)
      call getEnthalpies_RT(mixture,hi(:))
      hi(:)=hi(:)*rgspec(:)*tmp
      !
    else
      !
      it=1
      do is=1,num_species
        !
        if(is==1) then
          do jt=1,ntint(is)
            if(tmp>tinthi(ntint(is),is)) then
              print*,' !! ENTH Error - Temperature out of bound!!, T =',tmp
            endif
            if(tmp>tinthi(jt,is)) it=it+1
          enddo
        endif
        !
        fornow=amasch(ncpoly(it,is),it,is)
        do icp=ncpom1(it,is),1,-1
          fornow=amasch(icp,it,is)+fornow*tmp
        enddo
        hi(is)=amasch(ncenth(it,is),it,is)+fornow*tmp
        !
      enddo
      !
    endif
    !
#endif
    !
  end subroutine enthpy
  !+-------------------------------------------------------------------+
  !| The end of the subroutine enthpy.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine computes chemcal reaction rates                   |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 1-Nov-2020  | Created by Z.X. Chen @ Cambridge                   |
  !+-------------------------------------------------------------------+
  subroutine chemrate(den,tmp,spc,wi)
    !
    ! arguments
    real(8),intent(in) :: den,tmp,spc(:)
    real(8),optional,intent(out) :: wi(:)
    !
#ifdef COMB
    ! local data
    integer :: js,jr,j3b,jsspec,it,jt,icp
    real(8) :: &
      fornow,tbconc,rfwd,rfld,rfln,rbln,gibbs,rbwd,stoi,p_rdc,ftc, &
      const1,const2,rfsr,rftr,wi_molar(num_species),prs
    logical :: flag3by
    !
    if(ctrflag) then
      call setState_TRY(mixture,tmp,den,spc(:))
      call getNetProductionRates(mixture,wi_molar(:))
      !
      wirate(:)=wi_molar(:)*wmolar(:)

      if(present(wi)) wi(:)=wirate(:)
    endif
#endif
    !
  end subroutine chemrate
  !+-------------------------------------------------------------------+
  !| The end of the subroutine chemrate.                               |
  !+-------------------------------------------------------------------+
  !
  !
end module chem