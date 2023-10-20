!+---------------------------------------------------------------------+
!|This module is used to write field in tecplot formate.                |
!+---------------------------------------------------------------------+
module WriteTec
  !
  implicit none
  !
  real(4) :: eohmarker         = 357.0,                                &
             zonemarker        = 299.0,                                &
             geometrymarker    = 399.0,                                &
             textmarker        = 499.0,                                &
             customlabelmarker = 599.0,                                &
             userrecmarker     = 699.0,                                &
             datasetauxmarker  = 799.0,                                &
             varauxmarker      = 899.0
  logical :: tecinfout=.true.
  !
  Interface writetecbin
    !
    module procedure writetecbin3d3var
    module procedure writetecbin3d4var
    module procedure writetecbin3d5var
    module procedure writetecbin3d5var_r4
    module procedure writetecbin3d6var
    module procedure writetecbin3d7var
    module procedure writetecbin3d8var
    module procedure writetecbin3d9var
    module procedure writetecbin3d10var
    module procedure writetecbin3d11var
    module procedure writetecbin3d12var
    module procedure writetecbin3d13var
    module procedure writetecbin3d14var
    module procedure writetecbin3d15var
    module procedure writetecbin3d17var
    !
    module procedure writetecbin2d2var
    module procedure writetecbin2d3var
    module procedure writetecbin2d4var
    module procedure writetecbin2d5var
    module procedure writetecbin2d6var
    module procedure writetecbin2d7var
    module procedure writetecbin2d8var
    module procedure writetecbin2d9var
    module procedure writetecbin2d10var
    module procedure writetecbin2d11var
    module procedure writetecbin2d12var
    module procedure writetecbin2d13var
    module procedure writetecbin2d14var
    !
  end Interface writetecbin
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !|This subroutine is used to write bin file for 3d tecplot field.    |
  !|ifort compiler only                                                |
  !+-------------------------------------------------------------------+
  subroutine writetecbin3d3var(filename,var1,var1name,var2,var2name,   &
                                        var3,var3name,                 &
                                        imax,jmax,kmax,                &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax,kmax
    real(8),intent(in) :: var1(0:imax,0:jmax,0:kmax),                 &
                          var2(0:imax,0:jmax,0:kmax),                 &
                          var3(0:imax,0:jmax,0:kmax)
    character(len=*),intent(in) :: var1name,var2name,var3name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=3
	  !
	  allocate(var(0:imax,0:jmax,0:kmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,0:kmax,1)=real(var1(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,2)=real(var2(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,3)=real(var3(0:imax,0:jmax,0:kmax))
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)kmax+1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
  end subroutine writetecbin3d3var
  !!
  subroutine writetecbin3d4var(filename,var1,var1name,var2,var2name,   &
                                        var3,var3name,var4,var4name,   &
                                        imax,jmax,kmax,                &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax,kmax
    real(8),intent(in) :: var1(0:imax,0:jmax,0:kmax),                 &
                          var2(0:imax,0:jmax,0:kmax),                 &
                          var3(0:imax,0:jmax,0:kmax),                 &
                          var4(0:imax,0:jmax,0:kmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=4
	  !
	  allocate(var(0:imax,0:jmax,0:kmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,0:kmax,1)=real(var1(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,2)=real(var2(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,3)=real(var3(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,4)=real(var4(0:imax,0:jmax,0:kmax))
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)kmax+1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
  end subroutine writetecbin3d4var
  !!
  subroutine writetecbin3d5var(filename,var1,var1name,var2,var2name,   &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,                 &
                                        imax,jmax,kmax,                &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax,kmax
    real(8),intent(in) :: var1(0:imax,0:jmax,0:kmax),                 &
                          var2(0:imax,0:jmax,0:kmax),                 &
                          var3(0:imax,0:jmax,0:kmax),                 &
                          var4(0:imax,0:jmax,0:kmax),                 &
                          var5(0:imax,0:jmax,0:kmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=5
	  !
	  allocate(var(0:imax,0:jmax,0:kmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,0:kmax,1)=real(var1(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,2)=real(var2(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,3)=real(var3(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,4)=real(var4(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,5)=real(var5(0:imax,0:jmax,0:kmax))
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  vname(5)=var5name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)kmax+1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
  end subroutine writetecbin3d5var
  !!
  subroutine writetecbin3d5var_r4(filename,var1,var1name,var2,var2name,&
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,                 &
                                        imax,jmax,kmax,                &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax,kmax
    real(4),intent(in) :: var1(0:imax,0:jmax,0:kmax),                 &
                          var2(0:imax,0:jmax,0:kmax),                 &
                          var3(0:imax,0:jmax,0:kmax),                 &
                          var4(0:imax,0:jmax,0:kmax),                 &
                          var5(0:imax,0:jmax,0:kmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
    if( .not. present(solutiontime) ) then
      solutiontime1=0.d0
    else
      solutiontime1=solutiontime
    endif
    if( .not. present(zonenumber) )  then
      zonenumber1=1
    else
      zonenumber1=zonenumber
    endif
    if( .not. present(title) ) then
      title1="Bin field for tecplot"
    else
      title1=title
    endif
    !
    nbrvar=5
    !
    allocate(var(0:imax,0:jmax,0:kmax,nbrvar))
    allocate(vname(nbrvar))
    !
    var(0:imax,0:jmax,0:kmax,1)=var1(0:imax,0:jmax,0:kmax)
    var(0:imax,0:jmax,0:kmax,2)=var2(0:imax,0:jmax,0:kmax)
    var(0:imax,0:jmax,0:kmax,3)=var3(0:imax,0:jmax,0:kmax)
    var(0:imax,0:jmax,0:kmax,4)=var4(0:imax,0:jmax,0:kmax)
    var(0:imax,0:jmax,0:kmax,5)=var5(0:imax,0:jmax,0:kmax)
    !
    vname(1)=var1name
    vname(2)=var2name
    vname(3)=var3name
    vname(4)=var4name
    vname(5)=var5name
    !
    unitf=18
    !
    if(tecinfout) then
      write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
      write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
      write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
    endif
    !
    open(18,file=filename,form='unformatted',access='stream')
      !
      !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
        call ecrirebin(unitf,vname(n))
        if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
      if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)kmax+1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
      !
    close(18)
    !
    print*,' << ',filename
    !
    deallocate(var)
  end subroutine writetecbin3d5var_r4
  !!
  subroutine writetecbin3d6var(filename,var1,var1name,var2,var2name,   &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        imax,jmax,kmax,                &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax,kmax
    real(8),intent(in) :: var1(0:,0:,0:),                 &
                          var2(0:,0:,0:),                 &
                          var3(0:,0:,0:),                 &
                          var4(0:,0:,0:),                 &
                          var5(0:,0:,0:),                 &
                          var6(0:,0:,0:)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name,var6name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=6
	  !
	  allocate(var(0:imax,0:jmax,0:kmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,0:kmax,1)=real(var1(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,2)=real(var2(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,3)=real(var3(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,4)=real(var4(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,5)=real(var5(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,6)=real(var6(0:imax,0:jmax,0:kmax))
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  vname(5)=var5name
	  vname(6)=var6name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)kmax+1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
  end subroutine writetecbin3d6var
  !!
  subroutine writetecbin3d7var(filename,var1,var1name,var2,var2name,   &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        var7,var7name,                 &
                                        imax,jmax,kmax,                &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax,kmax
    real(8),intent(in) :: var1(0:imax,0:jmax,0:kmax),                 &
                          var2(0:imax,0:jmax,0:kmax),                 &
                          var3(0:imax,0:jmax,0:kmax),                 &
                          var4(0:imax,0:jmax,0:kmax),                 &
                          var5(0:imax,0:jmax,0:kmax),                 &
                          var6(0:imax,0:jmax,0:kmax),                 &
                          var7(0:imax,0:jmax,0:kmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name,var6name,var7name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=7
	  !
	  allocate(var(0:imax,0:jmax,0:kmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,0:kmax,1)=real(var1(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,2)=real(var2(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,3)=real(var3(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,4)=real(var4(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,5)=real(var5(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,6)=real(var6(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,7)=real(var7(0:imax,0:jmax,0:kmax))
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  vname(5)=var5name
	  vname(6)=var6name
	  vname(7)=var7name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)kmax+1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
  end subroutine writetecbin3d7var
  !!
  subroutine writetecbin3d8var(filename,var1,var1name,var2,var2name,   &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        var7,var7name,var8,var8name,   &
                                        imax,jmax,kmax,                &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax,kmax
    real(8),intent(in) :: var1(0:imax,0:jmax,0:kmax),                 &
                          var2(0:imax,0:jmax,0:kmax),                 &
                          var3(0:imax,0:jmax,0:kmax),                 &
                          var4(0:imax,0:jmax,0:kmax),                 &
                          var5(0:imax,0:jmax,0:kmax),                 &
                          var6(0:imax,0:jmax,0:kmax),                 &
                          var7(0:imax,0:jmax,0:kmax),                 &
                          var8(0:imax,0:jmax,0:kmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name,var6name,var7name,var8name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=8
	  !
	  allocate(var(0:imax,0:jmax,0:kmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,0:kmax,1)=real(var1(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,2)=real(var2(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,3)=real(var3(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,4)=real(var4(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,5)=real(var5(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,6)=real(var6(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,7)=real(var7(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,8)=real(var8(0:imax,0:jmax,0:kmax))
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  vname(5)=var5name
	  vname(6)=var6name
	  vname(7)=var7name
	  vname(8)=var8name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)kmax+1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
  end subroutine writetecbin3d8var
  !!
  subroutine writetecbin3d9var(filename,var1,var1name,var2,var2name,   &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        var7,var7name,var8,var8name,   &
                                        var9,var9name,                 &
                                        imax,jmax,kmax,                &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax,kmax
    real(8),intent(in) :: var1(0:imax,0:jmax,0:kmax),                 &
                          var2(0:imax,0:jmax,0:kmax),                 &
                          var3(0:imax,0:jmax,0:kmax),                 &
                          var4(0:imax,0:jmax,0:kmax),                 &
                          var5(0:imax,0:jmax,0:kmax),                 &
                          var6(0:imax,0:jmax,0:kmax),                 &
                          var7(0:imax,0:jmax,0:kmax),                 &
                          var8(0:imax,0:jmax,0:kmax),                 &
                          var9(0:imax,0:jmax,0:kmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name,var6name,var7name,var8name,&
                                   var9name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=9
	  !
	  allocate(var(0:imax,0:jmax,0:kmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,0:kmax,1)=real(var1(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,2)=real(var2(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,3)=real(var3(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,4)=real(var4(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,5)=real(var5(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,6)=real(var6(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,7)=real(var7(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,8)=real(var8(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,9)=real(var9(0:imax,0:jmax,0:kmax))
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  vname(5)=var5name
	  vname(6)=var6name
	  vname(7)=var7name
	  vname(8)=var8name
	  vname(9)=var9name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)kmax+1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
  end subroutine writetecbin3d9var
  !
  subroutine writetecbin3d10var(filename,var1,var1name,var2,var2name,  &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        var7,var7name,var8,var8name,   &
                                        var9,var9name,var10,var10name, &
                                        imax,jmax,kmax,                &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax,kmax
    real(8),intent(in) :: var1(0:imax,0:jmax,0:kmax),                 &
                          var2(0:imax,0:jmax,0:kmax),                 &
                          var3(0:imax,0:jmax,0:kmax),                 &
                          var4(0:imax,0:jmax,0:kmax),                 &
                          var5(0:imax,0:jmax,0:kmax),                 &
                          var6(0:imax,0:jmax,0:kmax),                 &
                          var7(0:imax,0:jmax,0:kmax),                 &
                          var8(0:imax,0:jmax,0:kmax),                 &
                          var9(0:imax,0:jmax,0:kmax),                 &
                          var10(0:imax,0:jmax,0:kmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name,var6name,var7name,var8name,&
                                   var9name,var10name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=10
	  !
	  allocate(var(0:imax,0:jmax,0:kmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,0:kmax,1)=real(var1(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,2)=real(var2(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,3)=real(var3(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,4)=real(var4(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,5)=real(var5(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,6)=real(var6(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,7)=real(var7(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,8)=real(var8(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,9)=real(var9(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,10)=real(var10(0:imax,0:jmax,0:kmax))
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  vname(5)=var5name
	  vname(6)=var6name
	  vname(7)=var7name
	  vname(8)=var8name
	  vname(9)=var9name
	  vname(10)=var10name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)kmax+1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
    !
  end subroutine writetecbin3d10var
  !
  subroutine writetecbin3d11var(filename,var1,var1name,var2,var2name,  &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        var7,var7name,var8,var8name,   &
                                        var9,var9name,var10,var10name, &
                                        var11,var11name,               &
                                        imax,jmax,kmax,                &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax,kmax
    real(8),intent(in) :: var1(0:imax,0:jmax,0:kmax),                 &
                          var2(0:imax,0:jmax,0:kmax),                 &
                          var3(0:imax,0:jmax,0:kmax),                 &
                          var4(0:imax,0:jmax,0:kmax),                 &
                          var5(0:imax,0:jmax,0:kmax),                 &
                          var6(0:imax,0:jmax,0:kmax),                 &
                          var7(0:imax,0:jmax,0:kmax),                 &
                          var8(0:imax,0:jmax,0:kmax),                 &
                          var9(0:imax,0:jmax,0:kmax),                 &
                          var10(0:imax,0:jmax,0:kmax),                &
                          var11(0:imax,0:jmax,0:kmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name,var6name,var7name,var8name,&
                                   var9name,var10name,var11name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=11
	  !
	  allocate(var(0:imax,0:jmax,0:kmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,0:kmax,1)=real(var1(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,2)=real(var2(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,3)=real(var3(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,4)=real(var4(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,5)=real(var5(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,6)=real(var6(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,7)=real(var7(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,8)=real(var8(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,9)=real(var9(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,10)=real(var10(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,11)=real(var11(0:imax,0:jmax,0:kmax))
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  vname(5)=var5name
	  vname(6)=var6name
	  vname(7)=var7name
	  vname(8)=var8name
	  vname(9)=var9name
	  vname(10)=var10name
	  vname(11)=var11name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)kmax+1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
    !
  end subroutine writetecbin3d11var
  !!
  subroutine writetecbin3d12var(filename,var1,var1name,var2,var2name,  &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        var7,var7name,var8,var8name,   &
                                        var9,var9name,var10,var10name, &
                                      var11,var11name,var12,var12name, &
                                        imax,jmax,kmax,                &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax,kmax
    real(8),intent(in) :: var1(0:imax,0:jmax,0:kmax),                 &
                          var2(0:imax,0:jmax,0:kmax),                 &
                          var3(0:imax,0:jmax,0:kmax),                 &
                          var4(0:imax,0:jmax,0:kmax),                 &
                          var5(0:imax,0:jmax,0:kmax),                 &
                          var6(0:imax,0:jmax,0:kmax),                 &
                          var7(0:imax,0:jmax,0:kmax),                 &
                          var8(0:imax,0:jmax,0:kmax),                 &
                          var9(0:imax,0:jmax,0:kmax),                 &
                          var10(0:imax,0:jmax,0:kmax),                &
                          var11(0:imax,0:jmax,0:kmax),                &
                          var12(0:imax,0:jmax,0:kmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name,var6name,var7name,var8name,&
                                   var9name,var10name,var11name,       &
                                   var12name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=12
	  !
	  allocate(var(0:imax,0:jmax,0:kmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,0:kmax,1)=real(var1(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,2)=real(var2(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,3)=real(var3(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,4)=real(var4(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,5)=real(var5(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,6)=real(var6(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,7)=real(var7(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,8)=real(var8(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,9)=real(var9(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,10)=real(var10(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,11)=real(var11(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,12)=real(var12(0:imax,0:jmax,0:kmax))
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  vname(5)=var5name
	  vname(6)=var6name
	  vname(7)=var7name
	  vname(8)=var8name
	  vname(9)=var9name
	  vname(10)=var10name
	  vname(11)=var11name
	  vname(12)=var12name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)kmax+1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
    !
  end subroutine writetecbin3d12var
  !
  subroutine writetecbin3d13var(filename,var1,var1name,var2,var2name,  &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        var7,var7name,var8,var8name,   &
                                        var9,var9name,var10,var10name, &
                                      var11,var11name,var12,var12name, &
                                      var13,var13name,                 &
                                        imax,jmax,kmax,                &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax,kmax
    real(8),intent(in) :: var1(0:imax,0:jmax,0:kmax),                 &
                          var2(0:imax,0:jmax,0:kmax),                 &
                          var3(0:imax,0:jmax,0:kmax),                 &
                          var4(0:imax,0:jmax,0:kmax),                 &
                          var5(0:imax,0:jmax,0:kmax),                 &
                          var6(0:imax,0:jmax,0:kmax),                 &
                          var7(0:imax,0:jmax,0:kmax),                 &
                          var8(0:imax,0:jmax,0:kmax),                 &
                          var9(0:imax,0:jmax,0:kmax),                 &
                          var10(0:imax,0:jmax,0:kmax),                &
                          var11(0:imax,0:jmax,0:kmax),                &
                          var12(0:imax,0:jmax,0:kmax),                &
                          var13(0:imax,0:jmax,0:kmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name,var6name,var7name,var8name,&
                                   var9name,var10name,var11name,       &
                                   var12name,var13name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=13
	  !
	  allocate(var(0:imax,0:jmax,0:kmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,0:kmax,1)=real(var1(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,2)=real(var2(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,3)=real(var3(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,4)=real(var4(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,5)=real(var5(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,6)=real(var6(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,7)=real(var7(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,8)=real(var8(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,9)=real(var9(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,10)=real(var10(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,11)=real(var11(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,12)=real(var12(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,13)=real(var13(0:imax,0:jmax,0:kmax))
	  !
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  vname(5)=var5name
	  vname(6)=var6name
	  vname(7)=var7name
	  vname(8)=var8name
	  vname(9)=var9name
	  vname(10)=var10name
	  vname(11)=var11name
	  vname(12)=var12name
	  vname(13)=var13name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)kmax+1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
    !
  end subroutine writetecbin3d13var
  !
  subroutine writetecbin3d14var(filename,var1,var1name,var2,var2name,  &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        var7,var7name,var8,var8name,   &
                                        var9,var9name,var10,var10name, &
                                      var11,var11name,var12,var12name, &
                                      var13,var13name,var14,var14name, &
                                        imax,jmax,kmax,                &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax,kmax
    real(8),intent(in) :: var1(0:imax,0:jmax,0:kmax),                 &
                          var2(0:imax,0:jmax,0:kmax),                 &
                          var3(0:imax,0:jmax,0:kmax),                 &
                          var4(0:imax,0:jmax,0:kmax),                 &
                          var5(0:imax,0:jmax,0:kmax),                 &
                          var6(0:imax,0:jmax,0:kmax),                 &
                          var7(0:imax,0:jmax,0:kmax),                 &
                          var8(0:imax,0:jmax,0:kmax),                 &
                          var9(0:imax,0:jmax,0:kmax),                 &
                          var10(0:imax,0:jmax,0:kmax),                &
                          var11(0:imax,0:jmax,0:kmax),                &
                          var12(0:imax,0:jmax,0:kmax),                &
                          var13(0:imax,0:jmax,0:kmax),                &
                          var14(0:imax,0:jmax,0:kmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name,var6name,var7name,var8name,&
                                   var9name,var10name,var11name,       &
                                   var12name,var13name,var14name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=14
	  !
	  allocate(var(0:imax,0:jmax,0:kmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,0:kmax,1)=real(var1(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,2)=real(var2(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,3)=real(var3(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,4)=real(var4(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,5)=real(var5(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,6)=real(var6(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,7)=real(var7(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,8)=real(var8(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,9)=real(var9(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,10)=real(var10(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,11)=real(var11(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,12)=real(var12(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,13)=real(var13(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,14)=real(var14(0:imax,0:jmax,0:kmax))
	  !
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  vname(5)=var5name
	  vname(6)=var6name
	  vname(7)=var7name
	  vname(8)=var8name
	  vname(9)=var9name
	  vname(10)=var10name
	  vname(11)=var11name
	  vname(12)=var12name
	  vname(13)=var13name
	  vname(14)=var14name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)kmax+1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
    !
  end subroutine writetecbin3d14var
  !
  subroutine writetecbin3d15var(filename,var1,var1name,var2,var2name,  &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        var7,var7name,var8,var8name,   &
                                        var9,var9name,var10,var10name, &
                                      var11,var11name,var12,var12name, &
                                      var13,var13name,var14,var14name, &
                                      var15,var15name,                 &
                                        imax,jmax,kmax,                &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax,kmax
    real(8),intent(in) :: var1(0:imax,0:jmax,0:kmax),                 &
                          var2(0:imax,0:jmax,0:kmax),                 &
                          var3(0:imax,0:jmax,0:kmax),                 &
                          var4(0:imax,0:jmax,0:kmax),                 &
                          var5(0:imax,0:jmax,0:kmax),                 &
                          var6(0:imax,0:jmax,0:kmax),                 &
                          var7(0:imax,0:jmax,0:kmax),                 &
                          var8(0:imax,0:jmax,0:kmax),                 &
                          var9(0:imax,0:jmax,0:kmax),                 &
                          var10(0:imax,0:jmax,0:kmax),                &
                          var11(0:imax,0:jmax,0:kmax),                &
                          var12(0:imax,0:jmax,0:kmax),                &
                          var13(0:imax,0:jmax,0:kmax),                &
                          var14(0:imax,0:jmax,0:kmax),                &
                          var15(0:imax,0:jmax,0:kmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name,var6name,var7name,var8name,&
                                   var9name,var10name,var11name,       &
                                   var12name,var13name,var14name,var15name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=15
	  !
	  allocate(var(0:imax,0:jmax,0:kmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,0:kmax,1)=real(var1(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,2)=real(var2(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,3)=real(var3(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,4)=real(var4(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,5)=real(var5(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,6)=real(var6(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,7)=real(var7(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,8)=real(var8(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,9)=real(var9(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,10)=real(var10(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,11)=real(var11(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,12)=real(var12(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,13)=real(var13(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,14)=real(var14(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,15)=real(var15(0:imax,0:jmax,0:kmax))
	  !
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  vname(5)=var5name
	  vname(6)=var6name
	  vname(7)=var7name
	  vname(8)=var8name
	  vname(9)=var9name
	  vname(10)=var10name
	  vname(11)=var11name
	  vname(12)=var12name
	  vname(13)=var13name
	  vname(14)=var14name
	  vname(15)=var15name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)kmax+1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
    !
  end subroutine writetecbin3d15var
  !
  subroutine writetecbin3d17var(filename,var1,var1name,var2,var2name,  &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        var7,var7name,var8,var8name,   &
                                        var9,var9name,var10,var10name, &
                                      var11,var11name,var12,var12name, &
                                      var13,var13name,var14,var14name, &
                                      var15,var15name,var16,var16name, &
                                      var17,var17name,                 &
                                        imax,jmax,kmax,                &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax,kmax
    real(8),intent(in) :: var1(0:imax,0:jmax,0:kmax),                 &
                          var2(0:imax,0:jmax,0:kmax),                 &
                          var3(0:imax,0:jmax,0:kmax),                 &
                          var4(0:imax,0:jmax,0:kmax),                 &
                          var5(0:imax,0:jmax,0:kmax),                 &
                          var6(0:imax,0:jmax,0:kmax),                 &
                          var7(0:imax,0:jmax,0:kmax),                 &
                          var8(0:imax,0:jmax,0:kmax),                 &
                          var9(0:imax,0:jmax,0:kmax),                 &
                          var10(0:imax,0:jmax,0:kmax),                &
                          var11(0:imax,0:jmax,0:kmax),                &
                          var12(0:imax,0:jmax,0:kmax),                &
                          var13(0:imax,0:jmax,0:kmax),                &
                          var14(0:imax,0:jmax,0:kmax),                &
                          var15(0:imax,0:jmax,0:kmax),                &
                          var16(0:imax,0:jmax,0:kmax),                &
                          var17(0:imax,0:jmax,0:kmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name,var6name,var7name,var8name,&
                                   var9name,var10name,var11name,       &
                                   var12name,var13name,var14name,      &
                                   var15name,var16name,var17name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=17
	  !
	  allocate(var(0:imax,0:jmax,0:kmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,0:kmax,1)=real(var1(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,2)=real(var2(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,3)=real(var3(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,4)=real(var4(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,5)=real(var5(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,6)=real(var6(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,7)=real(var7(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,8)=real(var8(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,9)=real(var9(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,10)=real(var10(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,11)=real(var11(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,12)=real(var12(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,13)=real(var13(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,14)=real(var14(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,15)=real(var15(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,17)=real(var16(0:imax,0:jmax,0:kmax))
	  var(0:imax,0:jmax,0:kmax,17)=real(var17(0:imax,0:jmax,0:kmax))
	  !
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  vname(5)=var5name
	  vname(6)=var6name
	  vname(7)=var7name
	  vname(8)=var8name
	  vname(9)=var9name
	  vname(10)=var10name
	  vname(11)=var11name
	  vname(12)=var12name
	  vname(13)=var13name
	  vname(14)=var14name
	  vname(15)=var15name
	  vname(16)=var16name
	  vname(17)=var17name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)kmax+1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,0:kmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
    !
  end subroutine writetecbin3d17var
  !
  !+-------------------------------------------------------------------+
  !|This subroutine is used to write writetecbin3d.                    |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !|This subroutine is used to write bin file for 2d tecplot field.    |
  !|ifort compiler only                                                |
  !+-------------------------------------------------------------------+
  subroutine writetecbin2d2var(filename,var1,var1name,var2,var2name,   &
                                        imax,jmax,                     &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax
    real(8),intent(in) :: var1(0:imax,0:jmax),                         &
                          var2(0:imax,0:jmax)
    character(len=*),intent(in) :: var1name,var2name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=2
	  !
	  allocate(var(0:imax,0:jmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,1)=real(var1(0:imax,0:jmax))
	  var(0:imax,0:jmax,2)=real(var2(0:imax,0:jmax))
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
    !
  end subroutine writetecbin2d2var
  !!
  subroutine writetecbin2d3var(filename,var1,var1name,var2,var2name,   &
                                        var3,var3name,                 &
                                        imax,jmax,                     &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax
    real(8),intent(in) :: var1(0:imax,0:jmax),                         &
                          var2(0:imax,0:jmax),                         &
                          var3(0:imax,0:jmax)
    character(len=*),intent(in) :: var1name,var2name,var3name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=3
	  !
	  allocate(var(0:imax,0:jmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,1)=real(var1(0:imax,0:jmax))
	  var(0:imax,0:jmax,2)=real(var2(0:imax,0:jmax))
	  var(0:imax,0:jmax,3)=real(var3(0:imax,0:jmax))
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
    !
  end subroutine writetecbin2d3var
  !!
  subroutine writetecbin2d4var(filename,var1,var1name,var2,var2name,   &
                                        var3,var3name,var4,var4name,   &
                                        imax,jmax,                     &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax
    real(8),intent(in) :: var1(0:imax,0:jmax),                         &
                          var2(0:imax,0:jmax),                         &
                          var3(0:imax,0:jmax),                         &
                          var4(0:imax,0:jmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=4
	  !
	  allocate(var(0:imax,0:jmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,1)=real(var1(0:imax,0:jmax))
	  var(0:imax,0:jmax,2)=real(var2(0:imax,0:jmax))
	  var(0:imax,0:jmax,3)=real(var3(0:imax,0:jmax))
	  var(0:imax,0:jmax,4)=real(var4(0:imax,0:jmax))
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
    !
  end subroutine writetecbin2d4var
  !!
  subroutine writetecbin2d5var(filename,var1,var1name,var2,var2name,   &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,                 &
                                        imax,jmax,                     &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax
    real(8),intent(in) :: var1(0:imax,0:jmax),                         &
                          var2(0:imax,0:jmax),                         &
                          var3(0:imax,0:jmax),                         &
                          var4(0:imax,0:jmax),                         &
                          var5(0:imax,0:jmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=5
	  !
	  allocate(var(0:imax,0:jmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,1)=real(var1(0:imax,0:jmax))
	  var(0:imax,0:jmax,2)=real(var2(0:imax,0:jmax))
	  var(0:imax,0:jmax,3)=real(var3(0:imax,0:jmax))
	  var(0:imax,0:jmax,4)=real(var4(0:imax,0:jmax))
	  var(0:imax,0:jmax,5)=real(var5(0:imax,0:jmax))
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  vname(5)=var5name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
    !
  end subroutine writetecbin2d5var
  !!
  subroutine writetecbin2d6var(filename,var1,var1name,var2,var2name,   &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        imax,jmax,                     &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax
    real(8),intent(in) :: var1(0:imax,0:jmax),                         &
                          var2(0:imax,0:jmax),                         &
                          var3(0:imax,0:jmax),                         &
                          var4(0:imax,0:jmax),                         &
                          var5(0:imax,0:jmax),                         &
                          var6(0:imax,0:jmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name,var6name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=6
	  !
	  allocate(var(0:imax,0:jmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,1)=real(var1(0:imax,0:jmax))
	  var(0:imax,0:jmax,2)=real(var2(0:imax,0:jmax))
	  var(0:imax,0:jmax,3)=real(var3(0:imax,0:jmax))
	  var(0:imax,0:jmax,4)=real(var4(0:imax,0:jmax))
	  var(0:imax,0:jmax,5)=real(var5(0:imax,0:jmax))
	  var(0:imax,0:jmax,6)=real(var6(0:imax,0:jmax))
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  vname(5)=var5name
	  vname(6)=var6name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
    !
  end subroutine writetecbin2d6var
  !!
  subroutine writetecbin2d7var(filename,var1,var1name,var2,var2name,   &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        var7,var7name,                 &
                                        imax,jmax,                     &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax
    real(8),intent(in) :: var1(0:imax,0:jmax),                         &
                          var2(0:imax,0:jmax),                         &
                          var3(0:imax,0:jmax),                         &
                          var4(0:imax,0:jmax),                         &
                          var5(0:imax,0:jmax),                         &
                          var6(0:imax,0:jmax),                         &
                          var7(0:imax,0:jmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name,var6name,var7name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=7
	  !
	  allocate(var(0:imax,0:jmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,1)=real(var1(0:imax,0:jmax))
	  var(0:imax,0:jmax,2)=real(var2(0:imax,0:jmax))
	  var(0:imax,0:jmax,3)=real(var3(0:imax,0:jmax))
	  var(0:imax,0:jmax,4)=real(var4(0:imax,0:jmax))
	  var(0:imax,0:jmax,5)=real(var5(0:imax,0:jmax))
	  var(0:imax,0:jmax,6)=real(var6(0:imax,0:jmax))
	  var(0:imax,0:jmax,7)=real(var7(0:imax,0:jmax))
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  vname(5)=var5name
	  vname(6)=var6name
	  vname(7)=var7name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
    !
  end subroutine writetecbin2d7var
  !!
  subroutine writetecbin2d8var(filename,var1,var1name,var2,var2name,   &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        var7,var7name,var8,var8name,   &
                                        imax,jmax,                     &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax
    real(8),intent(in) :: var1(0:imax,0:jmax),                         &
                          var2(0:imax,0:jmax),                         &
                          var3(0:imax,0:jmax),                         &
                          var4(0:imax,0:jmax),                         &
                          var5(0:imax,0:jmax),                         &
                          var6(0:imax,0:jmax),                         &
                          var7(0:imax,0:jmax),                         &
                          var8(0:imax,0:jmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name,var6name,var7name,var8name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=8
	  !
	  allocate(var(0:imax,0:jmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,1)=real(var1(0:imax,0:jmax))
	  var(0:imax,0:jmax,2)=real(var2(0:imax,0:jmax))
	  var(0:imax,0:jmax,3)=real(var3(0:imax,0:jmax))
	  var(0:imax,0:jmax,4)=real(var4(0:imax,0:jmax))
	  var(0:imax,0:jmax,5)=real(var5(0:imax,0:jmax))
	  var(0:imax,0:jmax,6)=real(var6(0:imax,0:jmax))
	  var(0:imax,0:jmax,7)=real(var7(0:imax,0:jmax))
	  var(0:imax,0:jmax,8)=real(var8(0:imax,0:jmax))
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  vname(5)=var5name
	  vname(6)=var6name
	  vname(7)=var7name
	  vname(8)=var8name
	  !
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
    !
  end subroutine writetecbin2d8var
  !!
  subroutine writetecbin2d9var(filename,var1,var1name,var2,var2name,   &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        var7,var7name,var8,var8name,   &
                                        var9,var9name,                 &
                                        imax,jmax,                     &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax
    real(8),intent(in) :: var1(0:,0:),                         &
                          var2(0:,0:),                         &
                          var3(0:,0:),                         &
                          var4(0:,0:),                         &
                          var5(0:,0:),                         &
                          var6(0:,0:),                         &
                          var7(0:,0:),                         &
                          var8(0:,0:),                         &
                          var9(0:,0:)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name,var6name,var7name,var8name,&
                                   var9name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=9
	  !
	  allocate(var(0:imax,0:jmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,1)=real(var1(0:imax,0:jmax))
	  var(0:imax,0:jmax,2)=real(var2(0:imax,0:jmax))
	  var(0:imax,0:jmax,3)=real(var3(0:imax,0:jmax))
	  var(0:imax,0:jmax,4)=real(var4(0:imax,0:jmax))
	  var(0:imax,0:jmax,5)=real(var5(0:imax,0:jmax))
	  var(0:imax,0:jmax,6)=real(var6(0:imax,0:jmax))
	  var(0:imax,0:jmax,7)=real(var7(0:imax,0:jmax))
	  var(0:imax,0:jmax,8)=real(var8(0:imax,0:jmax))
	  var(0:imax,0:jmax,9)=real(var9(0:imax,0:jmax))
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  vname(5)=var5name
	  vname(6)=var6name
	  vname(7)=var7name
	  vname(8)=var8name
	  vname(9)=var9name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
    !
  end subroutine writetecbin2d9var
  !!
  subroutine writetecbin2d10var(filename,var1,var1name,var2,var2name,  &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        var7,var7name,var8,var8name,   &
                                        var9,var9name,var10,var10name, &
                                        imax,jmax,                     &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax
    real(8),intent(in) :: var1(0:,0:),                         &
                          var2(0:,0:),                         &
                          var3(0:,0:),                         &
                          var4(0:,0:),                         &
                          var5(0:,0:),                         &
                          var6(0:,0:),                         &
                          var7(0:,0:),                         &
                          var8(0:,0:),                         &
                          var9(0:,0:),                         &
                         var10(0:,0:)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name,var6name,var7name,var8name,&
                                   var9name,var10name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=10
	  !
	  allocate(var(0:imax,0:jmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,1)=real(var1(0:imax,0:jmax))
	  var(0:imax,0:jmax,2)=real(var2(0:imax,0:jmax))
	  var(0:imax,0:jmax,3)=real(var3(0:imax,0:jmax))
	  var(0:imax,0:jmax,4)=real(var4(0:imax,0:jmax))
	  var(0:imax,0:jmax,5)=real(var5(0:imax,0:jmax))
	  var(0:imax,0:jmax,6)=real(var6(0:imax,0:jmax))
	  var(0:imax,0:jmax,7)=real(var7(0:imax,0:jmax))
	  var(0:imax,0:jmax,8)=real(var8(0:imax,0:jmax))
	  var(0:imax,0:jmax,9)=real(var9(0:imax,0:jmax))
	  var(0:imax,0:jmax,10)=real(var10(0:imax,0:jmax))
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  vname(5)=var5name
	  vname(6)=var6name
	  vname(7)=var7name
	  vname(8)=var8name
	  vname(9)=var9name
	  vname(10)=var10name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
    !
  end subroutine writetecbin2d10var
  !
  subroutine writetecbin2d11var(filename,var1,var1name,var2,var2name,  &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        var7,var7name,var8,var8name,   &
                                        var9,var9name,var10,var10name, &
                                        var11,var11name,               &
                                        imax,jmax,                     &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax
    real(8),intent(in) :: var1(0:imax,0:jmax),                         &
                          var2(0:imax,0:jmax),                         &
                          var3(0:imax,0:jmax),                         &
                          var4(0:imax,0:jmax),                         &
                          var5(0:imax,0:jmax),                         &
                          var6(0:imax,0:jmax),                         &
                          var7(0:imax,0:jmax),                         &
                          var8(0:imax,0:jmax),                         &
                          var9(0:imax,0:jmax),                         &
                          var10(0:imax,0:jmax),                         &
                          var11(0:imax,0:jmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name,var6name,var7name,var8name,&
                                   var9name,var10name,var11name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=11
	  !
	  allocate(var(0:imax,0:jmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,1)=real(var1(0:imax,0:jmax))
	  var(0:imax,0:jmax,2)=real(var2(0:imax,0:jmax))
	  var(0:imax,0:jmax,3)=real(var3(0:imax,0:jmax))
	  var(0:imax,0:jmax,4)=real(var4(0:imax,0:jmax))
	  var(0:imax,0:jmax,5)=real(var5(0:imax,0:jmax))
	  var(0:imax,0:jmax,6)=real(var6(0:imax,0:jmax))
	  var(0:imax,0:jmax,7)=real(var7(0:imax,0:jmax))
	  var(0:imax,0:jmax,8)=real(var8(0:imax,0:jmax))
	  var(0:imax,0:jmax,9)=real(var9(0:imax,0:jmax))
	  var(0:imax,0:jmax,10)=real(var10(0:imax,0:jmax))
	  var(0:imax,0:jmax,11)=real(var11(0:imax,0:jmax))
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  vname(5)=var5name
	  vname(6)=var6name
	  vname(7)=var7name
	  vname(8)=var8name
	  vname(9)=var9name
	  vname(10)=var10name
	  vname(11)=var11name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
    !
  end subroutine writetecbin2d11var
  !
  subroutine writetecbin2d12var(filename,var1,var1name,var2,var2name,  &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        var7,var7name,var8,var8name,   &
                                        var9,var9name,var10,var10name, &
                                      var11,var11name,var12,var12name, &
                                        imax,jmax,                     &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax
    real(8),intent(in) :: var1(0:imax,0:jmax),                         &
                          var2(0:imax,0:jmax),                         &
                          var3(0:imax,0:jmax),                         &
                          var4(0:imax,0:jmax),                         &
                          var5(0:imax,0:jmax),                         &
                          var6(0:imax,0:jmax),                         &
                          var7(0:imax,0:jmax),                         &
                          var8(0:imax,0:jmax),                         &
                          var9(0:imax,0:jmax),                         &
                          var10(0:imax,0:jmax),                        &
                          var11(0:imax,0:jmax),                        &
                          var12(0:imax,0:jmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name,var6name,var7name,var8name,&
                                   var9name,var10name,var11name,var12name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
	  if( .not. present(solutiontime) ) then
	  	solutiontime1=0.d0
	  else
	  	solutiontime1=solutiontime
	  endif
	  if( .not. present(zonenumber) )  then
	  	zonenumber1=1
	  else
	  	zonenumber1=zonenumber
	  endif
	  if( .not. present(title) ) then
	  	title1="Bin field for tecplot"
	  else
	  	title1=title
	  endif
	  !
	  nbrvar=12
	  !
	  allocate(var(0:imax,0:jmax,nbrvar))
	  allocate(vname(nbrvar))
	  !
	  var(0:imax,0:jmax,1)=real(var1(0:imax,0:jmax))
	  var(0:imax,0:jmax,2)=real(var2(0:imax,0:jmax))
	  var(0:imax,0:jmax,3)=real(var3(0:imax,0:jmax))
	  var(0:imax,0:jmax,4)=real(var4(0:imax,0:jmax))
	  var(0:imax,0:jmax,5)=real(var5(0:imax,0:jmax))
	  var(0:imax,0:jmax,6)=real(var6(0:imax,0:jmax))
	  var(0:imax,0:jmax,7)=real(var7(0:imax,0:jmax))
	  var(0:imax,0:jmax,8)=real(var8(0:imax,0:jmax))
	  var(0:imax,0:jmax,9)=real(var9(0:imax,0:jmax))
	  var(0:imax,0:jmax,10)=real(var10(0:imax,0:jmax))
	  var(0:imax,0:jmax,11)=real(var11(0:imax,0:jmax))
	  var(0:imax,0:jmax,12)=real(var12(0:imax,0:jmax))
	  !
	  vname(1)=var1name
	  vname(2)=var2name
	  vname(3)=var3name
	  vname(4)=var4name
	  vname(5)=var5name
	  vname(6)=var6name
	  vname(7)=var7name
	  vname(8)=var8name
	  vname(9)=var9name
	  vname(10)=var10name
	  vname(11)=var11name
	  vname(12)=var12name
	  !
	  unitf=18
	  !
	  if(tecinfout) then
	    write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
	    write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
	    write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
	  endif
	  !
	  open(18,file=filename,form='unformatted',access='stream')
      !
	    !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
      	call ecrirebin(unitf,vname(n))
	      if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
	    if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
	    !
	  close(18)
	  !
	  print*,' << ',filename
	  !
	  deallocate(var)
    !
  end subroutine writetecbin2d12var
  !!
  subroutine writetecbin2d13var(filename,var1,var1name,var2,var2name,  &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        var7,var7name,var8,var8name,   &
                                        var9,var9name,var10,var10name, &
                                      var11,var11name,var12,var12name, &
                                      var13,var13name,                 &
                                        imax,jmax,                     &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax
    real(8),intent(in) :: var1(0:imax,0:jmax),                         &
                          var2(0:imax,0:jmax),                         &
                          var3(0:imax,0:jmax),                         &
                          var4(0:imax,0:jmax),                         &
                          var5(0:imax,0:jmax),                         &
                          var6(0:imax,0:jmax),                         &
                          var7(0:imax,0:jmax),                         &
                          var8(0:imax,0:jmax),                         &
                          var9(0:imax,0:jmax),                         &
                          var10(0:imax,0:jmax),                        &
                          var11(0:imax,0:jmax),                        &
                          var12(0:imax,0:jmax),                        &
                          var13(0:imax,0:jmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name,var6name,var7name,var8name,&
                                   var9name,var10name,var11name,       &
                                   var12name,var13name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
    if( .not. present(solutiontime) ) then
      solutiontime1=0.d0
    else
      solutiontime1=solutiontime
    endif
    if( .not. present(zonenumber) )  then
      zonenumber1=1
    else
      zonenumber1=zonenumber
    endif
    if( .not. present(title) ) then
      title1="Bin field for tecplot"
    else
      title1=title
    endif
    !
    nbrvar=13
    !
    allocate(var(0:imax,0:jmax,nbrvar))
    allocate(vname(nbrvar))
    !
    var(0:imax,0:jmax,1)=real(var1(0:imax,0:jmax))
    var(0:imax,0:jmax,2)=real(var2(0:imax,0:jmax))
    var(0:imax,0:jmax,3)=real(var3(0:imax,0:jmax))
    var(0:imax,0:jmax,4)=real(var4(0:imax,0:jmax))
    var(0:imax,0:jmax,5)=real(var5(0:imax,0:jmax))
    var(0:imax,0:jmax,6)=real(var6(0:imax,0:jmax))
    var(0:imax,0:jmax,7)=real(var7(0:imax,0:jmax))
    var(0:imax,0:jmax,8)=real(var8(0:imax,0:jmax))
    var(0:imax,0:jmax,9)=real(var9(0:imax,0:jmax))
    var(0:imax,0:jmax,10)=real(var10(0:imax,0:jmax))
    var(0:imax,0:jmax,11)=real(var11(0:imax,0:jmax))
    var(0:imax,0:jmax,12)=real(var12(0:imax,0:jmax))
    var(0:imax,0:jmax,13)=real(var13(0:imax,0:jmax))
    !
    vname(1)=var1name
    vname(2)=var2name
    vname(3)=var3name
    vname(4)=var4name
    vname(5)=var5name
    vname(6)=var6name
    vname(7)=var7name
    vname(8)=var8name
    vname(9)=var9name
    vname(10)=var10name
    vname(11)=var11name
    vname(12)=var12name
    vname(13)=var13name
    !
    unitf=18
    !
    if(tecinfout) then
      write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
      write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
      write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
    endif
    !
    open(18,file=filename,form='unformatted',access='stream')
      !
      !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
        call ecrirebin(unitf,vname(n))
        if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
      if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
      !
    close(18)
    !
    print*,' << ',filename
    !
    deallocate(var)
    !
  end subroutine writetecbin2d13var
  !
  subroutine writetecbin2d14var(filename,var1,var1name,var2,var2name,  &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        var7,var7name,var8,var8name,   &
                                        var9,var9name,var10,var10name, &
                                      var11,var11name,var12,var12name, &
                                      var13,var13name,var14,var14name, &
                                        imax,jmax,                     &
                                        solutiontime,zonenumber,title)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax
    real(8),intent(in) :: var1(0:imax,0:jmax),                         &
                          var2(0:imax,0:jmax),                         &
                          var3(0:imax,0:jmax),                         &
                          var4(0:imax,0:jmax),                         &
                          var5(0:imax,0:jmax),                         &
                          var6(0:imax,0:jmax),                         &
                          var7(0:imax,0:jmax),                         &
                          var8(0:imax,0:jmax),                         &
                          var9(0:imax,0:jmax),                         &
                          var10(0:imax,0:jmax),                        &
                          var11(0:imax,0:jmax),                        &
                          var12(0:imax,0:jmax),                        &
                          var13(0:imax,0:jmax),                        &
                          var14(0:imax,0:jmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                   var5name,var6name,var7name,var8name,&
                                   var9name,var10name,var11name,       &
                                   var12name,var13name,var14name
    real(8),optional,intent(in) :: solutiontime
    integer,optional,intent(in) :: zonenumber
    character(len=*),optional,intent(in) :: title
    !
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
    if( .not. present(solutiontime) ) then
      solutiontime1=0.d0
    else
      solutiontime1=solutiontime
    endif
    if( .not. present(zonenumber) )  then
      zonenumber1=1
    else
      zonenumber1=zonenumber
    endif
    if( .not. present(title) ) then
      title1="Bin field for tecplot"
    else
      title1=title
    endif
    !
    nbrvar=14
    !
    allocate(var(0:imax,0:jmax,nbrvar))
    allocate(vname(nbrvar))
    !
    var(0:imax,0:jmax,1)=real(var1(0:imax,0:jmax))
    var(0:imax,0:jmax,2)=real(var2(0:imax,0:jmax))
    var(0:imax,0:jmax,3)=real(var3(0:imax,0:jmax))
    var(0:imax,0:jmax,4)=real(var4(0:imax,0:jmax))
    var(0:imax,0:jmax,5)=real(var5(0:imax,0:jmax))
    var(0:imax,0:jmax,6)=real(var6(0:imax,0:jmax))
    var(0:imax,0:jmax,7)=real(var7(0:imax,0:jmax))
    var(0:imax,0:jmax,8)=real(var8(0:imax,0:jmax))
    var(0:imax,0:jmax,9)=real(var9(0:imax,0:jmax))
    var(0:imax,0:jmax,10)=real(var10(0:imax,0:jmax))
    var(0:imax,0:jmax,11)=real(var11(0:imax,0:jmax))
    var(0:imax,0:jmax,12)=real(var12(0:imax,0:jmax))
    var(0:imax,0:jmax,13)=real(var13(0:imax,0:jmax))
    var(0:imax,0:jmax,14)=real(var14(0:imax,0:jmax))
    !
    vname(1)=var1name
    vname(2)=var2name
    vname(3)=var3name
    vname(4)=var4name
    vname(5)=var5name
    vname(6)=var6name
    vname(7)=var7name
    vname(8)=var8name
    vname(9)=var9name
    vname(10)=var10name
    vname(11)=var11name
    vname(12)=var12name
    vname(13)=var13name
    vname(14)=var14name
    !
    unitf=18
    !
    if(tecinfout) then
      write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
      write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
      write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
    endif
    !
    open(18,file=filename,form='unformatted',access='stream')
      !
      !i. header section
      ! i. magic number, version number
      ! +------------+
      ! | "#!tdv112" |
      ! +------------+
      write(unitf)"#!TDV112"
      ! ii. integer value of 1
      ! +------------+
      ! | int32      |
      ! +------------+
      int32=1
      write(unitf)int32
      ! iii. title and variable names
      ! +------------+
      ! | int32      | filetype: 0=full, 1=grid, 2=solution
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*n    | the title
      ! +------------+?
      call ecrirebin(unitf,title1)
      ! +------------+
      ! | int32      | number of variables in the datafile
      ! +------------+
      write(unitf)nbrvar
      ! +------------+
      ! | int32*n    | variable names
      ! +------------+
      do n=1,nbrvar
        call ecrirebin(unitf,vname(n))
        if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
      enddo
      ! iv. zones
      ! +------------+
      ! | float32    | zone marker. value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | zone name
      ! +------------+ 
      Ligne=""
      write(Ligne,"(A,I3.3)")"Zone",zonenumber1
      call EcrireBin(UnitF,Ligne)
      if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
      ! +------------+
      ! | int32      | parentzone
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | strandid
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | float64    | solution time
      ! +------------+
      write(unitf)solutiontime1
      ! +------------+
      ! | int32      | not used. set to -1
      ! +------------+
      int32=-1
      write(unitf)int32
      ! +------------+
      ! | int32      | zonetype 0=ordered,       1=felineseg,
      ! +------------+          2=fetriangle,    3=fequadrilateral,
      !                         4=fetetrahedron, 5=febrick,
      !                         6=fepolygon,     7=fepolyhedron
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | specify var location
      ! +------------+    0 = don't specify, 1 = specify
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | are raw local 1-to-1 face neighbors supplied?
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32      | number of miscellaneous user-defined face neighbor connections
      ! +------------+
      int32=0
      write(unitf)int32
      ! +------------+
      ! | int32*3    | imax,jmax,kmax
      ! +------------+
      write(unitf)imax+1
      write(unitf)jmax+1
      write(unitf)1
      ! +------------+
      ! | int32      | 1=auxiliary name/value pair to follow
      ! +------------+ 0=no more auxiliary name/value pairs
      int32=0
      write(unitf)int32
      ! +------------+
      ! | float32    | eohmarker, value = 357.0, end of header section
      ! +------------+
      write(unitf)eohmarker
      !ii. data section
      ! i. for both ordered and fe zones
      ! +------------+
      ! | float32    | zone marker value = 299.0
      ! +------------+
      write(unitf)zonemarker
      ! +------------+
      ! | int32*n    | variable data format, n=total number of vars
      ! +------------+     1=float,    2=double, 3=longint
      !                    4=shortint, 5=byte,   6=bit
      do n=1,nbrvar
        int32=1
        write(18)int32
      enddo
      ! +------------+
      ! | int32      | has passive variables: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | has variable sharing: 0=no, 1=yes
      ! +------------+
      int32=0
      write(18)int32
      ! +------------+
      ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
      ! +------------+
      int32=-1
      write(18)int32
      !
      do n=1,nbrvar
        !+------------+
        !| float64    | min value
        !+------------+
        float32=minval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
        !+------------+
        !| float64    | max value
        !+------------+
        float32=maxval(var(0:imax,0:jmax,n))
        float64=real(float32,8)
        write(18)float64
      enddo
      ! +------------+
      ! | xxxxxxxxxx | zone data
      ! +------------+
      write(18)var
      !
      !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
      !do n=1,nbrvar
      !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
        
      !end do
      !
    close(18)
    !
    print*,' << ',filename
    !
    deallocate(var)
    !
  end subroutine writetecbin2d14var
  !+-------------------------------------------------------------------+
  ! This subroutine is used to write writetecbin2d.                    |
  !+-------------------------------------------------------------------+
  !

  subroutine ecrirebin(unitf,string)
    !
    integer unitf
    integer i
    character*80 string
    !
    do i=1,len_trim(string)
      write(unitf)ichar(string(i:i))
    enddo
    write(unitf)0
    !
  end subroutine ecrirebin
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! this subroutine is used to write lay file to show TKE budget profiles using tecplot
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writebudgetlay(datfilename,layfilename)
    !
    character(len=*),intent(in) :: datfilename,layfilename
    !
    open(18,file=layfilename)
    write(18,'(A)')'#!MC 1410'
    write(18,'(6(A))')'$!VarSet |LFDSFN1| = ',"'",'"',datfilename,'"',"'"
    write(18,'(4(A))')'$!VarSet |LFDSVL1| = ',"'",'"V1" "V2" "V3" "V4" "V5" "V6" "V7" "V8" "V9" "V10" "V11" "V12"',"'"
    write(18,'(A)')"$!SETSTYLEBASE FACTORY"
    write(18,'(A)')"$!PLOTOPTIONS "
    write(18,'(A)')"  SUBDIVIDEALLCELLS = NO"
    write(18,'(A)')"$!GLOBALPAPER "
    write(18,'(A)')"  PAPERSIZEINFO"
    write(18,'(A)')"    {"
    write(18,'(A)')"    LETTER"
    write(18,'(A)')"      {"
    write(18,'(A)')"      WIDTH = 8.5"
    write(18,'(A)')"      HEIGHT = 11"
    write(18,'(A)')"      LEFTHARDCLIPOFFSET = 0.125"
    write(18,'(A)')"      RIGHTHARDCLIPOFFSET = 0.125"
    write(18,'(A)')"      TOPHARDCLIPOFFSET = 0.125"
    write(18,'(A)')"      BOTTOMHARDCLIPOFFSET = 0.125"
    write(18,'(A)')"      }"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!PAGE "
    write(18,'(A)')"  NAME = 'Untitled'"
    write(18,'(A)')"  PAPERATTRIBUTES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    BACKGROUNDCOLOR = WHITE"
    write(18,'(A)')"    ISTRANSPARENT = YES"
    write(18,'(A)')"    ORIENTPORTRAIT = NO"
    write(18,'(A)')"    SHOWGRID = YES"
    write(18,'(A)')"    SHOWRULER = NO"
    write(18,'(A)')"    SHOWPAPER = NO"
    write(18,'(A)')"    PAPERSIZE = LETTER"
    write(18,'(A)')"    RULERSPACING = ONEINCH"
    write(18,'(A)')"    PAPERGRIDSPACING = HALFINCH"
    write(18,'(A)')"    REGIONINWORKAREA"
    write(18,'(A)')"      {"
    write(18,'(A)')"      X1 = 1"
    write(18,'(A)')"      Y1 = 0.25"
    write(18,'(A)')"      X2 = 10"
    write(18,'(A)')"      Y2 = 8.25"
    write(18,'(A)')"      }"
    write(18,'(A)')"    }"
    write(18,'(A)')"### Frame Number 1 ###"
    write(18,'(A)')"$!READDATASET  '|LFDSFN1|'"
    write(18,'(A)')"  INITIALPLOTTYPE = XYLINE"
    write(18,'(A)')"  INCLUDETEXT = NO"
    write(18,'(A)')"  INCLUDEGEOM = NO"
    write(18,'(A)')"  ASSIGNSTRANDIDS = YES"
    write(18,'(A)')"  VARLOADMODE = BYNAME"
    write(18,'(A)')"  VARNAMELIST = '|LFDSVL1|'"
    write(18,'(A)')"$!REMOVEVAR |LFDSVL1|"
    write(18,'(A)')"$!REMOVEVAR |LFDSFN1|"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 3"
    write(18,'(A)')"  NAME = 'V3'"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 1"
    write(18,'(A)')"  NAME = 'y'"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 2"
    write(18,'(A)')"  NAME = 'yplus'"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 3"
    write(18,'(A)')"  NAME = 'convection'"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 4"
    write(18,'(A)')"  NAME = 'pressure_accelaration'"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 5"
    write(18,'(A)')"  NAME = 'pressure_dilatation'"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 6"
    write(18,'(A)')"  NAME = 'pressure_transport'"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 7"
    write(18,'(A)')"  NAME = 'production'"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 8"
    write(18,'(A)')"  NAME = 'turbulent_transport'"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 9"
    write(18,'(A)')"  NAME = 'viscous_accelaration'"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 10"
    write(18,'(A)')"  NAME = 'viscous_diffusion'"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 11"
    write(18,'(A)')"  NAME = 'dissipation'"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 12"
    write(18,'(A)')"  NAME = 'balance'"
    write(18,'(A)')"$!FRAMELAYOUT "
    write(18,'(A)')"  SHOWHEADER = NO"
    write(18,'(A)')"  HEADERCOLOR = RED"
    write(18,'(A)')"  XYPOS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    X = -1.3245"
    write(18,'(A)')"    Y = 0.51299"
    write(18,'(A)')"    }"
    write(18,'(A)')"  WIDTH = 12.02"
    write(18,'(A)')"  HEIGHT = 7.1601"
    write(18,'(A)')"$!THREEDAXIS "
    write(18,'(A)')"  ASPECTRATIOLIMIT = 25"
    write(18,'(A)')"  BOXASPECTRATIOLIMIT = 25"
    write(18,'(A)')"$!PLOTTYPE  = XYLINE"
    write(18,'(A)')"$!FRAMENAME  = 'Frame 001'"
    write(18,'(A)')"$!GLOBALTIME "
    write(18,'(A)')"  SOLUTIONTIME = 0"
    write(18,'(A)')"$!DELETELINEMAPS "
    write(18,'(A)')"$!ACTIVELINEMAPS  =  [2-11]"
    write(18,'(A)')"$!GLOBALLINEPLOT "
    write(18,'(A)')"  DATALABELS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    DISTANCESKIP = 5"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LEGEND"
    write(18,'(A)')"    {"
    write(18,'(A)')"    SHOW = YES"
    write(18,'(A)')"    BOX"
    write(18,'(A)')"      {"
    write(18,'(A)')"      BOXTYPE = NONE"
    write(18,'(A)')"      }"
    write(18,'(A)')"    XYPOS"
    write(18,'(A)')"      {"
    write(18,'(A)')"      X = 99.729"
    write(18,'(A)')"      Y = 100"
    write(18,'(A)')"      }"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [1]"
    write(18,'(A)')"  NAME = '&DV&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 1"
    write(18,'(A)')"    YAXISVAR = 2"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = RED"
    write(18,'(A)')"    LINETHICKNESS = 0.4"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = RED"
    write(18,'(A)')"    FILLCOLOR = RED"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = RED"
    write(18,'(A)')"    FILLCOLOR = RED"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = RED"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [2]"
    write(18,'(A)')"  NAME = '&DV&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 1"
    write(18,'(A)')"    YAXISVAR = 3"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = GREEN"
    write(18,'(A)')"    LINETHICKNESS = 0.4"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = GREEN"
    write(18,'(A)')"    FILLCOLOR = GREEN"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = GREEN"
    write(18,'(A)')"    FILLCOLOR = GREEN"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = GREEN"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [3]"
    write(18,'(A)')"  NAME = '&DV&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 1"
    write(18,'(A)')"    YAXISVAR = 4"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = BLUE"
    write(18,'(A)')"    LINETHICKNESS = 0.4"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = BLUE"
    write(18,'(A)')"    FILLCOLOR = BLUE"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = BLUE"
    write(18,'(A)')"    FILLCOLOR = BLUE"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = BLUE"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [4]"
    write(18,'(A)')"  NAME = '&DV&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 1"
    write(18,'(A)')"    YAXISVAR = 5"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM1"
    write(18,'(A)')"    LINETHICKNESS = 0.4"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM1"
    write(18,'(A)')"    FILLCOLOR = CUSTOM1"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM1"
    write(18,'(A)')"    FILLCOLOR = CUSTOM1"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM1"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [5]"
    write(18,'(A)')"  NAME = '&DV&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 1"
    write(18,'(A)')"    YAXISVAR = 6"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM3"
    write(18,'(A)')"    LINETHICKNESS = 0.4"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM3"
    write(18,'(A)')"    FILLCOLOR = CUSTOM3"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM3"
    write(18,'(A)')"    FILLCOLOR = CUSTOM3"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM3"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [6]"
    write(18,'(A)')"  NAME = '&DV&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 1"
    write(18,'(A)')"    YAXISVAR = 7"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = PURPLE"
    write(18,'(A)')"    LINETHICKNESS = 0.4"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = PURPLE"
    write(18,'(A)')"    FILLCOLOR = PURPLE"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = PURPLE"
    write(18,'(A)')"    FILLCOLOR = PURPLE"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = PURPLE"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [7]"
    write(18,'(A)')"  NAME = '&DV&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 1"
    write(18,'(A)')"    YAXISVAR = 8"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = RED"
    write(18,'(A)')"    LINETHICKNESS = 0.4"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = RED"
    write(18,'(A)')"    FILLCOLOR = RED"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = RED"
    write(18,'(A)')"    FILLCOLOR = RED"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = RED"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [8]"
    write(18,'(A)')"  NAME = '&DV&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 1"
    write(18,'(A)')"    YAXISVAR = 9"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM54"
    write(18,'(A)')"    LINETHICKNESS = 0.4"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = GREEN"
    write(18,'(A)')"    FILLCOLOR = GREEN"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = GREEN"
    write(18,'(A)')"    FILLCOLOR = GREEN"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = GREEN"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [9]"
    write(18,'(A)')"  NAME = '&DV&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 1"
    write(18,'(A)')"    YAXISVAR = 10"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM10"
    write(18,'(A)')"    LINETHICKNESS = 0.4"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = BLUE"
    write(18,'(A)')"    FILLCOLOR = BLUE"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = BLUE"
    write(18,'(A)')"    FILLCOLOR = BLUE"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = BLUE"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [10]"
    write(18,'(A)')"  NAME = '&DV&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 1"
    write(18,'(A)')"    YAXISVAR = 11"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CYAN"
    write(18,'(A)')"    LINETHICKNESS = 0.4"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM1"
    write(18,'(A)')"    FILLCOLOR = CUSTOM1"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM1"
    write(18,'(A)')"    FILLCOLOR = CUSTOM1"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM1"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [11]"
    write(18,'(A)')"  NAME = '&DV&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 1"
    write(18,'(A)')"    YAXISVAR = 12"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = BLACK"
    write(18,'(A)')"    LINETHICKNESS = 0.4"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM3"
    write(18,'(A)')"    FILLCOLOR = CUSTOM3"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM3"
    write(18,'(A)')"    FILLCOLOR = CUSTOM3"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM3"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!XYLINEAXIS "
    write(18,'(A)')"  DEPXTOYRATIO = 1"
    write(18,'(A)')"  VIEWPORTPOSITION"
    write(18,'(A)')"    {"
    write(18,'(A)')"    X1 = 8.7653"
    write(18,'(A)')"    Y1 = 10.171"
    write(18,'(A)')"    X2 = 98.657"
    write(18,'(A)')"    Y2 = 97.479"
    write(18,'(A)')"    }"
    write(18,'(A)')"  VIEWPORTTOPSNAPTARGET = 97.4786729858"
    write(18,'(A)')"$!XYLINEAXIS "
    write(18,'(A)')"  XDETAIL 1"
    write(18,'(A)')"    {"
    write(18,'(A)')"    RANGEMIN = -0.0014807059527370547"
    write(18,'(A)')"    RANGEMAX = 0.30148070595273707"
    write(18,'(A)')"    GRSPACING = 0.05"
    write(18,'(A)')"    TITLE"
    write(18,'(A)')"      {"
    write(18,'(A)')"      TITLEMODE = USETEXT"
    write(18,'(A)')"      TEXT = 'y'"
    write(18,'(A)')"      }"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!XYLINEAXIS "
    write(18,'(A)')"  YDETAIL 1"
    write(18,'(A)')"    {"
    write(18,'(A)')"    RANGEMIN = -0.31770759159266598"
    write(18,'(A)')"    RANGEMAX = 0.33621201984226351"
    write(18,'(A)')"    GRSPACING = 0.1"
    write(18,'(A)')"    TITLE"
    write(18,'(A)')"      {"
    write(18,'(A)')"      TITLEMODE = USETEXT"
    write(18,'(A)')"      TEXT = 'TKE budget'"
    write(18,'(A)')"      }"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!FRAMECONTROL ACTIVATEBYNUMBER"
    write(18,'(A)')"  FRAME = 1"
    write(18,'(A)')"$!SETSTYLEBASE CONFIG"
    close(18)
    !
    print*,' << ',layfilename,' link to file: ',datfilename
    !
  end subroutine writebudgetlay
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! end end of the subroutine writebudgetlay
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! this subroutine is used to write lay file to show uplus-yplus profile using tecplot
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writeupluslay(datfilename,layfilename)
    !
    character(len=*),intent(in) :: datfilename,layfilename
    !
    open(18,file=layfilename)
    write(18,'(A)')'#!MC 1410'
    write(18,'(6(A))')'$!VarSet |LFDSFN1| = ',"'",'"',datfilename,'"',"'"
    write(18,'(4(A))')'$!VarSet |LFDSVL1| = ',"'",'"V1" "V2"',"'"
    write(18,'(A)')"$!SETSTYLEBASE FACTORY"
    write(18,'(A)')"$!PLOTOPTIONS "
    write(18,'(A)')"  SUBDIVIDEALLCELLS = NO"
    write(18,'(A)')"$!GLOBALPAPER "
    write(18,'(A)')"  PAPERSIZEINFO"
    write(18,'(A)')"    {"
    write(18,'(A)')"    LETTER"
    write(18,'(A)')"      {"
    write(18,'(A)')"      WIDTH = 8.5"
    write(18,'(A)')"      HEIGHT = 11"
    write(18,'(A)')"      LEFTHARDCLIPOFFSET = 0.125"
    write(18,'(A)')"      RIGHTHARDCLIPOFFSET = 0.125"
    write(18,'(A)')"      TOPHARDCLIPOFFSET = 0.125"
    write(18,'(A)')"      BOTTOMHARDCLIPOFFSET = 0.125"
    write(18,'(A)')"      }"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!PAGE "
    write(18,'(A)')"  NAME = 'Untitled'"
    write(18,'(A)')"  PAPERATTRIBUTES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    BACKGROUNDCOLOR = WHITE"
    write(18,'(A)')"    ISTRANSPARENT = YES"
    write(18,'(A)')"    ORIENTPORTRAIT = NO"
    write(18,'(A)')"    SHOWGRID = YES"
    write(18,'(A)')"    SHOWRULER = NO"
    write(18,'(A)')"    SHOWPAPER = NO"
    write(18,'(A)')"    PAPERSIZE = LETTER"
    write(18,'(A)')"    RULERSPACING = ONEINCH"
    write(18,'(A)')"    PAPERGRIDSPACING = HALFINCH"
    write(18,'(A)')"    REGIONINWORKAREA"
    write(18,'(A)')"      {"
    write(18,'(A)')"      X1 = 1"
    write(18,'(A)')"      Y1 = 0.25"
    write(18,'(A)')"      X2 = 10"
    write(18,'(A)')"      Y2 = 8.25"
    write(18,'(A)')"      }"
    write(18,'(A)')"    }"
    write(18,'(A)')"### Frame Number 1 ###"
    write(18,'(A)')"$!READDATASET  '|LFDSFN1|'"
    write(18,'(A)')"  INITIALPLOTTYPE = XYLINE"
    write(18,'(A)')"  INCLUDETEXT = NO"
    write(18,'(A)')"  INCLUDEGEOM = NO"
    write(18,'(A)')"  ASSIGNSTRANDIDS = YES"
    write(18,'(A)')"  VARLOADMODE = BYNAME"
    write(18,'(A)')"  VARNAMELIST = '|LFDSVL1|'"
    write(18,'(A)')"$!REMOVEVAR |LFDSVL1|"
    write(18,'(A)')"$!REMOVEVAR |LFDSFN1|"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 1"
    write(18,'(A)')"  NAME = 'yplus'"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 2"
    write(18,'(A)')"  NAME = 'uplus'"
    write(18,'(A)')"$!ALTERDATA "
    write(18,'(A)')"  EQUATION = 'v3=v1'"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 3"
    write(18,'(A)')"  NAME = 'uplus'"
    write(18,'(A)')"$!ALTERDATA "
    write(18,'(A)')"  IRANGE"
    write(18,'(A)')"    {"
    write(18,'(A)')"    MIN = 2"
    write(18,'(A)')"    }"
    write(18,'(A)')"  EQUATION = 'v4=1/0.41*log(v1)+5.2'"
    write(18,'(A)')"$!FRAMELAYOUT "
    write(18,'(A)')"  SHOWHEADER = NO"
    write(18,'(A)')"  HEADERCOLOR = RED"
    write(18,'(A)')"  XYPOS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    X = -0.54401"
    write(18,'(A)')"    Y = 0.47906"
    write(18,'(A)')"    }"
    write(18,'(A)')"  WIDTH = 10.332"
    write(18,'(A)')"  HEIGHT = 7.245"
    write(18,'(A)')"$!THREEDAXIS "
    write(18,'(A)')"  ASPECTRATIOLIMIT = 25"
    write(18,'(A)')"  BOXASPECTRATIOLIMIT = 25"
    write(18,'(A)')"$!PLOTTYPE  = XYLINE"
    write(18,'(A)')"$!FRAMENAME  = 'Frame 001'"
    write(18,'(A)')"$!GLOBALTIME "
    write(18,'(A)')"  SOLUTIONTIME = 0"
    write(18,'(A)')"$!DELETELINEMAPS "
    write(18,'(A)')"$!ACTIVELINEMAPS  =  [1-3]"
    write(18,'(A)')"$!GLOBALLINEPLOT "
    write(18,'(A)')"  DATALABELS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    DISTANCESKIP = 5"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LEGEND"
    write(18,'(A)')"    {"
    write(18,'(A)')"    XYPOS"
    write(18,'(A)')"      {"
    write(18,'(A)')"      X = 95"
    write(18,'(A)')"      }"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [1]"
    write(18,'(A)')"  NAME = '&ZN&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 1"
    write(18,'(A)')"    YAXISVAR = 2"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = RED"
    write(18,'(A)')"    LINETHICKNESS = 0.8"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = RED"
    write(18,'(A)')"    FILLCOLOR = RED"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = RED"
    write(18,'(A)')"    FILLCOLOR = RED"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = RED"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [2]"
    write(18,'(A)')"  NAME = '&ZN&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 1"
    write(18,'(A)')"    YAXISVAR = 3"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = GREEN"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = GREEN"
    write(18,'(A)')"    FILLCOLOR = GREEN"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = GREEN"
    write(18,'(A)')"    FILLCOLOR = GREEN"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = GREEN"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [3]"
    write(18,'(A)')"  NAME = '&ZN&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 1"
    write(18,'(A)')"    YAXISVAR = 4"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = BLUE"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = BLUE"
    write(18,'(A)')"    FILLCOLOR = BLUE"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = BLUE"
    write(18,'(A)')"    FILLCOLOR = BLUE"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = BLUE"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!XYLINEAXIS "
    write(18,'(A)')"  DEPXTOYRATIO = 1"
    write(18,'(A)')"  VIEWPORTPOSITION"
    write(18,'(A)')"    {"
    write(18,'(A)')"    X1 = 9.6335"
    write(18,'(A)')"    Y1 = 9.3607"
    write(18,'(A)')"    X2 = 97.607"
    write(18,'(A)')"    Y2 = 94.323"
    write(18,'(A)')"    }"
    write(18,'(A)')"  VIEWPORTTOPSNAPTARGET = 94.3231850117"
    write(18,'(A)')"$!XYLINEAXIS "
    write(18,'(A)')"  XDETAIL 1"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COORDSCALE = LOG"
    write(18,'(A)')"    RANGEMIN = 0.50201039584694418"
    write(18,'(A)')"    RANGEMAX = 4008.8788624500326"
    write(18,'(A)')"    GRSPACING = 500"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!XYLINEAXIS "
    write(18,'(A)')"  YDETAIL 1"
    write(18,'(A)')"    {"
    write(18,'(A)')"    RANGEMIN = -0.22252724804847054"
    write(18,'(A)')"    RANGEMAX = 23.163585021279445"
    write(18,'(A)')"    GRSPACING = 5"
    write(18,'(A)')"    TITLE"
    write(18,'(A)')"      {"
    write(18,'(A)')"      TITLEMODE = USETEXT"
    write(18,'(A)')"      TEXT = 'uplus'"
    write(18,'(A)')"      }"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!ATTACHTEXT "
    write(18,'(A)')"  ANCHORPOS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    X = 46.83840534549869"
    write(18,'(A)')"    Y = 87.11955632545789"
    write(18,'(A)')"    }"
    write(18,'(A)')"  TEXTSHAPE"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ISBOLD = NO"
    write(18,'(A)')"    }"
    write(18,'(A)')"  TEXT = 'u+=y+'"
    write(18,'(A)')"$!ATTACHTEXT "
    write(18,'(A)')"  ANCHORPOS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    X = 67.28369387655397"
    write(18,'(A)')"    Y = 71.42877019804476"
    write(18,'(A)')"    }"
    write(18,'(A)')"  TEXTSHAPE"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ISBOLD = NO"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ANGLE = 30"
    write(18,'(A)')"  TEXT = 'u+=1/0.41*ln(y+)+5.2'"
    write(18,'(A)')"$!FRAMECONTROL ACTIVATEBYNUMBER"
    write(18,'(A)')"  FRAME = 1"
    write(18,'(A)')"$!SETSTYLEBASE CONFIG"
    close(18)
    !
    print*,' << ',layfilename,' link to file: ',datfilename
    !
  end subroutine writeupluslay
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! end end of the subroutine writebudgetlay
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! this subroutine is used to write .lay to show Reynolds stress profiles using tecplot
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writerestresslay(datfilename,layfilename)
    !
    character(len=*),intent(in) :: datfilename,layfilename
    !
    open(18,file=layfilename)
    write(18,'(A)')'#!MC 1410'
    write(18,'(6(A))')'$!VarSet |LFDSFN1| = ',"'",'"',datfilename,'"',"'"
    write(18,'(4(A))')'$!VarSet |LFDSVL1| = ',"'",'"V1" "V2" "V3" "V4" "V5" "V6" "V7" "V8"',"'"
    write(18,'(A)')"$!SETSTYLEBASE FACTORY"
    write(18,'(A)')"$!PLOTOPTIONS "
    write(18,'(A)')"  SUBDIVIDEALLCELLS = NO"
    write(18,'(A)')"$!GLOBALPAPER "
    write(18,'(A)')"  PAPERSIZEINFO"
    write(18,'(A)')"    {"
    write(18,'(A)')"    LETTER"
    write(18,'(A)')"      {"
    write(18,'(A)')"      WIDTH = 8.5"
    write(18,'(A)')"      HEIGHT = 11"
    write(18,'(A)')"      LEFTHARDCLIPOFFSET = 0.125"
    write(18,'(A)')"      RIGHTHARDCLIPOFFSET = 0.125"
    write(18,'(A)')"      TOPHARDCLIPOFFSET = 0.125"
    write(18,'(A)')"      BOTTOMHARDCLIPOFFSET = 0.125"
    write(18,'(A)')"      }"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!PAGE "
    write(18,'(A)')"  NAME = 'Untitled'"
    write(18,'(A)')"  PAPERATTRIBUTES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    BACKGROUNDCOLOR = WHITE"
    write(18,'(A)')"    ISTRANSPARENT = YES"
    write(18,'(A)')"    ORIENTPORTRAIT = NO"
    write(18,'(A)')"    SHOWGRID = YES"
    write(18,'(A)')"    SHOWRULER = NO"
    write(18,'(A)')"    SHOWPAPER = NO"
    write(18,'(A)')"    PAPERSIZE = LETTER"
    write(18,'(A)')"    RULERSPACING = ONEINCH"
    write(18,'(A)')"    PAPERGRIDSPACING = HALFINCH"
    write(18,'(A)')"    REGIONINWORKAREA"
    write(18,'(A)')"      {"
    write(18,'(A)')"      X1 = 1"
    write(18,'(A)')"      Y1 = 0.25"
    write(18,'(A)')"      X2 = 10"
    write(18,'(A)')"      Y2 = 8.25"
    write(18,'(A)')"      }"
    write(18,'(A)')"    }"
    write(18,'(A)')"### Frame Number 1 ###"
    write(18,'(A)')"$!READDATASET  '|LFDSFN1|'"
    write(18,'(A)')"  INITIALPLOTTYPE = XYLINE"
    write(18,'(A)')"  INCLUDETEXT = NO"
    write(18,'(A)')"  INCLUDEGEOM = NO"
    write(18,'(A)')"  ASSIGNSTRANDIDS = YES"
    write(18,'(A)')"  VARLOADMODE = BYNAME"
    write(18,'(A)')"  VARNAMELIST = '|LFDSVL1|'"
    write(18,'(A)')"$!REMOVEVAR |LFDSVL1|"
    write(18,'(A)')"$!REMOVEVAR |LFDSFN1|"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 1"
    write(18,'(A)')"  NAME = 'y'"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 2"
    write(18,'(A)')"  NAME = 'yplus'"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 3"
    write(18,'(A)')"  NAME = 'uu'"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 4"
    write(18,'(A)')"  NAME = 'vv'"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 5"
    write(18,'(A)')"  NAME = 'ww'"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 6"
    write(18,'(A)')"  NAME = 'TKE'"
    write(18,'(A)')"$!RENAMEDATASETVAR "
    write(18,'(A)')"  VAR = 7"
    write(18,'(A)')"  NAME = 'uv'"
    write(18,'(A)')"$!FRAMELAYOUT "
    write(18,'(A)')"  SHOWHEADER = NO"
    write(18,'(A)')"  HEADERCOLOR = RED"
    write(18,'(A)')"  XYPOS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    X = -1.333"
    write(18,'(A)')"    Y = 0.33484"
    write(18,'(A)')"    }"
    write(18,'(A)')"  WIDTH = 13.428"
    write(18,'(A)')"  HEIGHT = 7.8134"
    write(18,'(A)')"$!THREEDAXIS "
    write(18,'(A)')"  ASPECTRATIOLIMIT = 25"
    write(18,'(A)')"  BOXASPECTRATIOLIMIT = 25"
    write(18,'(A)')"$!PLOTTYPE  = XYLINE"
    write(18,'(A)')"$!FRAMENAME  = 'Frame 001'"
    write(18,'(A)')"$!GLOBALTIME "
    write(18,'(A)')"  SOLUTIONTIME = 0"
    write(18,'(A)')"$!DELETELINEMAPS "
    write(18,'(A)')"$!ACTIVELINEMAPS  =  [2-6]"
    write(18,'(A)')"$!GLOBALLINEPLOT "
    write(18,'(A)')"  DATALABELS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    DISTANCESKIP = 5"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LEGEND"
    write(18,'(A)')"    {"
    write(18,'(A)')"    SHOW = YES"
    write(18,'(A)')"    TEXTSHAPE"
    write(18,'(A)')"      {"
    write(18,'(A)')"      HEIGHT = 4"
    write(18,'(A)')"      }"
    write(18,'(A)')"    BOX"
    write(18,'(A)')"      {"
    write(18,'(A)')"      BOXTYPE = NONE"
    write(18,'(A)')"      }"
    write(18,'(A)')"    XYPOS"
    write(18,'(A)')"      {"
    write(18,'(A)')"      X = 98.691"
    write(18,'(A)')"      Y = 98.132"
    write(18,'(A)')"      }"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [1]"
    write(18,'(A)')"  NAME = '&DV&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 2"
    write(18,'(A)')"    YAXISVAR = 2"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = RED"
    write(18,'(A)')"    LINETHICKNESS = 0.8"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = RED"
    write(18,'(A)')"    FILLCOLOR = RED"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = RED"
    write(18,'(A)')"    FILLCOLOR = RED"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = RED"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [2]"
    write(18,'(A)')"  NAME = '&DV&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 2"
    write(18,'(A)')"    YAXISVAR = 3"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = GREEN"
    write(18,'(A)')"    LINETHICKNESS = 0.8"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = GREEN"
    write(18,'(A)')"    FILLCOLOR = GREEN"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = GREEN"
    write(18,'(A)')"    FILLCOLOR = GREEN"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = GREEN"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [3]"
    write(18,'(A)')"  NAME = '&DV&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 2"
    write(18,'(A)')"    YAXISVAR = 4"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = BLUE"
    write(18,'(A)')"    LINETHICKNESS = 0.8"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = BLUE"
    write(18,'(A)')"    FILLCOLOR = BLUE"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = BLUE"
    write(18,'(A)')"    FILLCOLOR = BLUE"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = BLUE"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [4]"
    write(18,'(A)')"  NAME = '&DV&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 2"
    write(18,'(A)')"    YAXISVAR = 5"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM1"
    write(18,'(A)')"    LINETHICKNESS = 0.8"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM1"
    write(18,'(A)')"    FILLCOLOR = CUSTOM1"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM1"
    write(18,'(A)')"    FILLCOLOR = CUSTOM1"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM1"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [5]"
    write(18,'(A)')"  NAME = '&DV&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 2"
    write(18,'(A)')"    YAXISVAR = 6"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM3"
    write(18,'(A)')"    LINETHICKNESS = 0.8"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM3"
    write(18,'(A)')"    FILLCOLOR = CUSTOM3"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM3"
    write(18,'(A)')"    FILLCOLOR = CUSTOM3"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = CUSTOM3"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [6]"
    write(18,'(A)')"  NAME = '&DV&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 2"
    write(18,'(A)')"    YAXISVAR = 7"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = PURPLE"
    write(18,'(A)')"    LINETHICKNESS = 0.8"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = PURPLE"
    write(18,'(A)')"    FILLCOLOR = PURPLE"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = PURPLE"
    write(18,'(A)')"    FILLCOLOR = PURPLE"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = PURPLE"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [7]"
    write(18,'(A)')"  NAME = '&DV&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 2"
    write(18,'(A)')"    YAXISVAR = 8"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = RED"
    write(18,'(A)')"    LINETHICKNESS = 0.8"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = RED"
    write(18,'(A)')"    FILLCOLOR = RED"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = RED"
    write(18,'(A)')"    FILLCOLOR = RED"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = RED"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!LINEMAP  [8]"
    write(18,'(A)')"  NAME = '&DV&'"
    write(18,'(A)')"  ASSIGN"
    write(18,'(A)')"    {"
    write(18,'(A)')"    ZONE = 1"
    write(18,'(A)')"    XAXISVAR = 2"
    write(18,'(A)')"    YAXISVAR = 9"
    write(18,'(A)')"    }"
    write(18,'(A)')"  LINES"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = GREEN"
    write(18,'(A)')"    LINETHICKNESS = 0.8"
    write(18,'(A)')"    }"
    write(18,'(A)')"  SYMBOLS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = GREEN"
    write(18,'(A)')"    FILLCOLOR = GREEN"
    write(18,'(A)')"    }"
    write(18,'(A)')"  BARCHARTS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = GREEN"
    write(18,'(A)')"    FILLCOLOR = GREEN"
    write(18,'(A)')"    }"
    write(18,'(A)')"  ERRORBARS"
    write(18,'(A)')"    {"
    write(18,'(A)')"    COLOR = GREEN"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!XYLINEAXIS "
    write(18,'(A)')"  DEPXTOYRATIO = 1"
    write(18,'(A)')"  VIEWPORTPOSITION"
    write(18,'(A)')"    {"
    write(18,'(A)')"    X1 = 10.423"
    write(18,'(A)')"    Y1 = 9.8056"
    write(18,'(A)')"    X2 = 97.681"
    write(18,'(A)')"    Y2 = 97.555"
    write(18,'(A)')"    }"
    write(18,'(A)')"  VIEWPORTTOPSNAPTARGET = 97.5548317047"
    write(18,'(A)')"$!XYLINEAXIS "
    write(18,'(A)')"  XDETAIL 1"
    write(18,'(A)')"    {"
    write(18,'(A)')"    RANGEMIN = 0"
    write(18,'(A)')"    RANGEMAX = 250"
    write(18,'(A)')"    GRSPACING = 50"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!XYLINEAXIS "
    write(18,'(A)')"  YDETAIL 1"
    write(18,'(A)')"    {"
    write(18,'(A)')"    RANGEMIN = -1.5140145566044976"
    write(18,'(A)')"    RANGEMAX = 8.5981075602747783"
    write(18,'(A)')"    GRSPACING = 2"
    write(18,'(A)')"    TITLE"
    write(18,'(A)')"      {"
    write(18,'(A)')"      TITLEMODE = USETEXT"
    write(18,'(A)')"      TEXT = 'Reynolds Stress'"
    write(18,'(A)')"      }"
    write(18,'(A)')"    }"
    write(18,'(A)')"$!FRAMECONTROL ACTIVATEBYNUMBER"
    write(18,'(A)')"  FRAME = 1"
    write(18,'(A)')"$!SETSTYLEBASE CONFIG"
    close(18)
    !
    print*,' << ',layfilename,' link to file: ',datfilename
    !
  end subroutine writerestresslay
  !
end module writetec
!+---------------------------------------------------------------------+
!|The end of the module WriteTec                                       |
!+---------------------------------------------------------------------+
