!+---------------------------------------------------------------------+
!| This module contains subroutines needed for read control files.     |
!+---------------------------------------------------------------------+
module controlfileread
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This is the subroutine that read commonds from keyboard           |
  !+-------------------------------------------------------------------+
  subroutine readkeyboard(npost,command,number1,number2,number3,      &
                                     filename1,filename2,nthread,inputfile)
    !
    integer,intent(out) :: npost
    character(len=*),intent(out) :: command
    integer,intent(out) :: nthread
    !
    integer,intent(out),optional :: number1,number2,number3
    character(len=*),intent(out),optional :: filename1,filename2,inputfile
    !
    integer :: i,nargu,ierr,cli_len,cli_len1,cli_len2,cli_len3,cli_len4
    character(len=255) :: cli_arg
    !
    nargu = command_argument_count()
    !
    npost=0
    !
    if(nargu>0) then
      !
      call get_command_argument(1,command,cli_len,ierr)
      !
      if(trim(command) == 'wing') then
        npost=-5
      elseif(trim(command) == 'scheme') then
        npost=-4
      elseif(trim(command) == 'test') then
        npost=-3
      elseif(trim(command) == 'gridgen') then
        npost=-1
      elseif(trim(command) == 'flowgen') then
        npost=-2
      elseif(trim(command) == 'notdefined') then
        npost=1
      elseif(trim(command) == 'tec3d') then
        npost=2
      elseif(trim(command) == 'tecxy') then
        npost=211
      elseif(trim(command) == 'tecxz') then
        npost=212
      elseif(trim(command) == 'tecyz') then
        npost=213
      elseif(trim(command) == 'teccs') then
        npost=210
      elseif(trim(command) == 'prof') then
        npost=214
      elseif(trim(command) == 'techan') then
        npost=215
      elseif(trim(command) == 'writetec2d') then
        npost=202
      elseif(trim(command) == 'monitor') then
        npost=204
      elseif(trim(command) == 'sepzone') then
        npost=205
      elseif(trim(command) == 'shock') then
        npost=206
      elseif(trim(command) == 'meanflow') then
        npost=3
      elseif(trim(command) == 'chanmean') then
        npost=301
      elseif(trim(command) == 'blmean') then
        npost=302
      elseif(trim(command) == 'diffmean') then
        npost=303
      elseif(trim(command) == 'cylinder') then
        npost=304
      elseif(trim(command) == 'tlv') then
        npost=305
      elseif(trim(command) == 'streamline') then
        npost=306
      elseif(trim(command) == 'cavity') then
        npost=307
      elseif(trim(command) == 'budget') then
        npost=4
      elseif(trim(command) == 'budgetana') then
        npost=401
      elseif(trim(command) == 'flowinterp') then
        npost=5
      elseif(trim(command) == 'inletgen') then
        npost=501
      elseif(trim(command) == 'flowrenorm') then
        npost=502
      elseif(trim(command) == '2pointcorre') then
        npost=6
      elseif(trim(command) == 'spectra') then
        npost=7
      elseif(trim(command) == 'specout') then
        npost=701
      elseif(trim(command) == 'specins') then
        npost=702
      elseif(trim(command) == 'stokes') then
        npost=8
      elseif(trim(command) == 'kslice') then
        npost=801
      elseif(trim(command) == 'islice') then
        npost=802
      elseif(trim(command) == 'jslice') then
        npost=803
      elseif(trim(command) == 'stasdatacon') then
        npost=9
      elseif(trim(command) == 'islicecon') then
        npost=901
      elseif(trim(command) == 'kslicecon') then
        npost=902
      elseif(trim(command) == 'meanflowcon') then
        npost=903
      elseif(trim(command) == 'flowfieldcon') then
        npost=904
      elseif(trim(command) == '3dto2d') then
        npost=905
      elseif(trim(command) == 'supplycon') then
        npost=906
      elseif(trim(command) == 'invcon') then
        npost=907
      elseif(trim(command) == 'jslicecon') then
        npost=908
      elseif(trim(command) == 'udf') then
        npost=1000
      elseif(trim(command) == 'hipstar') then
        npost=1001
      elseif(trim(command) == 'swtbli') then
        npost=1002
      elseif(trim(command) == 'inflowgen') then
        npost=1003
      elseif(trim(command) == 'circup') then
        npost=1004
      elseif(trim(command) == 'upm') then
        npost=1005
      elseif(trim(command) == 'tecscript') then
        npost=1006
      elseif(trim(command) == 'particle') then
        npost=1007
      elseif(trim(command) == 'incompact') then
        npost=1008
      else
        print*,'!! post command not recognized !!'
        call printcommand
      endif
      !
      number1=-1
      number2=-1
      number3=-1
      !
      filename1=' '
      filename2=' '
      !
      nthread=1
      inputfile='datin/input.3dp.dat'
      !
      print*,' ** nargu',nargu
      !
      if(nargu>=2) then
        !
        i=1
        do while(i<nargu)
          !
          i=i+1
          !
          call get_command_argument(i,cli_arg,cli_len,ierr)
          !
          if(isnum(trim(cli_arg))==1) then
            !
            ! the input is a number
            if(number1==-1) then
              read(cli_arg,*) number1
            elseif(number2==-1) then
              read(cli_arg,*) number2
            elseif(number3==-1) then
              read(cli_arg,*) number3
            else
              print*,' !! no more number input allowed !!'
            endif
            !
          elseif(trim(cli_arg)=='-nt') then
            ! the input is number of thread
            i=i+1
            call get_command_argument(i,cli_arg,cli_len,ierr)
            read(cli_arg,*) nthread
          elseif(trim(cli_arg)=='-input') then
            ! the input is name of input file
            i=i+1
            call get_command_argument(i,cli_arg,cli_len,ierr)
            inputfile=trim(cli_arg)
          else
            !
            ! the input is strings
            if(filename1==' ') then
              filename1=trim(cli_arg)
            elseif(filename2==' ') then
              filename2=trim(cli_arg)
            endif
            !
          endif
          !
          ! if(isnum(trim(cli_arg(i)))==1) then
          !   ! the input is a number
          !   if(number1==-1) then
          !     read(cli_arg(i),*) number1
          !   elseif(number2==-1) then
          !     read(cli_arg(i),*) number2
          !   elseif(number3==-1) then
          !     read(cli_arg(i),*) number3
          !   else
          !     stop ' !! error in get number1/2/3 !!'
          !   endif
          ! elseif(isnum(trim(cli_arg(i)))==0) then
          !   if(filename1==' ') then
          !     filename1=trim(cli_arg(i))
          !   elseif(filename2==' ') then
          !     filename2=trim(cli_arg(i))
          !   else
          !     print*,cli_arg,filename2
          !     stop ' !! error in get filename1/2 !!'
          !   endif
          !   if(trim(cli_arg(i))=='-nt') then
          !     i=i+1
          !     call get_command_argument(i,cli_arg(i),cli_len,ierr)
          !     read(cli_arg(i),*) nthread
          !   else
          !     nthread=1
          !   endif
          !   if(present(inputfile) .and. trim(cli_arg(i))=='-input') then
          !     i=i+1
          !     call get_command_argument(i,cli_arg(i),cli_len,ierr)
          !     inputfile=trim(cli_arg(i))
          !   endif
          ! else
          !   stop ' !! error in readkeyboard !! '
          ! endif
          !
        enddo
        !
      endif
      !
    else
      print*,'!! please input the command for postprocess !!'
    endif
    !
    print*,'=========================readkeyboard========================'
    write(*,'(A18,A15)')'  Repeat command:',trim(command)
    write(*,'(A18,I15)')'           npost:',npost
    !
    print*,'=========================readkeyboard========================'
    print*
    !
  end subroutine readkeyboard
  !+-------------------------------------------------------------------+
  !| The end of the subroutine readkeyboard                            |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to print commands.                        |
  !+-------------------------------------------------------------------+
  subroutine printcommand
    !
    print*,'==================================================='
    print*,'           supportted command line input           '
    print*,'------------+-----+-------------------------------+'
    print*,'     command|npost|                      arguments|'
    print*,'------------+-----+-------------------------------+'
    print*,'        wing|   -5|                No input needed|'
    print*,'------------+-----+-------------------------------+'
    print*,'      scheme|   -4|                No input needed|'
    print*,'------------+-----+-------------------------------+'
    print*,'        test|   -3|                No input needed|'
    print*,'------------+-----+-------------------------------+'
    print*,'     gridana|   -2|                No input needed|'
    print*,'------------+-----+-------------------------------+'
    print*,'     flowgen|   -1|                No input needed|'
    print*,'------------+-----+-------------------------------+'
    print*,'       tec3d|    2|     file input|    file output|'
    print*,'------------+-----+---------------+---------------+'
    print*,'       tecxy|  211|     file input|        k slice|'
    print*,'------------+-----+---------------+---------------+'
    print*,'       tecxz|  212|     file input|        j slice|'
    print*,'------------+-----+---------------+---------------+'
    print*,'       tecyz|  213|     file input|        i slice|'
    print*,'            |     +---------------+---------------+---------------+'
    print*,'            |     |     1st flow #|    last flow #|        i slice|'
    print*,'            +-----+---------------+---------------+---------------+'
    print*,'       teccs|  210|     file input|    file output|'
    print*,'            |     +---------------+---------------+'
    print*,'            |     |     1st flow #|    last flow #|'
    print*,'------------+-----+---------------+---------------+'
    print*,'        prof|  214|     file input|    file output|'
    print*,'------------+-----+---------------+---------------+'
    print*,'      techan|  215|     file input|    file output|'
    print*,'------------+-----+---------------+---------------+'
    print*,'  writetec2d|  202|     1st flow #|    last flow #|'
    print*,'------------+-----+---------------+---------------+'
    print*,' instflowsta|  203|                   input needed|'
    print*,'------------+-----+---------------+---------------+'
    print*,'     monitor|  204|     # monitors|  step of stat.|'
    print*,'------------+-----+---------------+---------------+'
    print*,'     sepzone|  205|     1st flow #|    last flow #|'
    print*,'------------+-----+---------------+---------------+'
    print*,'       shock|  206|                               |'
    print*,'------------+-----+---------------+---------------+'
    print*,'    meanflow|    3|     No input if single h5 file|'
    print*,'            |     +---------------+---------------+'
    print*,'            |     |     1st flow #|    last flow #|'
    print*,'------------+-----+---------------+---------------+'
    print*,'    chanmean|  301|                No input needed|'
    print*,'------------+-----+---------------+---------------+'
    print*,'      blmean|  302|             i slice of profile|'
    print*,'------------+-----+---------------+---------------+'
    print*,'    diffmean|  303|                No input needed|'
    print*,'------------+-----+---------------+---------------+'
    print*,'    cylinder|  304|                No input needed|'
    print*,'------------+-----+---------------+---------------+'
    print*,'         tlv|  305|                No input needed|'
    print*,'------------+-----+---------------+---------------+'
    print*,'  streamline|  306|  the first point of streamline|'
    print*,'------------+-----+---------------+---------------+'
    print*,'      cavity|  307|                               |'
    print*,'------------+-----+---------------+---------------+'
    print*,'      budget|    4|                No input needed|'
    print*,'------------+-----+---------------+---------------+'
    print*,'   budgetana|  401|                No input needed|'
    print*,'------------+-----+---------------+---------------+'
    print*,'  flowinterp|    5|                No input needed|'
    print*,'------------+-----+---------------+---------------+'
    print*,'    inletgen|  501|                No input needed|'
    print*,'------------+-----+---------------+---------------+'
    print*,'  flowrenorm|  502|                No input needed|'
    print*,'------------+-----+---------------+---------------+'
    print*,' 2pointcorre|    6|     1st flow #|    last flow #|'
    print*,'------------+-----+---------------+---------------+'
    print*,'     spectra|    7|     1st flow #|    last flow #|'
    print*,'------------+-----+---------------+---------------+'
    print*,'     specout|  701|                No input needed|'
    print*,'------------+-----+---------------+---------------+'
    print*,'     specins|  702|     1st flow #|    last flow #|'
    print*,'------------+-----+---------------+---------------+'
    print*,'      stokes|    8|     1st flow #|    last flow #|'
    print*,'------------+-----+---------------+---------------+---------------+'
    print*,'      kslice|  801|     1st flow #|    last flow #|         skip #|'
    print*,'------------+-----+---------------+---------------+---------------+'
    print*,'      islice|  802|     1st flow #|    last flow #|         skip #|'
    print*,'------------+-----+---------------+---------------+---------------+'
    print*,'      jslice|  803|     1st flow #|    last flow #|         skip #|'
    print*,'------------+-----+---------------+---------------+---------------+'
    print*,' stasdatacon|    9|                No input needed|'
    print*,'------------+-----+---------------+---------------+'
    print*,'   islicecon|  901|     1st flow #|    last flow #|'
    print*,'------------+-----+---------------+---------------+'
    print*,'   kslicecon|  902|     1st flow #|    last flow #|'
    print*,'------------+-----+---------------+---------------+'
    print*,' meanflowcon|  903|                No input needed|'
    print*,'------------+-----+---------------+---------------+'
    print*,'flowfieldcon|  904|     1st flow #|    last flow #|'
    print*,'------------+-----+---------------+---------------+'
    print*,'      3dto2d|  905|                No input needed|'
    print*,'------------+-----+---------------+---------------+'
    print*,'   supplycon|  906|                                '
    print*,'------------+-----+---------------+---------------+'
    print*,'      invcon|  907|       the number of processoer|'
    print*,'------------+-----+---------------+---------------+'
    print*,'   jslicecon|  908|     1st flow #|    last flow #|'
    print*,'------------+-----+---------------+---------------+'
    print*,'         udf| 1000|                    user define|'
    print*,'------------+-----+---------------+---------------+'
    print*,'     hipstar| 1001|                No input needed|'
    print*,'------------+-----+---------------+---------------+'
    print*,'     hipstar| 1001|                No input needed|'
    print*,'------------+-----+---------------+---------------+'
    print*,'      swtbli| 1002|                No input needed|'
    print*,'------------+-----+---------------+---------------+'
    print*,'   inflowgen| 1003|                No input needed|'
    print*,'------------+-----+---------------+---------------+'
    print*,'      circup| 1004|                No input needed|'
    print*,'------------+-----+---------------+---------------+'
    print*,'         upm| 1005|                No input needed|'
    print*,'------------+-----+---------------+---------------+---------------+'
    print*,'   tecscript| 1006|   1st number #|   2nd number #| in file name #|'
    print*,'------------+-----+---------------+---------------+---------------+'
    print*,'    particle| 1007|                No input needed|'
    print*,'------------+-----+---------------+---------------+'
    print*,'   incompact| 1008|                No input needed|'
    print*,'============+================================================'
    !
    stop
    !
  end subroutine printcommand
  !+-------------------------------------------------------------------+
  !| The end of the subroutine printcommand                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !|Description. Returns the entire command by which the program was   |
  !|invoked.                                                           |
  !|                                                                   |
  !|Class. Subroutine.                                                 |
  !|                                                                   |
  !|Arguments.                                                         |
  !|COMMAND (optional) shall be scalar and of type default character.  |
  !|  It is an INTENT(OUT) argument. It is assigned the entire command |
  !|  by which the program was invoked. If the command cannot be       |
  !|  determined, COMMAND is assigned all blanks.                      |
  !|LENGTH (optional) shall be scalar and of type default integer. It  |
  !|  is an INTENT(OUT) argument. It is assigned the significant length|
  !|  of the command by which the program was invoked. The significant |
  !|  length may include trailing blanks if the processor allows       |
  !|  commands with significant trailing blanks. This length does not  |
  !|  consider any possible truncation or padding in assigning the     |
  !|  command to the COMMAND argument; in fact the COMMAND argument    |
  !|  need not even be present. If the command length cannot be        |
  !|  determined, a length of 0 is assigned.                           |
  !|STATUS (optional) shall be scalar and of type default integer. It  |
  !|  is an INTENT(OUT) argument. It is assigned the value 0 if the    |
  !|  command retrieval is sucessful. It is assigned a                 |
  !|  processor-dependent non-zero value if the command retrieval fails|
  !+-------------------------------------------------------------------+
  subroutine get_command(command,length,status)
    !
    character(len=*), intent(out), optional :: command
    integer         , intent(out), optional :: length
    integer         , intent(out), optional :: status
    !
    integer                   :: iarg,narg,ipos
    integer            , save :: lenarg
    character(len=2000), save :: argstr
    logical            , save :: getcmd = .true.
    !
    integer :: iargc
    !
    ! Under Unix we must reconstruct the command line from its 
    ! constituent parts. This will not be the original command line.
    ! Rather it will be the expanded command line as generated by 
    ! the shell.
    if (getcmd) then
      narg = iargc()
      if (narg > 0) then
        ipos = 1
        do iarg = 1,narg
          call getarg(iarg,argstr(ipos:))
          lenarg = len_trim(argstr)
          ipos   = lenarg + 2
          if (ipos > len(argstr)) exit
        end do
      else
        argstr = ' '
        lenarg = 0
      endif
      getcmd = .false.
    endif
    if (present(command)) command = argstr
    if (present(length))  length  = lenarg
    if (present(status))  status  = 0
    return
  end subroutine get_command
  !+-------------------------------------------------------------------+
  !| The end of the subroutine get_command                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !|Description. Returns the number of command arguments.              |
  !|                                                                   |
  !|Class. Inquiry function                                            |
  !|                                                                   |
  !|Arguments. None.                                                   |
  !|                                                                   |
  !|Result Characteristics. Scalar default integer.                    |
  !|                                                                   |
  !|Result Value. The result value is equal to the number of command   |
  !|  arguments available. If there are no command arguments available |
  !|  or if the processor does not support command arguments, then     |
  !|  the result value is 0. If the processor has a concept of a       |
  !|  command name, the command name does not count as one of the      |
  !|  command arguments.                                               |
  !|                                                                   |
  !|The following INTEGER/EXTERNAL declarations of IARGC should not    |
  !|really be necessary. However, at least one compiler (PGI) comments |
  !|on their absence, so they are included for completeness.           |
  !+-------------------------------------------------------------------+
  integer function command_argument_count()
    !
    integer :: iargc
    !
    command_argument_count = iargc()
    return
    !
  end function command_argument_count
  !+-------------------------------------------------------------------+
  !| The end of the function command_argument_count                    |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !|Description. Returns a command argument.                           |
  !|                                                                   |
  !|Class. Subroutine.                                                 |
  !|                                                                   |
  !|Arguments.                                                         |
  !|NUMBER shall be scalar and of type default integer. It is an       |
  !|  INTENT(IN) argument. It specifies the number of the command      |
  !|  argument that the other arguments give information about. Useful |
  !|  values of NUMBER are those between 0 and the argument count      |
  !|  returned by the COMMAND_ARGUMENT_COUNT intrinsic.                |
  !|  Other values are allowed, but will result in error status return |
  !|  (see below).  Command argument 0 is defined to be the command    |
  !|  name by which the program was invoked if the processor has such  |
  !|  a concept. It is allowed to call the GET_COMMAND_ARGUMENT        |
  !|  procedure for command argument number 0, even if the processor   |
  !|  does not define command names or other command arguments.        |
  !|  The remaining command arguments are numbered consecutively from  |
  !|  1 to the argument count in an order determined by the processor. |
  !|VALUE (optional) shall be scalar and of type default character.    |
  !|  It is an INTENT(OUT) argument. It is assigned the value of the   |
  !|  command argument specified by NUMBER. If the command argument    |
  !|  value cannot be determined, VALUE is assigned all blanks.        |
  !|LENGTH (optional) shall be scalar and of type default integer.     |
  !|  It is an INTENT(OUT) argument. It is assigned the significant    |
  !|  length of the command argument specified by NUMBER. The          |
  !|  significant length may include trailing blanks if the processor  |
  !|   allows command arguments with significant trailing blanks. This |
  !|  length does not consider any possible truncation or padding in   |
  !|   assigning the command argument value to the VALUE argument; in  |
  !|  fact the VALUE argument need not even be present. If the command |
  !|  argument length cannot be determined, a length of 0 is assigned. |
  !|STATUS (optional) shall be scalar and of type default integer.     |
  !|  It is an INTENT(OUT) argument. It is assigned the value 0 if     |
  !|  the argument retrieval is sucessful. It is assigned a            |
  !|  processor-dependent non-zero value if the argument retrieval.    |
  !|  fails                                                            |
  !|NOTE                                                               |
  !|  One possible reason for failure is that NUMBER is negative or    |
  !|  greater than COMMAND_ARGUMENT_COUNT().                           |
  !+-------------------------------------------------------------------+
  subroutine get_command_argument(number,value,length,status)
    integer         , intent(in)            :: number
    character(len=*), intent(out), optional :: value
    integer         , intent(out), optional :: length
    integer         , intent(out), optional :: status
    !
    !  A temporary variable for the rare case case where LENGTH is
    !  specified but VALUE is not. An arbitrary maximum argument length
    !  of 1000 characters should cover virtually all situations.
    character(len=1000) :: tmpval
    !
    integer :: iargc
    ! Possible error codes:
    ! 1 = Argument number is less than minimum
    ! 2 = Argument number exceeds maximum
    if (number < 0) then
      if (present(value )) value  = ' '
      if (present(length)) length = 0
      if (present(status)) status = 1
      return
    else if (number > iargc()) then
      if (present(value )) value  = ' '
      if (present(length)) length = 0
      if (present(status)) status = 2
      return
    end if
    !
    ! Get the argument if VALUE is present
    !
    if (present(value)) call getarg(number,value)
    !
    ! The LENGTH option is fairly pointless under Unix.
    ! Trailing spaces can only be specified using quotes.
    ! Since the command line has already been processed by the
    ! shell before the application sees it, we have no way of
    ! knowing the true length of any quoted arguments. LEN_TRIM
    ! is used to ensure at least some sort of meaningful result.
    !
    if (present(length)) then
      if (present(value)) then
        length = len_trim(value)
      else
        call getarg(number,tmpval)
        length = len_trim(tmpval)
      end if
    end if
    !
    ! Since GETARG does not return a result code, assume success
    !
    if (present(status)) status = 0
    return
    !
  end subroutine get_command_argument
  !+-------------------------------------------------------------------+
  !| The end of the subroutine get_command_argument                    |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This function is to to determine if a string is numeric.          |
  !+-------------------------------------------------------------------+
  !| ref: http://fcode.cn/code_gen-115-1.html                          |
  !+-------------------------------------------------------------------+
  integer function isnum(zval)
    !
    ! Verify that a character string represents a numerical value
    ! 确定字符是否是数值类型：
    !   0-非数值的字符串
    !   1-整数(integer)
    !   2-小数(fixed point real)
    !   3-指数类型实数(exponent type real)
    !   4-双精度实数指数形式(exponent type double)
    !
    character(len=*),intent (in) :: zval
    integer :: num, nmts, nexp, kmts, ifexp, ichr
    integer, parameter :: kint = 1 ! integer
    integer, parameter :: kfix = 2 ! fixed point real
    integer, parameter :: kexp = 3 ! exponent type real
    integer, parameter :: kdbl = 4 ! exponent type double
    ! 
    ! initialise
    !
    num = 0  ! 数字的格式，最后传递给isnum返回
    nmts = 0 ! 整数或浮点数的数字个数
    nexp = 0 ! 指数形式的数字个数
    kmts = 0 ! 有+-号为1，否则为0
    ifexp = 0! 似乎没用
    !
    ! loop over characters
    !
    ichr = 0
    !
    do
      !
      if(ichr>=len(zval)) then
        !
        ! last check
        !
        if (nmts==0) exit
        !
        if (num>=kexp .and. nexp==0) exit
        !
        isnum = num
        !
        return
        !
      end if
      !
      ichr = ichr + 1
      !
      select case (zval(ichr:ichr))
        ! process blanks
      case (' ')
        continue
      ! process digits
      case ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
        if (num==0) num = kint
        if (num<kexp) then
          nmts = nmts + 1
          ! 整数或浮点数+1
        else
          nexp = nexp + 1
          ! 指数形式+1
        end if
      ! process signs
      case ('+', '-')
        !
        if (num==0) then
          if (kmts>0) exit
          ! 出现2个符号，非数字
          kmts = 1
          num = kint
        else
          if (num<kexp) exit
          if (ifexp>0) exit
          ifexp = 1
        end if
        ! process decimal point
      case ('.')
        if (num/=kint .and. ichr/=1) exit
        ! 前面不是整数，小数点也不是第一个字符，则非数字
        num = kfix
        ! process exponent
      case ('e', 'E')
        if (num>=kexp) exit
        if (nmts==0) exit
        num = kexp
      case ('d', 'D')
        if (num>=kexp) exit
        if (nmts==0) exit
        num = kdbl
        ! any other character means the string is non-numeric
      case default
        exit
      end select
      !
    end do
    !
    ! if this point is reached, the string is non-numeric
    !
    isnum = 0
    !  
    return
    !
  end function isnum
  !+-------------------------------------------------------------------+
  !| The end of the function isnum                                     |
  !+-------------------------------------------------------------------+
  !
  subroutine mongen
    !
    integer :: imon,jmon,kmon,i
    !
    open(18,file='monitor.dat')
    imon=480
    jmon=0
    kmon=256
    !
    write(18,'(A)')'# monitor dat'
    do while(imon<1700)
    write(18,'(4(I4,1x))')imon,jmon,kmon
    imon=imon+10
    enddo
    close(18)
    print*,' >> monitor.dat'
    !
  end subroutine mongen
  !
end module controlfileread
!+---------------------------------------------------------------------+
!| The end of the module controlfileread                               |
!+---------------------------------------------------------------------+
