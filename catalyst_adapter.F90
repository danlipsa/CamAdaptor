module catalyst

  use cam_logfile,     only: iulog
!
! !PUBLIC TYPES:
  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: catalyst_init
  public :: catalyst_coprocess
  public :: catalyst_finalize

!
!================================================================================
CONTAINS
!================================================================================

  subroutine catalyst_init()

    !-----------------------------------------------------------------------
    !
    ! Arguments
    !

    !
    ! Locals
    !
    write(iulog,'(a13)') "catalyst_init"
    call coprocessorinitializewithpython("catalyst_coprocess.py",21)
 end subroutine catalyst_init

!==============================================================================

  subroutine catalyst_coprocess()

    !-----------------------------------------------------------------------
    !
    ! Arguments
    !

    !
    ! Locals
    !
    write(iulog,'(a18)') "catalyst_coprocess"

 end subroutine catalyst_coprocess

!==============================================================================
  subroutine catalyst_finalize()

    !-----------------------------------------------------------------------
    !
    ! Arguments
    !

    !
    ! Locals
    !
    write(iulog,'(a17)') "catalyst_finalize"
    call coprocessorfinalize()
 end subroutine catalyst_finalize

!==============================================================================

end module catalyst
