! This routine uses modules and variables from cam5 which is an executable,
! so this routine must be linked in.

subroutine build_params_and_coprocess(phys_state)
  use time_manager, only: get_nstep, get_curr_time
  use dyn_grid,     only: get_horiz_grid_dim_d, get_dyn_grid_parm_real1d, &
       get_dyn_grid_parm
  use hycoef,       only: hyam, hybm, ps0
  use cam_history_support, only : registeredmdims, hist_mdims
  use ppgrid,       only: pver
  use cam_catalyst_adapter  , only: catalyst_coprocess, catalyst_create_grid, &
       catalyst_add_chunk
  ! for getting the MPI rank
  use cam_pio_utils, only: pio_subsystem

  ! Input/Output arguments
  !
  type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
  !-----------------------------------------------------------------------
  !
  ! Locals
  !
  integer :: nstep          ! current timestep number
  real(kind=8) :: time      ! current time
  integer :: ndcur          ! day component of current time
  integer :: nscur          ! seconds component of current time
  integer, dimension(1:3) :: dim    ! lon, lat and lev
  real(r8), pointer :: latdeg(:)    ! degrees gaussian latitudes 
  integer :: plon
  real(r8), allocatable :: alon(:)  ! longitude values (degrees)
  real(r8) :: alev(pver)    ! level values (pascals)
  integer :: i,f,c          ! indexes
  integer :: nPoints2D

  call t_startf ('catalyst_coprocess')

  ! current time step and time
  nstep = get_nstep()
  call get_curr_time(ndcur, nscur)
  time = ndcur + nscur/86400._r8

  ! lon, lat and lev
  call get_horiz_grid_dim_d(dim(1),dim(2))
  dim(3) = 1
  do i=1,registeredmdims
     if (trim(hist_mdims(i)%name) == 'lev') then
        dim(3) = hist_mdims(i)%value
        exit
     endif
  end do

  ! longitude
  plon = get_dyn_grid_parm('plon')
  allocate(alon(plon))
  do i=1,plon
     alon(i) = (i-1) * 360.0_r8 / plon
  end do

  ! latitude
  latdeg => get_dyn_grid_parm_real1d('latdeg')

  ! levels
  ! converts Pascals to millibars
  alev(:pver) = 0.01_r8*ps0*(hyam(:pver) + hybm(:pver))

  ! total number of points on a MPI node
  nPoints2D = 0
  do c=begchunk, endchunk
     nPoints2D = nPoints2D + get_ncols_p(c)
  end do
  if (catalyst_create_grid(nstep, time, dim, alon, latdeg, alev, &
       nPoints2D, pio_subsystem%comp_rank)) then
     do c=begchunk, endchunk
        call catalyst_add_chunk(nstep, time, phys_state(c), get_ncols_p(c))
     end do
     call catalyst_coprocess(nstep, time)
  end if
  deallocate(alon)
  call t_stopf ('catalyst_coprocess')
end subroutine build_params_and_coprocess
