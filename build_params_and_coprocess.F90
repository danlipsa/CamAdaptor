! This routine uses modules and variables from cam5 which is an executable,
! so this routine must be linked in.
subroutine build_params_and_coprocess()
  use time_manager, only: get_nstep, get_curr_time
  use dyn_grid,     only: get_horiz_grid_dim_d, get_dyn_grid_parm_real1d, get_dyn_grid_parm
  use hycoef,       only: hyam, hybm, ps0
  use cam_history_support, only : registeredmdims, hist_mdims
  use ppgrid,       only: pcols, pver
  use cam_catalyst_adapter  , only: catalyst_coprocess

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
  integer :: i,f            ! indexes
  integer, allocatable :: fieldIndex(:) ! indexes into array of filds for tape 0
                                        ! that have 4 dimensions
  integer :: t              ! tape number
  character(len=max_chars) :: fname_tmp ! local copy of field name

  call t_startf ('catalyst_coprocess')
  ! we look at fileds defined for tape 0
  t = 1

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

  ! londeg
  plon = get_dyn_grid_parm('plon')
  allocate(alon(plon))
  do i=1,plon
     alon(i) = (i-1) * 360.0_r8 / plon
  end do
  ! latdeg
  latdeg => get_dyn_grid_parm_real1d('latdeg')
  ! levdeg
  ! converts Pascals to millibars
  alev(:pver) = 0.01_r8*ps0*(hyam(:pver) + hybm(:pver))

  ! fields
  allocate(fieldIndex(nflds(t) + 1))
  i = 0
  do f=1,nflds(t)
     if (associated(tape(t)%hlist(f)%field%mdims) .and. size(tape(t)%hlist(f)%field%mdims) == 1) then
        i = i + 1
        fieldIndex(i) = f
        !fname_tmp = strip_suffix(tape(t)%hlist(f)%field%name)        
        !write(iulog, '(a20)') trim(fname_tmp)
     endif
  enddo
  call catalyst_coprocess(nstep, time, dim, alon, latdeg, alev, fieldIndex, i, tape(t)%hlist)
  deallocate(fieldIndex)
  deallocate(alon)
  call t_stopf ('catalyst_coprocess')
end subroutine build_params_and_coprocess
