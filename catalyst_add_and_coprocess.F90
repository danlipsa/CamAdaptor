! Creates the grid (if necessary), copies the data into Catalyst and
! calls coprocess
! This routine uses functions and variables from cam5 which is an executable,
! so this routine is included into a physpkg.F90, a file in cam5
subroutine catalyst_add_and_coprocess(phys_state)
  use time_manager, only: get_nstep, get_curr_time
  use physconst,     only: pi
  use dyn_grid,     only: get_horiz_grid_dim_d, get_horiz_grid_d, get_dyn_grid_parm_real1d, &
       get_dyn_grid_parm
  use dycore,       only: dycore_is
  use hycoef,       only: hyam, hybm, ps0
  use ppgrid,       only: pver
  use fv_catalyst_adapter  , only: fv_catalyst_coprocess, fv_catalyst_create_grid, &
       fv_catalyst_add_chunk
  use se_catalyst_adapter  , only: se_catalyst_coprocess, se_catalyst_create_grid, &
       se_catalyst_add_chunk
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

  integer             :: dim1s,dim2s             ! global size of the first and second horizontal dim.
  integer             :: ncol
  real(r8), pointer   :: alat(:)                 ! latitude values (degrees)
  integer             :: ne, np                  ! SE grid parameters

  call t_startf ('catalyst_add_and_coprocess')

  ! current time step and time
  nstep = get_nstep()
  call get_curr_time(ndcur, nscur)
  time = ndcur + nscur/86400._r8

  if (dycore_is('LR') ) then
     ! FV dynamic core

     ! lon, lat
     call get_horiz_grid_dim_d(dim(1),dim(2))
     ! lev
     dim(3) = pver

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
     if (fv_catalyst_create_grid(nstep, time, dim, alon, latdeg, alev, &
          nPoints2D, pio_subsystem%comp_rank)) then
        do c=begchunk, endchunk
           call fv_catalyst_add_chunk(nstep, time, phys_state(c), get_ncols_p(c))
        end do
        call fv_catalyst_coprocess(nstep, time)
     end if
     deallocate(alon)

  else if (dycore_is('UNSTRUCTURED')) then
     ! SE dynamic core
     call get_horiz_grid_dim_d(dim1s, dim2s)
     ncol = dim1s*dim2s

     ! read longitude and latitude
     allocate(alon(ncol))
     call get_horiz_grid_d(ncol, clon_d_out=alon)

     allocate(alat(ncol))
     call get_horiz_grid_d(ncol, clat_d_out=alat)

     ! levels
     ! converts Pascals to millibars
     alev(:pver) = 0.01_r8*ps0*(hyam(:pver) + hybm(:pver))
     
     ! total number of points on a MPI node
     nPoints2D = 0
     do c=begchunk, endchunk
        nPoints2D = nPoints2D + get_ncols_p(c)
     end do

     ! ne, np
     ne = get_dyn_grid_parm('ne')
     np = get_dyn_grid_parm('np')
     
     if (se_catalyst_create_grid(nstep, time, ne, np, ncol, alon, ncol, alat, pver, alev, &
          nPoints2D, pio_subsystem%comp_rank)) then
        do c=begchunk, endchunk
           call se_catalyst_add_chunk(nstep, time, phys_state(c), get_ncols_p(c))
        end do
        call se_catalyst_coprocess(nstep, time)
     end if
     deallocate(alon)
     deallocate(alat)
  endif

  call t_stopf ('catalyst_add_and_coprocess')
end subroutine catalyst_add_and_coprocess
