program test_loop_reproducer
   implicit none

   integer, parameter :: r8 = selected_real_kind(12, 300)
   integer, parameter :: above_canopy = 1
   integer, parameter :: below_canopy = 2

   type patch_type
      integer, allocatable :: column(:)
      integer, allocatable :: gridcell(:)
      integer, allocatable :: itype(:)
      logical, allocatable :: is_fates(:)
   end type patch_type

   type params_type
      real(r8) :: cv
      real(r8) :: a_coef
      real(r8) :: a_exp
      real(r8) :: csoilc
   end type params_type

   integer :: fn
   integer :: num_patches
   integer :: num_columns
   integer :: num_types
   integer :: f, p, c
   integer :: ierr
   character(len=64) :: arg

   type(patch_type) :: patch
   type(params_type) :: params_inst

   integer, allocatable :: filterp(:)

   real(r8), allocatable :: tlbef(:), t_veg(:), del2(:), del(:)
   real(r8), allocatable :: ram1(:), rah(:,:), raw(:,:)
   real(r8), allocatable :: ustar(:), um(:), temp1(:), temp2(:)
   real(r8), allocatable :: uaf(:), uuc(:)
   real(r8), allocatable :: dleaf_patch(:), dleaf(:)
   real(r8), allocatable :: rb(:), rb1(:)
   real(r8), allocatable :: elai(:), esai(:)
   real(r8), allocatable :: z0mg(:), htop(:), taf(:), t_grnd(:)
   real(r8), allocatable :: grnd_ch4_cond(:)
   real(r8), allocatable :: svpts(:), el(:), eah(:), forc_pbot(:), qaf(:)
   real(r8), allocatable :: rhaf(:)
   real(r8), allocatable :: rah1(:), raw1(:), rah2(:), raw2(:)
   real(r8), allocatable :: vpd(:)

   real(r8) :: vkc
   real(r8) :: nu_param
   real(r8) :: grav
   real(r8) :: ria
   real(r8) :: t0, t1
   real(r8) :: checksum

   logical :: use_undercanopy_stability
   logical :: use_biomass_heat_storage
   logical :: use_lch4

   ! ------------------------------------------------------------------
   ! Default scalable problem size.
   !
   ! Optional command line:
   !
   !   ./test
   !   ./test 1000000
   !   ./test 1000000 1000
   !   ./test 1000000 1000 16
   !
   ! Arguments are:
   !
   !   fn num_columns num_types
   !
   ! num_patches is set equal to fn.
   ! ------------------------------------------------------------------

   fn = 10000000
   num_columns = 10000
   num_types = 160

   call get_command_argument(1, arg)
   if (len_trim(arg) > 0) then
      read(arg, *, iostat=ierr) fn
      if (ierr /= 0) stop "Could not parse argument 1: fn"
   end if

   call get_command_argument(2, arg)
   if (len_trim(arg) > 0) then
      read(arg, *, iostat=ierr) num_columns
      if (ierr /= 0) stop "Could not parse argument 2: num_columns"
   end if

   call get_command_argument(3, arg)
   if (len_trim(arg) > 0) then
      read(arg, *, iostat=ierr) num_types
      if (ierr /= 0) stop "Could not parse argument 3: num_types"
   end if

   if (fn < 1) stop "fn must be positive"
   if (num_columns < 1) stop "num_columns must be positive"
   if (num_types < 1) stop "num_types must be positive"

   num_patches = fn

   print *, "Problem size:"
   print *, "  fn          = ", fn
   print *, "  num_patches = ", num_patches
   print *, "  num_columns = ", num_columns
   print *, "  num_types   = ", num_types

   ! ------------------------------------------------------------------
   ! Allocate.
   ! ------------------------------------------------------------------

   allocate(filterp(fn))

   allocate(patch%column(num_patches))
   allocate(patch%gridcell(num_patches))
   allocate(patch%itype(num_patches))
   allocate(patch%is_fates(num_patches))

   allocate(tlbef(num_patches), t_veg(num_patches))
   allocate(del2(num_patches), del(num_patches))
   allocate(ram1(num_patches))
   allocate(rah(num_patches, 2))
   allocate(raw(num_patches, 2))
   allocate(ustar(num_patches), um(num_patches))
   allocate(temp1(num_patches), temp2(num_patches))
   allocate(uaf(num_patches), uuc(num_patches))
   allocate(dleaf_patch(num_patches), dleaf(num_types))
   allocate(rb(num_patches), rb1(num_patches))
   allocate(elai(num_patches), esai(num_patches))
   allocate(z0mg(num_columns))
   allocate(htop(num_patches), taf(num_patches), t_grnd(num_columns))
   allocate(grnd_ch4_cond(num_patches))
   allocate(svpts(num_patches), el(num_patches))
   allocate(eah(num_patches), forc_pbot(num_columns), qaf(num_patches))
   allocate(rhaf(num_patches))
   allocate(rah1(num_patches), raw1(num_patches))
   allocate(rah2(num_patches), raw2(num_patches))
   allocate(vpd(num_patches))

   ! ------------------------------------------------------------------
   ! Constants and switches.
   ! ------------------------------------------------------------------

   params_inst%cv     = 0.01_r8
   params_inst%a_coef = 13.0_r8
   params_inst%a_exp  = 0.45_r8
   params_inst%csoilc = 0.004_r8

   vkc = 0.40_r8
   nu_param = 1.5e-5_r8
   grav = 9.80616_r8
   ria = 5.0_r8

   use_undercanopy_stability = .true.
   use_biomass_heat_storage = .false.
   use_lch4 = .true.

   ! ------------------------------------------------------------------
   ! Scalable initialization.
   !
   ! These loops preserve valid index ranges:
   !
   !   1 <= filterp(f)       <= num_patches
   !   1 <= patch%column(p)  <= num_columns
   !   1 <= patch%itype(p)   <= num_types
   !
   ! No hard-coded size-5 array constructors.
   ! ------------------------------------------------------------------

   do f = 1, fn
      filterp(f) = f
   end do

   do p = 1, num_patches
      patch%column(p)   = mod(p - 1, num_columns) + 1
      patch%gridcell(p) = patch%column(p)
      patch%itype(p)    = mod(p - 1, num_types) + 1
      patch%is_fates(p) = mod(p, 7) == 0
   end do

   do p = 1, num_types
      dleaf(p) = 0.03_r8 + 0.002_r8 * real(p, r8)
   end do

   do c = 1, num_columns
      z0mg(c)     = 0.05_r8 + 0.001_r8 * real(mod(c, 100), r8)
      t_grnd(c)   = 283.0_r8 + 0.01_r8 * real(mod(c, 300), r8)
      forc_pbot(c) = 95000.0_r8 + 10.0_r8 * real(mod(c, 1000), r8)
   end do

   do p = 1, num_patches
      t_veg(p) = 285.0_r8 + 0.01_r8 * real(mod(p, 100), r8)
      del(p)   = 0.10_r8 + 0.001_r8 * real(mod(p, 50), r8)

      ustar(p) = 0.20_r8 + 0.001_r8 * real(mod(p, 300), r8)
      um(p)    = 1.00_r8 + 0.010_r8 * real(mod(p, 500), r8)

      temp1(p) = 0.10_r8 + 0.001_r8 * real(mod(p, 50), r8)
      temp2(p) = 0.20_r8 + 0.001_r8 * real(mod(p, 50), r8)

      dleaf_patch(p) = 0.04_r8 + 0.001_r8 * real(mod(p, 20), r8)

      elai(p) = 1.0_r8 + 0.01_r8 * real(mod(p, 200), r8)
      esai(p) = 0.1_r8 + 0.01_r8 * real(mod(p, 20), r8)

      htop(p) = 5.0_r8 + 0.1_r8 * real(mod(p, 200), r8)
      taf(p)  = 285.0_r8 + 0.01_r8 * real(mod(p, 300), r8)

      el(p)  = 1500.0_r8 + 2.0_r8 * real(mod(p, 500), r8)
      qaf(p) = 0.005_r8 + 0.00001_r8 * real(mod(p, 500), r8)
   end do

   tlbef = 0.0_r8
   del2 = 0.0_r8
   ram1 = 0.0_r8
   rah = 0.0_r8
   raw = 0.0_r8
   uaf = 0.0_r8
   uuc = 0.0_r8
   rb = 0.0_r8
   rb1 = 0.0_r8
   grnd_ch4_cond = 0.0_r8
   svpts = 0.0_r8
   eah = 0.0_r8
   rhaf = 0.0_r8
   rah1 = 0.0_r8
   raw1 = 0.0_r8
   rah2 = 0.0_r8
   raw2 = 0.0_r8
   vpd = 0.0_r8

   call cpu_time(t0)

   call reproducer_loop( &
      fn, filterp, patch, params_inst, &
      tlbef, t_veg, del2, del, ram1, rah, raw, &
      ustar, um, temp1, temp2, uaf, uuc, &
      dleaf_patch, dleaf, rb, rb1, elai, esai, &
      z0mg, htop, taf, t_grnd, grnd_ch4_cond, &
      svpts, el, eah, forc_pbot, qaf, rhaf, &
      rah1, raw1, rah2, raw2, vpd, &
      vkc, nu_param, grav, ria, &
      use_undercanopy_stability, use_biomass_heat_storage, use_lch4)

   call cpu_time(t1)

   checksum = sum(tlbef) + sum(del2) + sum(ram1) + sum(rah) + sum(raw) + &
              sum(uaf) + sum(uuc) + sum(rb) + sum(rb1) + sum(grnd_ch4_cond) + &
              sum(svpts) + sum(eah) + sum(rhaf) + sum(rah1) + sum(raw1) + &
              sum(rah2) + sum(raw2) + sum(vpd)

   print *, "Loop completed."
   print *, "CPU time in loop section, seconds = ", t1 - t0
   print *, "Checksum = ", checksum

contains

   subroutine reproducer_loop( &
      fn, filterp, patch, params_inst, &
      tlbef, t_veg, del2, del, ram1, rah, raw, &
      ustar, um, temp1, temp2, uaf, uuc, &
      dleaf_patch, dleaf, rb, rb1, elai, esai, &
      z0mg, htop, taf, t_grnd, grnd_ch4_cond, &
      svpts, el, eah, forc_pbot, qaf, rhaf, &
      rah1, raw1, rah2, raw2, vpd, &
      vkc, nu_param, grav, ria, &
      use_undercanopy_stability, use_biomass_heat_storage, use_lch4)

      implicit none

      integer, intent(in) :: fn
      integer, intent(in) :: filterp(:)

      type(patch_type), intent(in) :: patch
      type(params_type), intent(in) :: params_inst

      real(r8), intent(inout) :: tlbef(:)
      real(r8), intent(in)    :: t_veg(:)

      real(r8), intent(inout) :: del2(:)
      real(r8), intent(in)    :: del(:)

      real(r8), intent(inout) :: ram1(:)
      real(r8), intent(inout) :: rah(:,:)
      real(r8), intent(inout) :: raw(:,:)

      real(r8), intent(in)    :: ustar(:)
      real(r8), intent(in)    :: um(:)
      real(r8), intent(in)    :: temp1(:)
      real(r8), intent(in)    :: temp2(:)

      real(r8), intent(inout) :: uaf(:)
      real(r8), intent(inout) :: uuc(:)

      real(r8), intent(inout) :: dleaf_patch(:)
      real(r8), intent(in)    :: dleaf(:)

      real(r8), intent(inout) :: rb(:)
      real(r8), intent(inout) :: rb1(:)

      real(r8), intent(in)    :: elai(:)
      real(r8), intent(in)    :: esai(:)
      real(r8), intent(in)    :: z0mg(:)
      real(r8), intent(in)    :: htop(:)
      real(r8), intent(in)    :: taf(:)
      real(r8), intent(in)    :: t_grnd(:)

      real(r8), intent(inout) :: grnd_ch4_cond(:)

      real(r8), intent(inout) :: svpts(:)
      real(r8), intent(in)    :: el(:)

      real(r8), intent(inout) :: eah(:)
      real(r8), intent(in)    :: forc_pbot(:)
      real(r8), intent(in)    :: qaf(:)

      real(r8), intent(inout) :: rhaf(:)
      real(r8), intent(inout) :: rah1(:)
      real(r8), intent(inout) :: raw1(:)
      real(r8), intent(inout) :: rah2(:)
      real(r8), intent(inout) :: raw2(:)
      real(r8), intent(inout) :: vpd(:)
      real(r8), allocatable :: hold(:)

      real(r8), intent(in) :: vkc
      real(r8), intent(in) :: nu_param
      real(r8), intent(in) :: grav
      real(r8), intent(in) :: ria

      logical, intent(in) :: use_undercanopy_stability
      logical, intent(in) :: use_biomass_heat_storage
      logical, intent(in) :: use_lch4

      integer :: f
      integer :: p
      integer :: c
      integer :: g

      real(r8) :: cf
      real(r8) :: w
      real(r8) :: csoilb
      real(r8) :: ri
      real(r8) :: ricsoilc
      real(r8) :: csoilcn
      real(r8) :: t0, t1

      allocate(hold(fn))
  
call cpu_time(t0)
      do concurrent (p = 1:fn)
      if (.not. patch%is_fates(p)) then
           hold(p) = dleaf(patch%itype(p))
      else
           hold(p) = dleaf_patch(p)
      end if
      enddo

      !PGI$ VECTOR ALWAYS
      !DIR$ IVDEP
      do f = 1,fn
         p = filterp(f)
         c = patch%column(p)
         g = patch%gridcell(p)
         tlbef(p) = t_veg(p)
         del2(p) = del(p)
         ram1(p) = 1._r8 / (ustar(p) * ustar(p) / um(p))
         rah(p,above_canopy) = 1._r8 / (temp1(p) * ustar(p))
         raw(p,above_canopy) = 1._r8 / (temp2(p) * ustar(p))
         uaf(p) = um(p) * sqrt(1._r8 / (ram1(p) * um(p)))
         uuc(p) = min(0.4_r8, 0.03_r8 * um(p) / ustar(p))
!         if (.not. patch%is_fates(p)) then
           dleaf_patch(p) = hold(p) 
!         end if
         cf = params_inst%cv / (sqrt(uaf(p)) * sqrt(dleaf_patch(p)))
         rb(p) = 1._r8 / (cf * uaf(p))
         rb1(p) = rb(p)
         w = exp(-(elai(p) + esai(p)))
         csoilb = vkc / &
            (params_inst%a_coef * (z0mg(c) * uaf(p) / nu_param)**params_inst%a_exp)
         ri = (grav * htop(p) * (taf(p) - t_grnd(c))) / &
            (taf(p) * uaf(p)**2.00_r8)
         if (use_undercanopy_stability .and. (taf(p) - t_grnd(c)) > 0._r8) then
            ricsoilc = params_inst%csoilc / (1.00_r8 + ria * min(ri, 10.0_r8))
            csoilcn = csoilb * w + ricsoilc * (1._r8 - w)
         else
            csoilcn = csoilb * w + params_inst%csoilc * (1._r8 - w)
         end if
         if (use_biomass_heat_storage) then
            rah(p,below_canopy) = 1._r8 / (csoilcn * uuc(p))
         else
            rah(p,below_canopy) = 1._r8 / (csoilcn * uaf(p))
         end if
         raw(p,below_canopy) = rah(p,below_canopy)
         if (use_lch4) then
            grnd_ch4_cond(p) = 1._r8 / &
               (raw(p,above_canopy) + raw(p,below_canopy))
         end if
         svpts(p) = el(p)
         eah(p) = forc_pbot(c) * qaf(p) / 0.622_r8
         rhaf(p) = eah(p) / svpts(p)
         rah1(p) = rah(p,above_canopy)
         raw1(p) = raw(p,above_canopy)
         rah2(p) = rah(p,below_canopy)
         raw2(p) = raw(p,below_canopy)
         vpd(p) = max((svpts(p) - eah(p)), 50._r8) * 0.001_r8
      end do
call cpu_time(t1)

   end subroutine reproducer_loop

end program test_loop_reproducer
