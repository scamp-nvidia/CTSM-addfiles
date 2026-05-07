! Standalone bucketed CPU reproducer for the CanopyFluxesMod.F90 can_iter region.
!
! The timed region mirrors the code bracketed by:
!   call t_startf('can_iter')
!   call t_stopf('can_iter')
!
! Large CTSM services that are not practical in a single-file reproducer
! (photosynthesis, FATES, CN fractionation, type infrastructure) are replaced
! with deterministic local kernels. FrictionVelocity and QSat are kept close to
! the original formulas because they directly feed the can_iter loop.
!
! CPU optimization strategy:
!   - Keep the shrinking filterp work list, which is faster on CPUs than a
!     full-length active mask for this workload.
!   - Dispatch once per can_iter call to a variant that matches the current
!     options/layout, falling back to the full general implementation.
!   - Inside the fast default variant, classify active patches into small
!     dynamic buckets and run branch-specialized loops over each bucket.
!   - Keep all state arrays and diagnostic calculations alive.
!
! Build examples:
!   gfortran  -O3 -march=native -ffree-line-length-none \
!             can_iter_cpu_bucketed_reproducer.f90 -o can_iter_cpu_bucketed.exe
!   ifx       -O3 /QxHost can_iter_cpu_bucketed_reproducer.f90 -o can_iter_cpu_bucketed.exe
!   nvfortran -O3 -Minfo=all can_iter_cpu_bucketed_reproducer.f90 -o can_iter_cpu_bucketed.exe
!
! Run:
!   ./can_iter_cpu_bucketed.exe [npts] [nrepeat] [nwarmup] [itmax]
!
! Example:
!   ./can_iter_cpu_bucketed.exe 65536 20 3 40

module can_iter_cpu_bucketed_reproducer_mod
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  private

  integer, parameter, public :: r8 = selected_real_kind(15, 307)

  integer, parameter :: above_canopy = 1
  integer, parameter :: below_canopy = 2

  real(r8), parameter :: pi = 3.141592653589793238462643383279502884_r8
  real(r8), parameter :: tkfrz = 273.15_r8
  real(r8), parameter :: sb = 5.670374419e-8_r8
  real(r8), parameter :: cpair = 1004.64_r8
  real(r8), parameter :: hvap = 2.501e6_r8
  real(r8), parameter :: vkc = 0.40_r8
  real(r8), parameter :: grav = 9.80616_r8
  real(r8), parameter :: nu_param = 1.5e-5_r8
  real(r8), parameter :: btran0 = 0.0_r8
  real(r8), parameter :: zii = 1000.0_r8
  real(r8), parameter :: beta = 1.0_r8
  real(r8), parameter :: delmax = 1.0_r8
  real(r8), parameter :: dlemin = 0.1_r8
  real(r8), parameter :: dtmin = 0.01_r8
  real(r8), parameter :: ria = 0.5_r8
  integer, parameter :: itmin = 2
  integer, parameter :: variant_general = 0
  integer, parameter :: variant_default_bucketed = 1

  type, public :: can_iter_params
    real(r8) :: lai_dl = 0.20_r8
    real(r8) :: z_dl = 0.05_r8
    real(r8) :: a_coef = 0.30_r8
    real(r8) :: a_exp = 0.50_r8
    real(r8) :: csoilc = 0.004_r8
    real(r8) :: cv = 0.010_r8
    real(r8) :: wind_min = 1.00_r8
    real(r8) :: zetamaxstable = 0.50_r8
  end type can_iter_params

  type, public :: can_iter_options
    logical :: use_lch4 = .true.
    logical :: use_hydrstress = .false.
    logical :: use_biomass_heat_storage = .true.
    logical :: use_undercanopy_stability = .true.
    logical :: use_fates = .false.
    logical :: use_cn = .false.
    logical :: use_c13 = .false.
    logical :: do_soilevap_beta = .false.
    logical :: do_soil_resistance_sl14 = .true.
    real(r8) :: dtime = 1800.0_r8
  end type can_iter_options

  type, public :: can_iter_state
    integer :: npts = 0
    integer :: fn0 = 0
    integer :: can_iter_variant = 0

    integer, allocatable :: filterp(:), fporig(:)
    integer, allocatable :: filter_under_stable(:), filter_under_neutral(:)
    integer, allocatable :: filter_delq_neg(:), filter_delq_nonneg(:)
    integer, allocatable :: column(:), gridcell(:), itype(:), snl(:), nmozsgn(:)
    logical, allocatable :: is_fates(:)

    real(r8), allocatable :: forc_u(:), forc_v(:), forc_q(:), forc_pbot(:)
    real(r8), allocatable :: forc_th(:), forc_rho(:), forc_lwrad(:)
    real(r8), allocatable :: forc_pco2(:), forc_po2(:)
    real(r8), allocatable :: forc_hgt_u_patch(:), forc_hgt_t_patch(:), forc_hgt_q_patch(:)

    real(r8), allocatable :: thv(:), thm(:), t_grnd(:), t_veg(:), t_stem(:)
    real(r8), allocatable :: t_soisno_snlp1(:), t_soisno_1(:), t_h2osfc(:)
    real(r8), allocatable :: emv(:), emg(:), htop(:), displa(:)
    real(r8), allocatable :: z0mv(:), z0hv(:), z0qv(:), z0mg(:)

    real(r8), allocatable :: elai(:), esai(:), laisun(:), laisha(:)
    real(r8), allocatable :: dleaf_patch(:), leafn_patch(:), dayl_factor(:)
    real(r8), allocatable :: fdry(:), fwet(:), liqcan(:), snocan(:), frac_veg_nosno(:)
    real(r8), allocatable :: qg(:), qg_snow(:), qg_soil(:), qg_h2osfc(:)
    real(r8), allocatable :: frac_sno(:), frac_h2osfc(:), soilbeta(:), soilresis(:)
    real(r8), allocatable :: snow_depth(:)

    real(r8), allocatable :: sabv(:), frac_rad_abs_by_stem(:)
    real(r8), allocatable :: sa_internal(:), sa_leaf(:), sa_stem(:)
    real(r8), allocatable :: cp_leaf(:), cp_stem(:), rstem(:), stem_biomass(:)
    real(r8), allocatable :: rssun(:), rssha(:), btran(:), bsun(:), bsha(:)

    real(r8), allocatable :: qflx_tran_veg(:), qflx_evap_veg(:)
    real(r8), allocatable :: eflx_sh_veg(:), eflx_sh_stem(:)
    real(r8), allocatable :: grnd_ch4_cond(:), canopy_cond(:), num_iter(:)

    real(r8), allocatable :: ram1(:), rb1(:), rah1(:), rah2(:), raw1(:), raw2(:)
    real(r8), allocatable :: ur(:), ustar(:), um(:), uaf(:), taf(:), qaf(:), obu(:), zeta(:), vpd(:)
    real(r8), allocatable :: rhaf(:)

    real(r8), allocatable :: air(:), bir(:), cir(:), co2(:), o2(:), svpts(:), eah(:)
    real(r8), allocatable :: el(:), qsatl(:), qsatldT(:)
    real(r8), allocatable :: del(:), del2(:), dele(:), det(:), efeb(:), efe(:), err(:)
    real(r8), allocatable :: dt_veg(:), tlbef(:), tl_ini(:), ts_ini(:), obuold(:)
    real(r8), allocatable :: dth(:), dthv(:), dqh(:), delq(:), zldis(:)
    real(r8), allocatable :: temp1(:), temp2(:), temp12m(:), temp22m(:), fm(:)

    real(r8), allocatable :: rah(:, :), raw(:, :), rb(:)
    real(r8), allocatable :: wtg(:), wta0(:), wtl0(:), wtstem0(:), wtal(:), wtga(:)
    real(r8), allocatable :: wtgq(:), wtaq0(:), wtlq0(:), wtalq(:)
    real(r8), allocatable :: lw_stem(:), lw_leaf(:), uuc(:)
  end type can_iter_state

  public :: allocate_state
  public :: reset_state
  public :: setup_can_iter
  public :: can_iter_cpu_bucketed
  public :: compute_checksum
  public :: print_summary

contains

  subroutine allocate_state(s, npts)
    type(can_iter_state), intent(inout) :: s
    integer, intent(in) :: npts

    s%npts = npts
    s%fn0 = npts
    s%can_iter_variant = variant_general

    allocate(s%filterp(npts), s%fporig(npts))
    allocate(s%filter_under_stable(npts), s%filter_under_neutral(npts))
    allocate(s%filter_delq_neg(npts), s%filter_delq_nonneg(npts))
    allocate(s%column(npts), s%gridcell(npts), s%itype(npts), s%snl(npts), s%nmozsgn(npts))
    allocate(s%is_fates(npts))

    allocate(s%forc_u(npts), s%forc_v(npts), s%forc_q(npts), s%forc_pbot(npts))
    allocate(s%forc_th(npts), s%forc_rho(npts), s%forc_lwrad(npts))
    allocate(s%forc_pco2(npts), s%forc_po2(npts))
    allocate(s%forc_hgt_u_patch(npts), s%forc_hgt_t_patch(npts), s%forc_hgt_q_patch(npts))

    allocate(s%thv(npts), s%thm(npts), s%t_grnd(npts), s%t_veg(npts), s%t_stem(npts))
    allocate(s%t_soisno_snlp1(npts), s%t_soisno_1(npts), s%t_h2osfc(npts))
    allocate(s%emv(npts), s%emg(npts), s%htop(npts), s%displa(npts))
    allocate(s%z0mv(npts), s%z0hv(npts), s%z0qv(npts), s%z0mg(npts))

    allocate(s%elai(npts), s%esai(npts), s%laisun(npts), s%laisha(npts))
    allocate(s%dleaf_patch(npts), s%leafn_patch(npts), s%dayl_factor(npts))
    allocate(s%fdry(npts), s%fwet(npts), s%liqcan(npts), s%snocan(npts))
    allocate(s%frac_veg_nosno(npts))
    allocate(s%qg(npts), s%qg_snow(npts), s%qg_soil(npts), s%qg_h2osfc(npts))
    allocate(s%frac_sno(npts), s%frac_h2osfc(npts), s%soilbeta(npts), s%soilresis(npts))
    allocate(s%snow_depth(npts))

    allocate(s%sabv(npts), s%frac_rad_abs_by_stem(npts))
    allocate(s%sa_internal(npts), s%sa_leaf(npts), s%sa_stem(npts))
    allocate(s%cp_leaf(npts), s%cp_stem(npts), s%rstem(npts), s%stem_biomass(npts))
    allocate(s%rssun(npts), s%rssha(npts), s%btran(npts), s%bsun(npts), s%bsha(npts))

    allocate(s%qflx_tran_veg(npts), s%qflx_evap_veg(npts))
    allocate(s%eflx_sh_veg(npts), s%eflx_sh_stem(npts))
    allocate(s%grnd_ch4_cond(npts), s%canopy_cond(npts), s%num_iter(npts))

    allocate(s%ram1(npts), s%rb1(npts), s%rah1(npts), s%rah2(npts), s%raw1(npts), s%raw2(npts))
    allocate(s%ur(npts), s%ustar(npts), s%um(npts), s%uaf(npts), s%taf(npts), s%qaf(npts))
    allocate(s%obu(npts), s%zeta(npts), s%vpd(npts), s%rhaf(npts))

    allocate(s%air(npts), s%bir(npts), s%cir(npts), s%co2(npts), s%o2(npts))
    allocate(s%svpts(npts), s%eah(npts), s%el(npts), s%qsatl(npts), s%qsatldT(npts))
    allocate(s%del(npts), s%del2(npts), s%dele(npts), s%det(npts), s%efeb(npts))
    allocate(s%efe(npts), s%err(npts), s%dt_veg(npts), s%tlbef(npts), s%tl_ini(npts))
    allocate(s%ts_ini(npts), s%obuold(npts))
    allocate(s%dth(npts), s%dthv(npts), s%dqh(npts), s%delq(npts), s%zldis(npts))
    allocate(s%temp1(npts), s%temp2(npts), s%temp12m(npts), s%temp22m(npts), s%fm(npts))

    allocate(s%rah(npts, 2), s%raw(npts, 2), s%rb(npts))
    allocate(s%wtg(npts), s%wta0(npts), s%wtl0(npts), s%wtstem0(npts))
    allocate(s%wtal(npts), s%wtga(npts))
    allocate(s%wtgq(npts), s%wtaq0(npts), s%wtlq0(npts), s%wtalq(npts))
    allocate(s%lw_stem(npts), s%lw_leaf(npts), s%uuc(npts))
  end subroutine allocate_state

  subroutine reset_state(s)
    type(can_iter_state), intent(inout) :: s
    integer :: p
    real(r8) :: x, wave1, wave2, lai_total, stem_area

    s%fn0 = s%npts
    s%can_iter_variant = variant_general

    do p = 1, s%npts
      x = real(mod(p - 1, 997), r8) / 997.0_r8
      wave1 = sin(0.017_r8 * real(p, r8))
      wave2 = cos(0.011_r8 * real(p, r8))

      s%filterp(p) = p
      s%fporig(p) = p
      s%column(p) = p
      s%gridcell(p) = p
      s%itype(p) = 1 + mod(p - 1, 4)
      s%snl(p) = 0
      s%is_fates(p) = .false.

      s%forc_u(p) = 1.4_r8 + 2.2_r8 * x + 0.30_r8 * wave1
      s%forc_v(p) = 0.4_r8 + 0.9_r8 * (1.0_r8 - x) + 0.20_r8 * wave2
      s%forc_q(p) = 0.0060_r8 + 0.0025_r8 * x
      s%forc_pbot(p) = 88000.0_r8 + 2500.0_r8 * x
      s%forc_th(p) = 292.0_r8 + 9.0_r8 * x + 0.8_r8 * wave2
      s%forc_rho(p) = 1.10_r8 + 0.08_r8 * (1.0_r8 - x)
      s%forc_lwrad(p) = 285.0_r8 + 45.0_r8 * x + 4.0_r8 * wave1
      s%forc_pco2(p) = 40.5_r8 + 0.2_r8 * wave1
      s%forc_po2(p) = 20900.0_r8

      s%htop(p) = 4.0_r8 + 18.0_r8 * x
      s%displa(p) = 0.65_r8 * s%htop(p)
      s%z0mv(p) = max(0.05_r8, 0.08_r8 * s%htop(p))
      s%z0hv(p) = s%z0mv(p)
      s%z0qv(p) = s%z0mv(p)
      s%z0mg(p) = 0.010_r8 + 0.015_r8 * x
      s%forc_hgt_u_patch(p) = 30.0_r8 + s%z0mv(p) + s%displa(p)
      s%forc_hgt_t_patch(p) = 30.0_r8 + s%z0hv(p) + s%displa(p)
      s%forc_hgt_q_patch(p) = 30.0_r8 + s%z0qv(p) + s%displa(p)

      s%t_grnd(p) = 291.0_r8 + 7.0_r8 * x + 0.5_r8 * wave1
      s%t_veg(p) = s%t_grnd(p) + 0.8_r8 - 0.4_r8 * wave2
      s%t_stem(p) = s%t_veg(p) - 0.3_r8
      s%t_soisno_snlp1(p) = s%t_grnd(p) - 1.0_r8
      s%t_soisno_1(p) = s%t_grnd(p) - 0.5_r8
      s%t_h2osfc(p) = s%t_grnd(p) + 0.2_r8
      s%thm(p) = s%forc_th(p) - 0.6_r8 + 0.3_r8 * wave1
      s%thv(p) = s%forc_th(p) * (1.0_r8 + 0.61_r8 * s%forc_q(p))
      s%emv(p) = 0.96_r8
      s%emg(p) = 0.94_r8

      s%elai(p) = 0.8_r8 + 4.2_r8 * x
      s%esai(p) = 0.15_r8 + 0.75_r8 * (1.0_r8 - x)
      lai_total = s%elai(p) + s%esai(p)
      s%laisun(p) = 0.45_r8 * s%elai(p) + 0.10_r8 * s%elai(p) * max(wave1, 0.0_r8)
      s%laisha(p) = max(0.05_r8, s%elai(p) - s%laisun(p))
      s%dleaf_patch(p) = 0.025_r8 + 0.010_r8 * real(s%itype(p), r8)
      s%leafn_patch(p) = 1.2_r8 + 1.8_r8 * x
      s%dayl_factor(p) = min(1.0_r8, max(0.01_r8, 0.25_r8 + 0.75_r8 * x))

      s%fdry(p) = 0.55_r8 + 0.25_r8 * x
      s%fwet(p) = 0.08_r8 + 0.15_r8 * (1.0_r8 - x)
      s%liqcan(p) = 0.020_r8 + 0.030_r8 * max(wave1, 0.0_r8)
      s%snocan(p) = 0.004_r8 * max(-wave2, 0.0_r8)
      s%frac_veg_nosno(p) = 1.0_r8

      s%qg(p) = s%forc_q(p) + 0.0015_r8 + 0.0004_r8 * wave2
      s%qg_snow(p) = s%qg(p) - 0.0008_r8
      s%qg_soil(p) = s%qg(p)
      s%qg_h2osfc(p) = s%qg(p) + 0.0005_r8
      s%frac_sno(p) = 0.04_r8 * max(-wave1, 0.0_r8)
      s%frac_h2osfc(p) = 0.02_r8 * max(wave2, 0.0_r8)
      s%soilbeta(p) = 0.65_r8 + 0.30_r8 * x
      s%soilresis(p) = 80.0_r8 + 160.0_r8 * (1.0_r8 - x)
      s%snow_depth(p) = 0.03_r8 * s%frac_sno(p)

      s%sabv(p) = 120.0_r8 + 260.0_r8 * x + 10.0_r8 * wave1
      s%frac_rad_abs_by_stem(p) = min(0.25_r8, 0.10_r8 * s%esai(p) / lai_total)
      s%sa_leaf(p) = 2.0_r8 * s%elai(p) + s%esai(p)
      stem_area = 0.25_r8 + 0.75_r8 * s%esai(p)
      s%sa_stem(p) = stem_area
      s%sa_internal(p) = 0.35_r8 * min(s%sa_leaf(p), s%sa_stem(p))
      s%stem_biomass(p) = 1.0_r8 + 4.0_r8 * x
      s%cp_leaf(p) = 12000.0_r8 + 7000.0_r8 * x
      s%cp_stem(p) = 25000.0_r8 + 25000.0_r8 * x
      s%rstem(p) = 20.0_r8 + 20.0_r8 * x

      s%rssun(p) = 120.0_r8
      s%rssha(p) = 180.0_r8
      s%btran(p) = 0.35_r8 + 0.60_r8 * x
      s%bsun(p) = s%btran(p)
      s%bsha(p) = s%btran(p)

      s%qflx_tran_veg(p) = 0.0_r8
      s%qflx_evap_veg(p) = 0.0_r8
      s%eflx_sh_veg(p) = 0.0_r8
      s%eflx_sh_stem(p) = 0.0_r8
      s%grnd_ch4_cond(p) = 0.0_r8
      s%canopy_cond(p) = 0.0_r8
      s%num_iter(p) = 0.0_r8

      s%ram1(p) = 0.0_r8
      s%rb1(p) = 0.0_r8
      s%rah1(p) = 0.0_r8
      s%rah2(p) = 0.0_r8
      s%raw1(p) = 0.0_r8
      s%raw2(p) = 0.0_r8
      s%ur(p) = 0.0_r8
      s%ustar(p) = 0.0_r8
      s%um(p) = 0.0_r8
      s%uaf(p) = 0.0_r8
      s%taf(p) = 0.0_r8
      s%qaf(p) = 0.0_r8
      s%obu(p) = 0.0_r8
      s%zeta(p) = 0.0_r8
      s%vpd(p) = 0.0_r8
      s%rhaf(p) = 0.0_r8

      s%air(p) = 0.0_r8
      s%bir(p) = 0.0_r8
      s%cir(p) = 0.0_r8
      s%co2(p) = 0.0_r8
      s%o2(p) = 0.0_r8
      s%svpts(p) = 0.0_r8
      s%eah(p) = 0.0_r8
      s%el(p) = 0.0_r8
      s%qsatl(p) = 0.0_r8
      s%qsatldT(p) = 0.0_r8

      s%del(p) = 0.0_r8
      s%del2(p) = 0.0_r8
      s%dele(p) = 0.0_r8
      s%det(p) = 0.0_r8
      s%efeb(p) = 0.0_r8
      s%efe(p) = 0.0_r8
      s%err(p) = 0.0_r8
      s%dt_veg(p) = 0.0_r8
      s%tlbef(p) = 0.0_r8
      s%tl_ini(p) = 0.0_r8
      s%ts_ini(p) = 0.0_r8
      s%obuold(p) = 0.0_r8

      s%dth(p) = 0.0_r8
      s%dthv(p) = 0.0_r8
      s%dqh(p) = 0.0_r8
      s%delq(p) = 0.0_r8
      s%zldis(p) = 0.0_r8
      s%temp1(p) = 0.0_r8
      s%temp2(p) = 0.0_r8
      s%temp12m(p) = 0.0_r8
      s%temp22m(p) = 0.0_r8
      s%fm(p) = 0.0_r8

      s%rah(p, above_canopy) = 0.0_r8
      s%rah(p, below_canopy) = 0.0_r8
      s%raw(p, above_canopy) = 0.0_r8
      s%raw(p, below_canopy) = 0.0_r8
      s%rb(p) = 0.0_r8
      s%wtg(p) = 0.0_r8
      s%wta0(p) = 0.0_r8
      s%wtl0(p) = 0.0_r8
      s%wtstem0(p) = 0.0_r8
      s%wtal(p) = 0.0_r8
      s%wtga(p) = 0.0_r8
      s%wtgq(p) = 0.0_r8
      s%wtaq0(p) = 0.0_r8
      s%wtlq0(p) = 0.0_r8
      s%wtalq(p) = 0.0_r8
      s%lw_stem(p) = 0.0_r8
      s%lw_leaf(p) = 0.0_r8
      s%uuc(p) = 0.0_r8
    end do
  end subroutine reset_state

  subroutine setup_can_iter(s, params, opts)
    type(can_iter_state), intent(inout) :: s
    type(can_iter_params), intent(in) :: params
    type(can_iter_options), intent(in) :: opts
    integer :: f, p, c, g

    do f = 1, s%fn0
      p = s%filterp(f)
      c = s%column(p)
      g = s%gridcell(p)

      s%del(p) = 0.0_r8
      s%efeb(p) = 0.0_r8
      s%wtlq0(p) = 0.0_r8
      s%wtalq(p) = 0.0_r8
      s%wtgq(p) = 0.0_r8
      s%wtaq0(p) = 0.0_r8
      s%obuold(p) = 0.0_r8
      s%eflx_sh_stem(p) = 0.0_r8

      if (.not. opts%use_biomass_heat_storage) then
        s%sa_leaf(p) = s%elai(p) + s%esai(p)
        s%frac_rad_abs_by_stem(p) = 0.0_r8
        s%sa_stem(p) = 0.0_r8
        s%sa_internal(p) = 0.0_r8
        s%cp_leaf(p) = 0.0_r8
        s%cp_stem(p) = 0.0_r8
        s%rstem(p) = 0.0_r8
      end if

      s%air(p) = s%emv(p) * (1.0_r8 + (1.0_r8 - s%emv(p)) * (1.0_r8 - s%emg(c))) &
           * s%forc_lwrad(c)
      s%bir(p) = -(2.0_r8 - s%emv(p) * (1.0_r8 - s%emg(c))) * s%emv(p) * sb
      s%cir(p) = s%emv(p) * s%emg(c) * sb

      call QSat(s%t_veg(p), s%forc_pbot(c), s%qsatl(p), es=s%el(p), qsdT=s%qsatldT(p))

      s%co2(p) = s%forc_pco2(g)
      s%o2(p) = s%forc_po2(g)
      s%nmozsgn(p) = 0

      s%taf(p) = (s%t_grnd(c) + s%thm(p)) / 2.0_r8
      s%qaf(p) = (s%forc_q(c) + s%qg(c)) / 2.0_r8

      s%ur(p) = max(params%wind_min, sqrt(s%forc_u(g) * s%forc_u(g) + s%forc_v(g) * s%forc_v(g)))
      s%dth(p) = s%thm(p) - s%taf(p)
      s%dqh(p) = s%forc_q(c) - s%qaf(p)
      s%delq(p) = s%qg(c) - s%qaf(p)
      s%dthv(p) = s%dth(p) * (1.0_r8 + 0.61_r8 * s%forc_q(c)) + 0.61_r8 * s%forc_th(c) * s%dqh(p)
      s%zldis(p) = s%forc_hgt_u_patch(p) - s%displa(p)

      call MoninObukIni(params, s%ur(p), s%thv(c), s%dthv(p), s%zldis(p), s%z0mv(p), &
           s%um(p), s%obu(p))

      s%num_iter(p) = 0.0_r8
      s%tl_ini(p) = s%t_veg(p)
      s%ts_ini(p) = s%t_stem(p)
    end do
    s%can_iter_variant = select_can_iter_variant(s, opts)
  end subroutine setup_can_iter

  subroutine can_iter_cpu_bucketed(s, params, opts, itmax_canopy_fluxes)
    type(can_iter_state), intent(inout) :: s
    type(can_iter_params), intent(in) :: params
    type(can_iter_options), intent(in) :: opts
    integer, intent(in) :: itmax_canopy_fluxes

    select case (s%can_iter_variant)
    case (variant_default_bucketed)
      call can_iter_cpu_default_bucketed(s, params, opts, itmax_canopy_fluxes)
    case default
      call can_iter_cpu_general(s, params, opts, itmax_canopy_fluxes)
    end select
  end subroutine can_iter_cpu_bucketed

  integer function select_can_iter_variant(s, opts) result(variant)
    type(can_iter_state), intent(in) :: s
    type(can_iter_options), intent(in) :: opts

    variant = variant_general

    if (opts%use_lch4 .and. .not. opts%use_hydrstress .and. opts%use_biomass_heat_storage &
         .and. opts%use_undercanopy_stability .and. .not. opts%use_fates &
         .and. .not. opts%use_cn .and. .not. opts%use_c13 &
         .and. .not. opts%do_soilevap_beta .and. opts%do_soil_resistance_sl14 &
         .and. opts%dtime == 1800.0_r8 .and. is_default_bucketed_state(s)) then
      variant = variant_default_bucketed
    end if
  end function select_can_iter_variant

  logical function is_default_bucketed_state(s) result(is_valid)
    type(can_iter_state), intent(in) :: s
    integer :: p

    is_valid = .true.
    do p = 1, s%fn0
      if (s%filterp(p) /= p .or. s%column(p) /= p .or. s%gridcell(p) /= p &
           .or. s%forc_hgt_q_patch(p) /= s%forc_hgt_t_patch(p) &
           .or. s%z0qv(p) /= s%z0hv(p) &
           .or. s%t_veg(p) <= tkfrz + 5.0_r8 &
           .or. s%fdry(p) <= 0.0_r8 .or. s%btran(p) <= btran0) then
        is_valid = .false.
        exit
      end if
    end do
  end function is_default_bucketed_state

  subroutine can_iter_cpu_general(s, params, opts, itmax_canopy_fluxes)
    type(can_iter_state), intent(inout) :: s
    type(can_iter_params), intent(in) :: params
    type(can_iter_options), intent(in) :: opts
    integer, intent(in) :: itmax_canopy_fluxes

    integer :: itlef, fn, fnorig, fnold
    integer :: f, p, c, g
    real(r8) :: cf, w, csoilb, ri, ricsoilc, csoilcn
    real(r8) :: wta, wtl, wtstem, wtshi, wtg0
    real(r8) :: rppdry, rpp, efpot, h2ocan
    real(r8) :: wtaq, wtlq, wtgq0, wtsqi, wtgaq
    real(r8) :: snow_depth_c, fsno_dl, elai_dl, rdl
    real(r8) :: dc1, dc2, efsh, erre, efeold, lw_grnd
    real(r8) :: dels, ecidif, tstar, qstar, thvstar, wc

    itlef = 0
    fn = s%fn0
    fnorig = fn
    s%fporig(1:fn) = s%filterp(1:fn)

    do while (itlef <= itmax_canopy_fluxes .and. fn > 0)
      call FrictionVelocityGeneral(s, params, fn, s%filterp(1:fn), itlef + 1)

      do f = 1, fn
        p = s%filterp(f)
        c = s%column(p)
        g = s%gridcell(p)

        s%tlbef(p) = s%t_veg(p)
        s%del2(p) = s%del(p)

        s%ram1(p) = 1.0_r8 / (s%ustar(p) * s%ustar(p) / s%um(p))
        s%rah(p, above_canopy) = 1.0_r8 / (s%temp1(p) * s%ustar(p))
        s%raw(p, above_canopy) = 1.0_r8 / (s%temp2(p) * s%ustar(p))

        s%uaf(p) = s%um(p) * sqrt(1.0_r8 / (s%ram1(p) * s%um(p)))
        s%uuc(p) = min(0.4_r8, 0.03_r8 * s%um(p) / s%ustar(p))

        if (.not. s%is_fates(p)) then
          s%dleaf_patch(p) = 0.025_r8 + 0.010_r8 * real(s%itype(p), r8)
        end if

        cf = params%cv / (sqrt(s%uaf(p)) * sqrt(s%dleaf_patch(p)))
        s%rb(p) = 1.0_r8 / (cf * s%uaf(p))
        s%rb1(p) = s%rb(p)

        w = exp(-(s%elai(p) + s%esai(p)))
        csoilb = vkc / (params%a_coef * (s%z0mg(c) * s%uaf(p) / nu_param)**params%a_exp)

        ri = (grav * s%htop(p) * (s%taf(p) - s%t_grnd(c))) / (s%taf(p) * s%uaf(p)**2.0_r8)
        if (opts%use_undercanopy_stability .and. (s%taf(p) - s%t_grnd(c)) > 0.0_r8) then
          ricsoilc = params%csoilc / (1.0_r8 + ria * min(ri, 10.0_r8))
          csoilcn = csoilb * w + ricsoilc * (1.0_r8 - w)
        else
          csoilcn = csoilb * w + params%csoilc * (1.0_r8 - w)
        end if

        if (opts%use_biomass_heat_storage) then
          s%rah(p, below_canopy) = 1.0_r8 / (csoilcn * s%uuc(p))
        else
          s%rah(p, below_canopy) = 1.0_r8 / (csoilcn * s%uaf(p))
        end if

        s%raw(p, below_canopy) = s%rah(p, below_canopy)
        if (opts%use_lch4) then
          s%grnd_ch4_cond(p) = 1.0_r8 / (s%raw(p, above_canopy) + s%raw(p, below_canopy))
        end if

        s%svpts(p) = s%el(p)
        s%eah(p) = s%forc_pbot(c) * s%qaf(p) / 0.622_r8
        s%rhaf(p) = s%eah(p) / s%svpts(p)
        s%rah1(p) = s%rah(p, above_canopy)
        s%raw1(p) = s%raw(p, above_canopy)
        s%rah2(p) = s%rah(p, below_canopy)
        s%raw2(p) = s%raw(p, below_canopy)
        s%vpd(p) = max((s%svpts(p) - s%eah(p)), 50.0_r8) * 0.001_r8
      end do

      if (opts%use_fates) then
        call PhotosynthesisStub(s, fn, s%filterp(1:fn), 'sun')
        call PhotosynthesisStub(s, fn, s%filterp(1:fn), 'sha')
      else
        if (opts%use_hydrstress) then
          call PhotosynthesisHydraulicStressStub(s, fn, s%filterp(1:fn))
        else
          call PhotosynthesisStub(s, fn, s%filterp(1:fn), 'sun')
        end if

        if (opts%use_cn .and. opts%use_c13) then
          call FractionationStub(s, fn, s%filterp(1:fn), 'sun')
        end if

        if (.not. opts%use_hydrstress) then
          call PhotosynthesisStub(s, fn, s%filterp(1:fn), 'sha')
        end if

        if (opts%use_cn .and. opts%use_c13) then
          call FractionationStub(s, fn, s%filterp(1:fn), 'sha')
        end if
      end if

      do f = 1, fn
        p = s%filterp(f)
        c = s%column(p)
        g = s%gridcell(p)

        wta = 1.0_r8 / s%rah(p, above_canopy)
        wtl = s%sa_leaf(p) / s%rb(p)
        s%wtg(p) = 1.0_r8 / s%rah(p, below_canopy)
        wtstem = s%sa_stem(p) / (s%rstem(p) + s%rb(p))

        wtshi = 1.0_r8 / (wta + wtl + wtstem + s%wtg(p))

        s%wtl0(p) = wtl * wtshi
        wtg0 = s%wtg(p) * wtshi
        s%wta0(p) = wta * wtshi

        s%wtstem0(p) = wtstem * wtshi
        s%wtga(p) = s%wta0(p) + wtg0 + s%wtstem0(p)
        s%wtal(p) = s%wta0(p) + s%wtl0(p) + s%wtstem0(p)

        s%lw_stem(p) = s%sa_internal(p) * s%emv(p) * sb * s%t_stem(p)**4
        s%lw_leaf(p) = s%sa_internal(p) * s%emv(p) * sb * s%t_veg(p)**4

        if (s%fdry(p) > 0.0_r8) then
          rppdry = s%fdry(p) * s%rb(p) * &
               (s%laisun(p) / (s%rb(p) + s%rssun(p)) + s%laisha(p) / (s%rb(p) + s%rssha(p))) &
               / s%elai(p)
        else
          rppdry = 0.0_r8
        end if

        if (opts%use_lch4) then
          s%canopy_cond(p) = (s%laisun(p) / (s%rb(p) + s%rssun(p)) &
               + s%laisha(p) / (s%rb(p) + s%rssha(p))) / max(s%elai(p), 0.01_r8)
        end if

        efpot = s%forc_rho(c) * ((s%elai(p) + s%esai(p)) / s%rb(p)) * (s%qsatl(p) - s%qaf(p))
        h2ocan = s%liqcan(p) + s%snocan(p)

        if (opts%use_hydrstress) then
          if (efpot > 0.0_r8) then
            if (s%btran(p) > btran0) then
              rpp = rppdry + s%fwet(p)
            else
              rpp = s%fwet(p)
            end if
            rpp = min(rpp, (s%qflx_tran_veg(p) + h2ocan / opts%dtime) / efpot)
          else
            rpp = 1.0_r8
          end if
        else
          if (efpot > 0.0_r8) then
            if (s%btran(p) > btran0) then
              s%qflx_tran_veg(p) = efpot * rppdry
              rpp = rppdry + s%fwet(p)
            else
              rpp = s%fwet(p)
              s%qflx_tran_veg(p) = 0.0_r8
            end if
            rpp = min(rpp, (s%qflx_tran_veg(p) + h2ocan / opts%dtime) / efpot)
          else
            rpp = 1.0_r8
            s%qflx_tran_veg(p) = 0.0_r8
          end if
        end if

        wtaq = s%frac_veg_nosno(p) / s%raw(p, above_canopy)
        wtlq = s%frac_veg_nosno(p) * (s%elai(p) + s%esai(p)) / s%rb(p) * rpp

        snow_depth_c = params%z_dl
        fsno_dl = s%snow_depth(c) / snow_depth_c
        elai_dl = params%lai_dl * (1.0_r8 - min(fsno_dl, 1.0_r8))
        rdl = (1.0_r8 - exp(-elai_dl)) / (0.004_r8 * s%uaf(p))

        if (s%delq(p) < 0.0_r8) then
          s%wtgq(p) = s%frac_veg_nosno(p) / (s%raw(p, below_canopy) + rdl)
        else
          if (opts%do_soilevap_beta) then
            s%wtgq(p) = s%soilbeta(c) * s%frac_veg_nosno(p) / (s%raw(p, below_canopy) + rdl)
          end if
          if (opts%do_soil_resistance_sl14) then
            s%wtgq(p) = s%frac_veg_nosno(p) / (s%raw(p, below_canopy) + s%soilresis(c))
          end if
        end if

        wtsqi = 1.0_r8 / (wtaq + wtlq + s%wtgq(p))

        wtgq0 = s%wtgq(p) * wtsqi
        s%wtlq0(p) = wtlq * wtsqi
        s%wtaq0(p) = wtaq * wtsqi

        wtgaq = s%wtaq0(p) + wtgq0
        s%wtalq(p) = s%wtaq0(p) + s%wtlq0(p)

        dc1 = s%forc_rho(c) * cpair * wtl
        dc2 = hvap * s%forc_rho(c) * wtlq

        efsh = dc1 * (s%wtga(p) * s%t_veg(p) - wtg0 * s%t_grnd(c) &
             - s%wta0(p) * s%thm(p) - s%wtstem0(p) * s%t_stem(p))
        s%eflx_sh_stem(p) = s%forc_rho(c) * cpair * wtstem &
             * ((s%wta0(p) + wtg0 + s%wtl0(p)) * s%t_stem(p) &
             - wtg0 * s%t_grnd(c) - s%wta0(p) * s%thm(p) - s%wtl0(p) * s%t_veg(p))
        s%efe(p) = dc2 * (wtgaq * s%qsatl(p) - wtgq0 * s%qg(c) - s%wtaq0(p) * s%forc_q(c))

        erre = 0.0_r8
        if (s%efe(p) * s%efeb(p) < 0.0_r8) then
          efeold = s%efe(p)
          s%efe(p) = 0.1_r8 * efeold
          erre = s%efe(p) - efeold
        end if

        lw_grnd = s%frac_sno(c) * s%t_soisno_snlp1(c)**4 &
             + (1.0_r8 - s%frac_sno(c) - s%frac_h2osfc(c)) * s%t_soisno_1(c)**4 &
             + s%frac_h2osfc(c) * s%t_h2osfc(c)**4

        s%dt_veg(p) = ((1.0_r8 - s%frac_rad_abs_by_stem(p)) &
             * (s%sabv(p) + s%air(p) + s%bir(p) * s%t_veg(p)**4 + s%cir(p) * lw_grnd) &
             - efsh - s%efe(p) - s%lw_leaf(p) + s%lw_stem(p) &
             - (s%cp_leaf(p) / opts%dtime) * (s%t_veg(p) - s%tl_ini(p))) &
             / ((1.0_r8 - s%frac_rad_abs_by_stem(p)) * (-4.0_r8 * s%bir(p) * s%t_veg(p)**3) &
             + 4.0_r8 * s%sa_internal(p) * s%emv(p) * sb * s%t_veg(p)**3 &
             + dc1 * s%wtga(p) + dc2 * wtgaq * s%qsatldT(p) + s%cp_leaf(p) / opts%dtime)

        s%t_veg(p) = s%tlbef(p) + s%dt_veg(p)

        dels = s%dt_veg(p)
        s%del(p) = abs(dels)
        s%err(p) = 0.0_r8
        if (s%del(p) > delmax) then
          s%dt_veg(p) = delmax * dels / s%del(p)
          s%t_veg(p) = s%tlbef(p) + s%dt_veg(p)
          s%err(p) = (1.0_r8 - s%frac_rad_abs_by_stem(p)) &
               * (s%sabv(p) + s%air(p) + s%bir(p) * s%tlbef(p)**3 &
               * (s%tlbef(p) + 4.0_r8 * s%dt_veg(p)) + s%cir(p) * lw_grnd) &
               - s%sa_internal(p) * s%emv(p) * sb * s%tlbef(p)**3 &
               * (s%tlbef(p) + 4.0_r8 * s%dt_veg(p)) + s%lw_stem(p) &
               - (efsh + dc1 * s%wtga(p) * s%dt_veg(p)) &
               - (s%efe(p) + dc2 * wtgaq * s%qsatldT(p) * s%dt_veg(p)) &
               - (s%cp_leaf(p) / opts%dtime) * (s%t_veg(p) - s%tl_ini(p))
        end if

        efpot = s%forc_rho(c) * ((s%elai(p) + s%esai(p)) / s%rb(p)) &
             * (wtgaq * (s%qsatl(p) + s%qsatldT(p) * s%dt_veg(p)) &
             - wtgq0 * s%qg(c) - s%wtaq0(p) * s%forc_q(c))
        s%qflx_evap_veg(p) = rpp * efpot

        if (opts%use_hydrstress) then
          ecidif = max(0.0_r8, s%qflx_evap_veg(p) - s%qflx_tran_veg(p) - h2ocan / opts%dtime)
          s%qflx_evap_veg(p) = min(s%qflx_evap_veg(p), s%qflx_tran_veg(p) + h2ocan / opts%dtime)
        else
          ecidif = 0.0_r8
          if (efpot > 0.0_r8 .and. s%btran(p) > btran0) then
            s%qflx_tran_veg(p) = efpot * rppdry
          else
            s%qflx_tran_veg(p) = 0.0_r8
          end if
          ecidif = max(0.0_r8, s%qflx_evap_veg(p) - s%qflx_tran_veg(p) - h2ocan / opts%dtime)
          s%qflx_evap_veg(p) = min(s%qflx_evap_veg(p), s%qflx_tran_veg(p) + h2ocan / opts%dtime)
        end if

        s%eflx_sh_veg(p) = efsh + dc1 * s%wtga(p) * s%dt_veg(p) + s%err(p) + erre + hvap * ecidif

        s%eflx_sh_stem(p) = s%eflx_sh_stem(p) &
             + s%forc_rho(c) * cpair * wtstem * (-s%wtl0(p) * s%dt_veg(p))
        s%lw_leaf(p) = s%sa_internal(p) * s%emv(p) * sb * s%tlbef(p)**3 &
             * (s%tlbef(p) + 4.0_r8 * s%dt_veg(p))

        call QSat(s%t_veg(p), s%forc_pbot(c), s%qsatl(p), es=s%el(p), qsdT=s%qsatldT(p))

        s%taf(p) = wtg0 * s%t_grnd(c) + s%wta0(p) * s%thm(p) &
             + s%wtl0(p) * s%t_veg(p) + s%wtstem0(p) * s%t_stem(p)
        s%qaf(p) = s%wtlq0(p) * s%qsatl(p) + wtgq0 * s%qg(c) + s%forc_q(c) * s%wtaq0(p)

        s%dth(p) = s%thm(p) - s%taf(p)
        s%dqh(p) = s%forc_q(c) - s%qaf(p)
        s%delq(p) = s%wtalq(p) * s%qg(c) - s%wtlq0(p) * s%qsatl(p) - s%wtaq0(p) * s%forc_q(c)

        tstar = s%temp1(p) * s%dth(p)
        qstar = s%temp2(p) * s%dqh(p)

        thvstar = tstar * (1.0_r8 + 0.61_r8 * s%forc_q(c)) + 0.61_r8 * s%forc_th(c) * qstar
        s%zeta(p) = s%zldis(p) * vkc * grav * thvstar / (s%ustar(p)**2 * s%thv(c))

        if (s%zeta(p) >= 0.0_r8) then
          s%zeta(p) = min(params%zetamaxstable, max(s%zeta(p), 0.01_r8))
          s%um(p) = max(s%ur(p), 0.1_r8)
        else
          s%zeta(p) = max(-100.0_r8, min(s%zeta(p), -0.01_r8))
          if (s%ustar(p) * thvstar > 0.0_r8) then
            wc = 0.0_r8
          else
            wc = beta * (-grav * s%ustar(p) * thvstar * zii / s%thv(c))**0.333_r8
          end if
          s%um(p) = sqrt(s%ur(p) * s%ur(p) + wc * wc)
        end if
        s%obu(p) = s%zldis(p) / s%zeta(p)

        if (s%obuold(p) * s%obu(p) < 0.0_r8) s%nmozsgn(p) = s%nmozsgn(p) + 1
        if (s%nmozsgn(p) >= 4) s%obu(p) = s%zldis(p) / (-0.01_r8)
        s%obuold(p) = s%obu(p)
      end do

      itlef = itlef + 1
      if (itlef > itmin) then
        do f = 1, fn
          p = s%filterp(f)
          s%dele(p) = abs(s%efe(p) - s%efeb(p))
          s%efeb(p) = s%efe(p)
          s%det(p) = max(s%del(p), s%del2(p))
          s%num_iter(p) = real(itlef, r8)
        end do

        fnold = fn
        fn = 0
        do f = 1, fnold
          p = s%filterp(f)
          if (.not. (s%det(p) < dtmin .and. s%dele(p) < dlemin)) then
            fn = fn + 1
            s%filterp(fn) = p
          end if
        end do
      end if
    end do

    s%filterp(1:fnorig) = s%fporig(1:fnorig)
  end subroutine can_iter_cpu_general

  subroutine can_iter_cpu_default_bucketed(s, params, opts, itmax_canopy_fluxes)
    type(can_iter_state), intent(inout) :: s
    type(can_iter_params), intent(in) :: params
    type(can_iter_options), intent(in) :: opts
    integer, intent(in) :: itmax_canopy_fluxes

    integer :: itlef, fn, fnorig, fnold
    integer :: f, p
    integer :: n_under_stable, n_under_neutral
    integer :: n_delq_neg, n_delq_nonneg

    itlef = 0
    fn = s%fn0
    fnorig = fn
    s%fporig(1:fn) = s%filterp(1:fn)

    do while (itlef <= itmax_canopy_fluxes .and. fn > 0)
      call FrictionVelocityDefault(s, params, fn, s%filterp(1:fn), itlef + 1)

      call build_understory_buckets(s, fn, s%filterp(1:fn), n_under_stable, n_under_neutral)
      if (n_under_stable > 0) then
        call can_iter_first_bucket(s, params, opts, n_under_stable, &
             s%filter_under_stable(1:n_under_stable), .true.)
      end if
      if (n_under_neutral > 0) then
        call can_iter_first_bucket(s, params, opts, n_under_neutral, &
             s%filter_under_neutral(1:n_under_neutral), .false.)
      end if

      call build_delq_buckets(s, fn, s%filterp(1:fn), n_delq_neg, n_delq_nonneg)
      if (n_delq_neg > 0) then
        call can_iter_wtgq_bucket(s, params, n_delq_neg, s%filter_delq_neg(1:n_delq_neg), .true.)
      end if
      if (n_delq_nonneg > 0) then
        call can_iter_wtgq_bucket(s, params, n_delq_nonneg, s%filter_delq_nonneg(1:n_delq_nonneg), .false.)
      end if
      call can_iter_main_after_wtgq(s, params, opts, fn, s%filterp(1:fn))

      itlef = itlef + 1
      if (itlef > itmin) then
        do f = 1, fn
          p = s%filterp(f)
          s%dele(p) = abs(s%efe(p) - s%efeb(p))
          s%efeb(p) = s%efe(p)
          s%det(p) = max(s%del(p), s%del2(p))
          s%num_iter(p) = real(itlef, r8)
        end do

        fnold = fn
        fn = 0
        do f = 1, fnold
          p = s%filterp(f)
          if (.not. (s%det(p) < dtmin .and. s%dele(p) < dlemin)) then
            fn = fn + 1
            s%filterp(fn) = p
          end if
        end do
      end if
    end do

    s%filterp(1:fnorig) = s%fporig(1:fnorig)
  end subroutine can_iter_cpu_default_bucketed

  subroutine build_understory_buckets(s, fn, filterp, n_stable, n_neutral)
    type(can_iter_state), intent(inout) :: s
    integer, intent(in) :: fn, filterp(fn)
    integer, intent(out) :: n_stable, n_neutral
    integer :: f, p

    n_stable = 0
    n_neutral = 0
    do f = 1, fn
      p = filterp(f)
      if ((s%taf(p) - s%t_grnd(p)) > 0.0_r8) then
        n_stable = n_stable + 1
        s%filter_under_stable(n_stable) = p
      else
        n_neutral = n_neutral + 1
        s%filter_under_neutral(n_neutral) = p
      end if
    end do
  end subroutine build_understory_buckets

  subroutine build_delq_buckets(s, fn, filterp, n_neg, n_nonneg)
    type(can_iter_state), intent(inout) :: s
    integer, intent(in) :: fn, filterp(fn)
    integer, intent(out) :: n_neg, n_nonneg
    integer :: f, p

    n_neg = 0
    n_nonneg = 0
    do f = 1, fn
      p = filterp(f)
      if (s%delq(p) < 0.0_r8) then
        n_neg = n_neg + 1
        s%filter_delq_neg(n_neg) = p
      else
        n_nonneg = n_nonneg + 1
        s%filter_delq_nonneg(n_nonneg) = p
      end if
    end do
  end subroutine build_delq_buckets

  subroutine can_iter_first_bucket(s, params, opts, n, filterp, under_stable)
    type(can_iter_state), intent(inout) :: s
    type(can_iter_params), intent(in) :: params
    type(can_iter_options), intent(in) :: opts
    integer, intent(in) :: n, filterp(n)
    logical, intent(in) :: under_stable

    integer :: f, p, c
    real(r8) :: cf, w, csoilb, ri, ricsoilc, csoilcn
    real(r8) :: vpd_kpa, rh_can, light, nitrogen, temp_stress, co2_fac, conductance

    if (under_stable) then
      !DIR$ IVDEP
      !DIR$ VECTOR ALWAYS
      do f = 1, n
        p = filterp(f)
        c = p

        s%tlbef(p) = s%t_veg(p)
        s%del2(p) = s%del(p)

        s%ram1(p) = 1.0_r8 / (s%ustar(p) * s%ustar(p) / s%um(p))
        s%rah(p, above_canopy) = 1.0_r8 / (s%temp1(p) * s%ustar(p))
        s%raw(p, above_canopy) = 1.0_r8 / (s%temp2(p) * s%ustar(p))

        s%uaf(p) = s%um(p) * sqrt(1.0_r8 / (s%ram1(p) * s%um(p)))
        s%uuc(p) = min(0.4_r8, 0.03_r8 * s%um(p) / s%ustar(p))
        s%dleaf_patch(p) = 0.025_r8 + 0.010_r8 * real(s%itype(p), r8)

        cf = params%cv / (sqrt(s%uaf(p)) * sqrt(s%dleaf_patch(p)))
        s%rb(p) = 1.0_r8 / (cf * s%uaf(p))
        s%rb1(p) = s%rb(p)

        w = exp(-(s%elai(p) + s%esai(p)))
        csoilb = vkc / (params%a_coef * (s%z0mg(c) * s%uaf(p) / nu_param)**params%a_exp)
        ri = (grav * s%htop(p) * (s%taf(p) - s%t_grnd(c))) / (s%taf(p) * s%uaf(p)**2.0_r8)
        ricsoilc = params%csoilc / (1.0_r8 + ria * min(ri, 10.0_r8))
        csoilcn = csoilb * w + ricsoilc * (1.0_r8 - w)

        s%rah(p, below_canopy) = 1.0_r8 / (csoilcn * s%uuc(p))
        s%raw(p, below_canopy) = s%rah(p, below_canopy)
        s%grnd_ch4_cond(p) = 1.0_r8 / (s%raw(p, above_canopy) + s%raw(p, below_canopy))

        s%svpts(p) = s%el(p)
        s%eah(p) = s%forc_pbot(c) * s%qaf(p) / 0.622_r8
        s%rhaf(p) = s%eah(p) / s%svpts(p)
        s%rah1(p) = s%rah(p, above_canopy)
        s%raw1(p) = s%raw(p, above_canopy)
        s%rah2(p) = s%rah(p, below_canopy)
        s%raw2(p) = s%raw(p, below_canopy)
        s%vpd(p) = max((s%svpts(p) - s%eah(p)), 50.0_r8) * 0.001_r8

        vpd_kpa = max((s%svpts(p) - s%eah(p)) * 0.001_r8, 0.05_r8)
        rh_can = min(1.0_r8, max(0.05_r8, s%eah(p) / max(s%svpts(p), 1.0_r8)))
        light = max(0.05_r8, s%dayl_factor(p))
        nitrogen = max(0.25_r8, min(2.0_r8, s%leafn_patch(p) / 2.0_r8))
        temp_stress = exp(-((s%t_veg(p) - 298.0_r8) / 18.0_r8)**2)
        co2_fac = max(0.5_r8, min(1.5_r8, s%co2(p) / 40.0_r8))

        conductance = 0.0030_r8 + 0.0180_r8 * light * nitrogen * temp_stress &
             * co2_fac * max(s%btran(p), 0.05_r8) * rh_can / (1.0_r8 + 0.10_r8 * vpd_kpa)
        s%rssun(p) = 1.0_r8 / max(conductance, 1.0e-4_r8)
        conductance = 0.0020_r8 + 0.0100_r8 * light * nitrogen * temp_stress &
             * co2_fac * max(s%btran(p), 0.05_r8) * rh_can / (1.0_r8 + 0.14_r8 * vpd_kpa)
        s%rssha(p) = 1.0_r8 / max(conductance, 1.0e-4_r8)
      end do
    else
      !DIR$ IVDEP
      !DIR$ VECTOR ALWAYS
      do f = 1, n
        p = filterp(f)
        c = p

        s%tlbef(p) = s%t_veg(p)
        s%del2(p) = s%del(p)

        s%ram1(p) = 1.0_r8 / (s%ustar(p) * s%ustar(p) / s%um(p))
        s%rah(p, above_canopy) = 1.0_r8 / (s%temp1(p) * s%ustar(p))
        s%raw(p, above_canopy) = 1.0_r8 / (s%temp2(p) * s%ustar(p))

        s%uaf(p) = s%um(p) * sqrt(1.0_r8 / (s%ram1(p) * s%um(p)))
        s%uuc(p) = min(0.4_r8, 0.03_r8 * s%um(p) / s%ustar(p))
        s%dleaf_patch(p) = 0.025_r8 + 0.010_r8 * real(s%itype(p), r8)

        cf = params%cv / (sqrt(s%uaf(p)) * sqrt(s%dleaf_patch(p)))
        s%rb(p) = 1.0_r8 / (cf * s%uaf(p))
        s%rb1(p) = s%rb(p)

        w = exp(-(s%elai(p) + s%esai(p)))
        csoilb = vkc / (params%a_coef * (s%z0mg(c) * s%uaf(p) / nu_param)**params%a_exp)
        csoilcn = csoilb * w + params%csoilc * (1.0_r8 - w)

        s%rah(p, below_canopy) = 1.0_r8 / (csoilcn * s%uuc(p))
        s%raw(p, below_canopy) = s%rah(p, below_canopy)
        s%grnd_ch4_cond(p) = 1.0_r8 / (s%raw(p, above_canopy) + s%raw(p, below_canopy))

        s%svpts(p) = s%el(p)
        s%eah(p) = s%forc_pbot(c) * s%qaf(p) / 0.622_r8
        s%rhaf(p) = s%eah(p) / s%svpts(p)
        s%rah1(p) = s%rah(p, above_canopy)
        s%raw1(p) = s%raw(p, above_canopy)
        s%rah2(p) = s%rah(p, below_canopy)
        s%raw2(p) = s%raw(p, below_canopy)
        s%vpd(p) = max((s%svpts(p) - s%eah(p)), 50.0_r8) * 0.001_r8

        vpd_kpa = max((s%svpts(p) - s%eah(p)) * 0.001_r8, 0.05_r8)
        rh_can = min(1.0_r8, max(0.05_r8, s%eah(p) / max(s%svpts(p), 1.0_r8)))
        light = max(0.05_r8, s%dayl_factor(p))
        nitrogen = max(0.25_r8, min(2.0_r8, s%leafn_patch(p) / 2.0_r8))
        temp_stress = exp(-((s%t_veg(p) - 298.0_r8) / 18.0_r8)**2)
        co2_fac = max(0.5_r8, min(1.5_r8, s%co2(p) / 40.0_r8))

        conductance = 0.0030_r8 + 0.0180_r8 * light * nitrogen * temp_stress &
             * co2_fac * max(s%btran(p), 0.05_r8) * rh_can / (1.0_r8 + 0.10_r8 * vpd_kpa)
        s%rssun(p) = 1.0_r8 / max(conductance, 1.0e-4_r8)
        conductance = 0.0020_r8 + 0.0100_r8 * light * nitrogen * temp_stress &
             * co2_fac * max(s%btran(p), 0.05_r8) * rh_can / (1.0_r8 + 0.14_r8 * vpd_kpa)
        s%rssha(p) = 1.0_r8 / max(conductance, 1.0e-4_r8)
      end do
    end if
  end subroutine can_iter_first_bucket

  subroutine can_iter_wtgq_bucket(s, params, n, filterp, delq_negative)
    type(can_iter_state), intent(inout) :: s
    type(can_iter_params), intent(in) :: params
    integer, intent(in) :: n, filterp(n)
    logical, intent(in) :: delq_negative

    integer :: f, p, c
    real(r8) :: snow_depth_c, fsno_dl, elai_dl, rdl

    if (delq_negative) then
      snow_depth_c = params%z_dl

      !DIR$ IVDEP
      !DIR$ VECTOR ALWAYS
      do f = 1, n
        p = filterp(f)
        c = p

        fsno_dl = s%snow_depth(c) / snow_depth_c
        elai_dl = params%lai_dl * (1.0_r8 - min(fsno_dl, 1.0_r8))
        rdl = (1.0_r8 - exp(-elai_dl)) / (0.004_r8 * s%uaf(p))
        s%wtgq(p) = s%frac_veg_nosno(p) / (s%raw(p, below_canopy) + rdl)
      end do
    else
      !DIR$ IVDEP
      !DIR$ VECTOR ALWAYS
      do f = 1, n
        p = filterp(f)
        c = p

        s%wtgq(p) = s%frac_veg_nosno(p) / (s%raw(p, below_canopy) + s%soilresis(c))
      end do
    end if
  end subroutine can_iter_wtgq_bucket

  subroutine can_iter_main_after_wtgq(s, params, opts, n, filterp)
    type(can_iter_state), intent(inout) :: s
    type(can_iter_params), intent(in) :: params
    type(can_iter_options), intent(in) :: opts
    integer, intent(in) :: n, filterp(n)

    integer :: f, p, c
    real(r8) :: wta, wtl, wtstem, wtshi, wtg0
    real(r8) :: rppdry, rpp, efpot, h2ocan
    real(r8) :: wtaq, wtlq, wtgq0, wtsqi, wtgaq
    real(r8) :: dc1, dc2, efsh, erre, efeold, lw_grnd
    real(r8) :: dels, ecidif, tstar, qstar, thvstar, wc
    real(r8) :: td, es_local, esdT_local, vp, vp1, vp2

    !DIR$ IVDEP
    !DIR$ VECTOR ALWAYS
    do f = 1, n
      p = filterp(f)
      c = p

      wta = 1.0_r8 / s%rah(p, above_canopy)
      wtl = s%sa_leaf(p) / s%rb(p)
      s%wtg(p) = 1.0_r8 / s%rah(p, below_canopy)
      wtstem = s%sa_stem(p) / (s%rstem(p) + s%rb(p))

      wtshi = 1.0_r8 / (wta + wtl + wtstem + s%wtg(p))

      s%wtl0(p) = wtl * wtshi
      wtg0 = s%wtg(p) * wtshi
      s%wta0(p) = wta * wtshi

      s%wtstem0(p) = wtstem * wtshi
      s%wtga(p) = s%wta0(p) + wtg0 + s%wtstem0(p)
      s%wtal(p) = s%wta0(p) + s%wtl0(p) + s%wtstem0(p)

      s%lw_stem(p) = s%sa_internal(p) * s%emv(p) * sb * s%t_stem(p)**4
      s%lw_leaf(p) = s%sa_internal(p) * s%emv(p) * sb * s%t_veg(p)**4

      rppdry = s%fdry(p) * s%rb(p) * &
           (s%laisun(p) / (s%rb(p) + s%rssun(p)) + s%laisha(p) / (s%rb(p) + s%rssha(p))) &
           / s%elai(p)

      s%canopy_cond(p) = (s%laisun(p) / (s%rb(p) + s%rssun(p)) &
           + s%laisha(p) / (s%rb(p) + s%rssha(p))) / max(s%elai(p), 0.01_r8)

      efpot = s%forc_rho(c) * ((s%elai(p) + s%esai(p)) / s%rb(p)) * (s%qsatl(p) - s%qaf(p))
      h2ocan = s%liqcan(p) + s%snocan(p)

      if (efpot > 0.0_r8) then
        s%qflx_tran_veg(p) = efpot * rppdry
        rpp = rppdry + s%fwet(p)
        rpp = min(rpp, (s%qflx_tran_veg(p) + h2ocan / opts%dtime) / efpot)
      else
        rpp = 1.0_r8
        s%qflx_tran_veg(p) = 0.0_r8
      end if

      wtaq = s%frac_veg_nosno(p) / s%raw(p, above_canopy)
      wtlq = s%frac_veg_nosno(p) * (s%elai(p) + s%esai(p)) / s%rb(p) * rpp
      wtsqi = 1.0_r8 / (wtaq + wtlq + s%wtgq(p))

      wtgq0 = s%wtgq(p) * wtsqi
      s%wtlq0(p) = wtlq * wtsqi
      s%wtaq0(p) = wtaq * wtsqi

      wtgaq = s%wtaq0(p) + wtgq0
      s%wtalq(p) = s%wtaq0(p) + s%wtlq0(p)

      dc1 = s%forc_rho(c) * cpair * wtl
      dc2 = hvap * s%forc_rho(c) * wtlq

      efsh = dc1 * (s%wtga(p) * s%t_veg(p) - wtg0 * s%t_grnd(c) &
           - s%wta0(p) * s%thm(p) - s%wtstem0(p) * s%t_stem(p))
      s%eflx_sh_stem(p) = s%forc_rho(c) * cpair * wtstem &
           * ((s%wta0(p) + wtg0 + s%wtl0(p)) * s%t_stem(p) &
           - wtg0 * s%t_grnd(c) - s%wta0(p) * s%thm(p) - s%wtl0(p) * s%t_veg(p))
      s%efe(p) = dc2 * (wtgaq * s%qsatl(p) - wtgq0 * s%qg(c) - s%wtaq0(p) * s%forc_q(c))

      erre = 0.0_r8
      if (s%efe(p) * s%efeb(p) < 0.0_r8) then
        efeold = s%efe(p)
        s%efe(p) = 0.1_r8 * efeold
        erre = s%efe(p) - efeold
      end if

      lw_grnd = s%frac_sno(c) * s%t_soisno_snlp1(c)**4 &
           + (1.0_r8 - s%frac_sno(c) - s%frac_h2osfc(c)) * s%t_soisno_1(c)**4 &
           + s%frac_h2osfc(c) * s%t_h2osfc(c)**4

      s%dt_veg(p) = ((1.0_r8 - s%frac_rad_abs_by_stem(p)) &
           * (s%sabv(p) + s%air(p) + s%bir(p) * s%t_veg(p)**4 + s%cir(p) * lw_grnd) &
           - efsh - s%efe(p) - s%lw_leaf(p) + s%lw_stem(p) &
           - (s%cp_leaf(p) / opts%dtime) * (s%t_veg(p) - s%tl_ini(p))) &
           / ((1.0_r8 - s%frac_rad_abs_by_stem(p)) * (-4.0_r8 * s%bir(p) * s%t_veg(p)**3) &
           + 4.0_r8 * s%sa_internal(p) * s%emv(p) * sb * s%t_veg(p)**3 &
           + dc1 * s%wtga(p) + dc2 * wtgaq * s%qsatldT(p) + s%cp_leaf(p) / opts%dtime)

      s%t_veg(p) = s%tlbef(p) + s%dt_veg(p)

      dels = s%dt_veg(p)
      s%del(p) = abs(dels)
      s%err(p) = 0.0_r8
      if (s%del(p) > delmax) then
        s%dt_veg(p) = delmax * dels / s%del(p)
        s%t_veg(p) = s%tlbef(p) + s%dt_veg(p)
        s%err(p) = (1.0_r8 - s%frac_rad_abs_by_stem(p)) &
             * (s%sabv(p) + s%air(p) + s%bir(p) * s%tlbef(p)**3 &
             * (s%tlbef(p) + 4.0_r8 * s%dt_veg(p)) + s%cir(p) * lw_grnd) &
             - s%sa_internal(p) * s%emv(p) * sb * s%tlbef(p)**3 &
             * (s%tlbef(p) + 4.0_r8 * s%dt_veg(p)) + s%lw_stem(p) &
             - (efsh + dc1 * s%wtga(p) * s%dt_veg(p)) &
             - (s%efe(p) + dc2 * wtgaq * s%qsatldT(p) * s%dt_veg(p)) &
             - (s%cp_leaf(p) / opts%dtime) * (s%t_veg(p) - s%tl_ini(p))
      end if

      efpot = s%forc_rho(c) * ((s%elai(p) + s%esai(p)) / s%rb(p)) &
           * (wtgaq * (s%qsatl(p) + s%qsatldT(p) * s%dt_veg(p)) &
           - wtgq0 * s%qg(c) - s%wtaq0(p) * s%forc_q(c))
      s%qflx_evap_veg(p) = rpp * efpot

      ecidif = 0.0_r8
      if (efpot > 0.0_r8) then
        s%qflx_tran_veg(p) = efpot * rppdry
      else
        s%qflx_tran_veg(p) = 0.0_r8
      end if
      ecidif = max(0.0_r8, s%qflx_evap_veg(p) - s%qflx_tran_veg(p) - h2ocan / opts%dtime)
      s%qflx_evap_veg(p) = min(s%qflx_evap_veg(p), s%qflx_tran_veg(p) + h2ocan / opts%dtime)

      s%eflx_sh_veg(p) = efsh + dc1 * s%wtga(p) * s%dt_veg(p) + s%err(p) + erre + hvap * ecidif

      s%eflx_sh_stem(p) = s%eflx_sh_stem(p) &
           + s%forc_rho(c) * cpair * wtstem * (-s%wtl0(p) * s%dt_veg(p))
      s%lw_leaf(p) = s%sa_internal(p) * s%emv(p) * sb * s%tlbef(p)**3 &
           * (s%tlbef(p) + 4.0_r8 * s%dt_veg(p))

      td = min(100.0_r8, max(-75.0_r8, s%t_veg(p) - tkfrz))

      es_local = 6.11213476_r8 + td * (0.444007856_r8 + td * (0.143064234e-1_r8 &
           + td * (0.264461437e-3_r8 + td * (0.305903558e-5_r8 + td * (0.196237241e-7_r8 &
           + td * (0.892344772e-10_r8 + td * (-0.373208410e-12_r8 &
           + td * 0.209339997e-15_r8)))))))
      esdT_local = 0.444017302_r8 + td * (0.286064092e-1_r8 + td * (0.794683137e-3_r8 &
           + td * (0.121211669e-4_r8 + td * (0.103354611e-6_r8 + td * (0.404125005e-9_r8 &
           + td * (-0.788037859e-12_r8 + td * (-0.114596802e-13_r8 &
           + td * 0.381294516e-16_r8)))))))

      s%el(p) = es_local * 100.0_r8
      esdT_local = esdT_local * 100.0_r8
      vp = 1.0_r8 / (s%forc_pbot(c) - 0.378_r8 * s%el(p))
      vp1 = 0.622_r8 * vp
      s%qsatl(p) = s%el(p) * vp1
      vp2 = vp1 * vp
      s%qsatldT(p) = esdT_local * vp2 * s%forc_pbot(c)

      s%taf(p) = wtg0 * s%t_grnd(c) + s%wta0(p) * s%thm(p) &
           + s%wtl0(p) * s%t_veg(p) + s%wtstem0(p) * s%t_stem(p)
      s%qaf(p) = s%wtlq0(p) * s%qsatl(p) + wtgq0 * s%qg(c) + s%forc_q(c) * s%wtaq0(p)

      s%dth(p) = s%thm(p) - s%taf(p)
      s%dqh(p) = s%forc_q(c) - s%qaf(p)
      s%delq(p) = s%wtalq(p) * s%qg(c) - s%wtlq0(p) * s%qsatl(p) - s%wtaq0(p) * s%forc_q(c)

      tstar = s%temp1(p) * s%dth(p)
      qstar = s%temp2(p) * s%dqh(p)

      thvstar = tstar * (1.0_r8 + 0.61_r8 * s%forc_q(c)) + 0.61_r8 * s%forc_th(c) * qstar
      s%zeta(p) = s%zldis(p) * vkc * grav * thvstar / (s%ustar(p)**2 * s%thv(c))

      if (s%zeta(p) >= 0.0_r8) then
        s%zeta(p) = min(params%zetamaxstable, max(s%zeta(p), 0.01_r8))
        s%um(p) = max(s%ur(p), 0.1_r8)
      else
        s%zeta(p) = max(-100.0_r8, min(s%zeta(p), -0.01_r8))
        if (s%ustar(p) * thvstar > 0.0_r8) then
          wc = 0.0_r8
        else
          wc = beta * (-grav * s%ustar(p) * thvstar * zii / s%thv(c))**0.333_r8
        end if
        s%um(p) = sqrt(s%ur(p) * s%ur(p) + wc * wc)
      end if
      s%obu(p) = s%zldis(p) / s%zeta(p)

      if (s%obuold(p) * s%obu(p) < 0.0_r8) s%nmozsgn(p) = s%nmozsgn(p) + 1
      if (s%nmozsgn(p) >= 4) s%obu(p) = s%zldis(p) / (-0.01_r8)
      s%obuold(p) = s%obu(p)
    end do
  end subroutine can_iter_main_after_wtgq

  ! Reference copy of the previous compact-inline fast path. The driver does not
  ! call this routine; it is left here to make the bucketed changes easy to audit.
  subroutine can_iter_cpu_default_inline(s, params, opts, itmax_canopy_fluxes)
    type(can_iter_state), intent(inout) :: s
    type(can_iter_params), intent(in) :: params
    type(can_iter_options), intent(in) :: opts
    integer, intent(in) :: itmax_canopy_fluxes

    integer :: itlef, fn, fnorig, fnold
    integer :: f, p, c
    real(r8) :: cf, w, csoilb, ri, ricsoilc, csoilcn
    real(r8) :: wta, wtl, wtstem, wtshi, wtg0
    real(r8) :: rppdry, rpp, efpot, h2ocan
    real(r8) :: wtaq, wtlq, wtgq0, wtsqi, wtgaq
    real(r8) :: snow_depth_c, fsno_dl, elai_dl, rdl
    real(r8) :: dc1, dc2, efsh, erre, efeold, lw_grnd
    real(r8) :: dels, ecidif, tstar, qstar, thvstar, wc
    real(r8) :: td, es_local, esdT_local, vp, vp1, vp2
    real(r8) :: vpd_kpa, rh_can, light, nitrogen, temp_stress, co2_fac, conductance

    itlef = 0
    fn = s%fn0
    fnorig = fn
    s%fporig(1:fn) = s%filterp(1:fn)

    do while (itlef <= itmax_canopy_fluxes .and. fn > 0)
      call FrictionVelocityDefault(s, params, fn, s%filterp(1:fn), itlef + 1)

      !DIR$ IVDEP
      !DIR$ VECTOR ALWAYS
      do f = 1, fn
        p = s%filterp(f)
        c = p

        s%tlbef(p) = s%t_veg(p)
        s%del2(p) = s%del(p)

        s%ram1(p) = 1.0_r8 / (s%ustar(p) * s%ustar(p) / s%um(p))
        s%rah(p, above_canopy) = 1.0_r8 / (s%temp1(p) * s%ustar(p))
        s%raw(p, above_canopy) = 1.0_r8 / (s%temp2(p) * s%ustar(p))

        s%uaf(p) = s%um(p) * sqrt(1.0_r8 / (s%ram1(p) * s%um(p)))
        s%uuc(p) = min(0.4_r8, 0.03_r8 * s%um(p) / s%ustar(p))

        s%dleaf_patch(p) = 0.025_r8 + 0.010_r8 * real(s%itype(p), r8)

        cf = params%cv / (sqrt(s%uaf(p)) * sqrt(s%dleaf_patch(p)))
        s%rb(p) = 1.0_r8 / (cf * s%uaf(p))
        s%rb1(p) = s%rb(p)

        w = exp(-(s%elai(p) + s%esai(p)))
        csoilb = vkc / (params%a_coef * (s%z0mg(c) * s%uaf(p) / nu_param)**params%a_exp)

        ri = (grav * s%htop(p) * (s%taf(p) - s%t_grnd(c))) / (s%taf(p) * s%uaf(p)**2.0_r8)
        if ((s%taf(p) - s%t_grnd(c)) > 0.0_r8) then
          ricsoilc = params%csoilc / (1.0_r8 + ria * min(ri, 10.0_r8))
          csoilcn = csoilb * w + ricsoilc * (1.0_r8 - w)
        else
          csoilcn = csoilb * w + params%csoilc * (1.0_r8 - w)
        end if

        s%rah(p, below_canopy) = 1.0_r8 / (csoilcn * s%uuc(p))

        s%raw(p, below_canopy) = s%rah(p, below_canopy)
        s%grnd_ch4_cond(p) = 1.0_r8 / (s%raw(p, above_canopy) + s%raw(p, below_canopy))

        s%svpts(p) = s%el(p)
        s%eah(p) = s%forc_pbot(c) * s%qaf(p) / 0.622_r8
        s%rhaf(p) = s%eah(p) / s%svpts(p)
        s%rah1(p) = s%rah(p, above_canopy)
        s%raw1(p) = s%raw(p, above_canopy)
        s%rah2(p) = s%rah(p, below_canopy)
        s%raw2(p) = s%raw(p, below_canopy)
        s%vpd(p) = max((s%svpts(p) - s%eah(p)), 50.0_r8) * 0.001_r8

        vpd_kpa = max((s%svpts(p) - s%eah(p)) * 0.001_r8, 0.05_r8)
        rh_can = min(1.0_r8, max(0.05_r8, s%eah(p) / max(s%svpts(p), 1.0_r8)))
        light = max(0.05_r8, s%dayl_factor(p))
        nitrogen = max(0.25_r8, min(2.0_r8, s%leafn_patch(p) / 2.0_r8))
        temp_stress = exp(-((s%t_veg(p) - 298.0_r8) / 18.0_r8)**2)
        co2_fac = max(0.5_r8, min(1.5_r8, s%co2(p) / 40.0_r8))

        conductance = 0.0030_r8 + 0.0180_r8 * light * nitrogen * temp_stress &
             * co2_fac * max(s%btran(p), 0.05_r8) * rh_can / (1.0_r8 + 0.10_r8 * vpd_kpa)
        s%rssun(p) = 1.0_r8 / max(conductance, 1.0e-4_r8)

        conductance = 0.0020_r8 + 0.0100_r8 * light * nitrogen * temp_stress &
             * co2_fac * max(s%btran(p), 0.05_r8) * rh_can / (1.0_r8 + 0.14_r8 * vpd_kpa)
        s%rssha(p) = 1.0_r8 / max(conductance, 1.0e-4_r8)
      end do

      !DIR$ IVDEP
      !DIR$ VECTOR ALWAYS
      do f = 1, fn
        p = s%filterp(f)
        c = p

        wta = 1.0_r8 / s%rah(p, above_canopy)
        wtl = s%sa_leaf(p) / s%rb(p)
        s%wtg(p) = 1.0_r8 / s%rah(p, below_canopy)
        wtstem = s%sa_stem(p) / (s%rstem(p) + s%rb(p))

        wtshi = 1.0_r8 / (wta + wtl + wtstem + s%wtg(p))

        s%wtl0(p) = wtl * wtshi
        wtg0 = s%wtg(p) * wtshi
        s%wta0(p) = wta * wtshi

        s%wtstem0(p) = wtstem * wtshi
        s%wtga(p) = s%wta0(p) + wtg0 + s%wtstem0(p)
        s%wtal(p) = s%wta0(p) + s%wtl0(p) + s%wtstem0(p)

        s%lw_stem(p) = s%sa_internal(p) * s%emv(p) * sb * s%t_stem(p)**4
        s%lw_leaf(p) = s%sa_internal(p) * s%emv(p) * sb * s%t_veg(p)**4

        rppdry = s%fdry(p) * s%rb(p) * &
             (s%laisun(p) / (s%rb(p) + s%rssun(p)) + s%laisha(p) / (s%rb(p) + s%rssha(p))) &
             / s%elai(p)

        s%canopy_cond(p) = (s%laisun(p) / (s%rb(p) + s%rssun(p)) &
             + s%laisha(p) / (s%rb(p) + s%rssha(p))) / max(s%elai(p), 0.01_r8)

        efpot = s%forc_rho(c) * ((s%elai(p) + s%esai(p)) / s%rb(p)) * (s%qsatl(p) - s%qaf(p))
        h2ocan = s%liqcan(p) + s%snocan(p)

        if (efpot > 0.0_r8) then
          s%qflx_tran_veg(p) = efpot * rppdry
          rpp = rppdry + s%fwet(p)
          rpp = min(rpp, (s%qflx_tran_veg(p) + h2ocan / opts%dtime) / efpot)
        else
          rpp = 1.0_r8
          s%qflx_tran_veg(p) = 0.0_r8
        end if

        wtaq = s%frac_veg_nosno(p) / s%raw(p, above_canopy)
        wtlq = s%frac_veg_nosno(p) * (s%elai(p) + s%esai(p)) / s%rb(p) * rpp

        snow_depth_c = params%z_dl
        fsno_dl = s%snow_depth(c) / snow_depth_c
        elai_dl = params%lai_dl * (1.0_r8 - min(fsno_dl, 1.0_r8))
        rdl = (1.0_r8 - exp(-elai_dl)) / (0.004_r8 * s%uaf(p))

        if (s%delq(p) < 0.0_r8) then
          s%wtgq(p) = s%frac_veg_nosno(p) / (s%raw(p, below_canopy) + rdl)
        else
          s%wtgq(p) = s%frac_veg_nosno(p) / (s%raw(p, below_canopy) + s%soilresis(c))
        end if

        wtsqi = 1.0_r8 / (wtaq + wtlq + s%wtgq(p))

        wtgq0 = s%wtgq(p) * wtsqi
        s%wtlq0(p) = wtlq * wtsqi
        s%wtaq0(p) = wtaq * wtsqi

        wtgaq = s%wtaq0(p) + wtgq0
        s%wtalq(p) = s%wtaq0(p) + s%wtlq0(p)

        dc1 = s%forc_rho(c) * cpair * wtl
        dc2 = hvap * s%forc_rho(c) * wtlq

        efsh = dc1 * (s%wtga(p) * s%t_veg(p) - wtg0 * s%t_grnd(c) &
             - s%wta0(p) * s%thm(p) - s%wtstem0(p) * s%t_stem(p))
        s%eflx_sh_stem(p) = s%forc_rho(c) * cpair * wtstem &
             * ((s%wta0(p) + wtg0 + s%wtl0(p)) * s%t_stem(p) &
             - wtg0 * s%t_grnd(c) - s%wta0(p) * s%thm(p) - s%wtl0(p) * s%t_veg(p))
        s%efe(p) = dc2 * (wtgaq * s%qsatl(p) - wtgq0 * s%qg(c) - s%wtaq0(p) * s%forc_q(c))

        erre = 0.0_r8
        if (s%efe(p) * s%efeb(p) < 0.0_r8) then
          efeold = s%efe(p)
          s%efe(p) = 0.1_r8 * efeold
          erre = s%efe(p) - efeold
        end if

        lw_grnd = s%frac_sno(c) * s%t_soisno_snlp1(c)**4 &
             + (1.0_r8 - s%frac_sno(c) - s%frac_h2osfc(c)) * s%t_soisno_1(c)**4 &
             + s%frac_h2osfc(c) * s%t_h2osfc(c)**4

        s%dt_veg(p) = ((1.0_r8 - s%frac_rad_abs_by_stem(p)) &
             * (s%sabv(p) + s%air(p) + s%bir(p) * s%t_veg(p)**4 + s%cir(p) * lw_grnd) &
             - efsh - s%efe(p) - s%lw_leaf(p) + s%lw_stem(p) &
             - (s%cp_leaf(p) / opts%dtime) * (s%t_veg(p) - s%tl_ini(p))) &
             / ((1.0_r8 - s%frac_rad_abs_by_stem(p)) * (-4.0_r8 * s%bir(p) * s%t_veg(p)**3) &
             + 4.0_r8 * s%sa_internal(p) * s%emv(p) * sb * s%t_veg(p)**3 &
             + dc1 * s%wtga(p) + dc2 * wtgaq * s%qsatldT(p) + s%cp_leaf(p) / opts%dtime)

        s%t_veg(p) = s%tlbef(p) + s%dt_veg(p)

        dels = s%dt_veg(p)
        s%del(p) = abs(dels)
        s%err(p) = 0.0_r8
        if (s%del(p) > delmax) then
          s%dt_veg(p) = delmax * dels / s%del(p)
          s%t_veg(p) = s%tlbef(p) + s%dt_veg(p)
          s%err(p) = (1.0_r8 - s%frac_rad_abs_by_stem(p)) &
               * (s%sabv(p) + s%air(p) + s%bir(p) * s%tlbef(p)**3 &
               * (s%tlbef(p) + 4.0_r8 * s%dt_veg(p)) + s%cir(p) * lw_grnd) &
               - s%sa_internal(p) * s%emv(p) * sb * s%tlbef(p)**3 &
               * (s%tlbef(p) + 4.0_r8 * s%dt_veg(p)) + s%lw_stem(p) &
               - (efsh + dc1 * s%wtga(p) * s%dt_veg(p)) &
               - (s%efe(p) + dc2 * wtgaq * s%qsatldT(p) * s%dt_veg(p)) &
               - (s%cp_leaf(p) / opts%dtime) * (s%t_veg(p) - s%tl_ini(p))
        end if

        efpot = s%forc_rho(c) * ((s%elai(p) + s%esai(p)) / s%rb(p)) &
             * (wtgaq * (s%qsatl(p) + s%qsatldT(p) * s%dt_veg(p)) &
             - wtgq0 * s%qg(c) - s%wtaq0(p) * s%forc_q(c))
        s%qflx_evap_veg(p) = rpp * efpot

        ecidif = 0.0_r8
        if (efpot > 0.0_r8) then
          s%qflx_tran_veg(p) = efpot * rppdry
        else
          s%qflx_tran_veg(p) = 0.0_r8
        end if
        ecidif = max(0.0_r8, s%qflx_evap_veg(p) - s%qflx_tran_veg(p) - h2ocan / opts%dtime)
        s%qflx_evap_veg(p) = min(s%qflx_evap_veg(p), s%qflx_tran_veg(p) + h2ocan / opts%dtime)

        s%eflx_sh_veg(p) = efsh + dc1 * s%wtga(p) * s%dt_veg(p) + s%err(p) + erre + hvap * ecidif

        s%eflx_sh_stem(p) = s%eflx_sh_stem(p) &
             + s%forc_rho(c) * cpair * wtstem * (-s%wtl0(p) * s%dt_veg(p))
        s%lw_leaf(p) = s%sa_internal(p) * s%emv(p) * sb * s%tlbef(p)**3 &
             * (s%tlbef(p) + 4.0_r8 * s%dt_veg(p))

        td = min(100.0_r8, max(-75.0_r8, s%t_veg(p) - tkfrz))

        es_local = 6.11213476_r8 + td * (0.444007856_r8 + td * (0.143064234e-1_r8 &
             + td * (0.264461437e-3_r8 + td * (0.305903558e-5_r8 + td * (0.196237241e-7_r8 &
             + td * (0.892344772e-10_r8 + td * (-0.373208410e-12_r8 &
             + td * 0.209339997e-15_r8)))))))
        esdT_local = 0.444017302_r8 + td * (0.286064092e-1_r8 + td * (0.794683137e-3_r8 &
             + td * (0.121211669e-4_r8 + td * (0.103354611e-6_r8 + td * (0.404125005e-9_r8 &
             + td * (-0.788037859e-12_r8 + td * (-0.114596802e-13_r8 &
             + td * 0.381294516e-16_r8)))))))

        s%el(p) = es_local * 100.0_r8
        esdT_local = esdT_local * 100.0_r8
        vp = 1.0_r8 / (s%forc_pbot(c) - 0.378_r8 * s%el(p))
        vp1 = 0.622_r8 * vp
        s%qsatl(p) = s%el(p) * vp1
        vp2 = vp1 * vp
        s%qsatldT(p) = esdT_local * vp2 * s%forc_pbot(c)

        s%taf(p) = wtg0 * s%t_grnd(c) + s%wta0(p) * s%thm(p) &
             + s%wtl0(p) * s%t_veg(p) + s%wtstem0(p) * s%t_stem(p)
        s%qaf(p) = s%wtlq0(p) * s%qsatl(p) + wtgq0 * s%qg(c) + s%forc_q(c) * s%wtaq0(p)

        s%dth(p) = s%thm(p) - s%taf(p)
        s%dqh(p) = s%forc_q(c) - s%qaf(p)
        s%delq(p) = s%wtalq(p) * s%qg(c) - s%wtlq0(p) * s%qsatl(p) - s%wtaq0(p) * s%forc_q(c)

        tstar = s%temp1(p) * s%dth(p)
        qstar = s%temp2(p) * s%dqh(p)

        thvstar = tstar * (1.0_r8 + 0.61_r8 * s%forc_q(c)) + 0.61_r8 * s%forc_th(c) * qstar
        s%zeta(p) = s%zldis(p) * vkc * grav * thvstar / (s%ustar(p)**2 * s%thv(c))

        if (s%zeta(p) >= 0.0_r8) then
          s%zeta(p) = min(params%zetamaxstable, max(s%zeta(p), 0.01_r8))
          s%um(p) = max(s%ur(p), 0.1_r8)
        else
          s%zeta(p) = max(-100.0_r8, min(s%zeta(p), -0.01_r8))
          if (s%ustar(p) * thvstar > 0.0_r8) then
            wc = 0.0_r8
          else
            wc = beta * (-grav * s%ustar(p) * thvstar * zii / s%thv(c))**0.333_r8
          end if
          s%um(p) = sqrt(s%ur(p) * s%ur(p) + wc * wc)
        end if
        s%obu(p) = s%zldis(p) / s%zeta(p)

        if (s%obuold(p) * s%obu(p) < 0.0_r8) s%nmozsgn(p) = s%nmozsgn(p) + 1
        if (s%nmozsgn(p) >= 4) s%obu(p) = s%zldis(p) / (-0.01_r8)
        s%obuold(p) = s%obu(p)
      end do

      itlef = itlef + 1
      if (itlef > itmin) then
        !DIR$ IVDEP
        !DIR$ VECTOR ALWAYS
        do f = 1, fn
          p = s%filterp(f)
          s%dele(p) = abs(s%efe(p) - s%efeb(p))
          s%efeb(p) = s%efe(p)
          s%det(p) = max(s%del(p), s%del2(p))
          s%num_iter(p) = real(itlef, r8)
        end do

        fnold = fn
        fn = 0
        do f = 1, fnold
          p = s%filterp(f)
          if (.not. (s%det(p) < dtmin .and. s%dele(p) < dlemin)) then
            fn = fn + 1
            s%filterp(fn) = p
          end if
        end do
      end if
    end do

    s%filterp(1:fnorig) = s%fporig(1:fnorig)
  end subroutine can_iter_cpu_default_inline

  subroutine FrictionVelocityGeneral(s, params, fn, filtern, iter)
    type(can_iter_state), intent(inout) :: s
    type(can_iter_params), intent(in) :: params
    integer, intent(in) :: fn, filtern(fn), iter

    real(r8), parameter :: zetam = 1.574_r8
    real(r8), parameter :: zetat = 0.465_r8
    integer :: f, n
    real(r8) :: zldis_local, zeta_local
    real(r8) :: tmp1, tmp2, tmp3, tmp4, fmnew, fm10, zeta10

    do f = 1, fn
      n = filtern(f)

      zldis_local = s%forc_hgt_u_patch(n) - s%displa(n)
      zeta_local = zldis_local / s%obu(n)
      if (zeta_local < -zetam) then
        s%ustar(n) = vkc * s%um(n) / (log(-zetam * s%obu(n) / s%z0mv(n)) &
             - StabilityFunc1(-zetam) + StabilityFunc1(s%z0mv(n) / s%obu(n)) &
             + 1.14_r8 * ((-zeta_local)**0.333_r8 - zetam**0.333_r8))
      else if (zeta_local < 0.0_r8) then
        s%ustar(n) = vkc * s%um(n) / (log(zldis_local / s%z0mv(n)) &
             - StabilityFunc1(zeta_local) + StabilityFunc1(s%z0mv(n) / s%obu(n)))
      else if (zeta_local <= 1.0_r8) then
        s%ustar(n) = vkc * s%um(n) / (log(zldis_local / s%z0mv(n)) &
             + 5.0_r8 * zeta_local - 5.0_r8 * s%z0mv(n) / s%obu(n))
      else
        s%ustar(n) = vkc * s%um(n) / (log(s%obu(n) / s%z0mv(n)) &
             + 5.0_r8 - 5.0_r8 * s%z0mv(n) / s%obu(n) &
             + (5.0_r8 * log(zeta_local) + zeta_local - 1.0_r8))
      end if

      zldis_local = s%forc_hgt_t_patch(n) - s%displa(n)
      zeta_local = zldis_local / s%obu(n)
      if (zeta_local < -zetat) then
        s%temp1(n) = vkc / (log(-zetat * s%obu(n) / s%z0hv(n)) &
             - StabilityFunc2(-zetat) + StabilityFunc2(s%z0hv(n) / s%obu(n)) &
             + 0.8_r8 * (zetat**(-0.333_r8) - (-zeta_local)**(-0.333_r8)))
      else if (zeta_local < 0.0_r8) then
        s%temp1(n) = vkc / (log(zldis_local / s%z0hv(n)) &
             - StabilityFunc2(zeta_local) + StabilityFunc2(s%z0hv(n) / s%obu(n)))
      else if (zeta_local <= 1.0_r8) then
        s%temp1(n) = vkc / (log(zldis_local / s%z0hv(n)) &
             + 5.0_r8 * zeta_local - 5.0_r8 * s%z0hv(n) / s%obu(n))
      else
        s%temp1(n) = vkc / (log(s%obu(n) / s%z0hv(n)) &
             + 5.0_r8 - 5.0_r8 * s%z0hv(n) / s%obu(n) &
             + (5.0_r8 * log(zeta_local) + zeta_local - 1.0_r8))
      end if

      if (s%forc_hgt_q_patch(n) == s%forc_hgt_t_patch(n) .and. s%z0qv(n) == s%z0hv(n)) then
        s%temp2(n) = s%temp1(n)
      else
        zldis_local = s%forc_hgt_q_patch(n) - s%displa(n)
        zeta_local = zldis_local / s%obu(n)
        if (zeta_local < -zetat) then
          s%temp2(n) = vkc / (log(-zetat * s%obu(n) / s%z0qv(n)) &
               - StabilityFunc2(-zetat) + StabilityFunc2(s%z0qv(n) / s%obu(n)) &
               + 0.8_r8 * (zetat**(-0.333_r8) - (-zeta_local)**(-0.333_r8)))
        else if (zeta_local < 0.0_r8) then
          s%temp2(n) = vkc / (log(zldis_local / s%z0qv(n)) &
               - StabilityFunc2(zeta_local) + StabilityFunc2(s%z0qv(n) / s%obu(n)))
        else if (zeta_local <= 1.0_r8) then
          s%temp2(n) = vkc / (log(zldis_local / s%z0qv(n)) &
               + 5.0_r8 * zeta_local - 5.0_r8 * s%z0qv(n) / s%obu(n))
        else
          s%temp2(n) = vkc / (log(s%obu(n) / s%z0qv(n)) &
               + 5.0_r8 - 5.0_r8 * s%z0qv(n) / s%obu(n) &
               + (5.0_r8 * log(zeta_local) + zeta_local - 1.0_r8))
        end if
      end if

      zldis_local = 2.0_r8 + s%z0hv(n)
      zeta_local = zldis_local / s%obu(n)
      if (zeta_local < -zetat) then
        s%temp12m(n) = vkc / (log(-zetat * s%obu(n) / s%z0hv(n)) &
             - StabilityFunc2(-zetat) + StabilityFunc2(s%z0hv(n) / s%obu(n)) &
             + 0.8_r8 * (zetat**(-0.333_r8) - (-zeta_local)**(-0.333_r8)))
      else if (zeta_local < 0.0_r8) then
        s%temp12m(n) = vkc / (log(zldis_local / s%z0hv(n)) &
             - StabilityFunc2(zeta_local) + StabilityFunc2(s%z0hv(n) / s%obu(n)))
      else if (zeta_local <= 1.0_r8) then
        s%temp12m(n) = vkc / (log(zldis_local / s%z0hv(n)) &
             + 5.0_r8 * zeta_local - 5.0_r8 * s%z0hv(n) / s%obu(n))
      else
        s%temp12m(n) = vkc / (log(s%obu(n) / s%z0hv(n)) &
             + 5.0_r8 - 5.0_r8 * s%z0hv(n) / s%obu(n) &
             + (5.0_r8 * log(zeta_local) + zeta_local - 1.0_r8))
      end if

      if (s%z0qv(n) == s%z0hv(n)) then
        s%temp22m(n) = s%temp12m(n)
      else
        zldis_local = 2.0_r8 + s%z0qv(n)
        zeta_local = zldis_local / s%obu(n)
        if (zeta_local < -zetat) then
          s%temp22m(n) = vkc / (log(-zetat * s%obu(n) / s%z0qv(n)) &
               - StabilityFunc2(-zetat) + StabilityFunc2(s%z0qv(n) / s%obu(n)) &
               + 0.8_r8 * (zetat**(-0.333_r8) - (-zeta_local)**(-0.333_r8)))
        else if (zeta_local < 0.0_r8) then
          s%temp22m(n) = vkc / (log(zldis_local / s%z0qv(n)) &
               - StabilityFunc2(zeta_local) + StabilityFunc2(s%z0qv(n) / s%obu(n)))
        else if (zeta_local <= 1.0_r8) then
          s%temp22m(n) = vkc / (log(zldis_local / s%z0qv(n)) &
               + 5.0_r8 * zeta_local - 5.0_r8 * s%z0qv(n) / s%obu(n))
        else
          s%temp22m(n) = vkc / (log(s%obu(n) / s%z0qv(n)) &
               + 5.0_r8 - 5.0_r8 * s%z0qv(n) / s%obu(n) &
               + (5.0_r8 * log(zeta_local) + zeta_local - 1.0_r8))
        end if
      end if

      zldis_local = s%forc_hgt_u_patch(n) - s%displa(n)
      zeta_local = zldis_local / s%obu(n)
      if (min(zeta_local, 1.0_r8) < 0.0_r8) then
        tmp1 = (1.0_r8 - 16.0_r8 * min(zeta_local, 1.0_r8))**0.25_r8
        tmp2 = log((1.0_r8 + tmp1 * tmp1) / 2.0_r8)
        tmp3 = log((1.0_r8 + tmp1) / 2.0_r8)
        fmnew = 2.0_r8 * tmp3 + tmp2 - 2.0_r8 * atan(tmp1) + pi * 0.5_r8
      else
        fmnew = -5.0_r8 * min(zeta_local, 1.0_r8)
      end if
      if (iter == 1) then
        s%fm(n) = fmnew
      else
        s%fm(n) = 0.5_r8 * (s%fm(n) + fmnew)
      end if

      zeta10 = min(10.0_r8 / s%obu(n), 1.0_r8)
      if (zeta_local == 0.0_r8) zeta10 = 0.0_r8
      if (zeta10 < 0.0_r8) then
        tmp1 = (1.0_r8 - 16.0_r8 * zeta10)**0.25_r8
        tmp2 = log((1.0_r8 + tmp1 * tmp1) / 2.0_r8)
        tmp3 = log((1.0_r8 + tmp1) / 2.0_r8)
        fm10 = 2.0_r8 * tmp3 + tmp2 - 2.0_r8 * atan(tmp1) + pi * 0.5_r8
      else
        fm10 = -5.0_r8 * zeta10
      end if
      tmp4 = log(max(1.0_r8, s%forc_hgt_u_patch(n) / 10.0_r8))
      s%fm(n) = s%ur(n) - s%ustar(n) / vkc * (tmp4 - s%fm(n) + fm10)
    end do
  end subroutine FrictionVelocityGeneral

  subroutine FrictionVelocityDefault(s, params, fn, filtern, iter)
    type(can_iter_state), intent(inout) :: s
    type(can_iter_params), intent(in) :: params
    integer, intent(in) :: fn, filtern(fn), iter

    real(r8), parameter :: zetam = 1.574_r8
    real(r8), parameter :: zetat = 0.465_r8
    integer :: f, n
    real(r8) :: zldis_local, zeta_local
    real(r8) :: tmp1, tmp2, tmp3, tmp4, fmnew, fm10, zeta10

    !DIR$ IVDEP
    !DIR$ VECTOR ALWAYS
    do f = 1, fn
      n = filtern(f)

      zldis_local = s%forc_hgt_u_patch(n) - s%displa(n)
      zeta_local = zldis_local / s%obu(n)
      if (zeta_local < -zetam) then
        s%ustar(n) = vkc * s%um(n) / (log(-zetam * s%obu(n) / s%z0mv(n)) &
             - StabilityFunc1(-zetam) + StabilityFunc1(s%z0mv(n) / s%obu(n)) &
             + 1.14_r8 * ((-zeta_local)**0.333_r8 - zetam**0.333_r8))
      else if (zeta_local < 0.0_r8) then
        s%ustar(n) = vkc * s%um(n) / (log(zldis_local / s%z0mv(n)) &
             - StabilityFunc1(zeta_local) + StabilityFunc1(s%z0mv(n) / s%obu(n)))
      else if (zeta_local <= 1.0_r8) then
        s%ustar(n) = vkc * s%um(n) / (log(zldis_local / s%z0mv(n)) &
             + 5.0_r8 * zeta_local - 5.0_r8 * s%z0mv(n) / s%obu(n))
      else
        s%ustar(n) = vkc * s%um(n) / (log(s%obu(n) / s%z0mv(n)) &
             + 5.0_r8 - 5.0_r8 * s%z0mv(n) / s%obu(n) &
             + (5.0_r8 * log(zeta_local) + zeta_local - 1.0_r8))
      end if

      zldis_local = s%forc_hgt_t_patch(n) - s%displa(n)
      zeta_local = zldis_local / s%obu(n)
      if (zeta_local < -zetat) then
        s%temp1(n) = vkc / (log(-zetat * s%obu(n) / s%z0hv(n)) &
             - StabilityFunc2(-zetat) + StabilityFunc2(s%z0hv(n) / s%obu(n)) &
             + 0.8_r8 * (zetat**(-0.333_r8) - (-zeta_local)**(-0.333_r8)))
      else if (zeta_local < 0.0_r8) then
        s%temp1(n) = vkc / (log(zldis_local / s%z0hv(n)) &
             - StabilityFunc2(zeta_local) + StabilityFunc2(s%z0hv(n) / s%obu(n)))
      else if (zeta_local <= 1.0_r8) then
        s%temp1(n) = vkc / (log(zldis_local / s%z0hv(n)) &
             + 5.0_r8 * zeta_local - 5.0_r8 * s%z0hv(n) / s%obu(n))
      else
        s%temp1(n) = vkc / (log(s%obu(n) / s%z0hv(n)) &
             + 5.0_r8 - 5.0_r8 * s%z0hv(n) / s%obu(n) &
             + (5.0_r8 * log(zeta_local) + zeta_local - 1.0_r8))
      end if

      s%temp2(n) = s%temp1(n)

      zldis_local = 2.0_r8 + s%z0hv(n)
      zeta_local = zldis_local / s%obu(n)
      if (zeta_local < -zetat) then
        s%temp12m(n) = vkc / (log(-zetat * s%obu(n) / s%z0hv(n)) &
             - StabilityFunc2(-zetat) + StabilityFunc2(s%z0hv(n) / s%obu(n)) &
             + 0.8_r8 * (zetat**(-0.333_r8) - (-zeta_local)**(-0.333_r8)))
      else if (zeta_local < 0.0_r8) then
        s%temp12m(n) = vkc / (log(zldis_local / s%z0hv(n)) &
             - StabilityFunc2(zeta_local) + StabilityFunc2(s%z0hv(n) / s%obu(n)))
      else if (zeta_local <= 1.0_r8) then
        s%temp12m(n) = vkc / (log(zldis_local / s%z0hv(n)) &
             + 5.0_r8 * zeta_local - 5.0_r8 * s%z0hv(n) / s%obu(n))
      else
        s%temp12m(n) = vkc / (log(s%obu(n) / s%z0hv(n)) &
             + 5.0_r8 - 5.0_r8 * s%z0hv(n) / s%obu(n) &
             + (5.0_r8 * log(zeta_local) + zeta_local - 1.0_r8))
      end if

      s%temp22m(n) = s%temp12m(n)

      zldis_local = s%forc_hgt_u_patch(n) - s%displa(n)
      zeta_local = zldis_local / s%obu(n)
      if (min(zeta_local, 1.0_r8) < 0.0_r8) then
        tmp1 = (1.0_r8 - 16.0_r8 * min(zeta_local, 1.0_r8))**0.25_r8
        tmp2 = log((1.0_r8 + tmp1 * tmp1) / 2.0_r8)
        tmp3 = log((1.0_r8 + tmp1) / 2.0_r8)
        fmnew = 2.0_r8 * tmp3 + tmp2 - 2.0_r8 * atan(tmp1) + pi * 0.5_r8
      else
        fmnew = -5.0_r8 * min(zeta_local, 1.0_r8)
      end if
      if (iter == 1) then
        s%fm(n) = fmnew
      else
        s%fm(n) = 0.5_r8 * (s%fm(n) + fmnew)
      end if

      zeta10 = min(10.0_r8 / s%obu(n), 1.0_r8)
      if (zeta_local == 0.0_r8) zeta10 = 0.0_r8
      if (zeta10 < 0.0_r8) then
        tmp1 = (1.0_r8 - 16.0_r8 * zeta10)**0.25_r8
        tmp2 = log((1.0_r8 + tmp1 * tmp1) / 2.0_r8)
        tmp3 = log((1.0_r8 + tmp1) / 2.0_r8)
        fm10 = 2.0_r8 * tmp3 + tmp2 - 2.0_r8 * atan(tmp1) + pi * 0.5_r8
      else
        fm10 = -5.0_r8 * zeta10
      end if
      tmp4 = log(max(1.0_r8, s%forc_hgt_u_patch(n) / 10.0_r8))
      s%fm(n) = s%ur(n) - s%ustar(n) / vkc * (tmp4 - s%fm(n) + fm10)
    end do
  end subroutine FrictionVelocityDefault

  pure subroutine MoninObukIni(params, ur, thv, dthv, zldis, z0m, um, obu)
    type(can_iter_params), intent(in) :: params
    real(r8), intent(in) :: ur, thv, dthv, zldis, z0m
    real(r8), intent(out) :: um, obu
    real(r8) :: wc, rib, zeta

    wc = 0.5_r8
    if (dthv >= 0.0_r8) then
      um = max(ur, 0.1_r8)
    else
      um = sqrt(ur * ur + wc * wc)
    end if

    rib = grav * zldis * dthv / (thv * um * um)

    if (rib >= 0.0_r8) then
      zeta = rib * log(zldis / z0m) / (1.0_r8 - 5.0_r8 * min(rib, 0.19_r8))
      zeta = min(params%zetamaxstable, max(zeta, 0.01_r8))
    else
      zeta = rib * log(zldis / z0m)
      zeta = max(-100.0_r8, min(zeta, -0.01_r8))
    end if

    obu = zldis / zeta
  end subroutine MoninObukIni

  elemental real(r8) function StabilityFunc1(zeta)
    real(r8), intent(in) :: zeta
    real(r8) :: chik, chik2

    chik2 = sqrt(1.0_r8 - 16.0_r8 * zeta)
    chik = sqrt(chik2)
    StabilityFunc1 = 2.0_r8 * log((1.0_r8 + chik) * 0.5_r8) &
         + log((1.0_r8 + chik2) * 0.5_r8) - 2.0_r8 * atan(chik) + pi * 0.5_r8
  end function StabilityFunc1

  elemental real(r8) function StabilityFunc2(zeta)
    real(r8), intent(in) :: zeta
    real(r8) :: chik2

    chik2 = sqrt(1.0_r8 - 16.0_r8 * zeta)
    StabilityFunc2 = 2.0_r8 * log((1.0_r8 + chik2) * 0.5_r8)
  end function StabilityFunc2

  subroutine PhotosynthesisStub(s, fn, filterp, phase)
    type(can_iter_state), intent(inout) :: s
    integer, intent(in) :: fn, filterp(fn)
    character(len=*), intent(in) :: phase
    integer :: f, p
    real(r8) :: vpd_kpa, rh_can, light, nitrogen, temp_stress, co2_fac, conductance

    do f = 1, fn
      p = filterp(f)
      vpd_kpa = max((s%svpts(p) - s%eah(p)) * 0.001_r8, 0.05_r8)
      rh_can = min(1.0_r8, max(0.05_r8, s%eah(p) / max(s%svpts(p), 1.0_r8)))
      light = max(0.05_r8, s%dayl_factor(p))
      nitrogen = max(0.25_r8, min(2.0_r8, s%leafn_patch(p) / 2.0_r8))
      temp_stress = exp(-((s%t_veg(p) - 298.0_r8) / 18.0_r8)**2)
      co2_fac = max(0.5_r8, min(1.5_r8, s%co2(p) / 40.0_r8))

      if (phase == 'sun') then
        conductance = 0.0030_r8 + 0.0180_r8 * light * nitrogen * temp_stress &
             * co2_fac * max(s%btran(p), 0.05_r8) * rh_can / (1.0_r8 + 0.10_r8 * vpd_kpa)
        s%rssun(p) = 1.0_r8 / max(conductance, 1.0e-4_r8)
      else
        conductance = 0.0020_r8 + 0.0100_r8 * light * nitrogen * temp_stress &
             * co2_fac * max(s%btran(p), 0.05_r8) * rh_can / (1.0_r8 + 0.14_r8 * vpd_kpa)
        s%rssha(p) = 1.0_r8 / max(conductance, 1.0e-4_r8)
      end if
    end do
  end subroutine PhotosynthesisStub

  subroutine PhotosynthesisHydraulicStressStub(s, fn, filterp)
    type(can_iter_state), intent(inout) :: s
    integer, intent(in) :: fn, filterp(fn)
    integer :: f, p
    real(r8) :: moisture

    do f = 1, fn
      p = filterp(f)
      moisture = max(0.05_r8, min(1.0_r8, s%btran(p)))
      s%bsun(p) = min(1.0_r8, moisture * (0.85_r8 + 0.10_r8 * s%dayl_factor(p)))
      s%bsha(p) = min(1.0_r8, moisture * (0.95_r8 + 0.05_r8 * s%dayl_factor(p)))
      s%btran(p) = 0.55_r8 * s%bsun(p) + 0.45_r8 * s%bsha(p)
    end do

    call PhotosynthesisStub(s, fn, filterp, 'sun')
    call PhotosynthesisStub(s, fn, filterp, 'sha')
  end subroutine PhotosynthesisHydraulicStressStub

  subroutine FractionationStub(s, fn, filterp, phase)
    type(can_iter_state), intent(inout) :: s
    integer, intent(in) :: fn, filterp(fn)
    character(len=*), intent(in) :: phase
    integer :: f, p
    real(r8) :: tiny_adjustment

    tiny_adjustment = 0.0_r8
    if (phase == 'sun') tiny_adjustment = 1.0e-12_r8
    do f = 1, fn
      p = filterp(f)
      s%canopy_cond(p) = s%canopy_cond(p) + tiny_adjustment * real(mod(p, 7), r8)
    end do
  end subroutine FractionationStub

  pure subroutine QSat(T, p, qs, es, qsdT)
    real(r8), intent(in) :: T, p
    real(r8), intent(out) :: qs
    real(r8), intent(out), optional :: es, qsdT

    real(r8), parameter :: a0 = 6.11213476_r8
    real(r8), parameter :: a1 = 0.444007856_r8
    real(r8), parameter :: a2 = 0.143064234e-1_r8
    real(r8), parameter :: a3 = 0.264461437e-3_r8
    real(r8), parameter :: a4 = 0.305903558e-5_r8
    real(r8), parameter :: a5 = 0.196237241e-7_r8
    real(r8), parameter :: a6 = 0.892344772e-10_r8
    real(r8), parameter :: a7 = -0.373208410e-12_r8
    real(r8), parameter :: a8 = 0.209339997e-15_r8
    real(r8), parameter :: b0 = 0.444017302_r8
    real(r8), parameter :: b1 = 0.286064092e-1_r8
    real(r8), parameter :: b2 = 0.794683137e-3_r8
    real(r8), parameter :: b3 = 0.121211669e-4_r8
    real(r8), parameter :: b4 = 0.103354611e-6_r8
    real(r8), parameter :: b5 = 0.404125005e-9_r8
    real(r8), parameter :: b6 = -0.788037859e-12_r8
    real(r8), parameter :: b7 = -0.114596802e-13_r8
    real(r8), parameter :: b8 = 0.381294516e-16_r8
    real(r8), parameter :: c0 = 6.11123516_r8
    real(r8), parameter :: c1 = 0.503109514_r8
    real(r8), parameter :: c2 = 0.188369801e-1_r8
    real(r8), parameter :: c3 = 0.420547422e-3_r8
    real(r8), parameter :: c4 = 0.614396778e-5_r8
    real(r8), parameter :: c5 = 0.602780717e-7_r8
    real(r8), parameter :: c6 = 0.387940929e-9_r8
    real(r8), parameter :: c7 = 0.149436277e-11_r8
    real(r8), parameter :: c8 = 0.262655803e-14_r8
    real(r8), parameter :: d0 = 0.503277922_r8
    real(r8), parameter :: d1 = 0.377289173e-1_r8
    real(r8), parameter :: d2 = 0.126801703e-2_r8
    real(r8), parameter :: d3 = 0.249468427e-4_r8
    real(r8), parameter :: d4 = 0.313703411e-6_r8
    real(r8), parameter :: d5 = 0.257180651e-8_r8
    real(r8), parameter :: d6 = 0.133268878e-10_r8
    real(r8), parameter :: d7 = 0.394116744e-13_r8
    real(r8), parameter :: d8 = 0.498070196e-16_r8

    real(r8) :: td, es_local, esdT_local, vp, vp1, vp2

    td = min(100.0_r8, max(-75.0_r8, T - tkfrz))

    if (td >= 0.0_r8) then
      es_local = a0 + td * (a1 + td * (a2 + td * (a3 + td * (a4 &
           + td * (a5 + td * (a6 + td * (a7 + td * a8)))))))
    else
      es_local = c0 + td * (c1 + td * (c2 + td * (c3 + td * (c4 &
           + td * (c5 + td * (c6 + td * (c7 + td * c8)))))))
    end if

    es_local = es_local * 100.0_r8
    vp = 1.0_r8 / (p - 0.378_r8 * es_local)
    vp1 = 0.622_r8 * vp
    qs = es_local * vp1
    if (present(es)) es = es_local

    if (present(qsdT)) then
      if (td >= 0.0_r8) then
        esdT_local = b0 + td * (b1 + td * (b2 + td * (b3 + td * (b4 &
             + td * (b5 + td * (b6 + td * (b7 + td * b8)))))))
      else
        esdT_local = d0 + td * (d1 + td * (d2 + td * (d3 + td * (d4 &
             + td * (d5 + td * (d6 + td * (d7 + td * d8)))))))
      end if
      esdT_local = esdT_local * 100.0_r8
      vp2 = vp1 * vp
      qsdT = esdT_local * vp2 * p
    end if
  end subroutine QSat

  real(r8) function compute_checksum(s) result(sumval)
    type(can_iter_state), intent(in) :: s
    integer :: p
    real(r8) :: weight

    sumval = 0.0_r8
    do p = 1, s%npts
      weight = 1.0_r8 + real(mod(p, 17), r8) * 0.001_r8
      sumval = sumval + weight * ( &
           0.50_r8 * s%t_veg(p) &
           + 0.03_r8 * s%t_stem(p) &
           + 0.20_r8 * s%obu(p) &
           + 0.10_r8 * s%eflx_sh_veg(p) &
           + 0.02_r8 * s%eflx_sh_stem(p) &
           + 1.0e5_r8 * s%qflx_evap_veg(p) &
           + 1.0e5_r8 * s%qflx_tran_veg(p) &
           + 10.0_r8 * s%canopy_cond(p) &
           + 10.0_r8 * s%grnd_ch4_cond(p) &
           + s%num_iter(p))
    end do
  end function compute_checksum

  subroutine print_summary(s)
    type(can_iter_state), intent(in) :: s
    real(r8) :: invn

    invn = 1.0_r8 / real(s%npts, r8)
    write(*,'(a,es20.12)') '  mean t_veg                  = ', sum(s%t_veg) * invn
    write(*,'(a,es20.12)') '  mean qflx_evap_veg          = ', sum(s%qflx_evap_veg) * invn
    write(*,'(a,es20.12)') '  mean eflx_sh_veg            = ', sum(s%eflx_sh_veg) * invn
    write(*,'(a,es20.12)') '  mean num_iter               = ', sum(s%num_iter) * invn
  end subroutine print_summary

end module can_iter_cpu_bucketed_reproducer_mod

program can_iter_cpu_bucketed_reproducer
  use, intrinsic :: iso_fortran_env, only : int64
  use can_iter_cpu_bucketed_reproducer_mod
  implicit none

  type(can_iter_state) :: state
  type(can_iter_params) :: params
  type(can_iter_options) :: opts
  integer :: npts, nrepeat, nwarmup, itmax_canopy_fluxes
  integer :: i
  integer(int64) :: count_rate, count_start, count_end
  real(r8) :: elapsed, avg_time, checksum, patch_iterations

  npts = 65536
  nrepeat = 20
  nwarmup = 3
  itmax_canopy_fluxes = 40

  call get_arg_int(1, npts)
  call get_arg_int(2, nrepeat)
  call get_arg_int(3, nwarmup)
  call get_arg_int(4, itmax_canopy_fluxes)

  if (npts < 1) stop 'npts must be positive'
  if (nrepeat < 1) stop 'nrepeat must be positive'
  if (nwarmup < 0) stop 'nwarmup must be nonnegative'
  if (itmax_canopy_fluxes < 0) stop 'itmax must be nonnegative'

  call allocate_state(state, npts)

  do i = 1, nwarmup
    call reset_state(state)
    call setup_can_iter(state, params, opts)
    call can_iter_cpu_bucketed(state, params, opts, itmax_canopy_fluxes)
  end do

  elapsed = 0.0_r8
  call system_clock(count_rate=count_rate)

  do i = 1, nrepeat
    call reset_state(state)
    call setup_can_iter(state, params, opts)
    call system_clock(count_start)
    call can_iter_cpu_bucketed(state, params, opts, itmax_canopy_fluxes)
    call system_clock(count_end)
    elapsed = elapsed + real(count_end - count_start, r8) / real(count_rate, r8)
  end do

  checksum = compute_checksum(state)
  avg_time = elapsed / real(nrepeat, r8)
  patch_iterations = real(nrepeat, r8) * sum(state%num_iter)

  write(*,'(a)') 'can_iter bucketed CPU reproducer'
  write(*,'(a,i0)') '  npts                         = ', npts
  write(*,'(a,i0)') '  nrepeat                      = ', nrepeat
  write(*,'(a,i0)') '  nwarmup                      = ', nwarmup
  write(*,'(a,i0)') '  itmax_canopy_fluxes          = ', itmax_canopy_fluxes
  write(*,'(a,es20.12)') '  timed can_iter seconds       = ', elapsed
  write(*,'(a,es20.12)') '  seconds per can_iter call    = ', avg_time
  write(*,'(a,es20.12)') '  patch iterations per second  = ', patch_iterations / max(elapsed, tiny(elapsed))
  write(*,'(a,es20.12)') '  checksum                     = ', checksum
  call print_summary(state)

contains

  subroutine get_arg_int(index, value)
    integer, intent(in) :: index
    integer, intent(inout) :: value
    character(len=64) :: arg
    integer :: status

    call get_command_argument(index, arg)
    if (len_trim(arg) == 0) return

    read(arg, *, iostat=status) value
    if (status /= 0) then
      write(*,'(a,i0,a,a)') 'Could not parse argument ', index, ': ', trim(arg)
      stop 2
    end if
  end subroutine get_arg_int

end program can_iter_cpu_bucketed_reproducer

