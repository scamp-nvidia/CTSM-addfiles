program driver_SoilBiogeochemCompetition

   use SoilBiogeochemCompetition_mod, only : SoilBiogeochemCompetition, r8

   implicit none

   integer :: begc, endc
   integer :: ncol, nlevdecomp, ndecomp_cascade_transitions
   integer :: num_bgc_soilc
   integer :: nrepeat, nwarmup
   integer :: i, j, k, c
   integer :: count_rate, count_start, count_end
   real(r8) :: elapsed, avg_time, calls_per_sec
   real(r8) :: checksum

   real(r8) :: dt, bdnr
   logical :: use_nitrif_denitrif, carbon_only
   integer :: decomp_method, mimics_decomp
   integer :: i_cop_mic, i_oli_mic
   real(r8) :: compet_plant_no3, compet_plant_nh4
   real(r8) :: compet_decomp_no3, compet_decomp_nh4
   real(r8) :: compet_denit, compet_nit

   integer, allocatable, target :: filter_bgc_soilc(:)
   real(r8), allocatable        :: dzsoi_decomp(:)

   integer, pointer :: cascade_receiver_pool(:)
integer, pointer :: landunit(:)

real(r8), pointer :: fpg(:)
real(r8), pointer :: fpi(:)
real(r8), pointer :: fpi_vr(:,:)
real(r8), pointer :: nfixation_prof(:,:)
real(r8), pointer :: plant_ndemand(:)

real(r8), pointer :: sminn_vr(:,:)
real(r8), pointer :: smin_nh4_vr(:,:)
real(r8), pointer :: smin_no3_vr(:,:)

real(r8), pointer :: c_overflow_vr(:,:,:)

real(r8), pointer :: pot_f_nit_vr(:,:)
real(r8), pointer :: pot_f_denit_vr(:,:)
real(r8), pointer :: f_nit_vr(:,:)
real(r8), pointer :: f_denit_vr(:,:)

real(r8), pointer :: potential_immob(:)
real(r8), pointer :: actual_immob(:)
real(r8), pointer :: sminn_to_plant(:)
real(r8), pointer :: sminn_to_denit_excess_vr(:,:)

real(r8), pointer :: actual_immob_no3_vr(:,:)
real(r8), pointer :: actual_immob_nh4_vr(:,:)
real(r8), pointer :: smin_no3_to_plant_vr(:,:)
real(r8), pointer :: smin_nh4_to_plant_vr(:,:)

real(r8), pointer :: n2_n2o_ratio_denit_vr(:,:)
real(r8), pointer :: f_n2o_denit_vr(:,:)
real(r8), pointer :: f_n2o_nit_vr(:,:)

real(r8), pointer :: supplement_to_sminn_vr(:,:)
real(r8), pointer :: sminn_to_plant_vr(:,:)
real(r8), pointer :: potential_immob_vr(:,:)
real(r8), pointer :: actual_immob_vr(:,:)

   
   real(r8), allocatable :: pmnf_decomp_cascade(:,:,:)
   real(r8), allocatable :: p_decomp_cn_gain(:,:,:)

   character(len=64) :: arg

   ! --------------------------------------------------------------------
   ! Defaults
   ! --------------------------------------------------------------------
   ncol                         = 50000
   nlevdecomp                   = 20
   ndecomp_cascade_transitions  = 8
   nrepeat                      = 100
   nwarmup                      = 10
   use_nitrif_denitrif          = .false.
   carbon_only                  = .false.
   mimics_decomp                = 2
   decomp_method                = 1

   ! Args:
   !   1 ncol
   !   2 nlev
   !   3 ntrans
   !   4 nrepeat
   !   5 nwarmup
   !   6 use_nitrif        0/1
   !   7 carbon_only       0/1
   !   8 use_mimics        0/1
   call get_arg_int(1, ncol)
   call get_arg_int(2, nlevdecomp)
   call get_arg_int(3, ndecomp_cascade_transitions)
   call get_arg_int(4, nrepeat)
   call get_arg_int(5, nwarmup)
   call get_arg_logical(6, use_nitrif_denitrif)
   call get_arg_logical(7, carbon_only)

   call get_command_argument(8, arg)
   if (len_trim(arg) > 0) then
      read(arg, *) i
      if (i /= 0) then
         decomp_method = mimics_decomp
      else
         decomp_method = 1
      end if
   end if

   begc = 1
   endc = ncol
   num_bgc_soilc = ncol

   dt   = 1800.0_r8
   bdnr = 1.0e-6_r8

   i_cop_mic = 3
   i_oli_mic = 4

   compet_plant_no3  = 1.0_r8
   compet_plant_nh4  = 1.0_r8
   compet_decomp_no3 = 1.0_r8
   compet_decomp_nh4 = 1.0_r8
   compet_denit      = 1.0_r8
   compet_nit        = 1.0_r8

   ! --------------------------------------------------------------------
   ! Allocate
   ! --------------------------------------------------------------------
   allocate(filter_bgc_soilc(num_bgc_soilc))
   allocate(dzsoi_decomp(nlevdecomp))

   allocate(cascade_receiver_pool(ndecomp_cascade_transitions))
   allocate(landunit(begc:endc))

   allocate(fpg(begc:endc))
   allocate(fpi(begc:endc))
   allocate(fpi_vr(begc:endc,nlevdecomp))
   allocate(nfixation_prof(begc:endc,nlevdecomp))
   allocate(plant_ndemand(begc:endc))

   allocate(sminn_vr(begc:endc,nlevdecomp))
   allocate(smin_nh4_vr(begc:endc,nlevdecomp))
   allocate(smin_no3_vr(begc:endc,nlevdecomp))

   allocate(c_overflow_vr(begc:endc,nlevdecomp,ndecomp_cascade_transitions))

   allocate(pot_f_nit_vr(begc:endc,nlevdecomp))
   allocate(pot_f_denit_vr(begc:endc,nlevdecomp))
   allocate(f_nit_vr(begc:endc,nlevdecomp))
   allocate(f_denit_vr(begc:endc,nlevdecomp))

   allocate(potential_immob(begc:endc))
   allocate(actual_immob(begc:endc))
   allocate(sminn_to_plant(begc:endc))
   allocate(sminn_to_denit_excess_vr(begc:endc,nlevdecomp))

   allocate(actual_immob_no3_vr(begc:endc,nlevdecomp))
   allocate(actual_immob_nh4_vr(begc:endc,nlevdecomp))
   allocate(smin_no3_to_plant_vr(begc:endc,nlevdecomp))
   allocate(smin_nh4_to_plant_vr(begc:endc,nlevdecomp))

   allocate(n2_n2o_ratio_denit_vr(begc:endc,nlevdecomp))
   allocate(f_n2o_denit_vr(begc:endc,nlevdecomp))
   allocate(f_n2o_nit_vr(begc:endc,nlevdecomp))

   allocate(supplement_to_sminn_vr(begc:endc,nlevdecomp))
   allocate(sminn_to_plant_vr(begc:endc,nlevdecomp))
   allocate(potential_immob_vr(begc:endc,nlevdecomp))
   allocate(actual_immob_vr(begc:endc,nlevdecomp))

   allocate(pmnf_decomp_cascade(begc:endc,nlevdecomp,ndecomp_cascade_transitions))
   allocate(p_decomp_cn_gain(begc:endc,nlevdecomp,ndecomp_cascade_transitions))

   call initialize_inputs()

   ! --------------------------------------------------------------------
   ! Warmup
   ! --------------------------------------------------------------------
!$acc data &
!$acc copy( &
!$acc   filter_bgc_soilc(:num_bgc_soilc), &
!$acc   dzsoi_decomp(:nlevdecomp), &
!$acc   cascade_receiver_pool(:ndecomp_cascade_transitions), &
!$acc   landunit(begc:endc), &
!$acc   nfixation_prof(begc:endc,1:nlevdecomp), &
!$acc   plant_ndemand(begc:endc), &
!$acc   sminn_vr(begc:endc,1:nlevdecomp), &
!$acc   smin_nh4_vr(begc:endc,1:nlevdecomp), &
!$acc   smin_no3_vr(begc:endc,1:nlevdecomp), &
!$acc   pot_f_nit_vr(begc:endc,1:nlevdecomp), &
!$acc   pot_f_denit_vr(begc:endc,1:nlevdecomp), &
!$acc   n2_n2o_ratio_denit_vr(begc:endc,1:nlevdecomp), &
!$acc   potential_immob_vr(begc:endc,1:nlevdecomp), &
!$acc   pmnf_decomp_cascade(begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions), &
!$acc   p_decomp_cn_gain(begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions) ) &
!$acc copy( &
!$acc   fpg(begc:endc), &
!$acc   fpi(begc:endc), &
!$acc   fpi_vr(begc:endc,1:nlevdecomp), &
!$acc   c_overflow_vr(begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions), &
!$acc   f_nit_vr(begc:endc,1:nlevdecomp), &
!$acc   f_denit_vr(begc:endc,1:nlevdecomp), &
!$acc   potential_immob(begc:endc), &
!$acc   actual_immob(begc:endc), &
!$acc   sminn_to_plant(begc:endc), &
!$acc   sminn_to_denit_excess_vr(begc:endc,1:nlevdecomp), &
!$acc   actual_immob_no3_vr(begc:endc,1:nlevdecomp), &
!$acc   actual_immob_nh4_vr(begc:endc,1:nlevdecomp), &
!$acc   smin_no3_to_plant_vr(begc:endc,1:nlevdecomp), &
!$acc   smin_nh4_to_plant_vr(begc:endc,1:nlevdecomp), &
!$acc   f_n2o_denit_vr(begc:endc,1:nlevdecomp), &
!$acc   f_n2o_nit_vr(begc:endc,1:nlevdecomp), &
!$acc   supplement_to_sminn_vr(begc:endc,1:nlevdecomp), &
!$acc   sminn_to_plant_vr(begc:endc,1:nlevdecomp), &
!$acc   actual_immob_vr(begc:endc,1:nlevdecomp) )
   do i = 1, nwarmup
      call reset_outputs()
      call run_kernel()
   end do

   ! --------------------------------------------------------------------
   ! Timed loop
   ! --------------------------------------------------------------------
   call system_clock(count_rate=count_rate)
   call system_clock(count_start)
  
   do i = 1, nrepeat
      call reset_outputs()
      call run_kernel()
   end do

   call system_clock(count_end)
!$ACC END DATA

   elapsed = real(count_end - count_start, r8) / real(count_rate, r8)
   avg_time = elapsed / real(nrepeat, r8)
   calls_per_sec = real(nrepeat, r8) / elapsed

   checksum = compute_checksum()

   write(*,'(a)')       'SoilBiogeochemCompetition benchmark'
   write(*,'(a,i0)')    '  ncol                         = ', ncol
   write(*,'(a,i0)')    '  nlevdecomp                   = ', nlevdecomp
   write(*,'(a,i0)')    '  ndecomp_cascade_transitions  = ', ndecomp_cascade_transitions
   write(*,'(a,i0)')    '  nrepeat                      = ', nrepeat
   write(*,'(a,i0)')    '  nwarmup                      = ', nwarmup
   write(*,'(a,l1)')    '  use_nitrif_denitrif          = ', use_nitrif_denitrif
   write(*,'(a,l1)')    '  carbon_only                  = ', carbon_only
   write(*,'(a,l1)')    '  mimics branch                = ', decomp_method == mimics_decomp
   write(*,'(a,es16.8)') '  elapsed seconds              = ', elapsed
   write(*,'(a,es16.8)') '  avg seconds / call           = ', avg_time
   write(*,'(a,es16.8)') '  calls / second               = ', calls_per_sec
   write(*,'(a,es16.8)') '  checksum                     = ', checksum

contains

   subroutine run_kernel()

      call SoilBiogeochemCompetition( &
           begc, endc, nlevdecomp, ndecomp_cascade_transitions, &
           num_bgc_soilc, filter_bgc_soilc,                     &
           dt, bdnr,                                            &
           use_nitrif_denitrif, carbon_only,                    &
           decomp_method, mimics_decomp, i_cop_mic, i_oli_mic,  &
           compet_plant_no3, compet_plant_nh4,                  &
           compet_decomp_no3, compet_decomp_nh4,                &
           compet_denit, compet_nit,                            &
           dzsoi_decomp, cascade_receiver_pool, landunit,       &
           fpg, fpi, fpi_vr, nfixation_prof, plant_ndemand,     &
           sminn_vr, smin_nh4_vr, smin_no3_vr,                  &
           c_overflow_vr,                                       &
           pot_f_nit_vr, pot_f_denit_vr, f_nit_vr, f_denit_vr,  &
           potential_immob, actual_immob, sminn_to_plant,       &
           sminn_to_denit_excess_vr,                            &
           actual_immob_no3_vr, actual_immob_nh4_vr,            &
           smin_no3_to_plant_vr, smin_nh4_to_plant_vr,          &
           n2_n2o_ratio_denit_vr, f_n2o_denit_vr, f_n2o_nit_vr, &
           supplement_to_sminn_vr, sminn_to_plant_vr,           &
           potential_immob_vr, actual_immob_vr,                 &
           pmnf_decomp_cascade, p_decomp_cn_gain)

   end subroutine run_kernel

   subroutine initialize_inputs()

      do c = begc, endc
         filter_bgc_soilc(c) = c
         landunit(c) = 1

         plant_ndemand(c) = 1.0e-7_r8 * (1.0_r8 + real(mod(c,17),r8) / 17.0_r8)

         fpg(c) = 0.0_r8
         fpi(c) = 0.0_r8
         potential_immob(c) = 0.0_r8
         actual_immob(c) = 0.0_r8
         sminn_to_plant(c) = 0.0_r8
      end do

      do j = 1, nlevdecomp
         dzsoi_decomp(j) = 0.05_r8 + 0.01_r8 * real(j, r8)
      end do

      do k = 1, ndecomp_cascade_transitions
         if (mod(k,4) == 1) then
            cascade_receiver_pool(k) = i_cop_mic
         else if (mod(k,4) == 2) then
            cascade_receiver_pool(k) = i_oli_mic
         else
            cascade_receiver_pool(k) = k
         end if
      end do

      do j = 1, nlevdecomp
         do c = begc, endc
            nfixation_prof(c,j) = 1.0_r8 / real(nlevdecomp, r8)

            smin_nh4_vr(c,j) = 2.0e-4_r8 * (1.0_r8 + real(mod(c+j,11),r8) / 20.0_r8)
            smin_no3_vr(c,j) = 3.0e-4_r8 * (1.0_r8 + real(mod(c+2*j,13),r8) / 20.0_r8)
            sminn_vr(c,j)    = smin_nh4_vr(c,j) + smin_no3_vr(c,j)

            potential_immob_vr(c,j) = 5.0e-8_r8 * (1.0_r8 + real(mod(c*j,19),r8) / 19.0_r8)

            pot_f_nit_vr(c,j)   = 1.0e-8_r8 * (1.0_r8 + real(mod(c+j,7),r8) / 7.0_r8)
            pot_f_denit_vr(c,j) = 1.0e-8_r8 * (1.0_r8 + real(mod(c+3*j,5),r8) / 5.0_r8)

            n2_n2o_ratio_denit_vr(c,j) = 10.0_r8

            fpi_vr(c,j) = 0.0_r8
            f_nit_vr(c,j) = 0.0_r8
            f_denit_vr(c,j) = 0.0_r8
            actual_immob_vr(c,j) = 0.0_r8
            sminn_to_plant_vr(c,j) = 0.0_r8
            supplement_to_sminn_vr(c,j) = 0.0_r8
            sminn_to_denit_excess_vr(c,j) = 0.0_r8

            actual_immob_no3_vr(c,j) = 0.0_r8
            actual_immob_nh4_vr(c,j) = 0.0_r8
            smin_no3_to_plant_vr(c,j) = 0.0_r8
            smin_nh4_to_plant_vr(c,j) = 0.0_r8
            f_n2o_denit_vr(c,j) = 0.0_r8
            f_n2o_nit_vr(c,j) = 0.0_r8

            do k = 1, ndecomp_cascade_transitions
               c_overflow_vr(c,j,k) = 0.0_r8
               pmnf_decomp_cascade(c,j,k) = 2.0e-8_r8 * &
                    (1.0_r8 + real(mod(c+j+k,23),r8) / 23.0_r8)
               p_decomp_cn_gain(c,j,k) = 25.0_r8 + real(mod(c+j+k,9),r8)
            end do
         end do
      end do

   end subroutine initialize_inputs

   subroutine reset_outputs()
!$ACC KERNELS
      fpg = 0.0_r8
      fpi = 0.0_r8
      fpi_vr = 0.0_r8

      c_overflow_vr = 0.0_r8

      f_nit_vr = 0.0_r8
      f_denit_vr = 0.0_r8

      potential_immob = 0.0_r8
      actual_immob = 0.0_r8
      sminn_to_plant = 0.0_r8
      sminn_to_denit_excess_vr = 0.0_r8

      actual_immob_no3_vr = 0.0_r8
      actual_immob_nh4_vr = 0.0_r8
      smin_no3_to_plant_vr = 0.0_r8
      smin_nh4_to_plant_vr = 0.0_r8

      f_n2o_denit_vr = 0.0_r8
      f_n2o_nit_vr = 0.0_r8

      supplement_to_sminn_vr = 0.0_r8
      sminn_to_plant_vr = 0.0_r8
      actual_immob_vr = 0.0_r8
!$ACC END KERNELS

      ! Keep potential_immob_vr as an input field.
      ! Keep sminn_vr, smin_nh4_vr, smin_no3_vr as input fields.

   end subroutine reset_outputs

   function compute_checksum() result(sumval)

      real(r8) :: sumval

      sumval = 0.0_r8

      sumval = sumval + sum(fpg)
      sumval = sumval + sum(fpi)
      sumval = sumval + sum(fpi_vr)
      sumval = sumval + sum(actual_immob)
      sumval = sumval + sum(potential_immob)
      sumval = sumval + sum(sminn_to_plant)
      sumval = sumval + sum(sminn_to_plant_vr)
      sumval = sumval + sum(actual_immob_vr)
      sumval = sumval + sum(f_nit_vr)
      sumval = sumval + sum(f_denit_vr)
      sumval = sumval + sum(f_n2o_denit_vr)
      sumval = sumval + sum(f_n2o_nit_vr)
      sumval = sumval + sum(c_overflow_vr)

   end function compute_checksum

   subroutine get_arg_int(pos, value)

      integer, intent(in)    :: pos
      integer, intent(inout) :: value
      character(len=64)      :: local_arg

      call get_command_argument(pos, local_arg)
      if (len_trim(local_arg) > 0) then
         read(local_arg, *) value
      end if

   end subroutine get_arg_int

   subroutine get_arg_logical(pos, value)

      integer, intent(in)    :: pos
      logical, intent(inout) :: value
      character(len=64)      :: local_arg
      integer                :: tmp

      call get_command_argument(pos, local_arg)
      if (len_trim(local_arg) > 0) then
         read(local_arg, *) tmp
         value = tmp /= 0
      end if

   end subroutine get_arg_logical

end program driver_SoilBiogeochemCompetition
