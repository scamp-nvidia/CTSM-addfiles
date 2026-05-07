module SoilBiogeochemCompetition_mod

  !-----------------------------------------------------------------------
  ! Standalone extraction of SoilBiogeochemCompetition from
  ! src/soilbiogeochem/SoilBiogeochemCompetitionMod.F90.
  !
  ! Stage 4: shr_kind_mod dependency removed; r8 is now defined locally
  ! via selected_real_kind(12), which matches CTSM's shr_kind_r8 in
  ! practice. The module has no non-intrinsic dependencies.
  !-----------------------------------------------------------------------

  implicit none
  private

  integer, parameter, public :: r8 = selected_real_kind(12)

  public :: SoilBiogeochemCompetition

contains

  !-----------------------------------------------------------------------
   subroutine SoilBiogeochemCompetition( &
       ! sizes / index ranges
       begc, endc, nlevdecomp, ndecomp_cascade_transitions, &
       num_bgc_soilc, filter_bgc_soilc,                     &
       ! scalar config (was module state / runtime flags / params_inst)
       dt, bdnr,                                            &
       use_nitrif_denitrif, carbon_only,                    &
       decomp_method, mimics_decomp, i_cop_mic, i_oli_mic,  &
       compet_plant_no3, compet_plant_nh4,                  &
       compet_decomp_no3, compet_decomp_nh4,                &
       compet_denit, compet_nit,                            &
       ! 1D arrays
       dzsoi_decomp, cascade_receiver_pool, landunit,       &
       ! state-block fields (was soilbiogeochem_state_inst%*_col)
       fpg, fpi, fpi_vr, nfixation_prof, plant_ndemand,     &
       ! n-state fields
       sminn_vr, smin_nh4_vr, smin_no3_vr,                  &
       ! c-flux fields
       c_overflow_vr,                                       &
       ! n-flux fields
       pot_f_nit_vr, pot_f_denit_vr, f_nit_vr, f_denit_vr,  &
       potential_immob, actual_immob, sminn_to_plant,       &
       sminn_to_denit_excess_vr,                            &
       actual_immob_no3_vr, actual_immob_nh4_vr,            &
       smin_no3_to_plant_vr, smin_nh4_to_plant_vr,          &
       n2_n2o_ratio_denit_vr, f_n2o_denit_vr, f_n2o_nit_vr, &
       supplement_to_sminn_vr, sminn_to_plant_vr,           &
       potential_immob_vr, actual_immob_vr,                 &
       ! 3D arrays
       pmnf_decomp_cascade, p_decomp_cn_gain)
    !
    ! !ARGUMENTS:
    integer , intent(in) :: begc, endc                                  ! column index range (was bounds%begc:bounds%endc)
    integer , intent(in) :: nlevdecomp                                  ! number of biogeochemically active soil layers
    integer , intent(in) :: ndecomp_cascade_transitions                 ! number of decomposition cascade transitions
    integer , intent(in) :: num_bgc_soilc                               ! number of soil columns in filter
    integer , intent(in) :: filter_bgc_soilc(:)                         ! filter for soil columns
    real(r8), intent(in) :: dt                                          ! decomp timestep (seconds)
    real(r8), intent(in) :: bdnr                                        ! bulk denitrification rate (1/s)
    logical , intent(in) :: use_nitrif_denitrif                         ! true => use nitrif/denitrif branch
    logical , intent(in) :: carbon_only                                 ! true => carbon-only mode (was allocate_carbon_only())
    integer , intent(in) :: decomp_method                               ! type of decomposition method
    integer , intent(in) :: mimics_decomp                               ! id value of MIMICS decomposition method
    integer , intent(in) :: i_cop_mic                                   ! copiotrophic microbial pool index
    integer , intent(in) :: i_oli_mic                                   ! oligotrophic microbial pool index
    real(r8), intent(in) :: compet_plant_no3                            ! relative competitiveness of plants for NO3
    real(r8), intent(in) :: compet_plant_nh4                            ! relative competitiveness of plants for NH4
    real(r8), intent(in) :: compet_decomp_no3                           ! relative competitiveness of immobilizers for NO3
    real(r8), intent(in) :: compet_decomp_nh4                           ! relative competitiveness of immobilizers for NH4
    real(r8), intent(in) :: compet_denit                                ! relative competitiveness of denitrifiers for NO3
    real(r8), intent(in) :: compet_nit                                  ! relative competitiveness of nitrifiers for NH4
    real(r8), intent(in) :: dzsoi_decomp(:)                             ! per-layer thickness for biogeochemical layers
    integer , pointer    :: cascade_receiver_pool(:)                    ! which pool is C added to for a given decomposition step
    integer , pointer    :: landunit(:)                                 ! landunit index per column (was col%landunit)
    real(r8), pointer    :: fpg(:)                                      ! fraction of potential gpp
    real(r8), pointer    :: fpi(:)                                      ! fraction of potential immobilization
    real(r8), pointer    :: fpi_vr(:,:)                                 ! fraction of potential immobilization (per layer)
    real(r8), pointer    :: nfixation_prof(:,:)
    real(r8), pointer    :: plant_ndemand(:)                            ! column-level plant N demand
    real(r8), pointer    :: sminn_vr(:,:)                               ! (gN/m3) soil mineral N
    real(r8), pointer    :: smin_nh4_vr(:,:)                            ! (gN/m3) soil mineral NH4
    real(r8), pointer    :: smin_no3_vr(:,:)                            ! (gN/m3) soil mineral NO3
    real(r8), pointer    :: c_overflow_vr(:,:,:)                        ! (gC/m3/s) C rejected by microbes that cannot process it
    real(r8), pointer    :: pot_f_nit_vr(:,:)                           ! (gN/m3/s) potential soil nitrification flux
    real(r8), pointer    :: pot_f_denit_vr(:,:)                         ! (gN/m3/s) potential soil denitrification flux
    real(r8), pointer    :: f_nit_vr(:,:)                               ! (gN/m3/s) soil nitrification flux
    real(r8), pointer    :: f_denit_vr(:,:)                             ! (gN/m3/s) soil denitrification flux
    real(r8), pointer    :: potential_immob(:)
    real(r8), pointer    :: actual_immob(:)
    real(r8), pointer    :: sminn_to_plant(:)
    real(r8), pointer    :: sminn_to_denit_excess_vr(:,:)
    real(r8), pointer    :: actual_immob_no3_vr(:,:)
    real(r8), pointer    :: actual_immob_nh4_vr(:,:)
    real(r8), pointer    :: smin_no3_to_plant_vr(:,:)
    real(r8), pointer    :: smin_nh4_to_plant_vr(:,:)
    real(r8), pointer    :: n2_n2o_ratio_denit_vr(:,:)                  ! ratio of N2 to N2O production by denitrification [gN/gN]
    real(r8), pointer    :: f_n2o_denit_vr(:,:)                         ! flux of N2O from denitrification [gN/m3/s]
    real(r8), pointer    :: f_n2o_nit_vr(:,:)                           ! flux of N2O from nitrification [gN/m3/s]
    real(r8), pointer    :: supplement_to_sminn_vr(:,:)
    real(r8), pointer    :: sminn_to_plant_vr(:,:)
    real(r8), pointer    :: potential_immob_vr(:,:)
    real(r8), pointer    :: actual_immob_vr(:,:)
    real(r8), intent(in) :: pmnf_decomp_cascade(begc:,1:,1:)            ! potential mineral N flux from one pool to another (gN/m3/s)
    real(r8), intent(in) :: p_decomp_cn_gain(begc:,1:,1:)               ! C:N ratio of the flux gained by the receiver pool
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: nitrif_n2o_loss_frac = 6.e-4_r8              ! fraction of N lost as N2O in nitrification (Li et al., 2000)
    integer  :: c,p,l,pi,j,k                                            ! indices
    integer  :: fc                                                      ! filter column index
    real(r8) :: amnf_immob_vr                                           ! actual mineral N flux from immobilization (gN/m3/s)
    real(r8) :: n_deficit_vr                                            ! microbial N deficit, vertically resolved (gN/m3/s)
    real(r8) :: fpi_no3_vr(begc:endc,1:nlevdecomp)                      ! fraction of potential immobilization supplied by no3
    real(r8) :: fpi_nh4_vr(begc:endc,1:nlevdecomp)                      ! fraction of potential immobilization supplied by nh4
    real(r8) :: sum_nh4_demand(begc:endc,1:nlevdecomp)
    real(r8) :: sum_nh4_demand_scaled(begc:endc,1:nlevdecomp)
    real(r8) :: sum_no3_demand(begc:endc,1:nlevdecomp)
    real(r8) :: sum_no3_demand_scaled(begc:endc,1:nlevdecomp)
    real(r8) :: sum_ndemand_vr(begc:endc, 1:nlevdecomp)                 ! total column N demand (gN/m3/s) at a given level
    real(r8) :: nuptake_prof(begc:endc, 1:nlevdecomp)
    real(r8) :: sminn_tot(begc:endc)
    integer  :: nlimit(begc:endc,0:nlevdecomp)                          ! flag for N limitation
    integer  :: nlimit_no3(begc:endc,0:nlevdecomp)                      ! flag for NO3 limitation
    integer  :: nlimit_nh4(begc:endc,0:nlevdecomp)                      ! flag for NH4 limitation
    real(r8) :: residual_sminn_vr(begc:endc, 1:nlevdecomp)
    real(r8) :: residual_sminn(begc:endc)
    real(r8) :: residual_smin_nh4_vr(begc:endc, 1:nlevdecomp)
    real(r8) :: residual_smin_no3_vr(begc:endc, 1:nlevdecomp)
    real(r8) :: residual_smin_nh4(begc:endc)
    real(r8) :: residual_smin_no3(begc:endc)
    real(r8) :: residual_plant_ndemand(begc:endc)

    ! Attempts at vectorization. ARRAY_filtered means a version of ARRAY excluding elements that the
    ! relevant filter would skip over.
    real(r8) :: sminn_tot_filtered(1:num_bgc_soilc)
    real(r8) :: sminn_vr_filtered(1:num_bgc_soilc,nlevdecomp)
    !-----------------------------------------------------------------------

      ! column loops to resolve plant/heterotroph competition for mineral N

         ! init sminn_tot(_filtered) and sminn_vr_filtered
         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            sminn_tot(c) = 0.
            do j = 1, nlevdecomp
               sminn_vr_filtered(fc,j) = sminn_vr(c,j)
            end do
         end do
         sminn_tot_filtered(:) = 0.

         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               sminn_tot_filtered(fc) = sminn_tot_filtered(fc) + sminn_vr_filtered(fc,j) * dzsoi_decomp(j)
            end do
         end do

         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            sminn_tot(c) = sminn_tot_filtered(fc)
         end do

         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               if (sminn_tot(c)  >  0.) then
                  nuptake_prof(c,j) = sminn_vr(c,j) / sminn_tot(c)
               else
                  nuptake_prof(c,j) = nfixation_prof(c,j)
               endif
            end do
         end do

         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               sum_ndemand_vr(c,j) = plant_ndemand(c) * nuptake_prof(c,j) + potential_immob_vr(c,j)
            end do
         end do

         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               l = landunit(c)
               if (sum_ndemand_vr(c,j)*dt < sminn_vr(c,j)) then
                  nlimit(c,j) = 0
                  fpi_vr(c,j) = 1.0_r8
                  actual_immob_vr(c,j) = potential_immob_vr(c,j)
                  sminn_to_plant_vr(c,j) = plant_ndemand(c) * nuptake_prof(c,j)
               else if ( carbon_only ) then !.or. &
                  nlimit(c,j) = 1
                  fpi_vr(c,j) = 1.0_r8
                  actual_immob_vr(c,j) = potential_immob_vr(c,j)
                  sminn_to_plant_vr(c,j) =  plant_ndemand(c) * nuptake_prof(c,j)
                  supplement_to_sminn_vr(c,j) = sum_ndemand_vr(c,j) - (sminn_vr(c,j)/dt)
               else
                  ! N availability can not satisfy the sum of immobilization and
                  ! plant growth demands, so these two demands compete for available
                  ! soil mineral N resource.

                  nlimit(c,j) = 1
                  if (sum_ndemand_vr(c,j) > 0.0_r8) then
                     actual_immob_vr(c,j) = (sminn_vr(c,j)/dt)*(potential_immob_vr(c,j) / sum_ndemand_vr(c,j))
                  else
                     actual_immob_vr(c,j) = 0.0_r8
                  end if

                  if (potential_immob_vr(c,j) > 0.0_r8) then
                     fpi_vr(c,j) = actual_immob_vr(c,j) / potential_immob_vr(c,j)
                  else
                     fpi_vr(c,j) = 0.0_r8
                  end if

                  sminn_to_plant_vr(c,j) = (sminn_vr(c,j)/dt) - actual_immob_vr(c,j)
               end if
            end do
         end do

         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               sminn_to_plant(c) = sminn_to_plant(c) + sminn_to_plant_vr(c,j) * dzsoi_decomp(j)
            end do
         end do

         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            residual_sminn(c) = 0._r8
         end do

         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            residual_plant_ndemand(c) = plant_ndemand(c) - sminn_to_plant(c)
         end do
         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               if (residual_plant_ndemand(c)  >  0._r8 ) then
                  if (nlimit(c,j) .eq. 0) then
                     residual_sminn_vr(c,j) = max(sminn_vr(c,j) - (actual_immob_vr(c,j) + sminn_to_plant_vr(c,j) ) * dt, 0._r8)
                     residual_sminn(c) = residual_sminn(c) + residual_sminn_vr(c,j) * dzsoi_decomp(j)
                  else
                     residual_sminn_vr(c,j)  = 0._r8
                  endif
               endif
            end do
         end do

         ! distribute residual N to plants
         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               if ( residual_plant_ndemand(c)  >  0._r8 .and. residual_sminn(c)  >  0._r8 .and. nlimit(c,j) .eq. 0) then
                  sminn_to_plant_vr(c,j) = sminn_to_plant_vr(c,j) + residual_sminn_vr(c,j) * &
                       min(( residual_plant_ndemand(c) *  dt ) / residual_sminn(c), 1._r8) / dt
               endif
            end do
         end do

         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            sminn_to_plant(c) = 0._r8
         end do
         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               sminn_to_plant(c) = sminn_to_plant(c) + sminn_to_plant_vr(c,j) * dzsoi_decomp(j)
               sum_ndemand_vr(c,j) = potential_immob_vr(c,j) + sminn_to_plant_vr(c,j)
            end do
         end do

         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               if ((sminn_to_plant_vr(c,j) + actual_immob_vr(c,j))*dt < sminn_vr(c,j)) then
                  sminn_to_denit_excess_vr(c,j) = max(bdnr*((sminn_vr(c,j)/dt) - sum_ndemand_vr(c,j)),0._r8)
               else
                  sminn_to_denit_excess_vr(c,j) = 0._r8
               endif
            end do
         end do

         do j = 1, nlevdecomp
            do fc=1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               actual_immob(c) = actual_immob(c) + actual_immob_vr(c,j) * dzsoi_decomp(j)
               potential_immob(c) = potential_immob(c) + potential_immob_vr(c,j) * dzsoi_decomp(j)
            end do
         end do

         do fc=1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            ! calculate the fraction of potential growth that can be
            ! acheived with the N available to plants
            if (plant_ndemand(c) > 0.0_r8) then
               fpg(c) = sminn_to_plant(c) / plant_ndemand(c)
            else
               fpg(c) = 1.0_r8
            end if

            ! calculate the fraction of immobilization realized (for diagnostic purposes)
            if (potential_immob(c) > 0.0_r8) then
               fpi(c) = actual_immob(c) / potential_immob(c)
            else
               fpi(c) = 1.0_r8
            end if
         end do

  end subroutine SoilBiogeochemCompetition

end module SoilBiogeochemCompetition_mod
