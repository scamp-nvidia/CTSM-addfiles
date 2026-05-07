module SoilBiogeochemCompetition_mod

  implicit none
  private

  integer, parameter, public :: r8 = selected_real_kind(12)

  public :: SoilBiogeochemCompetition

contains

   subroutine SoilBiogeochemCompetition( &
       ! sizes / index ranges
       begc, endc, nlevdecomp, ndecomp_cascade_transitions, &
       num_bgc_soilc, filter_bgc_soilc,                     &
       ! scalar config
       dt, bdnr,                                            &
       use_nitrif_denitrif, carbon_only,                    &
       decomp_method, mimics_decomp, i_cop_mic, i_oli_mic,  &
       compet_plant_no3, compet_plant_nh4,                  &
       compet_decomp_no3, compet_decomp_nh4,                &
       compet_denit, compet_nit,                            &
       ! 1D arrays
       dzsoi_decomp, cascade_receiver_pool, landunit,       &
       ! state-block fields
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

    integer , intent(in) :: begc, endc
    integer , intent(in) :: nlevdecomp
    integer , intent(in) :: ndecomp_cascade_transitions
    integer , intent(in) :: num_bgc_soilc
    integer , intent(in) :: filter_bgc_soilc(:)

    real(r8), intent(in) :: dt
    real(r8), intent(in) :: bdnr

    logical , intent(in) :: use_nitrif_denitrif
    logical , intent(in) :: carbon_only

    integer , intent(in) :: decomp_method
    integer , intent(in) :: mimics_decomp
    integer , intent(in) :: i_cop_mic
    integer , intent(in) :: i_oli_mic

    real(r8), intent(in) :: compet_plant_no3
    real(r8), intent(in) :: compet_plant_nh4
    real(r8), intent(in) :: compet_decomp_no3
    real(r8), intent(in) :: compet_decomp_nh4
    real(r8), intent(in) :: compet_denit
    real(r8), intent(in) :: compet_nit

    real(r8), intent(in) :: dzsoi_decomp(:)

    integer , pointer :: cascade_receiver_pool(:)
    integer , pointer :: landunit(:)

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

    real(r8), intent(in) :: pmnf_decomp_cascade(begc:,1:,1:)
    real(r8), intent(in) :: p_decomp_cn_gain(begc:,1:,1:)

    integer  :: c, fc, j
    real(r8) :: total_sminn
    real(r8) :: tmp_plant
    real(r8) :: tmp_residual
    real(r8) :: tmp_actual_immob
    real(r8) :: tmp_potential_immob

    real(r8) :: sum_ndemand_vr(begc:endc, 1:nlevdecomp)
    real(r8) :: nuptake_prof(begc:endc, 1:nlevdecomp)
    real(r8) :: sminn_tot(begc:endc)

    integer  :: nlimit(begc:endc, 0:nlevdecomp)

    real(r8) :: residual_sminn_vr(begc:endc, 1:nlevdecomp)
    real(r8) :: residual_sminn(begc:endc)
    real(r8) :: residual_plant_ndemand(begc:endc)

    ! ------------------------------------------------------------------
    ! This implementation intentionally executes only the old top branch:
    !
    !    if (.not. use_nitrif_denitrif) then
    !
    ! The nitrif/denitrif branch has been removed for this benchmark path.
    ! ------------------------------------------------------------------

    !$acc data &
    !$acc present( &
    !$acc   filter_bgc_soilc(1:num_bgc_soilc), &
    !$acc   dzsoi_decomp(1:nlevdecomp), &
    !$acc   fpg(begc:endc), &
    !$acc   fpi(begc:endc), &
    !$acc   fpi_vr(begc:endc,1:nlevdecomp), &
    !$acc   nfixation_prof(begc:endc,1:nlevdecomp), &
    !$acc   plant_ndemand(begc:endc), &
    !$acc   sminn_vr(begc:endc,1:nlevdecomp), &
    !$acc   potential_immob(begc:endc), &
    !$acc   actual_immob(begc:endc), &
    !$acc   sminn_to_plant(begc:endc), &
    !$acc   sminn_to_denit_excess_vr(begc:endc,1:nlevdecomp), &
    !$acc   supplement_to_sminn_vr(begc:endc,1:nlevdecomp), &
    !$acc   sminn_to_plant_vr(begc:endc,1:nlevdecomp), &
    !$acc   potential_immob_vr(begc:endc,1:nlevdecomp), &
    !$acc   actual_immob_vr(begc:endc,1:nlevdecomp) ) &
    !$acc create( &
    !$acc   sum_ndemand_vr(begc:endc,1:nlevdecomp), &
    !$acc   nuptake_prof(begc:endc,1:nlevdecomp), &
    !$acc   sminn_tot(begc:endc), &
    !$acc   nlimit(begc:endc,0:nlevdecomp), &
    !$acc   residual_sminn_vr(begc:endc,1:nlevdecomp), &
    !$acc   residual_sminn(begc:endc), &
    !$acc   residual_plant_ndemand(begc:endc) )

    ! ------------------------------------------------------------------
    ! Total soil mineral N per column.
    !
    ! Original code accumulated over j into sminn_tot_filtered(fc).
    ! That is a race if j is parallel. Here each fc owns its reduction.
    ! ------------------------------------------------------------------

    !$acc parallel loop gang vector private(c,j,total_sminn)
    do fc = 1, num_bgc_soilc
       c = filter_bgc_soilc(fc)

       total_sminn = 0.0_r8

       !$acc loop seq
       do j = 1, nlevdecomp
          total_sminn = total_sminn + sminn_vr(c,j) * dzsoi_decomp(j)
       end do

       sminn_tot(c) = total_sminn
    end do

    ! ------------------------------------------------------------------
    ! Define plant N uptake profile.
    ! ------------------------------------------------------------------

    !$acc parallel loop gang vector private(c,j)
    do fc = 1, num_bgc_soilc
       c = filter_bgc_soilc(fc)

       !$acc loop seq
       do j = 1, nlevdecomp
          if (sminn_tot(c) > 0.0_r8) then
             nuptake_prof(c,j) = sminn_vr(c,j) / sminn_tot(c)
          else
             nuptake_prof(c,j) = nfixation_prof(c,j)
          end if
       end do
    end do

    ! ------------------------------------------------------------------
    ! Total N demand per column/layer.
    ! ------------------------------------------------------------------

    !$acc parallel loop gang vector private(c,j)
    do fc = 1, num_bgc_soilc
       c = filter_bgc_soilc(fc)

       !$acc loop seq
       do j = 1, nlevdecomp
          sum_ndemand_vr(c,j) = plant_ndemand(c) * nuptake_prof(c,j) + &
                                potential_immob_vr(c,j)
       end do
    end do

    ! ------------------------------------------------------------------
    ! Resolve competition between plants and immobilizers.
    ! ------------------------------------------------------------------

    !$acc parallel loop gang vector private(c,j)
    do fc = 1, num_bgc_soilc
       c = filter_bgc_soilc(fc)

       !$acc loop seq
       do j = 1, nlevdecomp

          if (sum_ndemand_vr(c,j) * dt < sminn_vr(c,j)) then

             nlimit(c,j) = 0
             fpi_vr(c,j) = 1.0_r8
             actual_immob_vr(c,j) = potential_immob_vr(c,j)
             sminn_to_plant_vr(c,j) = plant_ndemand(c) * nuptake_prof(c,j)

          else if (carbon_only) then

             nlimit(c,j) = 1
             fpi_vr(c,j) = 1.0_r8
             actual_immob_vr(c,j) = potential_immob_vr(c,j)
             sminn_to_plant_vr(c,j) = plant_ndemand(c) * nuptake_prof(c,j)
             supplement_to_sminn_vr(c,j) = sum_ndemand_vr(c,j) - &
                                           (sminn_vr(c,j) / dt)

          else

             nlimit(c,j) = 1

             if (sum_ndemand_vr(c,j) > 0.0_r8) then
                actual_immob_vr(c,j) = (sminn_vr(c,j) / dt) * &
                     (potential_immob_vr(c,j) / sum_ndemand_vr(c,j))
             else
                actual_immob_vr(c,j) = 0.0_r8
             end if

             if (potential_immob_vr(c,j) > 0.0_r8) then
                fpi_vr(c,j) = actual_immob_vr(c,j) / potential_immob_vr(c,j)
             else
                fpi_vr(c,j) = 0.0_r8
             end if

             sminn_to_plant_vr(c,j) = (sminn_vr(c,j) / dt) - &
                                      actual_immob_vr(c,j)

          end if
       end do
    end do

    ! ------------------------------------------------------------------
    ! First sum of N fluxes to plant.
    !
    ! Original:
    !   do j
    !      do fc
    !         sminn_to_plant(c) = sminn_to_plant(c) + ...
    !
    ! That races across j. Here each fc owns the vertical accumulation.
    ! Preserve original behavior by starting from existing sminn_to_plant(c).
    ! ------------------------------------------------------------------

    !$acc parallel loop gang vector private(c,j,tmp_plant)
    do fc = 1, num_bgc_soilc
       c = filter_bgc_soilc(fc)

       tmp_plant = sminn_to_plant(c)

       !$acc loop seq
       do j = 1, nlevdecomp
          tmp_plant = tmp_plant + sminn_to_plant_vr(c,j) * dzsoi_decomp(j)
       end do

       sminn_to_plant(c) = tmp_plant
    end do

    ! ------------------------------------------------------------------
    ! Residual plant demand and residual soil mineral N initialization.
    ! ------------------------------------------------------------------

    !$acc parallel loop gang vector private(c)
    do fc = 1, num_bgc_soilc
       c = filter_bgc_soilc(fc)
       residual_sminn(c) = 0.0_r8
       residual_plant_ndemand(c) = plant_ndemand(c) - sminn_to_plant(c)
    end do

    ! ------------------------------------------------------------------
    ! Sum residual mineral N by column.
    !
    ! Original accumulated residual_sminn(c) over j, so j could not be
    ! parallel for a fixed c. This version makes the reduction private
    ! to each fc.
    ! ------------------------------------------------------------------

    !$acc parallel loop gang vector private(c,j,tmp_residual)
    do fc = 1, num_bgc_soilc
       c = filter_bgc_soilc(fc)

       tmp_residual = 0.0_r8

       !$acc loop seq
       do j = 1, nlevdecomp
          if (residual_plant_ndemand(c) > 0.0_r8) then
             if (nlimit(c,j) == 0) then
                residual_sminn_vr(c,j) = max(sminn_vr(c,j) - &
                     (actual_immob_vr(c,j) + sminn_to_plant_vr(c,j)) * dt, &
                     0.0_r8)

                tmp_residual = tmp_residual + &
                     residual_sminn_vr(c,j) * dzsoi_decomp(j)
             else
                residual_sminn_vr(c,j) = 0.0_r8
             end if
          else
             residual_sminn_vr(c,j) = 0.0_r8
          end if
       end do

       residual_sminn(c) = tmp_residual
    end do

    ! ------------------------------------------------------------------
    ! Distribute residual N to plants.
    ! ------------------------------------------------------------------

    !$acc parallel loop gang vector private(c,j)
    do fc = 1, num_bgc_soilc
       c = filter_bgc_soilc(fc)

       !$acc loop seq
       do j = 1, nlevdecomp
          if (residual_plant_ndemand(c) > 0.0_r8 .and. &
              residual_sminn(c) > 0.0_r8 .and. &
              nlimit(c,j) == 0) then

             sminn_to_plant_vr(c,j) = sminn_to_plant_vr(c,j) + &
                  residual_sminn_vr(c,j) * &
                  min((residual_plant_ndemand(c) * dt) / residual_sminn(c), &
                      1.0_r8) / dt
          end if
       end do
    end do

    ! ------------------------------------------------------------------
    ! Re-sum plant fluxes and update realized demand.
    !
    ! Original zeroed sminn_to_plant(c), then accumulated over j.
    ! This version preserves that behavior without a j-race.
    ! ------------------------------------------------------------------

    !$acc parallel loop gang vector private(c,j,tmp_plant)
    do fc = 1, num_bgc_soilc
       c = filter_bgc_soilc(fc)

       tmp_plant = 0.0_r8

       !$acc loop seq
       do j = 1, nlevdecomp
          tmp_plant = tmp_plant + sminn_to_plant_vr(c,j) * dzsoi_decomp(j)
          sum_ndemand_vr(c,j) = potential_immob_vr(c,j) + &
                                sminn_to_plant_vr(c,j)
       end do

       sminn_to_plant(c) = tmp_plant
    end do

    ! ------------------------------------------------------------------
    ! Excess N lost to denitrification.
    ! ------------------------------------------------------------------

    !$acc parallel loop gang vector private(c,j)
    do fc = 1, num_bgc_soilc
       c = filter_bgc_soilc(fc)

       !$acc loop seq
       do j = 1, nlevdecomp
          if ((sminn_to_plant_vr(c,j) + actual_immob_vr(c,j)) * dt < &
              sminn_vr(c,j)) then

             sminn_to_denit_excess_vr(c,j) = &
                  max(bdnr * ((sminn_vr(c,j) / dt) - sum_ndemand_vr(c,j)), &
                      0.0_r8)
          else
             sminn_to_denit_excess_vr(c,j) = 0.0_r8
          end if
       end do
    end do

    ! ------------------------------------------------------------------
    ! Sum N fluxes to immobilization.
    !
    ! Original accumulated actual_immob(c) and potential_immob(c) over j.
    ! This version makes that accumulation private per column.
    !
    ! Preserve original behavior by starting from existing values.
    ! ------------------------------------------------------------------

    !$acc parallel loop gang vector private(c,j,tmp_actual_immob,tmp_potential_immob)
    do fc = 1, num_bgc_soilc
       c = filter_bgc_soilc(fc)

       tmp_actual_immob = actual_immob(c)
       tmp_potential_immob = potential_immob(c)

       !$acc loop seq
       do j = 1, nlevdecomp
          tmp_actual_immob = tmp_actual_immob + &
               actual_immob_vr(c,j) * dzsoi_decomp(j)

          tmp_potential_immob = tmp_potential_immob + &
               potential_immob_vr(c,j) * dzsoi_decomp(j)
       end do

       actual_immob(c) = tmp_actual_immob
       potential_immob(c) = tmp_potential_immob
    end do

    ! ------------------------------------------------------------------
    ! Column-level diagnostics.
    ! ------------------------------------------------------------------

    !$acc parallel loop gang vector private(c)
    do fc = 1, num_bgc_soilc
       c = filter_bgc_soilc(fc)

       if (plant_ndemand(c) > 0.0_r8) then
          fpg(c) = sminn_to_plant(c) / plant_ndemand(c)
       else
          fpg(c) = 1.0_r8
       end if

       if (potential_immob(c) > 0.0_r8) then
          fpi(c) = actual_immob(c) / potential_immob(c)
       else
          fpi(c) = 1.0_r8
       end if
    end do

    !$acc end data

  end subroutine SoilBiogeochemCompetition

end module SoilBiogeochemCompetition_mod
