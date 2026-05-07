V100 - 73x speedup.
scamp@dev-sky5:/local/home/scamp/HACKATHON_CODE/CTSM-addfiles/perf_testing/Seth_Testing/SoilBiogeochemCompetition_mod$ ./make_and_run_original.sh
SoilBiogeochemCompetition benchmark
  ncol                         = 50000
  nlevdecomp                   = 20
  ndecomp_cascade_transitions  = 8
  nrepeat                      = 100
  nwarmup                      = 10
  use_nitrif_denitrif          = F
  carbon_only                  = F
  mimics branch                = F
  elapsed seconds              =   6.68786600E+00
  avg seconds / call           =   6.68786600E-02
  calls / second               =   1.49524527E+01
  checksum                     =   1.10000015E+06
scamp@dev-sky5:/local/home/scamp/HACKATHON_CODE/CTSM-addfiles/perf_testing/Seth_Testing/SoilBiogeochemCompetition_mod$ ./make_and_run_openacc.sh
SoilBiogeochemCompetition benchmark
  ncol                         = 50000
  nlevdecomp                   = 20
  ndecomp_cascade_transitions  = 8
  nrepeat                      = 100
  nwarmup                      = 10
  use_nitrif_denitrif          = F
  carbon_only                  = F
  mimics branch                = F
  elapsed seconds              =   9.12590000E-02
  avg seconds / call           =   9.12590000E-04
  calls / second               =   1.09578233E+03
  checksum                     =   1.10000015E+06
