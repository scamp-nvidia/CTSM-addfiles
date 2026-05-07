V100 - ~40x speedup in OpenACC GPU reproducer, and ~60x speedup in more optimized GPU reproducer.
scamp@dev-sky5:/local/home/scamp/HACKATHON_CODE/CTSM-addfiles/perf_testing/Seth_Testing/Canopy_Iteration$ ./original.sh ; ./openacc.sh ; ./openacc_opt.sh
can_iter CPU reproducer
  npts                         = 65536
  nrepeat                      = 20
  nwarmup                      = 3
  itmax_canopy_fluxes          = 40
  timed can_iter seconds       =   5.976871300000E+00
  seconds per can_iter call    =   2.988435650000E-01
  patch iterations per second  =   2.548316541465E+06
  checksum                     =   9.604064155120E+06
  mean t_veg                  =   2.933178742005E+02
  mean qflx_evap_veg          =   8.559098951719E-05
  mean eflx_sh_veg            =  -5.177479938736E+01
  mean num_iter               =   1.162030029297E+01
can_iter OpenACC GPU reproducer
  npts                         = 65536
  nrepeat                      = 20
  nwarmup                      = 3
  itmax_canopy_fluxes          = 40
  timed can_iter seconds       =   1.554968000000E-01
  seconds per can_iter call    =   7.774840000000E-03
  patch iterations per second  =   9.795031151766E+07
  checksum                     =   9.604064155120E+06
  mean t_veg                  =   2.933178742005E+02
  mean qflx_evap_veg          =   8.559098951719E-05
  mean eflx_sh_veg            =  -5.177479938736E+01
  mean num_iter               =   1.162030029297E+01
can_iter packed OpenACC GPU reproducer
  npts                         = 65536
  nrepeat                      = 20
  nwarmup                      = 3
  itmax_canopy_fluxes          = 40
  timed can_iter seconds       =   1.030824000000E-01
  seconds per can_iter call    =   5.154120000000E-03
  patch iterations per second  =   1.477551939031E+08
  checksum                     =   9.604064155120E+06
  mean t_veg                  =   2.933178742005E+02
  mean qflx_evap_veg          =   8.559098951719E-05
  mean eflx_sh_veg            =  -5.177479938736E+01
  mean num_iter               =   1.162030029297E+01
