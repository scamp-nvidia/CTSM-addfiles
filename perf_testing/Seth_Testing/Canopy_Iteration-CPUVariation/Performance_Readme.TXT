Sky5:
scamp@dev-sky5:/local/home/scamp/HACKATHON_CODE/CTSM-addfiles/perf_testing/Seth_Testing/Canopy_Iteration-CPUVariation$ ./test
can_iter CPU reproducer
  npts                         = 65536
  nrepeat                      = 20
  nwarmup                      = 3
  itmax_canopy_fluxes          = 40
  timed can_iter seconds       =   5.934502600000E+00
  seconds per can_iter call    =   2.967251300000E-01
  patch iterations per second  =   2.566509954853E+06
  checksum                     =   9.604064155120E+06
  mean t_veg                  =   2.933178742005E+02
  mean qflx_evap_veg          =   8.559098951719E-05
  mean eflx_sh_veg            =  -5.177479938736E+01
  mean num_iter               =   1.162030029297E+01
scamp@dev-sky5:/local/home/scamp/HACKATHON_CODE/CTSM-addfiles/perf_testing/Seth_Testing/Canopy_Iteration-CPUVariation$ nvfortran template_cut.F90 -O3 -o test
scamp@dev-sky5:/local/home/scamp/HACKATHON_CODE/CTSM-addfiles/perf_testing/Seth_Testing/Canopy_Iteration-CPUVariation$ ./test
can_iter bucketed CPU reproducer
  npts                         = 65536
  nrepeat                      = 20
  nwarmup                      = 3
  itmax_canopy_fluxes          = 40
  timed can_iter seconds       =   5.117982000000E+00
  seconds per can_iter call    =   2.558991000000E-01
  patch iterations per second  =   2.975969825607E+06
  checksum                     =   9.604064155120E+06
  mean t_veg                  =   2.933178742005E+02
  mean qflx_evap_veg          =   8.559098951719E-05
  mean eflx_sh_veg            =  -5.177479938736E+01
  mean num_iter               =   1.162030029297E+01
