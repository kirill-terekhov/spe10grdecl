cmake_minimum_required(VERSION 2.6)
project(SPE10)



add_executable(spe10grdecl spe10grdecl.cpp)




install(TARGETS spe10grdecl EXPORT SPE10-targets RUNTIME DESTINATION bin/spe10)
install(FILES spe_perm.dat spe_phi.dat DESTINATION bin/spe10)

