 oomph-convert -z soln*.dat; makePvd soln soln.pvd
 oomph-convert -z -p2 interface*.dat; makePvd interface interface.pvd
 oomph-convert -z fluid_coarse*.dat; makePvd fluid_coarse fluid_coarse.pvd
