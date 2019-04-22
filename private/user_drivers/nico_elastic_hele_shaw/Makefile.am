# Name of executables: The prefix "noinst" means the executables don't
# have to be installed anywhere.
noinst_PROGRAMS=elastic_hele_shaw_displ_steady_continuation elastic_hele_shaw_symmetric_displ_steady_continuation elastic_hele_shaw_displ_flow_rate elastic_hele_shaw_displ_unsteady membrane_collapse

#---------------------------------------------------------------------------

# Local sources 
elastic_hele_shaw_displ_steady_continuation_SOURCES = elastic_hele_shaw_displ_steady_continuation.cc

# Required libraries:
elastic_hele_shaw_displ_steady_continuation_LDADD = -L@libdir@  -lnico_elastic_hele_shaw -lconstitutive -lsolid  -lrigid_body -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# flags
elastic_hele_shaw_displ_steady_continuation_CXXFLAGS="-g"

#---------------------------------------------------------------------------

# Local sources 
elastic_hele_shaw_symmetric_displ_steady_continuation_SOURCES = elastic_hele_shaw_symmetric_displ_steady_continuation.cc

# Required libraries:
elastic_hele_shaw_symmetric_displ_steady_continuation_LDADD = -L@libdir@  -lnico_elastic_hele_shaw -lconstitutive -lsolid  -lrigid_body -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# flags
elastic_hele_shaw_symmetric_displ_steady_continuation_CXXFLAGS="-g"

#---------------------------------------------------------------------------

# Local sources 
elastic_hele_shaw_displ_flow_rate_SOURCES = elastic_hele_shaw_displ_flow_rate.cc

# Required libraries:
elastic_hele_shaw_displ_flow_rate_LDADD = -L@libdir@  -lnico_elastic_hele_shaw -lconstitutive -lsolid  -lrigid_body -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# flags
elastic_hele_shaw_displ_flow_rate_CXXFLAGS="-g"

#---------------------------------------------------------------------------

# Local sources 
elastic_hele_shaw_displ_unsteady_SOURCES = elastic_hele_shaw_displ_unsteady.cc

# Required libraries:
elastic_hele_shaw_displ_unsteady_LDADD = -L@libdir@  -lnico_elastic_hele_shaw -lconstitutive -lsolid  -lrigid_body -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# flags
#elastic_hele_shaw_displ_unsteady_CXXFLAGS="-g"

#---------------------------------------------------------------------------

# Local sources 
membrane_collapse_SOURCES = membrane_collapse.cc

# Required libraries:
membrane_collapse_LDADD = -L@libdir@  -lnico_elastic_hele_shaw -lconstitutive -lsolid  -lrigid_body -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# flags
membrane_collapse_CXXFLAGS="-g"

#---------------------------------------------------------------------------


# header files before splitting in to header and source 
EXTRA_DIST = post_process.bash run_steady.bash run_steady_restart.bash run_steady_symmetric.bash run_steady_symmetric_restart.bash deflection_with_pressure.pvsm \
hele_shaw_elements.h  hele_shaw_interface_elements.h  hele_shaw_flux_elements.h  hs_fvk_elements.h  hs_displ_fvk_elements.h  Thele_shaw_elements.h README documentation.pdf





# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
INCLUDES = -I@includedir@  
