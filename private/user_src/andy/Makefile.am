# A few file definitions
#-----------------------

# Define the sources
sources =  \
helmholtz_elements.cc spectral_helmholtz_elements.cc \
refineable_helmholtz_elements.cc Thelmholtz_elements.cc \
hello_world.cc

# Define the headers
headers =  \
helmholtz_elements.h helmholtz_flux_elements.h spectral_helmholtz_elements.h \
refineable_helmholtz_elements.h refineable_spectral_helmholtz_elements.h \
Thelmholtz_elements.h hello_world.h

# Define name of library
libname = andy

# Combine headers and sources
headers_and_sources = $(headers) $(sources)

# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = libandy.la

# Sources that the library depends on: 
#-------------------------------------
libandy_la_SOURCES = $(headers_and_sources)

# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below $(includedir)
include_HEADERS = $(libname).h 

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/andy

#These are the header files that are to be placed in subdirectory
library_include_HEADERS=$(headers) 


# Required libraries -- [assuming that we want to link against stuff in generic
#-------------------     add other oomph-lib libraries if you need them....]
# Include path for library headers -- need to refer to include files
# in their respective source directories as they will not have been 
# installed yet! 
INCLUDES = -I$(top_builddir)/src/generic

junk:
	echo $(top_builddir)


# Combined header file
#---------------------
# Rule for building combined header (dummy target never gets made so
# the combined header is remade every time)
$(libname).h: dummy_$(libname).h

dummy_$(libname).h:  $(headers)
	echo  $(libname) $(headers) > all_$(libname).aux
	$(AWK) -f $(top_builddir)/bin/headers.awk < \
	       all_$(libname).aux > $(libname).h
	rm all_$(libname).aux


# Extra hook for install: Optionally replace headers by symbolic links
#---------------------------------------------------------------------
if SYMBOLIC_LINKS_FOR_HEADERS
install-data-hook: 
	(cd $(library_includedir) && rm -f $(headers) )
	(echo "$(headers)" > include_files.list )
	($(top_builddir)/bin/change_headers_to_links.sh `pwd`)
	($(LN_S) `cat include_files.list.aux` $(library_includedir) )
	(rm -r include_files.list.aux include_files.list ) 
else
install-data-hook: 
endif


# Tidy up
#--------
clean-local:
	rm -f $(libname).h


