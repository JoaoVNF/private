# These are the directories we're actually going to build:
# They are defined by configure.ac which scans 
# config/configure.ac_scripts/private.dir_list
# to set the PRIVATE_SUBDIRECTORIES variable
SUBDIRS =  $(PRIVATE_SUBDIRECTORIES)

# Extra actions for clean and distclean
distclean-local:
	rm -f *~
clean-local:
	rm -f *~
doc_dirs:
	echo $(PRIVATE_SUBDIRECTORIES)

EXTRA_DIST=external_src
