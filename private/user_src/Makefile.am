# These are the src directories we're actually going to build:
# They are defined by configure.ac which scans 
# config/configure.ac_scripts/private_user_src.dir_list
# to set the PRIVATE_USER_SRC_SUBDIRECTORIES variable
SUBDIRS =  $(PRIVATE_USER_SRC_SUBDIRECTORIES)

# Extra actions for clean and distclean
distclean-local:
	rm -f *~
clean-local:
	rm -f *~
doc_dirs:
	echo $(PRIVATE_USER_SRC_SUBDIRECTORIES)