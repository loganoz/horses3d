#############################
# Add the different add-ons #
#############################

# Create horsesf90 compiling script
CONFIGURE_FILE(horsesf90.in.cmake horsesf90)

# ADD_SUBDIRECTORY(horses2tecplot)

###########
# Install #
###########

INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/horsesf90 DESTINATION ${BINDIR})
