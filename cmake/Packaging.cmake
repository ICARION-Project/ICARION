# =====================================================================
# ICARION - Install and package configuration
# =====================================================================

set(ICARION_INSTALL_DATADIR "${CMAKE_INSTALL_DATADIR}/icarion"
    CACHE PATH "Install directory for ICARION examples, schemas, and runtime data")

install(DIRECTORY "${CMAKE_SOURCE_DIR}/data"
        DESTINATION "${ICARION_INSTALL_DATADIR}"
        COMPONENT runtime
        PATTERN "jobs" EXCLUDE
        PATTERN "precomputed_lr1264" EXCLUDE
        PATTERN "*.Zone.Identifier" EXCLUDE
        PATTERN "*:Zone.Identifier" EXCLUDE)

install(DIRECTORY "${CMAKE_SOURCE_DIR}/examples"
        DESTINATION "${ICARION_INSTALL_DATADIR}"
        COMPONENT examples
        PATTERN "results" EXCLUDE
        PATTERN "*.h5" EXCLUDE
        PATTERN "__pycache__" EXCLUDE)

install(DIRECTORY "${CMAKE_SOURCE_DIR}/schema"
        DESTINATION "${ICARION_INSTALL_DATADIR}"
        COMPONENT runtime
        PATTERN "__pycache__" EXCLUDE)

install(DIRECTORY "${CMAKE_SOURCE_DIR}/analysis"
        DESTINATION "${ICARION_INSTALL_DATADIR}"
        COMPONENT runtime
        PATTERN "__pycache__" EXCLUDE
        PATTERN "output" EXCLUDE
        PATTERN "ccs_scaling" EXCLUDE
        PATTERN "initial_velocity_sensitivity" EXCLUDE
        PATTERN "initial_velocity_sensitivity_denseEN" EXCLUDE
        PATTERN "length_scaling" EXCLUDE
        PATTERN "mass_scaling" EXCLUDE
        PATTERN "variance_diagnostics" EXCLUDE)

install(FILES "${CMAKE_SOURCE_DIR}/requirements-analysis.txt"
        DESTINATION "${ICARION_INSTALL_DATADIR}"
        COMPONENT runtime)

install(DIRECTORY "${CMAKE_SOURCE_DIR}/docs"
        DESTINATION "${CMAKE_INSTALL_DOCDIR}"
        COMPONENT docs)

install(FILES
        "${CMAKE_SOURCE_DIR}/README.md"
        "${CMAKE_SOURCE_DIR}/CHANGELOG.md"
        "${CMAKE_SOURCE_DIR}/LICENSE"
        DESTINATION "${CMAKE_INSTALL_DOCDIR}"
        COMPONENT docs)

if(UNIX AND NOT APPLE)
    install(PROGRAMS "${CMAKE_SOURCE_DIR}/packaging/linux/ICARION-Launcher.sh"
            DESTINATION "${CMAKE_INSTALL_BINDIR}"
            RENAME "icarion-launcher"
            COMPONENT runtime)
endif()

set(CPACK_PACKAGE_NAME "icarion")
set(CPACK_PACKAGE_VENDOR "ICARION Project")
set(CPACK_PACKAGE_CONTACT "ICARION Project")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "Ion Collision And Reaction IntegratiON simulation framework")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://github.com/ICARION-Project/ICARION")
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README.md")
set(CPACK_PACKAGE_CHECKSUM SHA256)
set(CPACK_STRIP_FILES TRUE)

set(CPACK_GENERATOR "TGZ")
if(UNIX AND NOT APPLE)
    list(APPEND CPACK_GENERATOR "DEB")
    set(CPACK_DEBIAN_PACKAGE_SECTION "science")
    set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
    set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
    set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)
endif()

include(CPack)
