# build KalTest library

SET( lib_input_dirs geomlib kallib kaltracklib utils detector )

FOREACH( lib_input_dir ${lib_input_dirs} )
    LIST( APPEND ROOT_DICT_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/${lib_input_dir} )
ENDFOREACH()

#MESSAGE( STATUS "ROOT_DICT_INCLUDE_DIRS: ${ROOT_DICT_INCLUDE_DIRS}" )


FOREACH( lib_input_dir ${lib_input_dirs} )
    AUX_SOURCE_DIRECTORY( ${lib_input_dir} lib_sources )

    PREPARE_ROOT_DICT_HEADERS( ${lib_input_dir} )

    INSTALL_DIRECTORY( ${lib_input_dir}/ DESTINATION "include/kaltest"
        FILES_MATCHING PATTERN "*.h" PATTERN "LinkDef.h" EXCLUDE
    )

    GEN_ROOT_DICT_SOURCES( ${lib_input_dir}Dict.cxx )
    
    LIST( APPEND lib_sources ${ROOT_DICT_OUTPUT_SOURCES} )

ENDFOREACH()

INCLUDE_DIRECTORIES( ${ROOT_DICT_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )

########
AUX_SOURCE_DIRECTORY( processors lib_sources)
LIST( APPEND lib_sources ${ROOT_DICT_OUTPUT_SOURCES} )


#MESSAGE( STATUS "KalTest lib sources: ${calo_procesor_sources}" )

ADD_SHARED_LIBRARY( KalTestCalo ${lib_sources} )
INSTALL_SHARED_LIBRARY( KalTestCalo DESTINATION lib )
TARGET_LINK_LIBRARIES( KalTestCalo ${ROOT_LIBRARIES} )

ADD_EXECUTABLE( testcalo detector/EXVTXKalTest.cxx )
TARGET_LINK_LIBRARIES( testcalo KalTestCalo )
INSTALL( TARGETS testcalo DESTINATION bin )

IF( APPLE ) #---- need special linker flags for ROOT dictionary on MacOS
    SET_TARGET_PROPERTIES( KalTestCalo PROPERTIES
        LINK_FLAGS "-single_module -undefined dynamic_lookup -bind_at_load"
    )
ENDIF( APPLE )


