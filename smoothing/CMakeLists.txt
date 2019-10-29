get_filename_component(EXERCISENAME ${CMAKE_CURRENT_LIST_DIR} NAME)

file(GLOB_RECURSE HEADERS "*.h")
file(GLOB_RECURSE SHADERS "*.glsl")

#================================================================
# Smoothing
#================================================================

add_executable(${EXERCISENAME}
    main.cpp
    mesh_processing.cpp
    viewer.cpp
    ${HEADERS}
    ${SHADERS}
)
target_link_libraries(${EXERCISENAME} surface_mesh)
# Lastly, additional libraries may have been built for you.  In addition to linking
# against NanoGUI, we need to link against those as well.
target_link_libraries(${EXERCISENAME} nanogui ${NANOGUI_EXTRA_LIBS})

#================================================================
# Unit Test
#================================================================

add_executable(smoothing_test
    smoothing_test.cpp
    mesh_processing.cpp
    viewer.cpp
    ${HEADERS}
    ${SHADERS})
target_link_libraries(smoothing_test surface_mesh)
# Lastly, additional libraries may have been built for you.  In addition to linking
# against NanoGUI, we need to link against those as well.
target_link_libraries(smoothing_test nanogui ${NANOGUI_EXTRA_LIBS})