add_subdirectory( basictools )
add_subdirectory( defaulttool )
add_subdirectory( selectiontools )
add_subdirectory( tool_crop )
add_subdirectory( tool_polygon )
add_subdirectory( tool_polyline )
add_subdirectory( tool_transform2 )
add_subdirectory( tool_dyna )
add_subdirectory( karbonplugins )
add_subdirectory( tool_lazybrush )
add_subdirectory( tool_smart_patch )
add_subdirectory( svgtexttool )
add_subdirectory( tool_enclose_and_fill )

install( FILES
    tools.action
DESTINATION  ${KDE_INSTALL_DATADIR}/krita/actions)
