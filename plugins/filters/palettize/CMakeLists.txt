set(kritapalettize_SOURCES palettize.cpp)
ki18n_wrap_ui(kritapalettize_SOURCES palettize.ui)
kis_add_library(kritapalettize MODULE ${kritapalettize_SOURCES})
target_link_libraries(kritapalettize kritaui)
install(TARGETS kritapalettize  DESTINATION ${KRITA_PLUGIN_INSTALL_DIR})
install( FILES
    palettize.action
DESTINATION  ${KDE_INSTALL_DATADIR}/krita/actions)
