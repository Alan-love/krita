set(kritamultiarch_SRCS
    KisSupportedArchitectures.cpp
    KoMultiArchBuildSupport.cpp
)

kis_add_library(kritamultiarch SHARED ${kritamultiarch_SRCS})

generate_export_header(kritamultiarch)

target_link_libraries(kritamultiarch PUBLIC xsimd)

target_link_libraries(kritamultiarch
    PUBLIC
        Qt${QT_MAJOR_VERSION}::Core
    PRIVATE
        KF${KF_MAJOR}::ConfigCore
        kritaglobal
)

set_target_properties(kritamultiarch PROPERTIES
    VERSION ${GENERIC_KRITA_LIB_VERSION} SOVERSION ${GENERIC_KRITA_LIB_SOVERSION}
)

target_include_directories(kritamultiarch
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
)

install(TARGETS kritamultiarch ${INSTALL_TARGETS_DEFAULT_ARGS})
