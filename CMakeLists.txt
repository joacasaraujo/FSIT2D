cmake_minimum_required(VERSION 3.10)

project(FSIT2D)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

set(COMMON_SOURCES
    lbm/cell.cpp
    dem/body.cpp
    dem/interaction.cpp
    engine/engine.cpp
    property/shape.cpp
    scene/scene.cpp
    scene/output.cpp
)

file(GLOB SCENARIO_SOURCES
    CONFIGURE_DEPENDS
    scenarios/*.cpp
)

foreach(SCENARIO_FILE ${SCENARIO_SOURCES})

    get_filename_component(SCENARIO_NAME ${SCENARIO_FILE} NAME_WE)

    add_executable(${SCENARIO_NAME}
        ${SCENARIO_FILE}
        ${COMMON_SOURCES}
    )

    target_include_directories(${SCENARIO_NAME}
        PUBLIC
            ${PROJECT_SOURCE_DIR}
    )

    target_link_libraries(${SCENARIO_NAME}
        Eigen3::Eigen
    )

endforeach()
