cmake_minimum_required(VERSION 2.8)
project(composite_viewer_and_modeller)

set(CMAKE_BUILD_TYPE Debug)

add_definitions(-std=c++11)

aux_source_directory(./src SOURCES)
aux_source_directory(./src/image/gui SOURCES)
aux_source_directory(./src/image/algorithms SOURCES)
aux_source_directory(./src/geometry SOURCES)
aux_source_directory(./src/modeller SOURCES)
aux_source_directory(./src/modeller/components SOURCES)
aux_source_directory(./src/modeller/gui SOURCES)
aux_source_directory(./src/modeller/optimization SOURCES)
aux_source_directory(./src/osg SOURCES)
aux_source_directory(./src/utility SOURCES)
aux_source_directory(./src/wx SOURCES)

set(SOURCES ${SOURCES}
            src/geometry/Primitives.hpp
            src/modeller/optimization/ExtractPlaneNormals.hpp
            src/modeller/optimization/OptimizationUtility.hpp
            src/wx/WxGuiId.hpp)

include_directories(/usr/local/include/eigen3)

find_package(CGAL COMPONENTS Core REQUIRED)
if(CGAL_FOUND )
    include(${CGAL_USE_FILE})
    set(LINK_LIBS ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
    message(STATUS "CGAL_LIBRARIES: " ${CGAL_LIBRARIES})
    message(STATUS "CGAL_3RD_PARTY_LIBRARIES: " ${CGAL_3RD_PARTY_LIBRARIES})
else()
    message(FATAL_ERROR "CGAL is not found!")
endif()

find_package(wxWidgets REQUIRED core gl)
if(wxWidgets_FOUND)
    include(${wxWidgets_USE_FILE})
    message(STATUS "wxWidgets_USE_FILE: " ${wxWidgets_USE_FILE})
    set(LINK_LIBS ${LINK_LIBS} ${wxWidgets_LIBRARIES})
    message(STATUS "wxWidgets_LIBRARIES: " ${wxWidgets_LIBRARIES})
else()
    message(FATAL_ERROR "wxWidgets is not found!")
endif()


#find_package(OpenSceneGraph REQUIRED osgDB osgUtil osgGA osgViewer)
#if(OPENSCENEGRAPH_FOUND)
#    include_directories(${OPENSCENEGRAPH_INCLUDE_DIRS})
#    set(LINK_LIBS ${LINK_LIBS} ${OPENSCENEGRAPH_LIBRARIES})
#    message(STATUS "OPENSCENEGRAPH_LIBRARIES: " ${OPENSCENEGRAPH_LIBRARIES})
#else(OPENSCENEGRAPH_FOUND)
#    message(FATAL_ERROR "OSG is not found!")
#endif(OPENSCENEGRAPH_FOUND)

# If the find_package OpenSceneGraph fails use the following
set(LINK_LIBS ${LINK_LIBS} libosg.so libosgViewer.so libosgDB.so libosgGA.so libosgUtil.so)

# Names of the OpenSceneGraph Libraries
# libosg.so
# libosgViewer.so
# libOpenThreads.so
# libosgAnimation.so
# libosgDB.so
# libosgFX.so
# libosgGA.so
# libosgManipulator.so
# libosgParticle.so
# libosgPresentation.so
# libosgQt.so
# libosgShadow.so
# libosgSim.so
# libosgTerrain.so
# libosgText.so
# libosgUtil.so
# libosgVolume.so
# libosgWidget.so

find_package(OTB REQUIRED)
if(OTB_FOUND)
    include(${OTB_USE_FILE})
    set(LINK_LIBS ${LINK_LIBS} ${OTB_LIBRARIES})
else(OTB_FOUND)
    message(FATAL_ERROR "OTB is not found!")
endif(OTB_FOUND)

find_package(OriLight REQUIRED)
if(OriLight_FOUND)
    include_directories(${OriLight_INCLUDE_DIRS})
    add_definitions(${OriLight_DEFINITIONS})
    set(LINK_LIBS ${LINK_LIBS} ${OriLight_LIBRARIES})
    message(STATUS "OriLight_LIBRARIES: " ${OriLight_LIBRARIES})
else()
    message(FATAL_ERROR "OriLight is not found!")
endif()

find_package(Ceres 1.9.0 REQUIRED)
if(Ceres_FOUND)
    include_directories(${CERES_INCLUDE_DIRS})
    set(LINK_LIBS ${LINK_LIBS} ${CERES_LIBRARIES})
    message(STATUS "CERES_LIBRARIES: " ${CERES_LIBRARIES})
else()
    message(FATAL_ERROR "Ceres solver is not found!")
endif()

add_executable(cvm ${SOURCES})
target_link_libraries(cvm ${LINK_LIBS})
