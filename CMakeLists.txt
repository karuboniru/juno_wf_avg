cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
project(WfAverage)

set(CMAKE_CXX_STANDARD 23)

find_package(junosw)

add_library(WfAverage MODULE WfAverage.cxx)
target_link_libraries(WfAverage
  junosw::EDMUtil
  junosw::EvtNavigator
  junosw::ElecEvent
  junosw::Identifier
  junosw::Geometry
  junosw::RootIOSvc
  SNiPER::RootWriterLib
)
