# Copyright (c) 2018-2023 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause

add_library(hardspheres_common OBJECT hardspheres.cc)
# include library headers as SYSTEM headers to suppress compiler warnings
# https://gcc.gnu.org/onlinedocs/cpp/System-Headers.html
target_include_directories(hardspheres_common SYSTEM PUBLIC ${Boost_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
target_link_libraries(hardspheres_common PUBLIC ${Boost_LIBRARIES} ${HDF5_LIBRARIES})

add_executable(hardspheres main.cc solvent.cc config.cc snapshot.cc)
target_link_libraries(hardspheres PRIVATE hardspheres_common)

add_executable(hardspheres_test hardspheres_test.cc config.cc)
target_link_libraries(hardspheres_test PRIVATE hardspheres_common)
add_test(NAME hardspheres_test COMMAND hardspheres_test)

set_property(SOURCE main.cc APPEND PROPERTY COMPILE_DEFINITIONS VERSION="${VERSION}")
