# pthread
find_package(Threads)

# gmock/gtest
set(GMOCK_RootDir    ${PacBioConsensus_RootDir}/third-party/gmock-1.7.0)
set(GMOCK_IncludeDir ${GMOCK_RootDir})

file(GLOB           GMOCK_CC    "${GMOCK_RootDir}/*.cc")
file(GLOB_RECURSE   GMOCK_H     "${GMOCK_RootDir}/*.h")

include_directories(SYSTEM
    ${GMOCK_IncludeDir}
)

file(GLOB PacBioConsensusTests_CPP "*.cpp")

add_executable(test_pbconsensus EXCLUDE_FROM_ALL
    ${PacBioConsensusTests_CPP}
    ${GMOCK_CC}
    ${GMOCK_H}
)

target_link_libraries(test_pbconsensus
    ${PBCONSENSUS_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${CMAKE_DL_LIBS}
)

if (${ROOT_PROJECT_NAME} STREQUAL "PacBioConsensus")
    add_custom_target(check test_pbconsensus)
endif()
