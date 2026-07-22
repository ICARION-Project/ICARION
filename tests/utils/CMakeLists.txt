# SHA256 hash utility tests

add_executable(test_hash test_hash.cpp)
target_link_libraries(test_hash PRIVATE icarion_core Catch2::Catch2WithMain)
target_include_directories(test_hash PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_test(NAME test_hash COMMAND test_hash)

add_executable(test_user_annotation test_user_annotation.cpp)
target_link_libraries(test_user_annotation PRIVATE icarion_core Catch2::Catch2WithMain)
target_include_directories(test_user_annotation PRIVATE ${CMAKE_SOURCE_DIR}/src)
add_test(NAME UserAnnotation COMMAND test_user_annotation)
