add_executable(kdtree_test main.cpp load_ply.cpp)
target_include_directories(kdtree_test PRIVATE ../include . /usr/local/include)
target_link_libraries(kdtree_test kdtree)