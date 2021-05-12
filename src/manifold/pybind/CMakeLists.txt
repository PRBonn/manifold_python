pybind11_add_module(manifold_pybind manifold_pybind.cpp)
target_link_libraries(manifold_pybind PRIVATE manifold ${Open3D_LIBRARIES})
