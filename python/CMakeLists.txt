define_module(pymaxflow LIBRARY LINKS maxflow)
add_custom_target(rename_pymaxflow_lib ALL COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/python/libpymaxflow.so ${CMAKE_BINARY_DIR}/python/pymaxflow.so)
add_dependencies(rename_pymaxflow_lib pymaxflow)
