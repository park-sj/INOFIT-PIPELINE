# Dependencies: CGAL 5.0.2, Eigen 3.3.7, Boost 1.71.0, jsoncpp 1.7.4, and CMake 3.16.3
# Optional:
# (1) CUDA 10.2
# (2) OpenCL 1.2, and ViennaCL 1.7.1
# (3) CGAL 5.1, OpenGR, and libpointmatcher 1.1.0
# Development environment:
# (1) macOS 10.14.5, Visual Studio Code, packages installed via homebrew, with OpenGR and libpointmatcher installed manually
# (2) Ubuntu 20.04, Visual Studio Code, packages installed via apt, with CUDA installed manually

cmake_minimum_required(VERSION 3.16.3)

project(Data_Driven_Nonrigid_ICP CXX)

if(NOT CMAKE_BUILD_TYPE)
    message("-- No build type specified. Default is CMAKE_BUILD_TYPE=RelWithDebInfo")
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif(NOT CMAKE_BUILD_TYPE)

if(WIN32)
    set(CMAKE_CXX_COMPILER MSVC)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_FLAGS "")
    set(CMAKE_CXX_FLAGS_DEBUG "/Zi /EHsc")
    set(CMAKE_CXX_FLAGS_RELEASE "/EHsc")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/EHsc")
    set(CMAKE_CXX_FLAGS_MINSIZEREL "/EHsc")
else()
    set(CMAKE_CXX_COMPILER g++)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_FLAGS "-Wall")
    set(CMAKE_CXX_FLAGS_DEBUG "-g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O3")
    set(CMAKE_CXX_FLAGS_MINSIZEREL "")
endif()

#set(CMAKE_CUDA_HOST_COMPILER g++)
#set(CMAKE_CUDA_STANDARD 11)

# default
include_directories(include)
if(APPLE)
    include_directories(/usr/local/include)
    link_directories(/usr/local/lib)
elseif(UNIX AND NOT APPLE)
    include_directories(/usr/include)
    link_directories(/usr/lib/x86_64-linux-gnu)
endif()

if(WIN32)
    # CGAL 5.0.2
    set(CGAL_ROOT "C:\\vcpkg\\installed\\x64-windows")
    set(CGAL_INCLUDEDIR "C:\\vcpkg\\installed\\x64-windows\\include")
    set(CGAL_LIBRARYDIR "C:\\vcpkg\\installed\\x64-windows\\lib")

    # Eigen 3.3.7
    set(BOOST_ROOT "C:\\vcpkg\\installed\\x64-windows")
    set(BOOST_INCLUDEDIR "C:\\vcpkg\\installed\\x64-windows\\include")
    set(BOOST_LIBRARYDIR "C:\\vcpkg\\installed\\x64-windows\\lib")

    # Boost 1.71.0
    set(Eigen3_ROOT "C:\\vcpkg\\installed\\x64-windows")
    set(Eigen3_INCLUDEDIR "C:\\vcpkg\\installed\\x64-windows\\include\\eigen3")
    set(Eigen3_LIBRARYDIR "C:\\vcpkg\\installed\\x64-windows\\lib")

    # jsoncpp 1.7.4
    set(jsoncpp_ROOT "C:\\vcpkg\\installed\\x64-windows")
    set(jsoncpp_INCLUDEDIR "C:\\vcpkg\\installed\\x64-windows\\include")
    set(jsoncpp_LIBRARYDIR "C:\\vcpkg\\installed\\x64-windows\\lib")


    find_package(CGAL REQUIRED)
    find_package(Eigen3 3.3.7 REQUIRED NO_MODULE)
    find_package(Boost REQUIRED COMPONENTS thread system date_time chrono)
    find_package(jsoncpp REQUIRED)


    
else()
    # CGAL 5.0.2
    find_package(CGAL REQUIRED)
    include_directories(${CGAL_INCLUDE_DIR})

    # Eigen 3.3.7
    find_package(Eigen3 3.3.7 REQUIRED NO_MODULE)
    include_directories(${Eigen3_INCLUDE_DIR})

    # Boost 1.71.0
    find_package(Boost REQUIRED COMPONENTS thread system date_time chrono)
    include_directories(${Boost_INCLUDE_DIR})

    # jsoncpp 1.7.4
    #if(NOT jsoncpp_INCLUDE_DIR)
    #    find_package(jsoncpp REQUIRED)
    #    include_directories(${jsoncpp_INCLUDE_DIR})
    #    if(UNIX AND NOT APPLE)
    #        include_directories(/usr/include/jsoncpp)
    #    endif()
    #else()
    #    set(jsoncpp_FOUND TRUE)
    #    include_directories(${jsoncpp_INCLUDE_DIR})
    #endif()

    #find_package(jsoncpp REQUIRED)
    set(jsoncpp_ROOT "/home/shkim/anaconda3/envs/PCA")
    set(jsoncpp_INCLUDEDIR "/home/shkim/anaconda3/envs/PCA/include")
    set(jsoncpp_LIBRARYDIR "/home/shkim/anaconda3/envs/PCA/lib")
    find_package(jsoncpp REQUIRED)

    # CUDA 10.2
    #find_package(CUDA)
    #if(CUDA_FOUND)
    #    include_directories(${CUDA_INCLUDE_DIR})
    #    include_directories(/usr/local/cuda-10.2/include)
    #    enable_language(CUDA)
    #endif()

    # OpenCL 1.2
    find_package(OpenCL 1.2)
    if(OpenCL_FOUND)
        include_directories(${OpenCL_INCLUDE_DIR})
    endif()

    # ViennaCL 1.7.1
    if(NOT ViennaCL_INCLUDE_DIR)
        find_package(ViennaCL 1.7.1)
        if(ViennaCL_FOUND)
            add_definitions(-DVIENNACL_HAVE_EIGEN)
            if(OpenCL_FOUND)
                add_definitions(-DVIENNACL_WITH_OPENCL)
            endif()
            include_directories(${ViennaCL_INCLUDE_DIR})
        endif()
    else()
        set(ViennaCL_FOUND TRUE)
        add_definitions(-DVIENNACL_HAVE_EIGEN)
        if(OpenCL_FOUND)
            add_definitions(-DVIENNACL_WITH_OPENCL)
        endif()
        include_directories(${ViennaCL_INCLUDE_DIR})
    endif()

    # OpenGR
    if(CGAL_VERSION VERSION_GREATER_EQUAL "5.1")
        set(OpenGR_DIR /usr/local/manually_installed/OpenGR/build/install/lib/cmake)
        find_package(OpenGR)
        if(OpenGR_FOUND)
            include(/usr/local/lib/cmake/CGAL/CGAL_OpenGR_support.cmake)
        endif()
    endif()

    # libpointmatcher 1.1.0
    if(CGAL_VERSION VERSION_GREATER_EQUAL "5.1")
        find_package(libpointmatcher 1.1.0)
        if(libpointmatcher_FOUND)
            #include(/usr/local/lib/cmake/CGAL/CGAL_pointmatcher_support.cmake)
            include(/home/shkim/Libraries/CGAL-5.1.1/cmake/modules/CGAL_pointmatcher_support.cmake)
        endif()
    endif()

endif()

link_libraries(CGAL::CGAL)
link_libraries(Eigen3::Eigen)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin/test)

add_executable(generate_distorted_mesh
    src/test/generate_distorted_mesh.cpp
    src/lib/mesh_distortion.cpp
    src/lib/tools.cpp
)

add_executable(find_correspondence
    src/test/find_correspondence.cpp
    src/lib/mesh_correspondence.cpp
    src/lib/tools.cpp
)

add_executable(find_closest_point
    src/test/find_closest_point.cpp
    src/lib/mesh_correspondence.cpp
    src/lib/tools.cpp
)

add_executable(compute_pca_on_meshes
    src/test/compute_pca_on_meshes.cpp
    src/lib/mesh_pca.cpp
    src/lib/mesh_distortion.cpp
    src/lib/mesh_tools.cpp
    src/lib/tools.cpp
)

add_executable(test_pca_transpose_trick
    src/test/test_pca_transpose_trick.cpp
)

add_executable(read_regions
    src/test/read_regions.cpp
    src/lib/region_tools.cpp
)

add_executable(merge_regions
    src/test/merge_regions.cpp
    src/lib/region_tools.cpp
)

add_executable(remove_regions
    src/test/remove_regions.cpp
    src/lib/region_tools.cpp
)

add_executable(test_widen_region
    src/test/test_widen_region.cpp
    src/lib/region_tools.cpp
)

#if(CUDA_FOUND)
#    add_executable(test_cusolver
#        src/test/test_cusolver.cpp
#        src/lib/cuda_solver.cpp
#        src/lib/cuda_solver_internal.cu
#        src/lib/tools.cpp
#    )
#    target_link_directories(test_cusolver PRIVATE /usr/local/cuda-10.2/lib64)
#    target_link_libraries(test_cusolver cusolver cusparse cublas)
#    set_target_properties(test_cusolver PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#endif()

if(OpenCL_FOUND AND ViennaCL_FOUND)
    add_executable(test_viennacl
        src/test/test_viennacl.cpp
        src/lib/tools.cpp
    )
    if(APPLE)
        target_link_libraries(test_viennacl "-framework OpenCL")
    elseif(UNIX AND NOT APPLE)
        target_link_libraries(test_viennacl OpenCL)
    endif()
endif()

if(CGAL_VERSION VERSION_GREATER_EQUAL "5.1")
    if(OpenGR_FOUND AND libpointmatcher_FOUND)
        add_executable(test_3rd_party_mesh_alignment
            src/test/test_3rd_party_mesh_alignment.cpp
            src/lib/mesh_alignment_3rd_party.cpp
            src/lib/tools.cpp
        )
        target_link_libraries(test_3rd_party_mesh_alignment
            CGAL::OpenGR_support
            CGAL::pointmatcher_support
        )
    endif()
endif()


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin/app)

add_executable(preprocess
    src/app/preprocess.cpp
    src/lib/parameter_reader.cpp
    src/lib/template_manager.cpp
    src/lib/region_pca_manager.cpp
    src/lib/mesh_pca.cpp
    src/lib/mesh_tools.cpp
    src/lib/region_tools.cpp
    src/lib/eigen_tools.cpp
    src/lib/tools.cpp
)
target_link_libraries(preprocess /home/shkim/anaconda3/envs/PCA/lib/libjsoncpp.so)

add_executable(registration
    src/app/registration.cpp
    src/lib/parameter_reader.cpp
    src/lib/target_manager.cpp
    src/lib/region_pca_manager.cpp
    src/lib/data_driven_nonrigid_icp.cpp
    src/lib/mesh_correspondence.cpp
    src/lib/mesh_pca.cpp
    src/lib/mesh_tools.cpp
    src/lib/region_tools.cpp
    src/lib/eigen_tools.cpp
    src/lib/tools.cpp
)
target_link_libraries(registration /home/shkim/anaconda3/envs/PCA/lib/libjsoncpp.so)

add_executable(transfer_color_to_all_off
    src/app/transfer_color_to_all_off.cpp
    src/lib/region_tools.cpp
    src/lib/tools.cpp
)

add_executable(nose_tip_align
    src/app/nose_tip_align.cpp
    src/lib/mesh_tools.cpp
    src/lib/region_tools.cpp
    src/lib/tools.cpp
)

