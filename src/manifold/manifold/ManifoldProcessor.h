#ifndef ManifoldProcessor_H_
#define ManifoldProcessor_H_

#include <Eigen/Core>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "manifold/Octree.h"

class ManifoldProcessor {
public:
    /// @brief Constructor, takes as input vertices and trianlges of the
    /// original model we want to convert
    ManifoldProcessor(const std::vector<Eigen::Vector3d>& vertices,
                      const std::vector<Eigen::Vector3i>& triangles);

    /// @brief Process the input geometry and convert it to a manifold, returns
    /// true if the resulting mesh actually turns to be a manifold surface.
    /// @out get the vertices and triangles of the watertight mesh
    std::tuple<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3i>>
    GetManifoldMesh(int depth);

    /// @brief Check if a triangle mesh is a manifold or not
    static bool is_manifold(const std::vector<Eigen::Vector3d>& overtices,
                            const std::vector<Eigen::Vector3i>& otriangles);

private:
    void Calc_Bounding_Box();
    void Build_Tree(int depth);
    void Construct_Manifold();
    void Project_Manifold();
    glm::dvec3 Find_Closest(int i);
    bool Split_Grid(std::map<Grid_Index, int>& vcolor,
                    std::vector<glm::dvec3>& nvertices,
                    std::vector<glm::ivec4>& nface_indices,
                    std::vector<std::set<int>>& v_faces_,
                    std::vector<glm::ivec3>& triangles);

    // expose to set to 1
    int g_sharp_ = 0;

    std::vector<glm::dvec3> vertices_;
    std::vector<glm::ivec3> face_indices_;
    std::vector<glm::dvec3> vertices_buf_;
    std::vector<glm::dvec3> colors_;
    std::vector<glm::ivec3> face_indices_buf_;
    std::vector<glm::dvec3> face_normals_;

    std::vector<std::set<int>> v_faces_;
    std::vector<Grid_Index> v_info_;

    glm::dvec3 min_corner_, max_corner_;
    std::unique_ptr<Octree> tree_;
};
#endif
