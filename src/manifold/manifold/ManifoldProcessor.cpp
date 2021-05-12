#include "manifold/ManifoldProcessor.h"

#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include <queue>
#include <string>
#include <tuple>
#include <vector>

namespace {

glm::dvec3 Closest_Point(const glm::dvec3* triangle,
                         const glm::dvec3& sourcePosition) {
    glm::dvec3 edge0 = triangle[1] - triangle[0];
    glm::dvec3 edge1 = triangle[2] - triangle[0];
    glm::dvec3 v0 = triangle[0] - sourcePosition;

    double a = glm::dot(edge0, edge0);
    double b = glm::dot(edge0, edge1);
    double c = glm::dot(edge1, edge1);
    double d = glm::dot(edge0, v0);
    double e = glm::dot(edge1, v0);

    double det = a * c - b * b;
    double s = b * e - c * d;
    double t = b * d - a * e;

    const double kMax = 1.0F;
    const double kMin = 0.0F;

    if (s + t < det) {
        if (s < 0.f) {
            if (t < 0.f) {
                if (d < 0.f) {
                    s = std::clamp(-d / a, kMin, kMax);
                    t = 0.f;
                } else {
                    s = 0.f;
                    t = std::clamp(-e / c, kMin, kMax);
                }
            } else {
                s = 0.f;
                t = std::clamp(-e / c, kMin, kMax);
            }
        } else if (t < 0.f) {
            s = std::clamp(-d / a, kMin, kMax);
            t = 0.f;
        } else {
            double invDet = 1.f / det;
            s *= invDet;
            t *= invDet;
        }
    } else {
        if (s < 0.f) {
            double tmp0 = b + d;
            double tmp1 = c + e;
            if (tmp1 > tmp0) {
                double numer = tmp1 - tmp0;
                double denom = a - 2 * b + c;
                s = std::clamp(numer / denom, kMin, kMax);
                t = 1 - s;
            } else {
                t = std::clamp(-e / c, kMin, kMax);
                s = 0.f;
            }
        } else if (t < 0.f) {
            if (a + d > b + e) {
                double numer = c + e - b - d;
                double denom = a - 2 * b + c;
                s = std::clamp(numer / denom, kMin, kMax);
                t = 1 - s;
            } else {
                s = std::clamp(-e / c, kMin, kMax);
                t = 0.f;
            }
        } else {
            double numer = c + e - b - d;
            double denom = a - 2 * b + c;
            s = std::clamp(numer / denom, kMin, kMax);
            t = 1.f - s;
        }
    }

    return triangle[0] + s * edge0 + t * edge1;
}
}  // namespace

ManifoldProcessor::ManifoldProcessor(
        const std::vector<Eigen::Vector3d>& vertices,
        const std::vector<Eigen::Vector3i>& triangles) {
    vertices_.clear();
    vertices_.reserve(vertices.size());
    for (const auto& v : vertices) {
        vertices_.emplace_back(glm::dvec3(v[0], v[1], v[2]));
    }
    face_indices_.clear();
    face_indices_.reserve(triangles.size());
    for (const auto& f : triangles) {
        face_indices_.emplace_back(glm::ivec3(f[0], f[1], f[2]));
    }
}

void ManifoldProcessor::Calc_Bounding_Box() {
    min_corner_ = glm::dvec3(1e30, 1e30, 1e30);
    max_corner_ = -min_corner_;
    for (auto& vertice : vertices_) {
        for (int j = 0; j < 3; ++j) {
            if (vertice[j] < min_corner_[j]) {
                min_corner_[j] = vertice[j];
            }
            if (vertice[j] > max_corner_[j]) {
                max_corner_[j] = vertice[j];
            }
        }
    }
    glm::dvec3 length = max_corner_ - min_corner_;
    min_corner_ -= length * 0.2;
    max_corner_ += length * 0.2;
}

void ManifoldProcessor::Build_Tree(int depth) {
    Calc_Bounding_Box();
    tree_ = std::make_unique<Octree>(min_corner_, max_corner_, face_indices_,
                                     0.01);
    for (int iter = 0; iter < depth; ++iter) {
        tree_->Split(vertices_);
    }

    tree_->BuildConnection();
    tree_->BuildEmptyConnection();

    std::list<Octree*> empty_list;
    std::set<Octree*> empty_set;
    for (int i = 0; i < 6; ++i) {
        tree_->ExpandEmpty(empty_list, empty_set, i);
    }

    while ((int)empty_list.size() > 0) {
        Octree* empty = empty_list.front();
        empty->exterior = 1;
        for (auto& empty_neighbor : empty->empty_neighbors) {
            if (empty_set.find(empty_neighbor) == empty_set.end()) {
                empty_list.push_back(empty_neighbor);
                empty_set.insert(empty_neighbor);
            }
        }
        empty_list.pop_front();
    }
}

void ManifoldProcessor::Construct_Manifold() {
    std::map<Grid_Index, int> vcolor;
    std::vector<glm::dvec3> nvertices;
    std::vector<glm::ivec4> nface_indices;
    std::vector<glm::ivec3> triangles;
    tree_->ConstructFace(vcolor, glm::ivec3(0, 0, 0), nvertices, nface_indices,
                         v_faces_);
    Split_Grid(vcolor, nvertices, nface_indices, v_faces_, triangles);
    std::vector<int> hash_v(nvertices.size(), 0);
    for (auto& triangle : triangles) {
        for (int j = 0; j < 3; ++j) {
            hash_v[triangle[j]] = 1;
        }
    }
    vertices_.clear();
    for (int i = 0; i < (int)hash_v.size(); ++i) {
        if (hash_v[i] != 0) {
            hash_v[i] = (int)vertices_.size();
            v_faces_[vertices_.size()] = v_faces_[i];
            v_info_[vertices_.size()] = v_info_[i];
            vertices_.push_back(nvertices[i]);
            colors_.emplace_back(1, 1, 1);
        }
    }
    for (auto& triangle : triangles) {
        for (int j = 0; j < 3; ++j) {
            triangle[j] = hash_v[triangle[j]];
        }
    }
    face_indices_ = triangles;
}

glm::dvec3 ManifoldProcessor::Find_Closest(int i) {
    glm::dvec3 cpoint = glm::dvec3(1e20, 1e20, 1e20);
    glm::dvec3 tris[3];
    glm::dvec3 normal;
    for (auto it = v_faces_[i].begin(); it != v_faces_[i].end(); ++it) {
        int face_ind = *it;
        for (int j = 0; j < 3; ++j) {
            tris[j] = vertices_buf_[face_indices_buf_[face_ind][j]];
        }
        glm::dvec3 p = Closest_Point(tris, vertices_[i]);
        if (glm::length(p - vertices_[i]) <
            glm::length(cpoint - vertices_[i])) {
            normal = glm::normalize(
                    glm::cross(tris[1] - tris[0], tris[2] - tris[0]));
            if (glm::dot(normal, vertices_[i] - cpoint) < 0) {
                normal = -normal;
            }
            cpoint = p;
        }
    }
    return cpoint + normal * 5e-4;
}

void ManifoldProcessor::Project_Manifold() {
    const int ITER_NUM = 20;
    double len = glm::length(vertices_[face_indices_[0][1]] -
                             vertices_[face_indices_[0][0]]);
    double min_len = glm::length(vertices_[face_indices_[0][2]] -
                                 vertices_[face_indices_[0][0]]);
    if (min_len < len) {
        len = min_len;
    }
    colors_.clear();
    colors_.resize(vertices_.size(), glm::dvec3(1, 1, 1));

    std::vector<std::vector<int>> vertex_faces(vertices_.size());
    face_normals_.resize(face_indices_.size());
    for (int i = 0; i < (int)face_indices_.size(); ++i) {
        int id[3];
        id[0] = face_indices_[i][0];
        id[1] = face_indices_[i][1];
        id[2] = face_indices_[i][2];
        for (int j : id) {
            vertex_faces[j].push_back(i);
        }
    }
    std::vector<int> vertices_hash(vertices_.size(), 0);
    double min_step = 2.0 / ITER_NUM;

    for (int iter = 0; iter < ITER_NUM; ++iter) {
        for (int i = 0; i < (int)face_indices_.size(); ++i) {
            int id[3];
            id[0] = face_indices_[i][0];
            id[1] = face_indices_[i][1];
            id[2] = face_indices_[i][2];
            face_normals_[i] = glm::normalize(
                    glm::cross(vertices_[id[1]] - vertices_[id[0]],
                               vertices_[id[2]] - vertices_[id[0]]));
        }

        std::vector<int> invalid_vertices;
        std::vector<int> invalid_indices(vertices_.size(), -1);

        for (int i = 0; i < (int)vertices_.size(); ++i) {
            if (vertices_hash[i] != 0) {
                continue;
            }
            glm::dvec3 cpoint = Find_Closest(i);
            glm::dvec3 move_dir = cpoint - vertices_[i];
            double orig_step = glm::length(move_dir);
            move_dir /= orig_step;
            double step = orig_step;
            bool flag = step < 1e15;
            if (g_sharp_ != 0) {
                vertices_[i] = cpoint;
                continue;
            }
            glm::dvec3 normal(0, 0, 0);
            for (int j : vertex_faces[i]) {
                normal += face_normals_[j];
            }
            normal = glm::normalize(normal);
            if (flag) {
                bool convex = true;
                for (int j = 0; j < (int)vertex_faces[i].size(); ++j) {
                    for (int k = 0; k < 3; ++k) {
                        if (glm::dot(vertices_[face_indices_[vertex_faces[i][j]]
                                                            [k]] -
                                             vertices_[i],
                                     normal) > 0) {
                            convex = false;
                        }
                    }
                    if (!convex) {
                        break;
                    }
                }
                if (convex) {
                    for (int j : vertex_faces[i]) {
                        if (glm::dot(face_normals_[j], move_dir) > 0) {
                            flag = false;
                            break;
                        }
                    }
                } else {
                    flag = false;
                    for (int j : vertex_faces[i]) {
                        if (glm::dot(face_normals_[j], move_dir) < 0) {
                            flag = true;
                            break;
                        }
                    }
                }
            }
            if (flag) {
                if (step > min_step * len) {
                    step = min_step * len;
                }
                for (int j = 0; j < (int)vertex_faces[i].size(); ++j) {
                    glm::ivec3& face_index = face_indices_[vertex_faces[i][j]];
                    int t = 0;
                    while (face_index[t] != i) {
                        t += 1;
                    }
                    glm::dvec3 dir =
                            glm::normalize(vertices_[face_index[(t + 2) % 3]] -
                                           vertices_[face_index[(t + 1) % 3]]);
                    glm::dvec3 h =
                            vertices_[face_index[(t + 1) % 3]] +
                            glm::dot(vertices_[i] -
                                             vertices_[face_index[(t + 1) % 3]],
                                     dir) *
                                    dir -
                            vertices_[i];
                    double h_len = glm::length(h);
                    h /= h_len;
                    double h_step = glm::dot(h, move_dir) * step;
                    if (h_step > h_len * 0.7) {
                        step *= (h_len * 0.7) / h_step;
                        invalid_indices[i] = (int)invalid_vertices.size();
                        invalid_vertices.push_back(i);
                        colors_[i] = glm::dvec3(0, 0, 1);
                    }
                }
                if (fabs(step - orig_step) < 1e-6) {
                    vertices_[i] = cpoint + len * normal;
                    if (step > 1e-4) {
                        step -= 1e-4;
                    } else {
                        step = 0;
                    }
                    vertices_hash[i] = 1;
                } else {
                    vertices_[i] +=
                            step *
                            move_dir;  //* 0.1 * glm::normalize(move_dir);
                }
                for (int face : vertex_faces[i]) {
                    face_normals_[face] = glm::normalize(glm::cross(
                            vertices_[face_indices_[face][1]] -
                                    vertices_[face_indices_[face][0]],
                            vertices_[face_indices_[face][2]] -
                                    vertices_[face_indices_[face][0]]));
                }
            } else {
                invalid_indices[i] = (int)invalid_vertices.size();
                invalid_vertices.push_back(i);
                colors_[i] = glm::dvec3(0, 0, 1);
            }
        }
        //	cout << "Invalid " << invalid_vertices.size() << "\n";
        std::vector<int> invalid_colors(invalid_vertices.size(), -1);
        int c = 0;
        for (int i = 0; i < (int)invalid_vertices.size(); ++i) {
            if (invalid_colors[i] == -1) {
                //			colors_[invalid_vertices[i]] =
                // glm::dvec3(0,0,1);
                invalid_colors[i] = c;
                std::vector<int> queue;
                int f = 0;
                queue.push_back(i);
                while (f != (int)queue.size()) {
                    int id = invalid_vertices[queue[f]];
                    for (int j : vertex_faces[id]) {
                        for (int k = 0; k < 3; ++k) {
                            int index = invalid_indices[face_indices_[j][k]];
                            if (index != -1 && invalid_colors[index] == -1) {
                                invalid_colors[index] = c;
                                queue.push_back(index);
                            }
                        }
                    }
                    f++;
                }
                for (auto it = queue.rbegin(); it != queue.rend(); ++it) {
                    glm::dvec3 midpoint(0, 0, 0);
                    int count = 0;
                    int id = invalid_vertices[*it];
                    for (int j : vertex_faces[id]) {
                        for (int k = 0; k < 3; ++k) {
                            int vind = face_indices_[j][k];
                            if (invalid_indices[vind] == -1 ||
                                invalid_colors[invalid_indices[vind]] == -1) {
                                midpoint += vertices_[vind];
                                count += 1;
                            }
                        }
                    }
                    glm::dvec3 move_dir =
                            midpoint / (double)count - vertices_[id];
                    invalid_colors[*it] = -1;
                    double l = glm::length(move_dir);
                    if (l == 0 || count == 0) {
                        continue;
                    }
                    move_dir /= l;
                    for (int j = 0; j < (int)vertex_faces[id].size(); ++j) {
                        glm::ivec3& face_index =
                                face_indices_[vertex_faces[id][j]];
                        int t = 0;
                        while (face_index[t] != id) {
                            t += 1;
                        }
                        glm::dvec3 dir = glm::normalize(
                                vertices_[face_index[(t + 2) % 3]] -
                                vertices_[face_index[(t + 1) % 3]]);
                        glm::dvec3 h =
                                vertices_[face_index[(t + 1) % 3]] +
                                glm::dot(vertices_[id] -
                                                 vertices_[face_index[(t + 1) %
                                                                      3]],
                                         dir) *
                                        dir -
                                vertices_[id];
                        double h_len = glm::length(h);
                        h /= h_len;
                        double h_step = glm::dot(h, move_dir) * l;
                        if (h_step > h_len * 0.7) {
                            l *= (h_len * 0.7) / h_step;
                        }
                    }
                    move_dir *= l;
                    vertices_[id] += move_dir;
                }
                for (int i : queue) {
                    invalid_colors[i] = c;
                }
            }
            c++;
        }
    }
}

bool ManifoldProcessor::Split_Grid(std::map<Grid_Index, int>& vcolor,
                                   std::vector<glm::dvec3>& nvertices,
                                   std::vector<glm::ivec4>& nface_indices,
                                   std::vector<std::set<int>>& v_faces_,
                                   std::vector<glm::ivec3>& triangles) {
    double unit_len = 0;
    v_info_.resize(vcolor.size());
    for (auto& it : vcolor) {
        v_info_[it.second] = it.first;
    }
    std::set<int> marked_v;
    std::map<std::pair<int, int>, std::list<std::pair<int, int>>> edge_info;
    for (int i = 0; i < (int)nface_indices.size(); ++i) {
        for (int j = 0; j < 4; ++j) {
            int x = nface_indices[i][j];
            int y = nface_indices[i][(j + 1) % 4];
            if (x > y) {
                int temp = x;
                x = y;
                y = temp;
            }
            std::pair<int, int> edge = std::make_pair(x, y);
            auto it = edge_info.find(edge);
            if (it != edge_info.end()) {
                it->second.emplace_back(i, j);
            } else {
                std::list<std::pair<int, int>> buf;
                buf.emplace_back(i, j);
                edge_info.insert(std::make_pair(edge, buf));
            }
        }
    }
    for (auto& it : edge_info) {
        if (it.second.size() > 2) {
            marked_v.insert(it.first.first);
            marked_v.insert(it.first.second);
        }
    }
    triangles.clear();
    double half_len = glm::length(nvertices[nface_indices[0][1]] -
                                  nvertices[nface_indices[0][0]]) *
                      0.5;
    for (auto& nface_indice : nface_indices) {
        int t = 0;
        while (t < 4 && marked_v.find(nface_indice[t]) == marked_v.end()) {
            ++t;
        }
        if (t == 4) {
            triangles.emplace_back(nface_indice[0], nface_indice[2],
                                   nface_indice[1]);
            triangles.emplace_back(nface_indice[0], nface_indice[3],
                                   nface_indice[2]);
            continue;
        }
        int ind[4];
        for (int j = 0; j < 4; ++j) {
            ind[j] = nface_indice[(t + j) % 4];
        }
        bool flag1 = marked_v.find(ind[1]) != marked_v.end();
        bool flag2 = marked_v.find(ind[2]) != marked_v.end();
        bool flag3 = marked_v.find(ind[3]) != marked_v.end();
        Grid_Index pt1 = (v_info_[ind[0]] + v_info_[ind[1]]) / 2;
        Grid_Index pt2 = (v_info_[ind[0]] + v_info_[ind[3]]) / 2;
        Grid_Index pt3 = (v_info_[ind[2]] + v_info_[ind[3]]) / 2;
        Grid_Index pt4 = (v_info_[ind[1]] + v_info_[ind[2]]) / 2;
        int ind1 = 0;
        int ind2 = 0;
        int ind3 = 0;
        int ind4 = 0;
        auto it = vcolor.find(pt1);
        if (it == vcolor.end()) {
            vcolor.insert(std::make_pair(pt1, nvertices.size()));
            v_info_.push_back(pt1);
            ind1 = (int)nvertices.size();
            nvertices.push_back((nvertices[ind[0]] + nvertices[ind[1]]) * 0.5);
            v_faces_.push_back(v_faces_[ind[0]]);
        } else {
            ind1 = it->second;
        }
        it = vcolor.find(pt2);
        if (it == vcolor.end()) {
            vcolor.insert(std::make_pair(pt2, nvertices.size()));
            v_info_.push_back(pt2);
            ind2 = (int)nvertices.size();
            v_faces_.push_back(v_faces_[ind[0]]);
            nvertices.push_back((nvertices[ind[0]] + nvertices[ind[3]]) * 0.5);
        } else {
            ind2 = it->second;
        }
        if (flag1 || flag2) {
            it = vcolor.find(pt4);
            if (it == vcolor.end()) {
                vcolor.insert(std::make_pair(pt4, nvertices.size()));
                v_info_.push_back(pt4);
                ind4 = (int)nvertices.size();
                nvertices.push_back((nvertices[ind[1]] + nvertices[ind[2]]) *
                                    0.5);
                if (flag1) {
                    v_faces_.push_back(v_faces_[ind[1]]);
                } else {
                    v_faces_.push_back(v_faces_[ind[2]]);
                }
            } else {
                ind4 = it->second;
            }
        }
        if (flag2 || flag3) {
            it = vcolor.find(pt3);
            if (it == vcolor.end()) {
                vcolor.insert(std::make_pair(pt3, nvertices.size()));
                v_info_.push_back(pt3);
                ind3 = (int)nvertices.size();
                nvertices.push_back((nvertices[ind[2]] + nvertices[ind[3]]) *
                                    0.5);
                if (flag2) {
                    v_faces_.push_back(v_faces_[ind[2]]);
                } else {
                    v_faces_.push_back(v_faces_[ind[3]]);
                }
            } else {
                ind3 = it->second;
            }
        }
        if (!flag1 && !flag2 && !flag3) {
            triangles.emplace_back(ind1, ind[2], ind[1]);
            triangles.emplace_back(ind2, ind[2], ind1);
            triangles.emplace_back(ind[3], ind[2], ind2);
        } else if (!flag1 && !flag2 && flag3) {
            triangles.emplace_back(ind1, ind2, ind3);
            triangles.emplace_back(ind1, ind3, ind[2]);
            triangles.emplace_back(ind1, ind[2], ind[1]);
        } else if (!flag1 && flag2 && !flag3) {
            triangles.emplace_back(ind1, ind4, ind[1]);
            triangles.emplace_back(ind1, ind2, ind4);
            triangles.emplace_back(ind2, ind[3], ind3);
            triangles.emplace_back(ind2, ind3, ind4);
        } else if (!flag1 && flag2 && flag3) {
            triangles.emplace_back(ind1, ind4, ind[1]);
            triangles.emplace_back(ind1, ind2, ind4);
            triangles.emplace_back(ind2, ind3, ind4);
        } else if (flag1 && !flag2 && !flag3) {
            triangles.emplace_back(ind1, ind2, ind4);
            triangles.emplace_back(ind4, ind2, ind[3]);
            triangles.emplace_back(ind4, ind[3], ind[2]);
        } else if (flag1 && !flag2 && flag3) {
            triangles.emplace_back(ind1, ind2, ind4);
            triangles.emplace_back(ind4, ind2, ind3);
            triangles.emplace_back(ind4, ind3, ind[2]);
        } else if (flag1 && flag2 && !flag3) {
            triangles.emplace_back(ind1, ind2, ind4);
            triangles.emplace_back(ind2, ind3, ind4);
            triangles.emplace_back(ind2, ind[3], ind3);
        } else if (flag1 && flag2 && flag3) {
            triangles.emplace_back(ind1, ind2, ind3);
            triangles.emplace_back(ind1, ind3, ind4);
        }
    }
    for (int it : marked_v) {
        glm::dvec3 p = nvertices[it];
        for (int dimx = -1; dimx < 2; dimx += 2) {
            for (int dimy = -1; dimy < 2; dimy += 2) {
                for (int dimz = -1; dimz < 2; dimz += 2) {
                    glm::dvec3 p1 =
                            p + glm::dvec3(dimx * half_len, dimy * half_len,
                                           dimz * half_len);
                    if (tree_->Is_Exterior(p1)) {
                        Grid_Index ind = v_info_[it];
                        Grid_Index ind1 = ind;
                        Grid_Index ind2 = ind;
                        Grid_Index ind3 = ind;
                        ind1.id[0] += dimx;
                        ind2.id[1] += dimy;
                        ind3.id[2] += dimz;
                        if (vcolor.find(ind1) == vcolor.end()) {
                            vcolor.insert(
                                    std::make_pair(ind1, nvertices.size()));
                            v_info_.push_back(ind1);

                            nvertices.emplace_back(p[0] + half_len * dimx, p[1],
                                                   p[2]);
                            v_faces_.push_back(v_faces_[it]);
                        }
                        if (vcolor.find(ind2) == vcolor.end()) {
                            vcolor.insert(
                                    std::make_pair(ind2, nvertices.size()));
                            v_info_.push_back(ind2);

                            nvertices.emplace_back(p[0], p[1] + half_len * dimy,
                                                   p[2]);
                            v_faces_.push_back(v_faces_[it]);
                        }
                        if (vcolor.find(ind3) == vcolor.end()) {
                            vcolor.insert(
                                    std::make_pair(ind3, nvertices.size()));
                            v_info_.push_back(ind3);

                            nvertices.emplace_back(p[0], p[1],
                                                   p[2] + half_len * dimz);
                            v_faces_.push_back(v_faces_[it]);
                        }
                        int id1 = vcolor[ind1];
                        int id2 = vcolor[ind2];
                        int id3 = vcolor[ind3];
                        glm::dvec3 norm =
                                glm::cross(nvertices[id2] - nvertices[id1],
                                           nvertices[id3] - nvertices[id1]);
                        if (glm::dot(norm, glm::dvec3(dimx, dimy, dimz)) < 0) {
                            triangles.emplace_back(id1, id3, id2);
                        } else {
                            triangles.emplace_back(id1, id2, id3);
                        }
                    }
                }
            }
        }
    }
    std::map<int, std::set<std::pair<int, int>>> ocs;
    std::map<int, std::set<std::pair<int, int>>> ecs;
    std::set<int> odds;
    std::set<int> evens;
    for (int i = 0; i < (int)nvertices.size(); ++i) {
        bool flag = false;
        for (int k = 0; k < 3; ++k) {
            if (v_info_[i].id[k] % 2 == 1) {
                flag = true;
            }
        }
        if (flag) {
            odds.insert(i);
            ocs.insert(std::make_pair(i, std::set<std::pair<int, int>>()));
        }
    }
    for (int i = 0; i < (int)nvertices.size(); ++i) {
        Grid_Index ind = v_info_[i];
        int flag = 0;
        while (flag < 3 && ind.id[flag] % 2 == 0) {
            flag++;
        }
        if (flag < 3) {
            continue;
        }
        for (int j = -2; j < 5; j += 4) {
            if (flag < 3) {
                break;
            }
            for (int k = 0; k < 3; ++k) {
                Grid_Index ind1 = ind;
                ind1.id[k] += j;
                auto it = vcolor.find(ind1);
                if (it == vcolor.end()) {
                    flag = 0;
                    break;
                }
                int y = it->second;
                unit_len = glm::length(nvertices[y] - nvertices[i]);
                std::pair<int, int> edge_id;
                if (i < y) {
                    edge_id = std::make_pair(i, y);
                } else {
                    edge_id = std::make_pair(y, i);
                }
                if (edge_info.find(edge_id) == edge_info.end()) {
                    flag = 0;
                    break;
                }
            }
        }
        if (flag < 3) {
            continue;
        }
        evens.insert(i);
        ecs.insert(std::make_pair(i, std::set<std::pair<int, int>>()));
    }
    for (int i = 0; i < (int)triangles.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            int x = triangles[i][j];
            if (odds.find(x) != odds.end()) {
                ocs[x].insert(std::make_pair(i, j));
            }
            if (evens.find(x) != evens.end()) {
                ecs[x].insert(std::make_pair(i, j));
            }
        }
    }
    for (int i : evens) {
        glm::dvec3 dir;
        int count = 0;
        for (int j = 0; j < 8; ++j) {
            glm::dvec3 d((j & 0x04) > 0, (j & 0x02) > 0, (j & 0x01) > 0);
            d = d * 2.0 - glm::dvec3(1, 1, 1);
            d = glm::normalize(d) * (unit_len * 0.5);
            if (!tree_->Is_Exterior(nvertices[i] + d)) {
                dir = glm::normalize(d);
                count += 1;
            }
        }
        if (count > 2) {
            continue;
        }
        std::set<std::pair<int, int>>& p = ecs[i];
        for (const auto& it1 : p) {
            assert(triangles[it1->first][it1->second] == i);
            if (glm::dot(nvertices[triangles[it1.first][(it1.second + 1) % 3]] -
                                 nvertices[i],
                         dir) < 0) {
                triangles[it1.first][it1.second] = (int)nvertices.size();
            }
        }
        nvertices[i] += dir * (0.5 * unit_len);
        v_faces_.push_back(v_faces_[i]);
        nvertices.push_back(nvertices[i]);
        nvertices.back() -= unit_len * dir;
    }
    for (int i : odds) {
        int k = 0;
        while (v_info_[i].id[k] % 2 == 0) {
            k += 1;
        }
        Grid_Index id1;
        Grid_Index id2;
        id1 = v_info_[i];
        id2 = v_info_[i];
        id1.id[k] -= 1;
        id2.id[k] += 1;
        int x = vcolor[id1];
        int y = vcolor[id2];
        if (x > y) {
            int temp = x;
            x = y;
            y = temp;
        }
        if (edge_info[std::make_pair(x, y)].size() > 2) {
            glm::dvec3 vert = nvertices[x] - nvertices[y];
            double len = glm::length(vert);
            vert /= len;
            glm::dvec3 dir(len * 0.5, len * 0.5, len * 0.5);
            dir = dir - glm::dot(dir, vert) * vert;
            if (!tree_->Is_Exterior(nvertices[i] + dir)) {
                dir = glm::cross(vert, dir);
            }
            dir = glm::normalize(dir);
            std::set<std::pair<int, int>>& p = ocs[i];
            for (const auto& it1 : p) {
                assert(triangles[it1->first][it1->second] == i);
                if (glm::dot(nvertices[triangles[it1.first]
                                                [(it1.second + 1) % 3]] -
                                     nvertices[i],
                             dir) < 0) {
                    triangles[it1.first][it1.second] = (int)nvertices.size();
                }
            }
            nvertices[i] += dir * (0.5 * len);
            v_faces_.push_back(v_faces_[i]);
            nvertices.push_back(nvertices[i]);
            nvertices.back() -= len * dir;
        }
    }
    return true;
}

bool ManifoldProcessor::is_manifold(
        const std::vector<Eigen::Vector3d>& overtices,
        const std::vector<Eigen::Vector3i>& otriangles) {
    std::vector<glm::dvec3> vertices;
    std::vector<glm::ivec3> triangles;
    vertices.reserve(overtices.size());
    for (const auto& v : overtices) {
        vertices.emplace_back(glm::dvec3(v[0], v[1], v[2]));
    }
    triangles.reserve(otriangles.size());
    for (const auto& f : otriangles) {
        triangles.emplace_back(glm::ivec3(f[0], f[1], f[2]));
    }

    std::map<std::pair<int, int>, std::list<glm::dvec3>> edges;
    std::vector<std::set<int>> graph(vertices.size());
    int flag = 0;
    for (auto& face_indice : triangles) {
        for (int j = 0; j < 3; ++j) {
            int x = face_indice[j];
            int y = face_indice[(j + 1) % 3];
            if (x > y) {
                int temp = x;
                x = y;
                y = temp;
            }
            if (x >= (int)vertices.size() || y >= (int)vertices.size()) {
                flag = -1;
            }
            std::pair<int, int> edge_id = std::make_pair(x, y);
            auto it = edges.find(edge_id);
            if (it == edges.end()) {
                std::list<glm::dvec3> l;
                l.push_back(vertices[face_indice[(j + 2) % 3]]);
                edges.insert(std::make_pair(edge_id, l));
            } else {
                if (it->second.size() == 2) {
                    flag = -1;
                } else {
                    it->second.push_back(vertices[face_indice[(j + 2) % 3]]);
                    glm::dvec3 p1 = glm::cross(it->second.front() - vertices[x],
                                               vertices[y] - vertices[x]);
                    glm::dvec3 p2 = glm::cross(it->second.back() - vertices[x],
                                               vertices[y] - vertices[x]);
                    glm::dvec3 norm1 = glm::normalize(p1);
                    glm::dvec3 norm2 = glm::normalize(p2);
                    if (glm::dot(norm1, norm2) > 1 - 1e-7 &&
                        glm::length(p1) > 1e-4 && glm::length(p2) > 1e-4) {
                        flag = 0;
                    }
                }
            }
        }
    }
    for (auto& edge : edges) {
        if (edge.second.size() == 1) {
            flag = 1;
        }
    }
    return flag == 0;
}

std::tuple<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3i>>
ManifoldProcessor::GetManifoldMesh(int depth) {
    vertices_buf_ = vertices_;
    face_indices_buf_ = face_indices_;
    Build_Tree(depth);
    Construct_Manifold();
    Project_Manifold();

    for (int iter = 0; iter < 1; ++iter) {
        std::vector<glm::dvec3> dis(vertices_.size());
        std::vector<int> dis_weight(vertices_.size());
        for (auto& face_indice : face_indices_) {
            for (int j = 0; j < 3; ++j) {
                int x = face_indice[j];
                int y = face_indice[(j + 1) % 3];
                dis[x] += vertices_[y];
                dis[y] += vertices_[x];
                dis_weight[x] += 1;
                dis_weight[y] += 1;
            }
        }
        for (int i = 0; i < (int)vertices_.size(); ++i) {
            if (dis_weight[i] > 0) {
                vertices_[i] = dis[i] * (1.0 / dis_weight[i]);
            }
        }
    }

    // Convert to std::vector
    std::vector<Eigen::Vector3d> vertices;
    vertices.reserve(vertices_.size());
    for (const auto& v : vertices_) {
        vertices.emplace_back(Eigen::Vector3d(v[0], v[1], v[2]));
    }
    std::vector<Eigen::Vector3i> triangles;
    triangles.reserve(face_indices_.size());
    for (const auto& f : face_indices_) {
        triangles.emplace_back(Eigen::Vector3i(f[0], f[1], f[2]));
    }
    return std::make_tuple(vertices, triangles);
}
