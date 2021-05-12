#ifndef OCTREE_H_
#define OCTREE_H_

#include <list>
#include <map>
#include <set>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "manifold/Intersection.h"

class Grid_Index {
public:
    Grid_Index() = default;
    Grid_Index(int x, int y, int z) : id(x, y, z) {}
    bool operator<(const Grid_Index& ind) const {
        int i = 0;
        while (i < 3 && id[i] == ind.id[i]) {
            i++;
        }
        return (i < 3 && id[i] < ind.id[i]);
    }
    Grid_Index operator+(const Grid_Index& ind) const {
        Grid_Index grid(*this);
        grid.id += ind.id;
        return grid;
    }
    Grid_Index operator/(int x) const {
        return Grid_Index(id[0] / x, id[1] / x, id[2] / x);
    }
    glm::ivec3 id;
};

class Octree {
public:
    Octree(glm::dvec3& min_c,
           glm::dvec3& max_c,
           std::vector<glm::ivec3>& faces,
           float thickness) {
        min_corner = min_c;
        length = max_c - min_corner;
        int ind = 0;
        for (int i = 1; i < 3; ++i) {
            if (length[i] > length[ind]) {
                ind = i;
            }
        }
        for (int i = 0; i < 3; ++i) {
            min_corner[i] -= (length[ind] - length[i]) * 0.5 + thickness * 0.5;
        }
        length = glm::dvec3(length[ind] + thickness, length[ind] + thickness,
                            length[ind] + thickness);

        face_indices = faces;
        face_ind.resize(faces.size());
        for (int i = 0; i < (int)faces.size(); ++i) {
            face_ind[i] = i;
        }
    }

    Octree(glm::dvec3& min_c, glm::dvec3& length_) {
        min_corner = min_c;
        length = length_;
    }

    ~Octree() {
        for (int i = 0; i < 8; ++i) {
            if (children[i] != nullptr) {
                delete children[i];
            }
        }
    }

    bool Is_Exterior(const glm::dvec3& p) {
        for (int i = 0; i < 3; ++i) {
            if (p[i] < min_corner[i] || p[i] > min_corner[i] + length[i]) {
                return true;
            }
        }
        if (occupied == 0) {
            return exterior != 0;
        }
        if (level == 0) {
            return false;
        }
        int index = 0;
        for (int i = 0; i < 3; ++i) {
            index *= 2;
            if (p[i] > min_corner[i] + length[i] / 2) {
                index += 1;
            }
        }
        return children[index]->Is_Exterior(p);
    }

    bool Intersection(int face_index,
                      glm::dvec3& min_corner,
                      glm::dvec3& length,
                      std::vector<glm::dvec3>& vertices) {
        float boxcenter[3];
        float boxhalfsize[3];
        float triverts[3][3];
        for (int i = 0; i < 3; ++i) {
            boxhalfsize[i] = length[i] * 0.5;
            boxcenter[i] = min_corner[i] + boxhalfsize[i];
        }
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                triverts[i][j] = vertices[face_indices[face_index][i]][j];
            }
        }
        return triBoxOverlap(boxcenter, boxhalfsize, triverts) != 0;
    }

    void Split(std::vector<glm::dvec3>& vertices) {
        level += 1;
        number = 0;
        if (level > 1) {
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    for (int k = 0; k < 2; ++k) {
                        int ind = i * 4 + j * 2 + k;
                        if ((children[ind] != nullptr) &&
                            (children[ind]->occupied != 0)) {
                            children[ind]->Split(vertices);
                            number += children[ind]->number;
                        }
                    }
                }
            }
            face_indices.clear();
            face_ind.clear();
            return;
        }
        glm::dvec3 halfsize = length * glm::dvec3(0.5, 0.5, 0.5);
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 2; ++k) {
                    int ind = i * 4 + j * 2 + k;

                    glm::dvec3 startpoint = min_corner;
                    startpoint[0] += i * halfsize[0];
                    startpoint[1] += j * halfsize[1];
                    startpoint[2] += k * halfsize[2];

                    children[ind] = new Octree(startpoint, halfsize);
                    children[ind]->occupied = 0;
                    children[ind]->number = 0;

                    for (int face = 0; face < (int)face_indices.size();
                         ++face) {
                        if (Intersection(face, startpoint, halfsize,
                                         vertices)) {
                            children[ind]->face_indices.push_back(
                                    face_indices[face]);
                            children[ind]->face_ind.push_back(face_ind[face]);
                            if (children[ind]->occupied == 0) {
                                children[ind]->occupied = 1;
                                number += 1;
                                children[ind]->number = 1;
                            }
                        }
                    }
                }
            }
        }
        face_indices.clear();
        face_ind.clear();
    }

    void BuildConnection() {
        if (level == 0) {
            return;
        }
        for (int i = 0; i < 8; ++i) {
            if (children[i] != nullptr) {
                children[i]->BuildConnection();
            }
        }
        int y_index[] = {0, 1, 4, 5};
        for (int i = 0; i < 4; ++i) {
            if ((children[i * 2] != nullptr) &&
                (children[i * 2 + 1] != nullptr)) {
                ConnectTree(children[i * 2], children[i * 2 + 1], 2);
            }
            if ((children[y_index[i]] != nullptr) &&
                (children[y_index[i] + 2] != nullptr)) {
                ConnectTree(children[y_index[i]], children[y_index[i] + 2], 1);
            }
            if ((children[i] != nullptr) && (children[i + 4] != nullptr)) {
                ConnectTree(children[i], children[i + 4], 0);
            }
        }
    }

    void ConnectTree(Octree* l, Octree* r, int dim) {
        int y_index[] = {0, 1, 4, 5};
        if (dim == 2) {
            l->connection[2] = r;
            r->connection[5] = l;
            for (int i = 0; i < 4; ++i) {
                if ((l->children[i * 2 + 1] != nullptr) &&
                    (r->children[i * 2] != nullptr)) {
                    ConnectTree(l->children[i * 2 + 1], r->children[i * 2],
                                dim);
                }
            }
        } else if (dim == 1) {
            l->connection[1] = r;
            r->connection[4] = l;
            for (int i = 0; i < 4; ++i) {
                if ((l->children[y_index[i] + 2] != nullptr) &&
                    (r->children[y_index[i]] != nullptr)) {
                    ConnectTree(l->children[y_index[i] + 2],
                                r->children[y_index[i]], dim);
                }
            }
        } else if (dim == 0) {
            l->connection[0] = r;
            r->connection[3] = l;
            for (int i = 0; i < 4; ++i) {
                if ((l->children[i + 4] != nullptr) &&
                    (r->children[i] != nullptr)) {
                    ConnectTree(l->children[i + 4], r->children[i], dim);
                }
            }
        }
    }

    void ConnectEmptyTree(Octree* l, Octree* r, int dim) {
        int y_index[] = {0, 1, 4, 5};
        if ((l->occupied != 0) && (r->occupied != 0)) {
            if (l->level == 0) {
                return;
            }
            if (dim == 2) {
                for (int i = 0; i < 4; ++i) {
                    ConnectEmptyTree(l->children[i * 2 + 1], r->children[i * 2],
                                     dim);
                }
            } else if (dim == 1) {
                for (int i = 0; i < 4; ++i) {
                    ConnectEmptyTree(l->children[y_index[i] + 2],
                                     r->children[y_index[i]], dim);
                }
            } else if (dim == 0) {
                for (int i = 0; i < 4; ++i) {
                    ConnectEmptyTree(l->children[i + 4], r->children[i], dim);
                }
            }
            return;
        }
        if (!((l->occupied != 0) || (r->occupied != 0))) {
            l->empty_neighbors.push_back(r);
            r->empty_neighbors.push_back(l);
            return;
        }
        if (l->occupied == 0) {
            if (dim == 2) {
                r->empty_connection[5] = l;
                if (r->level > 0) {
                    for (int i = 0; i < 4; ++i) {
                        ConnectEmptyTree(l, r->children[i * 2], dim);
                    }
                }
            } else if (dim == 1) {
                r->empty_connection[4] = l;
                if (r->level > 0) {
                    for (int i = 0; i < 4; ++i) {
                        ConnectEmptyTree(l, r->children[y_index[i]], dim);
                    }
                }
            } else if (dim == 0) {
                r->empty_connection[3] = l;
                if (r->level > 0) {
                    for (int i = 0; i < 4; ++i) {
                        ConnectEmptyTree(l, r->children[i], dim);
                    }
                }
            }
            return;
        }
        if (r->occupied == 0) {
            if (dim == 2) {
                l->empty_connection[2] = r;
                if (l->level > 0) {
                    for (int i = 0; i < 4; ++i) {
                        ConnectEmptyTree(l->children[i * 2 + 1], r, dim);
                    }
                }
            } else if (dim == 1) {
                l->empty_connection[1] = r;
                if (l->level > 0) {
                    for (int i = 0; i < 4; ++i) {
                        ConnectEmptyTree(l->children[y_index[i] + 2], r, dim);
                    }
                }
            } else if (dim == 0) {
                l->empty_connection[0] = r;
                if (l->level > 0) {
                    for (int i = 0; i < 4; ++i) {
                        ConnectEmptyTree(l->children[i + 4], r, dim);
                    }
                }
            }
        }
    }

    void ExpandEmpty(std::list<Octree*>& empty_list,
                     std::set<Octree*>& empty_set,
                     int dim) {
        if (occupied == 0) {
            if (empty_set.find(this) == empty_set.end()) {
                empty_set.insert(this);
                empty_list.push_back(this);
            }
            return;
        }
        if (level == 0) {
            return;
        }
        int y_index[] = {0, 1, 4, 5};
        if (dim == 2 || dim == 5) {
            for (int i = 0; i < 4; ++i) {
                children[i * 2 + static_cast<int>(dim == 5)]->ExpandEmpty(
                        empty_list, empty_set, dim);
            }
            return;
        }
        if (dim == 1 || dim == 4) {
            for (int i = 0; i < 4; ++i) {
                children[y_index[i] + 2 * static_cast<int>(dim == 4)]
                        ->ExpandEmpty(empty_list, empty_set, dim);
            }
            return;
        }
        for (int i = 0; i < 4; ++i) {
            children[i + 4 * static_cast<int>(dim == 3)]->ExpandEmpty(
                    empty_list, empty_set, dim);
        }
    }

    void BuildEmptyConnection() {
        if (level == 0) {
            return;
        }

        for (int i = 0; i < 8; ++i) {
            if (children[i]->occupied != 0) {
                children[i]->BuildEmptyConnection();
            }
        }
        int pair_x[] = {0, 2, 4, 6, 0, 1, 4, 5, 0, 1, 2, 3};
        int pair_y[] = {1, 3, 5, 7, 2, 3, 6, 7, 4, 5, 6, 7};
        int dim[] = {2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0};
        for (int i = 0; i < 12; ++i) {
            ConnectEmptyTree(children[pair_x[i]], children[pair_y[i]], dim[i]);
        }
    }

    void ConstructFace(std::map<Grid_Index, int>& vcolor,
                       const glm::ivec3& start,
                       std::vector<glm::dvec3>& vertices,
                       std::vector<glm::ivec4>& faces,
                       std::vector<std::set<int>>& v_faces) {
        if (level == 0) {
            if (occupied == 0) {
                return;
            }
            glm::ivec3 offset[6][4] = {
                    {glm::ivec3(1, 0, 0), glm::ivec3(1, 0, 1),
                     glm::ivec3(1, 1, 1), glm::ivec3(1, 1, 0)},
                    {glm::ivec3(0, 1, 0), glm::ivec3(1, 1, 0),
                     glm::ivec3(1, 1, 1), glm::ivec3(0, 1, 1)},
                    {glm::ivec3(0, 0, 1), glm::ivec3(0, 1, 1),
                     glm::ivec3(1, 1, 1), glm::ivec3(1, 0, 1)},
                    {glm::ivec3(0, 0, 0), glm::ivec3(0, 1, 0),
                     glm::ivec3(0, 1, 1), glm::ivec3(0, 0, 1)},
                    {glm::ivec3(0, 0, 0), glm::ivec3(0, 0, 1),
                     glm::ivec3(1, 0, 1), glm::ivec3(1, 0, 0)},
                    {glm::ivec3(0, 0, 0), glm::ivec3(1, 0, 0),
                     glm::ivec3(1, 1, 0), glm::ivec3(0, 1, 0)}};
            for (int i = 0; i < 6; ++i) {
                if ((empty_connection[i] != nullptr) &&
                    (empty_connection[i]->exterior != 0)) {
                    if ((connection[i] != nullptr) &&
                        (connection[i]->occupied != 0)) {
                        std::cerr << "Error!\n";
                        exit(0);
                    }
                    int id[4];
                    for (int j = 0; j < 4; ++j) {
                        glm::ivec3 vind = start + offset[i][j];
                        Grid_Index v_id;
                        v_id.id = vind * 2;
                        auto it = vcolor.find(v_id);
                        if (it == vcolor.end()) {
                            glm::dvec3 d = min_corner;
                            for (int k = 0; k < 3; ++k) {
                                d[k] += offset[i][j][k] * length[k];
                            }
                            vcolor.insert(
                                    std::make_pair(v_id, vertices.size()));
                            id[j] = vertices.size();
                            vertices.push_back(d);
                            v_faces.emplace_back();
                            for (auto it1 = face_ind.begin();
                                 it1 != face_ind.end(); ++it1) {
                                v_faces[id[j]].insert(*it1);
                            }
                        } else {
                            id[j] = it->second;
                            for (auto it1 = face_ind.begin();
                                 it1 != face_ind.end(); ++it1) {
                                v_faces[it->second].insert(*it1);
                            }
                        }
                    }
                    faces.emplace_back(id[0], id[1], id[2], id[3]);
                }
            }
        } else {
            for (int i = 0; i < 8; ++i) {
                if ((children[i] != nullptr) && (children[i]->occupied != 0)) {
                    int x = i / 4;
                    int y = (i - x * 4) / 2;
                    int z = i - x * 4 - y * 2;
                    glm::ivec3 nstart = start * 2 + glm::ivec3(x, y, z);
                    children[i]->ConstructFace(vcolor, nstart, vertices, faces,
                                               v_faces);
                }
            }
        }
    }

    glm::dvec3 min_corner, length;
    int level{0};
    int number{1};
    int occupied{1};
    int exterior{0};

    std::vector<Octree*> children = std::vector<Octree*>(8);
    std::vector<Octree*> connection = std::vector<Octree*>(6);
    std::vector<Octree*> empty_connection = std::vector<Octree*>(6);
    std::list<Octree*> empty_neighbors;
    std::vector<glm::ivec3> face_indices;
    std::vector<int> face_ind;
};

#endif
