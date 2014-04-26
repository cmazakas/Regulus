/*

    It is of great aid in understanding this file to first explain how the 
    tetrahedral structure itself works.

    It contains the constructor which initializes the structure. Regulus is
    very similar to Arepo in the sense that we use a point and neighbour
    structure. We contain 4 pointers to points and 4 pointers to neighbouring
    tetrahedra. We then include storage for each face and the index of the
    point opposite the shared face, from the perspective of the opposite
    tetrahedron.

    struct tetra {

        array p = vertices (type point*)
        array ngb = neighbours (type tetra*)
        array nid = index of points opposite face (type unsigned int)
        2d array f = faces (type point**)
    }
    
    We use a specific ordering scheme to simplify various tasks. Namely, we use
    
    0 - 321
    1 - 023
    2 - 031
    3 - 012

    We use this to signify that ngb[0] (neighbour 0) is the neighbour through
    the face 321 of the tetrahedron. ngb[1] is through face 023, etc.
    nid[0] is the index of the point opposite face 321 contained by t->ngb[0].
    f[0] is direct storage of the addresses of the points 3, 2 and 1, in order.
    p[0] is naturally point 0.

*/

#include "structures.hpp"

unsigned int fast_inclusion(const tetra *t, const point *p);
unsigned int exact_inclusion(const tetra *t, const point *p);
tetra* step(const tetra *t, const point *p, std::vector<tetra*> &branches);
unsigned int fast_triangle_ray_inclusion(const std::array<point*, 3> &f, const point *a, const point *b);
void edgewalk(std::vector<std::pair<tetra*, unsigned int>> &leaves, const point *p);

/*

    Constructor

*/

tetra::tetra(point *a, point *b, point *c, point *d) {

    p[0] = a;
    p[1] = b;
    p[2] = c;
    p[3] = d;

    ngb.fill(nullptr);
    nid.fill(-1);

    f = { { { p[3], p[2], p[1] }, // 0 - 321
            { p[0], p[2], p[3] }, // 1 - 023
            { p[0], p[3], p[1] }, // 2 - 031
            { p[0], p[1], p[2] }  // 3 - 012
          } };

    #ifndef NDEBUG

        if (a) {

            Eigen::Matrix<double, 4, 4> ort;

            ort << 1, a->x, a->y, a->z,
                   1, b->x, b->y, b->z,
                   1, c->x, c->y, c->z,
                   1, d->x, d->y, d->z;

            assert(ort.determinant() > 0); // assert positive orientation
        }

    #endif

    return;
}

/*

    This routine performs the necessary fractures when a point is inserted into
    the current triangulation. Each fracture is positively oriented and then
    linked with each other, including ancestral neighbours from before the 
    fracture.

    To determine which faces are used in the fracture, we examine the first 4
    bits of the unsigned integer passed along with the tetrahedral address.
    If the bit is set, we construct a new tetrahedron in the buffer otherwise,
    we simply do nothing and continue iterating.

*/

void fracture(mesh &m, std::vector<std::pair<tetra*, unsigned int>> &leaves, std::vector<tetra*> &pile, point *p) {

    assert(pile.empty());

    if (verbose >= 1) {

        std::cout << "\nFracturing :\n" << std::endl;
        for (auto &l : leaves)
            std::cout << *(l.first) << std::endl;
    }

    std::unordered_set<tetra*> ancestors;
    ancestors.reserve(4 * leaves.size());

    /*
        Collect ancestors for future linking
    */

    for (auto it = leaves.begin(); it < leaves.end(); ++it) {

        for (int i = 0; i < 4; ++i)
            if ((*it).first->ngb[i])
                ancestors.emplace((*it).first->ngb[i]);
    }

    for (auto &leaf : leaves)
        ancestors.erase(leaf.first);

    /*
        Begin fracutring procedure...
    */

    for (auto it = leaves.begin(); it < leaves.end(); ++it) {

        tetra *t = (*it).first;
        unsigned int pos = (*it).second;

        tetra cpy = *t;

        unsigned int i = 0;

        for ( ; i < 4; ++i) {
            if (pos & (1 << i)) { // if bit is set...

                pile.emplace_back(new(t) tetra(cpy.f[i][0], cpy.f[i][1], cpy.f[i][2], p));
                ++i;
                break;
            }
        }

        for ( ; i < 4; ++i) {
            if (pos & (1 << i)) { // if bit is set...

                pile.emplace_back(new(m.t_position) tetra(cpy.f[i][0], cpy.f[i][1], cpy.f[i][2], p)); 

                ++m.t_position;
                ++m.num_tetra;
            }
        }
    }

    /*
        Link neighbouring tetrahedra...
    */

    for (auto it1 = pile.begin(); it1 < pile.end(); ++it1) {

        for (auto it2 = it1 + 1; it2 < pile.end(); ++it2) {
            link_tetra(*it1, *it2);
        }

        for (auto it2 = ancestors.begin(); it2 != ancestors.end(); ++it2) {
            link_tetra(*it1, *it2);
        }
    }

    if (verbose >= 1) {

        std::cout << "\nResuling fracture set : \n" << std::endl;

        for (auto &pl : pile)
            std::cout << *pl << std::endl;
    }

    m.m_position = m.t_position - 1;

    return;
}


/*

    We use a vector to collect relevant tetrahedra, a point to 
    search for and a starting point for our walking procedure

    Due to complexities that are entirely beyond the author, the walking
    procedure was modified to account for instances of "branching" in the walk.
    By branching, we mean that Regulus navigates by drawing a line between the
    centroid of a tetrahedron and the point to be inserted. We then examine
    which faces of the tetrahedron are intersected by this line. Sometimes this
    line intersects the edge or vertex of a tetrahedron so there are "breaks" in
    the path as an edge is the intersection of two faces and a vertex is the
    intersection of three. We store these alternate instances.

    And if that should fail (sometimes it does) we have a backup plan which is
    just starting at a random tetrahedron in the buffer and beginning the walk
    from there. Note that this is horribly ineffecient but a long walk that 
    terminates is far more efficient against walking in cirlces for all time.
    
*/

void walk(std::vector<std::pair<tetra*, unsigned>> &leaves, const point *p, tetra *place, tetra *backup) {

    assert(leaves.empty() && p && place);
    auto t = place;

    std::unordered_set<tetra*> history;
    history.reserve(2048); 

    std::vector<tetra*> branches;
    branches.reserve(512);

    while (t) {

        if (verbose >= 1)
            std::cout << "Walking... (" << t << ")" << std::endl;

        for (int i = 0; i < 4; ++i) {
            if (t->ngb[i]) {
                if (fast_inclusion(t->ngb[i], p) != 0) {

                    t = t->ngb[i];
                    break;
                }
            }
        }

        if (!history.insert(t).second) {

            if (branches.size() == 0) {

                while (!history.insert(backup).second || !backup->p[0])
                    ++backup;

                if (!backup)
                    throw std::runtime_error("Back-up plan failed horribly!");

                t = backup;

                continue;

            } else {

                t = branches.back();
                branches.pop_back();
                continue;

            }
        }

        unsigned int pos = fast_inclusion(t, p);

        assert(exact_inclusion(t, p) == pos);

        if (pos == 0) {

            t = step(t, p, branches);
            continue;
        }

        leaves.emplace_back(std::pair<tetra*, unsigned int>(t, pos));

        switch (pos) {

            case 15 :
                return;

            case 7 : // 4 + 2 + 1
            case 11 : // 8 + 2 + 1
            case 13 : // 8 + 4 + 1
            case 14 : // 8 + 4 + 2

                for (unsigned int n = 0; n < 4; ++n) {

                    if (!(pos & (1<< n))) // bit is clear...
                        if (t->ngb[n])
                            leaves.emplace_back(std::pair<tetra*, unsigned int>(t->ngb[n], fast_inclusion(t->ngb[n], p)));
                        
                }
                                
                return;

            case 3 : // 2 + 1
            case 5 : // 4 + 1
            case 6 : // 4 + 2
            case 9 : // 8 + 1
            case 10 : // 8 + 2
            case 12 : // 8 + 4

                edgewalk(leaves, p);
                return;

            case 1 :
            case 2 :
            case 4 :
            case 8 :

                std::cout << "duplicate vertex. not inserting..." << std::endl;
                return;

            default :
                throw std::runtime_error("Inclusion test errors!");
        }
    }

    return;
}

/*

    fast_inclusion takes a tetrahedron and a point p and will return an
    unsigned integer that has all necessary information encoded in it.
    exact_inclusion uses this exact method for its return as well.

    A tetrahedron has 4 faces so we use this to our advantage. We let the
    first 4 bits of an integer represent which faces the point intersects.

    A set bit (1) means that the point is above the face and a clear bit (0)
    means that the point is on the face.

    fast_inclusion uses 3 dimensional barycentric coordinates to evaluate
    point location.

*/

unsigned int fast_inclusion(const tetra *t, const point *p) {

    auto &a = *t->p[0],
         &b = *t->p[1],
         &c = *t->p[2],
         &d = *t->p[3];

    Eigen::Matrix<double, 3, 3> A;

    A.col(0) << b.x - a.x, b.y - a.y, b.z - a.z;
    A.col(1) << c.x - a.x, c.y - a.y, c.z - a.z;
    A.col(2) << d.x - a.x, d.y - a.y, d.z - a.z;

    if (verbose == 3)
        std::cout << A << std::endl;

    Eigen::Matrix<double, 3, 1> x, B(p->x - a.x, p->y - a.y, p->z - a.z);

    x = A.inverse() * B;

    double sum = 0;

    for (unsigned int i = 0; i < 3; ++i) {

        if (std::abs(x[i]) < 1e-10)
            x[i] = 0;

        if (x[i] < 0)
            return 0; // outside
        else
            sum += x[i];
    }


    if (std::abs(sum - 1) < 1e-10)
        sum = 1;//return exact_inclusion(t, p);

    if (std::abs(sum) < 1e-10)
        sum = 0;

    if (sum > 1)
        return 0; // outside

    if (sum == 0)
        return 1; // vertex 0

    double u(x[0]), v(x[1]), w(x[2]);

    if (u == 1) {

        return 2; // vertex 1
    }

    else if (u > 0) {

        if (v > 0) {
            
            if (w > 0) {

                if (sum == 1)
                    return 14; // surface 321
                else
                    return 15; // inside
            }

            else {

                if (sum == 1)
                    return 6; // edge 21
                else
                    return 7; // surface 012
            }
        }

        else {
            
            if (w > 0) {

                if (sum == 1)
                    return 10; // edge 31
                else
                    return 11; // surface 031
            }

            else {

                return 3; // edge 10
            }
        } 
    } else {

        if (v == 1)
            return 4; // vertex 2

        else if (v > 0) {

            if (w > 0) {

                if (sum == 1)
                    return 12; // edge 32
                else
                    return 13; // surface 023
            }

            else {

                return 5; // edge 20
            }
        }

        else {

            if (w == 1)
                return 8; // vertex 3
            else
                return 9; // edge 30
        }        
    }

    return 0;
}

/*

    exact_inclusion is, in theory, exact because it uses determinants to 
    evaluate whether or not a point is one a face. If the point is, the
    orientation of the proposed tetrahedron is 0 otherwise it is positive.

    We use the same bit encoding as fast_inclusion and the results agree

*/

unsigned int exact_inclusion(const tetra *t, const point *p) {

    Eigen::Matrix<double, 4, 4> A;

    std::array<int, 4> ort;

    for (int i = 0; i < 4; ++i) {

        point &a = *t->f[i][0];
        point &b = *t->f[i][1];
        point &c = *t->f[i][2];

        A << 1, a.x, a.y, a.z,
             1, b.x, b.y, b.z,
             1, c.x, c.y, c.z,
             1, p->x, p->y, p->z;

        auto tmp = A.determinant();

        if (tmp < 0)
            return 0;
        else if (tmp == 0)
            ort[i] = 0;
        else
            ort[i] = 1;
    }    

    return ort[0] + 2 * ort[1] + 4 * ort[2] + 8 * ort[3];
}

/*

    This is our step function and is used to determine which direction we are
    to head during our walk. We calculate the centroid of the tetrahedron and
    then determine whether or not the line between the center and the point to
    be inserted intersects a face.

*/

tetra* step(const tetra *t, const  point *p, std::vector<tetra*> &branches) {

    tetra *step = nullptr;

    point centroid((t->p[0]->x + t->p[1]->x + t->p[2]->x + t->p[3]->x)/4,
                   (t->p[0]->y + t->p[1]->y + t->p[2]->y + t->p[3]->y)/4,
                   (t->p[0]->z + t->p[1]->z + t->p[2]->z + t->p[3]->z)/4
                  );

    for (unsigned int i = 0; i < 4; ++i) {

        if (!t->ngb[i])
            continue;

        unsigned int included = fast_triangle_ray_inclusion(t->f[i], &centroid, p);

        if (included != 0) {

            if (!step)
                step = t->ngb[i];
            else
                branches.push_back(t->ngb[i]);
        }
    }

    return step;
}

/*

    This function takes two points (a and b) and forms a line. It then
    determines if the line intersects the triangle spanned by f. We
    continue to use barycentric coordinates.

*/

unsigned int fast_triangle_ray_inclusion(const std::array<point*, 3> &f, const point *a, const point *b) {

    Eigen::Matrix<double, 3, 3> A;

    A.col(0) << a->x - b->x,
                a->y - b->y,
                a->z - b->z;

    A.col(1) << f[1]->x - f[0]->x,
                f[1]->y - f[0]->y,
                f[1]->z - f[0]->z;

    A.col(2) << f[2]->x - f[0]->x,
                f[2]->y - f[0]->y,
                f[2]->z - f[0]->z;

    Eigen::Matrix<double, 3, 1> B;

    B << a->x - f[0]->x, 
         a->y - f[0]->y, 
         a->z - f[0]->z;

    Eigen::Matrix<double, 3, 3> inverse;
    bool invertible = false;
    
    A.computeInverseWithCheck(inverse, invertible);

    if (!invertible) {

        //if (verbose >= 0)
            //std::cout << "Uninvertible matrix. Skipping..." << std::endl;
        return 0;
    }

    Eigen::Matrix<double, 3, 1> x = inverse * B;

    for (int i = 0; i < 3; ++i)
        if (std::abs(x[i]) < 1e-10)
            x[i] = 0;

    double t(x[0]), u(x[1]), v(x[2]);

    if (verbose >= 3) {

        std::cout << "\nplane being examined : " << std::endl;
        for (int i = 0; i < 3; ++i)
            std::cout << *f[i] << std::endl;

        std::cout << "\ntwo points : " << std::endl;
        std::cout << *a << " <------> " << *b << std::endl;

        std::cout << "A : " << std::endl << A << std::endl;
        std::cout << "B : " << std::endl << B << std::endl;

        std::cout << "x : " << std::endl << x << std::endl;
    }

    assert(t == t && u == u && v == v); // NaN checks

    assert(t != 0);

    if (t == 1) {

        assert(u < 0 || v < 0 || u + v > 1);
        return 0; // outside triangle's interior
    }

    if (t > 0 && t < 1) {

        double sum = 0;

        for (int i = 1; i < 3; ++i) {

            if (x[i] < 0)
                return 0; // outside

            sum += x[i];
        }

        if (std::abs(sum - 1) < 1e-10)
            sum = 1;

        if (sum > 1)
            return 0; // outside

        if (u == 1) {

            return 2; // vertex 1
        }

        else if (u > 0) {

            if (v > 0) {

                if (sum == 1) {

                    return 6; // edge 21
                } else {

                    return 7; // inside
                }
            }

            else {

                return 3; // edge 10
            }
        }

        else {

            if (v == 1)
                return 4; // vertex 2

            else if (v > 0) {

                return 5; // edge 02

            }else
                return 1; // vertex 0
        }
    }

    return 0;
}

/*
    
    This is a much more precise function in the sense that fast_triangle_ray_inclusion
    doesn't determine which half-space the intersection point is in. This 
    function does exactly that.
    
*/

unsigned int halfspace_locator(const std::array<point*, 3> &f, const point *a, const point *b) {

    Eigen::Matrix<double, 3, 3> A;

    A.col(0) << a->x - b->x,
                a->y - b->y,
                a->z - b->z;

    A.col(1) << f[1]->x - f[0]->x,
                f[1]->y - f[0]->y,
                f[1]->z - f[0]->z;

    A.col(2) << f[2]->x - f[0]->x,
                f[2]->y - f[0]->y,
                f[2]->z - f[0]->z;


    bool invertible = false;
    Eigen::Matrix<double, 3, 3> inverse;

    A.computeInverseWithCheck(inverse, invertible);    

    if (!invertible)
        return -1;

    Eigen::Matrix<double, 3, 1> B;

    B << a->x - f[0]->x,
         a->y - f[0]->y,
         a->z - f[0]->z;

    Eigen::Matrix<double, 3, 1> x = inverse * B;

    double t(x[0]), u(x[1]), v(x[2]);

    assert(t == t && u == u && v == v);

    if (std::abs(u) < 1e-10)
        u = 0;

    if (std::abs(v) < 1e-10)
        v = 0;

    double sum = u + v;

    if (std::abs(1 - sum) < 1e-10)
        sum = 1;

    unsigned int location = 0;

    /*
        half-space 21
    */

    if (sum > 1) {

        location |= 1 << 0;
        location |= 1 << 1;
    }

    else if (sum < 1) {

        location |= 1 << 0;
    }

    /*
        half-space 02
    */

    if (u > 0) {

            location |= 1 << 2;
    }

    else if (u < 0) {

        location |= 1 << 2;
        location |= 1 << 3;
    }

    /*
        half-space 10
    */

    if (v > 0) {

        location |= 1 << 4;
    }

    else if (v < 0) {

        location |= 1 << 4;
        location |= 1 << 5;
    }

    return location;
}

/*

    When a point is on an edge, an arbitrary number of tetrahedra can share
    that edge. It is for this reason that a separate edge-walking procedure was
    created. This routine walks around the edge of a tetrahedra, assuming the 
    point p is on one of the edges of the tetrahedra.

*/

void edgewalk(std::vector<std::pair<tetra*, unsigned int>> &leaves, const point *p) {

    assert(leaves.size() == 1);

    auto pivot = leaves[0].first;
    auto pivot_pos = leaves[0].second;

    assert(pivot);

    tetra *left(nullptr), *right(nullptr);
    tetra *left_prev(pivot), *right_prev(pivot);

    for (unsigned int i = 0; i < 4; ++i) {
        if (!(pivot_pos & (1 << i))) {
            if (pivot->ngb[i]) {
                if (!left)
                    left = pivot->ngb[i];
                else
                    right = pivot->ngb[i];
            }
        }
    }

    if (!left) // no walk is possible
        return;

    while (left || right) {

        auto tmp = left;

        if (left && left == right) {

            leaves.emplace_back(std::pair<tetra*, unsigned int>(left, fast_inclusion(left, p)));
            return;
        }

        if (left) {

            auto left_pos = fast_inclusion(left, p);

            leaves.emplace_back(std::pair<tetra*, unsigned int>(left, left_pos));

            int i = 0;

            for ( ; i < 4; ++i) {

                if (!(left_pos & (1 << i))) {
                    if (left->ngb[i]) {
                        if (left->ngb[i] != left_prev && left->ngb[i] != right) {

                            left_prev = left;
                            left = left->ngb[i];
                            break;
                        }    
                    }
                }
            }

            if (i == 4)
                left = nullptr;
        }

        if (right) {

            auto right_pos = fast_inclusion(right, p);

            leaves.emplace_back(std::pair<tetra*, unsigned int>(right, right_pos));

            int i = 0;

            for ( ; i < 4; ++i) {

                if (!(right_pos & (1 << i))) {
                    if (right->ngb[i]) {
                        if (right->ngb[i] != right_prev && right->ngb[i] != tmp) {

                            right_prev = right;
                            right = right->ngb[i];
                            break;
                        }    
                    }
                }
            }

            if (i == 4)
                right = nullptr;
        }
    }

    return;
}
