/*

    This is the collection of our Delaunay refinement routines

*/

#include "structures.hpp"
#include "mpreal.h"

/*

    This is our insphere routine. There is potential support for arbitrary
    precision but at this time, it is not needed as our point sets are entirely
    integer in nature. Any results being < 1 in our case mean that they arise
    due to numerical methods found in the Eigen library. Extending the precision
    to 200 digits will show that the determinant gets progressively smaller and
    smaller, showing that it trends towards 0, e.g. 7.53e-15 as a double vs 
    1.23e-199 using mpfr::mpreal. This is just an example.

*/

bool insphere(const tetra &t, const point &p) {

    Eigen::Matrix<double, 5, 5> insphere;

    auto &a = *t.p[0], &b = *t.p[1], &c = *t.p[2], &d = *t.p[3];

    insphere << 1, a.x, a.y, a.z, a.q,
                1, b.x, b.y, b.z, b.q,
                1, c.x, c.y, c.z, c.q,
                1, d.x, d.y, d.z, d.q,
                1, p.x, p.y, p.z, p.q;

    double dist = insphere.determinant();

    if (verbose >= 2) {

        std::cout << "\nInsphere test matrix :\n" << std::endl;
        std::cout << insphere << std::endl;

        std::cout << "\nresult : " << dist << std::endl;
    }

    if (std::abs(dist) < 1)
        return true;
/*
    if (std::abs(dist) < 1) {

        //std::cout << insphere << std::endl;

        const int digits = 200;
        mpfr::mpreal::set_default_prec(mpfr::digits2bits(digits));

        Eigen::Matrix<mpfr::mpreal, 5, 5> exact_insphere;

        exact_insphere = insphere.cast<mpfr::mpreal>();

        mpfr::mpreal exact_dist = exact_insphere.determinant();
        //std::cout << "inexact : " << dist << std::endl;
        //std::cout << "exact : " << std::setprecision(digits) << exact_dist << std::endl << std::endl;

        if (exact_dist < 0)
            return false;
        else
            return true;
    }

*/

    if (dist < 0) // inside circumsphere
        return false;
    else
        return true; // on or outside sphere
}

/*

    This is our 2-to-3 flip routine and works rather simply.

    Let's say we have two tetrahedra made up of the points abcd and abce and
    they are not Delaunay. They obviously share the face abc and we find that
    the line de passes through abc. This is the criterion for a 2-to-3 flip
    and the resulting tetrahedra are as follows :

    abcd, abce => abde, bcde, cade

*/

void two_to_three_flip(mesh &m, tetra *t, tetra *v, std::vector<tetra*> &pile) {

    std::unordered_set<tetra*> ancestors;
    ancestors.reserve(16);

    for (int j = 0; j < 4; ++j) {

        if (t->ngb[j])
            ancestors.emplace(t->ngb[j]);

        if (v->ngb[j])
            ancestors.emplace(v->ngb[j]);
    }

    point *x(t->p[3]), *y(t->ngb[3]->p[t->nid[3]]);

    point *a(t->f[3][0]), *b(t->f[3][1]), *c(t->f[3][2]);

    #ifndef NDEBUG
    
        Eigen::Matrix<double, 4, 4> o1, o2, o3;

        o1 << 1, a->x, a->y, a->z,
              1, b->x, b->y, b->z,
              1, y->x, y->y, y->z,
              1, x->x, x->y, x->z;

        o2 << 1, b->x, b->y, b->z,
              1, c->x, c->y, c->z,
              1, y->x, y->y, y->z,
              1, x->x, x->y, x->z;

        o3 << 1, c->x, c->y, c->z,
              1, a->x, a->y, a->z,
              1, y->x, y->y, y->z,
              1, x->x, x->y, x->z;

        assert(o1.determinant() > 0);
        assert(o2.determinant() > 0);
        assert(o3.determinant() > 0);

    #endif

    auto u = m.t_position;

    if (verbose >= 2) {

        std::cout << "Potential 2-to-3 flip candidates : " << std::endl;
        std::cout << *t << std::endl << *v << std::endl;
    }

    new(t) tetra(a, b, y, x);
    new(v) tetra(b, c, y, x);
    new(u) tetra(c, a, y, x);

    ancestors.emplace(t);
    ancestors.emplace(v);
    ancestors.emplace(u);

    for (auto it = ancestors.begin(); it != ancestors.end(); ++it)
        for (auto it2 = std::next(it, 1); it2 != ancestors.end(); ++it2)
            link_tetra(*it, *it2);

    ++m.t_position;
    ++m.num_tetra;

    //pile.emplace_back(t);
    pile.emplace_back(v);
    pile.emplace_back(u);

    if (verbose >= 1) {

        std::cout << "\n2-to-3 mesh repair : " << std::endl;
        std::cout << *t << std::endl << *v << std::endl << *u << std::endl;
    }

    m.m_position = t;

    return;
}

/*

    A 4-to-4 flip is more interesting than a simple 2-to-3 one.

    Here when we have the tetrahedra abcd, abce, the line de intersects
    one of the edges of the triangle abc. Hence we have to examine for an
    appropriate number of neighbours and replace all 4!

    In our example, let's say de intersects the edge ab. For the sake of
    example, we have the following 4 tetrahdra :

    abcd, abce, abdf, abef

    We then replace these tetrahedra with :

    afde, fbde, bcde, cade

    Note that there must be exactly 4 tetrahedra in this scenario.

*/

bool four_to_four_flip(tetra *t, unsigned int x, std::vector<tetra*> &pile, mesh &m) {

    assert(t);

    tetra *v = nullptr;

    if (!t->ngb[3]) // no 4-to-4 flip is possible
        return false;
    else
        v = t->ngb[3];

    tetra *w = nullptr;

    switch (x) {

        case 3 : // edge 10

            if (!t->ngb[2])
                return false;
            else
                w = t->ngb[2];
            break;

        case 5 : // edge 02

            if (!t->ngb[1])
                return false;
            else
                w = t->ngb[1];
            break;

        case 6 : // edge 21

            if (!t->ngb[0])
                return false;
            else
                w = t->ngb[0];
            break;
    }

    tetra *r = nullptr;

    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            if (v->ngb[i] && v->ngb[i] != t && v->ngb[i] == w->ngb[j])
                r = v->ngb[i];

    if (!r)
        return false;


    point *p(t->p[3]), *q(v->p[t->nid[3]]); // p is insertion point

    std::array<point*, 4> vtx;
    vtx.fill(nullptr);

    switch (x) {

        case 3 :

            vtx[0] = t->f[3][0];
            vtx[1] = t->ngb[2]->p[t->nid[2]];
            vtx[2] = t->f[3][1];
            vtx[3] = t->f[3][2];

            break;

        case 5 :

            vtx[0] = t->f[3][0];
            vtx[1] = t->f[3][1];
            vtx[2] = t->f[3][2];
            vtx[3] = t->ngb[1]->p[t->nid[1]];

            break;

        case 6 :

            vtx[0] = t->f[3][0];
            vtx[1] = t->f[3][1];
            vtx[2] = t->ngb[0]->p[t->nid[0]];
            vtx[3] = t->f[3][2];

            break;
    }

    #ifndef NDEBUG

        Eigen::Matrix<double, 4, 4> o1, o2, o3, o4;

        o1 << 1, vtx[0]->x, vtx[0]->y, vtx[0]->z,
              1, vtx[1]->x, vtx[1]->y, vtx[1]->z,
              1, q->x, q->y, q->z,
              1, p->x, p->y, p->z;

        o2 << 1, vtx[1]->x, vtx[1]->y, vtx[1]->z,
              1, vtx[2]->x, vtx[2]->y, vtx[2]->z,
              1, q->x, q->y, q->z,
              1, p->x, p->y, p->z;

        o3 << 1, vtx[2]->x, vtx[2]->y, vtx[2]->z,
              1, vtx[3]->x, vtx[3]->y, vtx[3]->z,
              1, q->x, q->y, q->z,
              1, p->x, p->y, p->z;

        o4 << 1, vtx[3]->x, vtx[3]->y, vtx[3]->z,
              1, vtx[0]->x, vtx[0]->y, vtx[0]->z,
              1, q->x, q->y, q->z,
              1, p->x, p->y, p->z;

        assert(o1.determinant() > 0);
        assert(o2.determinant() > 0);
        assert(o3.determinant() > 0);
        assert(o4.determinant() > 0);

    #endif

    std::unordered_set<tetra*> ancestors;
    ancestors.reserve(32);
    
    for (int i = 0; i < 4; ++i) {

        if (t->ngb[i])
            ancestors.insert(t->ngb[i]);

        if (v->ngb[i])
            ancestors.insert(v->ngb[i]);

        if (w->ngb[i])
            ancestors.insert(w->ngb[i]);

        if (r->ngb[i])
            ancestors.insert(r->ngb[i]);
    }

    if (verbose >= 2) {

        std::cout << "Potential 4-to-4 flip candidates : " << std::endl;
        std::cout << *t << std::endl << *v << std::endl << *w << std::endl << *r << std::endl;
    }

    new(t) tetra(vtx[0], vtx[1], q, p);
    new(v) tetra(vtx[1], vtx[2], q, p);
    new(w) tetra(vtx[2], vtx[3], q, p);
    new(r) tetra(vtx[3], vtx[0], q, p);

    for (auto it = ancestors.begin(); it != ancestors.end(); ++it)
        for (auto it2 = std::next(it, 1); it2 != ancestors.end(); ++it2)
            link_tetra(*it, *it2);

    //pile.emplace_back(t);
    pile.emplace_back(v);
    pile.emplace_back(w);
    pile.emplace_back(r);

    if (verbose >= 1) {

        std::cout << "\n4-to-4 mesh repair : " << std::endl;
        std::cout << *t << std::endl << *v << std::endl << *w << std::endl << *r << std::endl;
    }

    m.m_position = t;

    return true;
}

/*

    When I first began this project, oh so long ago, a 3-to-2 flip was one of
    the most conceptually difficult things for me to visualize. First off,
    the flip requires that exactly 3 tetrahedra share an edge. Second, 
    the two tetradehra that violate Delaunay conditions must have points
    opposite the shared face that create a line intersecting only one
    outer half-space of the shared face.

    As complicated and cumbersome as that sounds, consider the following :

    abcd and abce violate the insphere criterion.

    We find that the line de intersects the plane spanned by abc outside the
    triangle of abc. We then use a routine to examine exactly which half-space
    the intersection point is outside. If it is outside two half-spaces, there
    is no repair possible so we can skip the refinement. However, if the point
    is outside exactly one half-space, we can do perform the flip.

    We examine for an appropriate neighbour. The criteria for
    this is that the third tetrahedra is neighbours with both our abcd and abce
    tetrahedra.

    So let's say we now have the following : abcd, abce, abde

    With all the proper conditions in place, the flip will yield the following,

    cdea, cdeb

    This may not seem like much but the flip is difficult to visualize because
    the repair is essentially a rotation of the tetrahedra. The shared edge
    suddently becomes the two points above and below the plane.
    
*/

bool three_to_two_flip(std::vector<tetra*> &pile, tetra *t, const std::array<point*, 3> &f, const point *a, const point *b, mesh &m) {

    tetra *w = nullptr;

    unsigned int hl = halfspace_locator(f, a, b);

    switch (hl) {

        case 13 :
        case 49 :
        case 61 :
        case 7 :
        case 52 :
        case 55 :
        case 19 :
        case 28 :
        case 31 :

            if (verbose >= 1)
                std::cout << "Outside two half-spaces. 3-to-2 repair impossible...\n" << std::endl;                
            return false;

        case 53 :
            w = t->ngb[2];
            break;

        case 23 :
            w = t->ngb[0];
            break;

        case 29 :
            w = t->ngb[1];
            break;

        default :
            throw std::runtime_error("Error in halfspace_locator routine...");    
    }

    tetra *v = t->ngb[3];

    bool t1(false), t2(false);

    for (int i = 0; i < 4; ++i) {
        if (w->ngb[i] == t)
            t1 = true;

        if (w->ngb[i] == v)
            t2 = true;
    }

    if (!t1 || !t2) {

        if (verbose >= 1)
            std::cout << "Too many tetrahedra for a 3-to-2 flip...\n" << std::endl;
        return false;
    }

    for (auto it = pile.begin(); it < pile.end(); ++it)
        if (*it == w) {
            
            *it = pile.back();
            pile.pop_back();
        }

    point *p0(nullptr), *p1(nullptr), *p2(nullptr), *x(nullptr), *y(nullptr);

    switch (hl) {

        case 53 :

            x = t->p[0];
            y = t->p[1];

            p0 = t->p[2];
            p1 = t->p[3];
            p2 = v->p[t->nid[3]];            
            break;

        case 23 :

            x = t->p[1];
            y = t->p[2];

            p0 = t->p[0];
            p1 = t->p[3];
            p2 = v->p[t->nid[3]];
            break;

        case 29 :

            x = t->p[2];
            y = t->p[0];

            p0 = t->p[1];
            p1 = t->p[3];
            p2 = v->p[t->nid[3]];
            break;
    }

    #ifndef NDEBUG

        Eigen::Matrix<double, 4, 4> o1, o2;

        o1 << 1, x->x, x->y, x->z,
              1, p2->x, p2->y, p2->z,
              1, p0->x, p0->y, p0->z,
              1, p1->x, p1->y, p1->z;

        o2 << 1, y->x, y->y, y->z,
              1, p0->x, p0->y, p0->z,
              1, p2->x, p2->y, p2->z,
              1, p1->x, p1->y, p1->z;

        assert(o1.determinant() > 0);
        assert(o2.determinant() > 0);

    #endif

    std::unordered_set<tetra*> ancestors;
    ancestors.reserve(16);

    for (int i = 0; i < 4; ++i) {

        if (t->ngb[i])
            ancestors.insert(t->ngb[i]);
        
        if (v->ngb[i])
            ancestors.insert(v->ngb[i]);

        if (w->ngb[i]) {
            ancestors.insert(w->ngb[i]);

            for (int j = 0; j < 4; ++j)
                if (w->ngb[i]->ngb[j] == w)
                    w->ngb[i]->ngb[j] = nullptr;
        }
    }

    if (verbose >= 2) {

        std::cout << "Potential 3-to-2 flip candidates : " << std::endl;
        std::cout << *t << std::endl << *v << std::endl << *w << std::endl;
    }

    new(w) tetra(nullptr, nullptr, nullptr, nullptr);
    new(t) tetra(x, p2, p0, p1);
    new(v) tetra(y, p0, p2, p1);

    ancestors.erase(w);

    for (auto it = ancestors.begin(); it != ancestors.end(); ++it)
        for (auto it2 = std::next(it, 1); it2 != ancestors.end(); ++it2)
            link_tetra(*it, *it2);

    //pile.emplace_back(t);
    pile.emplace_back(v);

    if (verbose >= 1) {

        std::cout << "\n3-to-2 mesh repair :\n" << std::endl;
        std::cout << *t << std::endl << *v << std::endl;
    }

    m.m_position = t;
    --m.num_tetra;

    return true;
}

/*

    This is the almighty refinement routine. It controls the stack of
    tetrahedra to be tested for Delaunayhood.

    The only face we actually need to test is the face opposite the inserted
    point which is fortunate for us because the inserted point is always the
    4th in our newly created tetrahedra and the face opposite that point is
    indexed by that as well. So t->p[3] is the inserted point and t->f[3] is
    the face opposite t->p[3]. t->ngb[3] is the neighbour that shares t->f[3].

    This assumes that every sibling fracture is Delaunay.

*/

void delaunay(mesh &m, std::vector<tetra*> &pile) {

    while (!pile.empty()) {

        tetra *t = pile.back();

        /*
            Perform mesh repairs, if necessary
        */

        if (t->ngb[3] && !insphere(*t, *(t->ngb[3]->p[t->nid[3]]))) {

            tetra *v = t->ngb[3];

            unsigned int x = fast_triangle_ray_inclusion(t->f[3], t->p[3], t->ngb[3]->p[t->nid[3]]);

            switch (x) {

                case 7 :
                    two_to_three_flip(m, t, v, pile);
                    continue;

                case 3 :
                case 5 :
                case 6 :
                    if (!four_to_four_flip(t, x, pile, m))
                        pile.pop_back();
                    continue;

                case 0 :
                    if (!three_to_two_flip(pile, t, t->f[3], t->p[3], t->ngb[3]->p[t->nid[3]], m))
                        pile.pop_back();
                    continue;

                default :
                    throw std::runtime_error("Region not accounted for in Delaunay refinement routine...");
            }
        }

        else
            pile.pop_back();
    }

    return;
}
