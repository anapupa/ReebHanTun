// (c) 2011 Oleksiy Busaryev

#ifndef CGAL_POLYHEDRON_ANNOTATOR_3_H
#define CGAL_POLYHEDRON_ANNOTATOR_3_H

#include "Annotated_polyhedron_3.h"

namespace CGAL {
    template<typename Polyhedron_>
    class Polyhedron_annotator_3 {
        typedef typename Polyhedron_::Traits::Point_3 Point;
        typedef typename Polyhedron_::Traits::Vector_3 Vector;

        typedef typename Polyhedron_::Vertex_handle Vertex_handle;
        typedef typename Polyhedron_::Halfedge_handle Halfedge_handle;
        typedef typename Polyhedron_::Face_handle Face_handle;

        typedef typename Polyhedron_::Vertex_iterator Vertex_iterator;
        typedef typename Polyhedron_::Halfedge_iterator Halfedge_iterator;
        typedef typename Polyhedron_::Edge_iterator Edge_iterator;
        typedef typename Polyhedron_::Facet_iterator Facet_iterator;

        typedef typename Polyhedron_::Halfedge_around_vertex_circulator
                Halfedge_around_vertex_circulator;
        typedef typename Polyhedron_::Halfedge_around_facet_circulator
                Halfedge_around_facet_circulator;

    public:

        Polyhedron_annotator_3(Polyhedron_ &polyhedron_);
        //~Polyhedron_annotator_3()
        //{

        //    std::cout << "deconstrutor" << std::endl;
        //}

        void fill_holes();

        void compute_dual_tree();

        void compute_dual_cotree();

        void compute_annotations();

    private:

        Polyhedron_ &m_polyhedron;
    };

    template<typename Polyhedron_>
    inline
    Polyhedron_annotator_3<Polyhedron_>::Polyhedron_annotator_3(
            Polyhedron_ &polyhedron_)
            : m_polyhedron(polyhedron_) {
    }

    template<typename Polyhedron_>
    inline void
    Polyhedron_annotator_3<Polyhedron_>::fill_holes() {
        Facet_iterator it_f(m_polyhedron.facets_begin());
        for (; it_f != m_polyhedron.facets_end(); ++it_f)
            it_f->is_hole = false;

        Halfedge_iterator it_he(m_polyhedron.halfedges_begin());
        for (; it_he != m_polyhedron.halfedges_end(); ++it_he) {
            if (it_he->is_border()) {
                m_polyhedron.fill_hole(it_he);
                it_he->face()->is_hole = true;
            }
        }
    }

    template<typename Polyhedron_>
    inline void
    Polyhedron_annotator_3<Polyhedron_>::compute_dual_tree() {
        Halfedge_iterator it_he(m_polyhedron.halfedges_begin());
        for (; it_he != m_polyhedron.halfedges_end(); ++it_he)
            it_he->is_in_tree = false;

        Facet_iterator it_f(m_polyhedron.facets_begin());
        for (; it_f != m_polyhedron.facets_end(); ++it_f) {
            it_f->parent = Face_handle();
            it_f->to_parent = Halfedge_handle();
        }

        std::vector<Face_handle> face_stack;
        face_stack.reserve(m_polyhedron.size_of_facets());

        it_f = m_polyhedron.facets_begin();
        for (; it_f != m_polyhedron.facets_end(); ++it_f) {
            if (!it_f->is_hole)
                continue;

            it_f->parent = it_f;
            face_stack.push_back(it_f);

            while (!face_stack.empty()) {
                Face_handle h_f(face_stack.back());
                face_stack.pop_back();

                Halfedge_around_facet_circulator ci_begin(h_f->facet_begin()),
                        ci(ci_begin);
                do {
                    Face_handle h_g(ci->opposite()->face());
                    if (h_g->is_hole)
                        continue;

                    if (h_g->parent == Face_handle()) {
                        h_g->parent = h_f;
                        h_g->to_parent = ci;
                        face_stack.push_back(h_g);
                        ci->is_in_tree = true;
                        ci->opposite()->is_in_tree = true;
                    }
                } while (++ci != ci_begin);
            }
        }

        it_f = m_polyhedron.facets_begin();
        for (; it_f != m_polyhedron.facets_end(); ++it_f) {
            if (it_f->parent != Face_handle())
                continue;

            it_f->parent = it_f;
            face_stack.push_back(it_f);

            while (!face_stack.empty()) {
                Face_handle h_f(face_stack.back());
                face_stack.pop_back();

                Halfedge_around_facet_circulator ci_begin(h_f->facet_begin()),
                        ci(ci_begin);
                do {
                    Face_handle h_g(ci->opposite()->face());
                    if (h_g->parent == Face_handle()) {
                        h_g->parent = h_f;
                        h_g->to_parent = ci;
                        face_stack.push_back(h_g);
                        ci->is_in_tree = true;
                        ci->opposite()->is_in_tree = true;
                    }
                } while (++ci != ci_begin);
            }
        }
    }

    template<typename Polyhedron_>
    inline void
    Polyhedron_annotator_3<Polyhedron_>::compute_dual_cotree() {
        Halfedge_iterator it_he(m_polyhedron.halfedges_begin());
        for (; it_he != m_polyhedron.halfedges_end(); ++it_he)
            it_he->is_in_cotree = false;

        Vertex_iterator it_v(m_polyhedron.vertices_begin());
        for (; it_v != m_polyhedron.vertices_end(); ++it_v) {
            it_v->parent = Vertex_handle();
            it_v->to_parent = Halfedge_handle();
        }

        std::vector<Vertex_handle> vertex_stack;
        vertex_stack.reserve(m_polyhedron.size_of_vertices());

        it_v = m_polyhedron.vertices_begin();
        for (; it_v != m_polyhedron.vertices_end(); ++it_v) {
            if (it_v->parent != Vertex_handle())
                continue;

            it_v->parent = it_v;
            vertex_stack.push_back(it_v);

            while (!vertex_stack.empty()) {
                Vertex_handle h_v(vertex_stack.back());
                vertex_stack.pop_back();

                Halfedge_around_vertex_circulator ci_begin(h_v->vertex_begin()),
                        ci(ci_begin);
                do {
                    if (ci->is_in_tree)
                        continue;

                    Vertex_handle h_w(ci->opposite()->vertex());
                    if (h_w->parent == Vertex_handle()) {
                        h_w->parent = h_v;
                        h_w->to_parent = ci;
                        vertex_stack.push_back(h_w);
                        ci->is_in_cotree = true;
                        ci->opposite()->is_in_cotree = true;
                    }
                } while (++ci != ci_begin);
            }
        }
    }

    template<typename Polyhedron_>
    inline void
    Polyhedron_annotator_3<Polyhedron_>::compute_annotations() {
        size_t rank(0);
        Edge_iterator it_e(m_polyhedron.edges_begin());
        for (; it_e != m_polyhedron.edges_end(); ++it_e) {
            if (!it_e->is_in_tree && !it_e->is_in_cotree)
                ++rank;
        }

        Halfedge_iterator it_he(m_polyhedron.halfedges_begin());
        for (; it_he != m_polyhedron.halfedges_end(); ++it_he)
            it_he->annotation.resize(rank);

        size_t index(0);
        it_e = m_polyhedron.edges_begin();
        for (; it_e != m_polyhedron.edges_end(); ++it_e) {
            if (it_e->is_in_tree || it_e->is_in_cotree)
                continue;

            it_e->annotation.flip(index);
            it_e->opposite()->annotation.flip(index);

            Face_handle h_f(it_e->face());
            Face_handle h_g(it_e->opposite()->face());

            while (h_f->parent != h_f) {
                h_f->to_parent->annotation.flip(index);
                h_f->to_parent->opposite()->annotation.flip(index);
                h_f = h_f->parent;
            }

            while (h_g->parent != h_g) {
                h_g->to_parent->annotation.flip(index);
                h_g->to_parent->opposite()->annotation.flip(index);
                h_g = h_g->parent;
            }

            ++index;
        }
    }

    template<typename Traits_, typename Items_>
    inline void
    annotate_edges(Annotated_polyhedron_3 <Traits_,
                   Items_> &polyhedron_) {
        typedef Annotated_polyhedron_3<Traits_, Items_> Polyhedron;
        Polyhedron_annotator_3<Polyhedron> annotator(polyhedron_);
        annotator.fill_holes();
        annotator.compute_dual_tree();
        annotator.compute_dual_cotree();
        annotator.compute_annotations();
    }
}

#endif // CGAL_POLYHEDRON_ANNOTATOR_3_H
