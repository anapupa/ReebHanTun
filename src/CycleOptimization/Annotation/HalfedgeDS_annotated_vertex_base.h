// (c) 2011 Oleksiy Busaryev

#ifndef CGAL_HALFEDGEDS_ANNOTATED_VERTEX_BASE_H
#define CGAL_HALFEDGEDS_ANNOTATED_VERTEX_BASE_H

#include <CGAL/HalfedgeDS_vertex_base.h>

namespace CGAL {
    template<class Refs_, class Point_>
    struct HalfedgeDS_annotated_vertex_base
            : public HalfedgeDS_vertex_base<Refs_, Tag_true, Point_> {
        typedef HalfedgeDS_vertex_base <Refs_, Tag_true, Point_> Base;
        typedef typename Refs_::Vertex_handle Vertex_handle;
        typedef typename Refs_::Halfedge_handle Halfedge_handle;

        HalfedgeDS_annotated_vertex_base()
                : Base() {
        }

        HalfedgeDS_annotated_vertex_base(Point_ p_)
                : Base(p_) {
        }

        size_t index;
        Vertex_handle parent;
        Halfedge_handle to_parent;
    };
}

#endif // CGAL_HALFEDGEDS_ANNOTATED_VERTEX_BASE_H
