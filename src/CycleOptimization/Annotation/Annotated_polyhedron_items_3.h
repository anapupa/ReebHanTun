// (c) 2011 Oleksiy Busaryev

#ifndef CGAL_ANNOTATED_POLYHEDRON_ITEMS_3_H
#define CGAL_ANNOTATED_POLYHEDRON_ITEMS_3_H

#include <CGAL/Polyhedron_items_3.h>

#include "HalfedgeDS_annotated_vertex_base.h"
#include "HalfedgeDS_annotated_halfedge_base.h"
#include "HalfedgeDS_annotated_face_base.h"

namespace CGAL {
    struct Annotated_polyhedron_items_3
            : public Polyhedron_items_3 {
        template<class Refs_, class Traits_>
        struct Vertex_wrapper {
            typedef HalfedgeDS_annotated_vertex_base<Refs_,
                    typename Traits_::Point_3> Vertex;
        };

        template<class Refs_, class Traits_>
        struct Halfedge_wrapper {
            typedef HalfedgeDS_annotated_halfedge_base<Refs_,
                    Traits_> Halfedge;
        };

        template<class Refs_, class Traits_>
        struct Face_wrapper {
            typedef HalfedgeDS_annotated_face_base<Refs_,
                    Traits_> Face;
        };
    };
}

#endif // CGAL_ANNOTATED_POLYHEDRON_ITEMS_3_H
