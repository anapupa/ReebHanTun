// (c) 2011 Oleksiy Busaryev

#ifndef CGAL_HALFEDGEDS_ANNOTATED_FACE_BASE_H
#define CGAL_HALFEDGEDS_ANNOTATED_FACE_BASE_H

#include <CGAL/HalfedgeDS_face_base.h>

#include <boost/dynamic_bitset.hpp>

namespace CGAL {
    template<class Refs_, class Traits_>
    struct HalfedgeDS_annotated_face_base
            : public HalfedgeDS_face_base<Refs_> {
        typedef typename Refs_::Face_handle Face_handle;
        typedef typename Refs_::Halfedge_handle Halfedge_handle;
        bool is_hole;
        Face_handle parent;
        Halfedge_handle to_parent;
    };
}

#endif // CGAL_HALFEDGEDS_ANNOTATED_FACE_BASE_H
