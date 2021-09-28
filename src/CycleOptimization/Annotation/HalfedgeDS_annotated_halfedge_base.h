// (c) 2011 Oleksiy Busaryev

#ifndef CGAL_HALFEDGEDS_ANNOTATED_HALFEDGE_BASE_H
#define CGAL_HALFEDGEDS_ANNOTATED_HALFEDGE_BASE_H

#include <CGAL/HalfedgeDS_halfedge_base.h>

#include <boost/dynamic_bitset.hpp>

namespace CGAL {
    template<class Refs_, class Traits_>
    struct HalfedgeDS_annotated_halfedge_base
            : public HalfedgeDS_halfedge_base<Refs_, Tag_false> {
        bool is_in_tree, is_in_cotree;
        boost::dynamic_bitset<> annotation;
    };
}

#endif // CGAL_HALFEDGEDS_ANNOTATED_HALFEDGE_BASE_H
