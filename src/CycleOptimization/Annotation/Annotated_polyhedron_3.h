// (c) 2011 Oleksiy Busaryev

#ifndef CGAL_ANNOTATED_POLYHEDRON_3_H
#define CGAL_ANNOTATED_POLYHEDRON_3_H

#include "Annotated_polyhedron_items_3.h"

#include <CGAL/Polyhedron_3.h>

namespace CGAL {
    template<typename Traits_, typename Items_ =
    Annotated_polyhedron_items_3>
    class Annotated_polyhedron_3
            : public Polyhedron_3<Traits_, Items_> {
    };
}

#endif // CGAL_ANNOTATED_POLYHEDRON_3_H
