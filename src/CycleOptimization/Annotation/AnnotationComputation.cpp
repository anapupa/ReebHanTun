/*
(c) 2012 Fengtao Fan
*/
#include "AnnotationComputation.h"


#include "Annotated_polyhedron_3.h"
#include "Polyhedron_annotator_3.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include<CGAL/Polyhedron_incremental_builder_3.h>

#include <boost/filesystem.hpp>

#include <fstream>

#include <boost/progress.hpp>

#include "SimpleMesh.h"

/********************************/
typedef CGAL::Simple_cartesian<float> Kernel;
typedef CGAL::Annotated_polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::HalfedgeDS HalfedgeDS;
/*Modifier*/
// A modifier creating a triangle with the incremental builder.
template<class HDS>
class polyhedron_builder : public CGAL::Modifier_base<HDS> {
public:
    std::vector<float> &coords;
    std::vector<int> &tris;

    polyhedron_builder(std::vector<float> &_coords, std::vector<int> &_tris) : coords(_coords), tris(_tris) {}

    void operator()(HDS &hds) {
        typedef typename HDS::Vertex Vertex;
        typedef typename Vertex::Point Point;

        // create a cgal incremental builder
        CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
        B.begin_surface(coords.size() / 3, tris.size() / 3); //100

        // add the polyhedron vertices
        for (int i = 0; i < (int) coords.size(); i += 3) {
            B.add_vertex(Point(coords[i + 0], coords[i + 1], coords[i + 2]));
        }

        // add the polyhedron triangles
        for (int i = 0; i < (int) tris.size(); i += 3) {
            B.begin_facet();
            B.add_vertex_to_facet(tris[i + 0]);
            B.add_vertex_to_facet(tris[i + 1]);
            B.add_vertex_to_facet(tris[i + 2]);
            B.end_facet();
        }

        // finish up the surface
        B.end_surface();
    }
};

/********************************/
int ComputeAnnotation(_SimpleMesh &inMesh, std::vector<std::vector<char> > &vecEdgeAnno, std::vector<int> &tris) {
    std::vector<float> coords;
    //std::vector<int>	tris;
    //
    for (unsigned int i = 0; i < inMesh.vecVertex.size(); i++) {
        coords.push_back(inMesh.vecVertex[i].x);
        coords.push_back(inMesh.vecVertex[i].y);
        coords.push_back(inMesh.vecVertex[i].z);
    }
    //
    //std::cout << coords.size() << std::endl;
    //
    //for (unsigned int i = 0; i < orientTri.size(); i++)
    //{
    //	tris.push_back(orientTri[i][0]);
    //	tris.push_back(orientTri[i][1]);
    //	tris.push_back(orientTri[i][2]);
    //}
    //std::cout << tris.size() << std::endl;
    /*******************************************/
    using namespace std;
    using namespace CGAL;
    using namespace boost;
    //using namespace filesystem;

    typedef Polyhedron::Vertex_handle Vertex_handle;
    typedef Polyhedron::Edge_iterator Edge_iterator;

//	string input_filename( "eight.off" );
//	ifstream input( input_filename.c_str() );
//	if ( !input.is_open() )
//	{
//		cerr << "Cannot open " << input_filename << " for reading" << endl;
//		return EXIT_FAILURE;
//	}
    Polyhedron polyhedron;
    //input >> polyhedron;
    //input.close();
    polyhedron_builder<HalfedgeDS> builder(coords, tris);
    polyhedron.delegate(builder);



    //cout << polyhedron.size_of_vertices() << " vertices" << endl;
    //cout << polyhedron.size_of_facets() << " facets" << endl;

    //CGAL_assertion( polyhedron.is_triangle( polyhedron.halfedges_begin()));

    std::cout << "Time for computing edge annotation : " << std::endl;
    {
        boost::progress_timer t;
        annotate_edges(polyhedron);
    }

    size_t index(0);
    Polyhedron::Vertex_iterator it_v(polyhedron.vertices_begin());
    for (; it_v != polyhedron.vertices_end(); ++it_v)
        it_v->index = index++;

    /////////////
    //string output_filename( "model_annotations.txt" );
    //ofstream output( output_filename );
    //if ( !output.is_open() )
    //{
    //	cerr << "Cannot open " << output_filename << " for writing" << endl;
    //	return EXIT_FAILURE;
    //}
    Edge_iterator it_e(polyhedron.edges_begin());
    //output << it_e->annotation.size() << endl;
    //
    vecEdgeAnno.resize(inMesh.vecEdge.size());
    std::vector<char> tempAnnotation(it_e->annotation.size());
    //
    for (; it_e != polyhedron.edges_end(); ++it_e) {
        Vertex_handle h_a(it_e->vertex());
        Vertex_handle h_b(it_e->opposite()->vertex());

        //output << h_a->index << ' ' << h_b->index << ' ' <<
        //	it_e->annotation << endl;
        for (boost::dynamic_bitset<>::size_type i = 0; i < it_e->annotation.size(); i++) {
            tempAnnotation[i] = it_e->annotation[i];
        }
        //
        int edge_index = -1;
        for (int i = 0; i < inMesh.vecVertex[h_a->index].adjEdges.size(); i++) {
            int loc_edge_index = inMesh.vecVertex[h_a->index].adjEdges[i];
            if (inMesh.vecEdge[loc_edge_index].v0 == h_b->index ||
                inMesh.vecEdge[loc_edge_index].v1 == h_b->index) {
                edge_index = loc_edge_index;
                break;
            }
        }
        if (edge_index < 0) {
            std::cout << "EDGE NOT MATCHED in ANNOTATION COMPUTATION" << std::endl;
            exit(9);
        }
        //
        vecEdgeAnno[edge_index] = tempAnnotation;
        //

    }
    //std::cout << vecEdgeAnno[0].size() << std::endl;
    //output.close();
    //cout << "Annotations written to " << output_filename << endl;

    return 0;
}
