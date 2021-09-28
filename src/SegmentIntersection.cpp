/*
(c) 2012 Fengtao Fan
*/
#include <CGAL/intersections.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <vector>


#include "SegmentIntersection.h"

namespace IntersectonComputing {
    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel_2;
    typedef Kernel_2::Vector_2 Vector_2;
    typedef Kernel_2::Point_2 Point_2;
    typedef Kernel_2::Segment_2 Segment_2;

    void NewValueForLinkNumber(const int crossing_type,
                               const int heightFlag,
                               const int incrementalValue,
                               int &_link_number) {
        if (crossing_type > 0) {// uv is clockwise crossing ab
            // amplify the crossing number 4 times
            // at the end it needs to divide 8 to get the link number
            _link_number += (incrementalValue * heightFlag);
        } else {// uv is counterclockwise crossing ab
            _link_number -= (incrementalValue * heightFlag);
        }
    }

    int LinkNumberOfTow3DLines(std::vector<Vector3> &inPly_a, std::vector<Vector3> &inPly_b) {
        CGAL::Cartesian_converter<Kernel, Kernel_2> converter;
        std::vector<Point_2> ply_a;
        std::vector<Point_2> ply_b;
        Segment_2 ab, uv;
        //
        CGAL::Object resIntersection;
        Point_2 resPoint;
        Segment_2 resSeg;
        //

        for (unsigned int i = 0; i < inPly_a.size(); i++) {
            ply_a.push_back(Point_2(inPly_a[i][0], inPly_a[i][1]));
        }

        for (unsigned int i = 0; i < inPly_b.size(); i++) {
            ply_b.push_back(Point_2(inPly_b[i][0], inPly_b[i][1]));
        }
        //
        const int cross_from_right_to_left = 1;
        const int cross_from_left_to_right = -1;
        //
        int cross_type_by_fixed_uv = 0;
        //
        double height_ab = 0.0;
        double height_uv = 0.0;
        //
        int _link_number = 0;
        //
        double ab_x_max = 0.0;
        double ab_y_max = 0.0;
        double ab_x_min = 0.0;
        double ab_y_min = 0.0;
        double uv_x_max = 0.0;
        double uv_y_max = 0.0;
        double uv_x_min = 0.0;
        double uv_y_min = 0.0;
        for (unsigned int i = 0; i < ply_a.size() - 1; i++) {
            ab = Segment_2(ply_a[i], ply_a[i + 1]);
            ab_x_max = std::max(ply_a[i].x(), ply_a[i + 1].x());
            ab_x_min = std::min(ply_a[i].x(), ply_a[i + 1].x());
            ab_y_max = std::max(ply_a[i].y(), ply_a[i + 1].y());
            ab_y_min = std::min(ply_a[i].y(), ply_a[i + 1].y());
            for (unsigned int j = 0; j < ply_b.size() - 1; j++) {
                uv = Segment_2(ply_b[j], ply_b[j + 1]);

                uv_x_max = std::max(ply_b[j].x(), ply_b[j + 1].x());
                uv_x_min = std::min(ply_b[j].x(), ply_b[j + 1].x());
                uv_y_max = std::max(ply_b[j].y(), ply_b[j + 1].y());
                uv_y_min = std::min(ply_b[j].y(), ply_b[j + 1].y());
                //
                if (!((uv_x_max < ab_x_min) || (ab_x_max < uv_x_min) ||
                      (uv_y_max < ab_y_min) || (ab_y_max < uv_y_min))
                        ) {
                    resIntersection = CGAL::intersection(ab, uv);
                    //
                    // check the intersection results
                    if (CGAL::assign(resPoint, resIntersection)) {// point intersection
                        if (resPoint == ab.source() || resPoint == ab.target() || resPoint == uv.source() ||
                            resPoint == uv.target()) {// end point intersection
                            // check the crossing type
                            std::cout << "point in segment" << std::endl;
                            exit(1);
                            if (resPoint == ab.source()) {
                                if (CGAL::left_turn(uv.source(), uv.target(), ab.target())) {// cross from right to left
                                    cross_type_by_fixed_uv = cross_from_right_to_left;
                                } else {
                                    cross_type_by_fixed_uv = cross_from_left_to_right;
                                    //std::cout << "cross from left to right" << std::endl;
                                }
                            } else {
                                if (CGAL::left_turn(uv.source(), uv.target(), ab.source())) {// cross from left to left
                                    cross_type_by_fixed_uv = cross_from_left_to_right;
                                    //std::cout << "cross from left to right" << std::endl;
                                } else {
                                    cross_type_by_fixed_uv = cross_from_right_to_left;
                                    //std::cout << "cross from right to left" << std::endl;
                                }
                            }
                            if (resPoint == uv.source() || resPoint == uv.target()) {// end point in both segment
                                if (resPoint == ab.source()) {
                                    height_ab = inPly_a[i][2];
                                } else {
                                    height_ab = inPly_a[i + 1][2];
                                }
                                //
                                if (resPoint == uv.source()) {
                                    height_uv = inPly_b[j][2];
                                } else {
                                    height_uv = inPly_b[j + 1][2];
                                }
                                if (height_uv > height_ab) {//
                                    if (CGAL::left_turn(uv.source(), uv.target(), ab.target())) {
                                        _link_number += 4;
                                        //NewValueForLinkNumber(cross_type_by_fixed_uv, 1, 1, _link_number);
                                    } else {
                                        _link_number -= 4;
                                    }
                                } else {
                                    if (CGAL::left_turn(ab.source(), ab.target(), uv.target())) {
                                        _link_number += 4;
                                        //NewValueForLinkNumber(cross_type_by_fixed_uv, 1, 1, _link_number);
                                    } else {
                                        _link_number -= 4;
                                    }
                                    //NewValueForLinkNumber(cross_type_by_fixed_uv, -1, 1, _link_number);
                                }
                                //std::cout << " both " << std::endl;
                            } else {// endpoint in ab and intermediate point for uv
                                // case
                                std::cout << " one end and the other in middle " << std::endl;
                                if (resPoint == ab.source()) {
                                    height_ab = inPly_a[i][2];
                                } else {
                                    height_ab = inPly_a[i + 1][2];
                                }
                                Vector_2 upDiff = resPoint - uv.source();
                                Kernel::FT para_t = upDiff.squared_length() / uv.squared_length();
                                double uv_t = sqrt(converter(para_t));
                                height_uv = inPly_b[j][2] + uv_t * (inPly_b[j + 1][2] - inPly_b[j][2]);
                                //
                                if (height_uv > height_ab) {//
                                    if (CGAL::left_turn(uv.source(), uv.target(), ab.source())) {
                                        _link_number += 4;
                                        //NewValueForLinkNumber(cross_type_by_fixed_uv, 1, 1, _link_number);
                                    } else {
                                        _link_number -= 4;
                                    }
                                    std::cout << "new" << std::endl;
                                    //NewValueForLinkNumber(cross_type_by_fixed_uv, 1, 2, _link_number);
                                } else {
                                    if (CGAL::left_turn(ab.source(), ab.target(), uv.source())) {
                                        _link_number += 4;
                                        //NewValueForLinkNumber(cross_type_by_fixed_uv, 1, 1, _link_number);
                                    } else {
                                        _link_number -= 4;
                                    }
                                    //NewValueForLinkNumber(cross_type_by_fixed_uv, -1, 2, _link_number);
                                }
                                //std::cout << "t " << uv_t << std::endl;
                            }
                        } else {// regular point
                            //std::cout << "regular " << std::endl;
                            if (CGAL::left_turn(uv.source(), uv.target(), ab.source())) {// cross from left to right
                                cross_type_by_fixed_uv = cross_from_left_to_right;
                                //std::cout << " cross from left to right "  << std::endl;
                            } else {
                                cross_type_by_fixed_uv = cross_from_right_to_left;
                                //std::cout << " cross from right to left "  << std::endl;
                            }
                            {//
                                Vector_2 upDiff = resPoint - uv.source();
                                Kernel::FT para_t = upDiff.squared_length() / uv.squared_length();
                                double uv_t = sqrt(converter(para_t));
                                height_uv = inPly_b[j][2] + uv_t * (inPly_b[j + 1][2] - inPly_b[j][2]);
                            }
                            {
                                Vector_2 upDiff = resPoint - ab.source();
                                Kernel::FT para_t = upDiff.squared_length() / ab.squared_length();
                                double ab_t = sqrt(converter(para_t));
                                height_ab = inPly_a[i][2] + ab_t * (inPly_a[i + 1][2] - inPly_a[i][2]);
                            }
                            //
                            if (height_uv > height_ab) {//
                                //if (CGAL::left_turn(uv.source(), uv.target(), ab.target()))
                                //	{
                                //		_link_number += 4;
                                //	}
                                //	else
                                //	{
                                //		_link_number -= 4;
                                //	}
                                //std::cout << "new"<< std::endl;
                                NewValueForLinkNumber(cross_type_by_fixed_uv, 1, 4, _link_number);
                            } else {
                                /*if (CGAL::left_turn(ab.source(), ab.target(), uv.target()))
                                    {
                                        _link_number += 4;
                                    }
                                    else
                                    {
                                        _link_number -= 4;
                                    }*/
                                NewValueForLinkNumber(cross_type_by_fixed_uv, -1, 4, _link_number);
                            }
                        }
                    } else {
                        if (CGAL::assign(resSeg, resIntersection)) {// segment intersection
                            // colinear , error
                            std::cout << "colinear " << std::endl;
                            exit(0);
                        }
                        //else
                        //{// no intersection
                        //	std::cout << "no intersection " << std::endl;
                        //}
                    }
                }
            } // j
        }// i
        //
        if (abs(_link_number) % 8 != 0) {
            std::cout << "WRONG LINK NUMBER COMPUTATION" << std::endl;
            exit(0);
        } else {
            _link_number = _link_number / 8;
        }
        return _link_number;
    }
}
