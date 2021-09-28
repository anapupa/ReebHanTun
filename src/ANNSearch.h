/*
(c) 2012 Fengtao Fan
*/
#ifndef _ANN_SEARCH_H_
#define _ANN_SEARCH_H_

#include <ANN/ANN.h>                    // ANN declarations
#include <ANN/ANNx.h>                    // more ANN declarations
#include <ANN/ANNperf.h>                // performance evaluation

#include <vector>
//
namespace ANNSearch {
    void GetMaxDistancePoint(std::vector<double> &xCoord,
                             std::vector<double> &yCoord,
                             const double xRange,
                             const double yRange,
                             const int xSample,
                             const int ySample,
                             double &xMaxCoord,
                             double &yMaxCoord);
};

#endif //_ANN_SEARCH_H_