/*
(c) 2012 Fengtao Fan
*/
#include "ANNSearch.h"

void ANNSearch::GetMaxDistancePoint(std::vector<double> &xCoord,
                                    std::vector<double> &yCoord,
                                    const double xRange,
                                    const double yRange,
                                    const int xSample,
                                    const int ySample,
                                    double &xMaxCoord,
                                    double &yMaxCoord) {
    int nPts = (int) xCoord.size();
    int dim = 2;
    //

    /*************/
    ANNpoint queryPt;
    ANNidxArray nnIdx = NULL;
    ANNdistArray dists = NULL;

    queryPt = annAllocPt(dim); // need to delloc

    nnIdx = new ANNidx[nPts];
    dists = new ANNdist[nPts];
    /**************************/

    ANNpointArray dataPts;
    dataPts = annAllocPts(nPts, dim);
    // adding vertices

    for (int i = 0; i < nPts; i++) {
        dataPts[i][0] = xCoord[i];
        dataPts[i][1] = yCoord[i];
    }
    ///////////////////
    ANNkd_tree *kdTree;
    kdTree = new ANNkd_tree(dataPts,
                            nPts,
                            dim);
    //
    int NeighborCounter = 2;
    //
    double xStep = xRange / xSample;
    double yStep = yRange / ySample;
    double curX = 0.0;
    double curY = 0.0;
    xMaxCoord = curX;
    yMaxCoord = curY;
    double maxDistance = 0.0;
    for (int i = 0; i < xSample; i++) {
        curX = i * xStep;
        for (int j = 0; j < ySample; j++) {
            curY = j * yStep;
            //
            queryPt[0] = curX;
            queryPt[1] = curY;
            //
            kdTree->annkSearch(queryPt,
                               1,
                               nnIdx,
                               dists);
            //
            if (dists[0] > maxDistance) {
                maxDistance = dists[0];
                xMaxCoord = curX;
                yMaxCoord = curY;
            }
        }
    }
    // clean data
    //annDeallocPts(dataPts);

    delete[] nnIdx;
    delete[] dists;
    delete kdTree;
    if (dataPts)
        annDeallocPts(dataPts);
    if (queryPt)
        annDeallocPt(queryPt);
    annClose();
    return;
}