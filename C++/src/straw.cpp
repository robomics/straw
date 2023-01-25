/*
  The MIT License (MIT)

  Copyright (c) 2017-2021 Aiden Lab, Rice University, Baylor College of Medicine

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
*/
#include "straw/straw.h"

#include <zlib.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <streambuf>
#include <utility>
#include <vector>

using namespace std;

/*
  Straw: fast C++ implementation of dump. Not as fully featured as the
  Java version. Reads the .hic file, finds the appropriate matrix and slice
  of data, and outputs as text in sparse upper triangular format.

  Currently only supporting matrices.

  Usage: straw [observed/oe/expected] <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2]
  <chr2>[:y1:y2] <BP/FRAG> <binsize>
 */

static void parsePositions(const string &chrLoc, string &chrom, int64_t &pos1, int64_t &pos2,
                           const HiCFile::ChromosomeMap &map) {
    string x, y;
    stringstream ss(chrLoc);
    getline(ss, chrom, ':');
    if (map.count(chrom) == 0) {
        throw std::runtime_error(chrom + " not found in the file");
    }

    if (getline(ss, x, ':') && getline(ss, y, ':')) {
        pos1 = stol(x);
        pos2 = stol(y);
    } else {
        pos1 = 0LL;
        pos2 = map.at(chrom).length;
    }
}

vector<contactRecord> straw(const string &matrixType, const string &norm, const string &fileName,
                            const string &chr1loc, const string &chr2loc, const string &unit,
                            int32_t binsize) {
    try {
        if (!(unit == "BP" || unit == "FRAG")) {
            throw std::runtime_error("Norm specified incorrectly, must be one of <BP/FRAG>");
        }

        HiCFile hiCFile(fileName);
        const auto &chroms = hiCFile.getChromosomeMap();
        string chr1, chr2;
        int64_t origRegionIndices[4] = {-100LL, -100LL, -100LL, -100LL};
        parsePositions(chr1loc, chr1, origRegionIndices[0], origRegionIndices[1], chroms);
        parsePositions((chr2loc), chr2, origRegionIndices[2], origRegionIndices[3], chroms);

        if (chroms.at(chr1).index > chroms.at(chr2).index) {
            auto mzd = hiCFile.getMatrixZoomData(chr2, chr1, matrixType, norm, unit, binsize);
            return mzd.getRecords(origRegionIndices[2], origRegionIndices[3], origRegionIndices[0],
                                  origRegionIndices[1]);
        } else {
            auto mzd = hiCFile.getMatrixZoomData(chr1, chr2, matrixType, norm, unit, binsize);
            return mzd.getRecords(origRegionIndices[0], origRegionIndices[1], origRegionIndices[2],
                                  origRegionIndices[3]);
        }
    } catch (const std::exception &e) {
        throw std::runtime_error(std::string("straw encountered the following error: ") + e.what());
    }
}

vector<vector<float> > strawAsMatrix(const string &matrixType, const string &norm,
                                     const string &fileName, const string &chr1loc,
                                     const string &chr2loc, const string &unit, int32_t binsize) {
    try {
        if (!(unit == "BP" || unit == "FRAG")) {
            throw std::runtime_error("Norm specified incorrectly, must be one of <BP/FRAG>");
        }

        HiCFile hiCFile(fileName);
        const auto &chroms = hiCFile.getChromosomeMap();
        string chr1, chr2;
        int64_t origRegionIndices[4] = {-100LL, -100LL, -100LL, -100LL};
        parsePositions(chr1loc, chr1, origRegionIndices[0], origRegionIndices[1], chroms);
        parsePositions((chr2loc), chr2, origRegionIndices[2], origRegionIndices[3], chroms);

        if (chroms.at(chr1).index > chroms.at(chr2).index) {
            auto mzd = hiCFile.getMatrixZoomData(chr2, chr1, matrixType, norm, unit, binsize);
            return mzd.getRecordsAsMatrix(origRegionIndices[2], origRegionIndices[3],
                                          origRegionIndices[0], origRegionIndices[1]);

        } else {
            auto mzd = hiCFile.getMatrixZoomData(chr1, chr2, matrixType, norm, unit, binsize);
            return mzd.getRecordsAsMatrix(origRegionIndices[0], origRegionIndices[1],
                                          origRegionIndices[2], origRegionIndices[3]);
        }
    } catch (const std::exception &e) {
        throw std::runtime_error(std::string("strawAsMatrix encountered the following error: ") +
                                 e.what());
    }
}

int64_t getNumRecordsForFile(const string &fileName, int32_t binsize, bool interOnly) {
    try {
        HiCFile hiCFile(fileName);
        int64_t totalNumRecords = 0;

        int32_t indexOffset = 0;
        if (interOnly) {
            indexOffset = 1;
        }

        vector<chromosome> chromosomes = hiCFile.getChromosomes();
        for (int32_t i = 0; i < chromosomes.size(); i++) {
            if (chromosomes[i].index <= 0) continue;
            for (int32_t j = i + indexOffset; j < chromosomes.size(); j++) {
                if (chromosomes[j].index <= 0) continue;
                const auto idx = std::minmax({chromosomes[i].index, chromosomes[j].index});
                const auto &chrom1 = chromosomes[idx.first].name;
                const auto &chrom2 = chromosomes[idx.second].name;
                totalNumRecords +=
                    hiCFile.getMatrixZoomData(chrom1, chrom2, "observed", "NONE", "BP", binsize)
                        .getNumberOfTotalRecords();
            }
        }

        return totalNumRecords;

    } catch (const std::exception &e) {
        throw std::runtime_error(
            std::string("getNumRecordsForFile encountered the following error: ") + e.what());
    }
}
