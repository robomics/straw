/*
  The MIT License (MIT)

  Copyright (c) 2011-2016 Broad Institute, Aiden Lab

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
#ifndef STRAW_H
#define STRAW_H

#include <curl/curl.h>

#include <cstdint>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <vector>

#include "straw/internal/common.h"

// pointer structure for reading blocks or matrices, holds the size and position
struct indexEntry {
    std::int64_t size;
    std::int64_t position;
};

// sparse matrixType entry
struct contactRecord {
    std::int32_t binX;
    std::int32_t binY;
    float counts;
};

// chromosome
struct chromosome {
    std::string name;
    std::int32_t index;
    std::int64_t length;
};

namespace internal {

class HiCFileStream {
    std::ifstream fin_{};
    CURL_ptr curl_{nullptr, &curl_easy_cleanup};

   public:
    explicit HiCFileStream(const std::string &fileName);
    bool isLocal() const noexcept;
    bool isRemote() const noexcept;

    // TODO: this is a bit of a code smell.
    //       Code requiring access to curl_ should probably
    //       be a member function of this class
    CURL_ptr &curl() noexcept;
    const CURL_ptr &curl() const noexcept;

    // TODO: same as above
    std::ifstream &fin() noexcept;
    const std::ifstream &fin() const noexcept;

    void readCompressedBytes(indexEntry idx, std::string &buffer);
    std::string readCompressedBytes(indexEntry idx);

   private:
    static CURL_ptr initRemoteFile(const std::string &url);
    static std::ifstream initRegularFile(const std::string &path);
};

class MatrixZoomData {
    bool isIntra;
    std::string fileName;
    std::int64_t myFilePos = 0LL;
    std::vector<double> expectedValues;
    std::vector<double> c1Norm;
    std::vector<double> c2Norm;
    std::int32_t c1 = 0;
    std::int32_t c2 = 0;
    std::string matrixType;
    std::string norm;
    std::int32_t version = 0;
    std::int32_t resolution = 0;
    std::int32_t numBins1 = 0;
    std::int32_t numBins2 = 0;
    float sumCounts;
    std::int32_t blockBinCount, blockColumnCount;
    std::map<std::int32_t, indexEntry> blockMap;
    double avgCount;

   public:
    MatrixZoomData(const chromosome &chrom1, const chromosome &chrom2,
                   const std::string &matrixType, const std::string &norm, const std::string &unit,
                   std::int32_t resolution, std::int32_t &version, std::int64_t &master, std::int64_t &totalFileSize,
                   const std::string &fileName);

    std::vector<contactRecord> getRecords(std::int64_t gx0, std::int64_t gx1, std::int64_t gy0, std::int64_t gy1);

    std::vector<std::vector<float>> getRecordsAsMatrix(std::int64_t gx0, std::int64_t gx1, std::int64_t gy0,
                                                       std::int64_t gy1);

    std::int64_t getNumberOfTotalRecords();

   private:
    static std::vector<double> readNormalizationVectorFromFooter(indexEntry cNormEntry,
                                                                 std::int32_t &version,
                                                                 const std::string &fileName);

    static bool isInRange(std::int32_t r, std::int32_t c, std::int32_t numRows, std::int32_t numCols);

    std::set<std::int32_t> getBlockNumbers(std::int64_t *regionIndices) const;

    std::vector<double> getNormVector(std::int32_t index);

    std::vector<double> getExpectedValues();
};

void readFooter(std::istream &fin, std::int64_t master, std::int32_t version, std::int32_t c1, std::int32_t c2,
                const std::string &matrixType, const std::string &norm, const std::string &unit,
                std::int32_t resolution, std::int64_t &myFilePos, indexEntry &c1NormEntry,
                indexEntry &c2NormEntry, std::vector<double> &expectedValues);

// reads the footer from the master pointer location. takes in the chromosomes,
// norm, unit (BP or FRAG) and resolution or binsize, and sets the file
// position of the matrix and the normalization vectors for those chromosomes
// at the given normalization and resolution
void readFooterURL(CURL_ptr &curl, std::int64_t master, std::int32_t version, std::int32_t c1, std::int32_t c2,
                   const std::string &matrixType, const std::string &norm, const std::string &unit,
                   std::int32_t resolution, std::int64_t &myFilePos, indexEntry &c1NormEntry,
                   indexEntry &c2NormEntry, std::vector<double> &expectedValues);

// TODO remove me!
inline std::string readCompressedBytesFromFile(const std::string &fileName, indexEntry idx) {
    return internal::HiCFileStream(fileName).readCompressedBytes(idx);
}

}  // namespace internal

class HiCFile {
   public:
    using ChromosomeMap = std::map<std::string, chromosome>;

   private:
    std::string fileName;
    std::int64_t master = 0LL;
    ChromosomeMap chromosomes;
    std::string genomeID;
    std::int32_t numChromosomes = 0;
    std::int32_t version = 0;
    std::int64_t nviPosition = 0LL;
    std::int64_t nviLength = 0LL;
    std::vector<std::int32_t> resolutions;
    std::int64_t totalFileSize;

   public:
    explicit HiCFile(std::string fileName_);

    const std::string &getGenomeID() const noexcept;

    const std::vector<std::int32_t> &getResolutions() const noexcept;

    std::vector<chromosome> getChromosomes() const;
    auto getChromosomeMap() const noexcept -> const ChromosomeMap &;

    internal::MatrixZoomData getMatrixZoomData(const std::string &chr1, const std::string &chr2,
                                               const std::string &matrixType,
                                               const std::string &norm, const std::string &unit,
                                               std::int32_t resolution);

   private:
    static std::int64_t readTotalFileSize(const std::string &url);

    auto readHeader(std::istream &fin, std::int64_t &masterIndexPosition, std::string &genomeID,
                    std::int32_t &numChromosomes, std::int32_t &version, std::int64_t &nviPosition,
                    std::int64_t &nviLength) -> ChromosomeMap;
};

std::map<std::int32_t, indexEntry> readMatrixZoomData(std::istream &fin, const std::string &myunit,
                                                 std::int32_t mybinsize, float &mySumCounts,
                                                 std::int32_t &myBlockBinCount,
                                                 std::int32_t &myBlockColumnCount, bool &found);

std::vector<double> readNormalizationVector(std::istream &fin, indexEntry entry);

std::vector<contactRecord> straw(const std::string &matrixType, const std::string &norm,
                                 const std::string &fname, const std::string &chr1loc,
                                 const std::string &chr2loc, const std::string &unit,
                                 std::int32_t binsize);

std::int64_t getNumRecordsForFile(const std::string &filename, std::int32_t binsize, bool interOnly);

#endif
