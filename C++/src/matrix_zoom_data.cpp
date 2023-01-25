/*
The MIT License (MIT)

Copyright (c) 2011-2016 Broad Institute, Aiden Lab

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
                                                              copies of the Software, and to permit
    persons to whom the Software is furnished to do so, subject to the following conditions:

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

#include <zlib.h>

#include <cmath>
#include <cstdint>
#include <ios>
#include <istream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "straw/internal/common.h"
#include "straw/straw.h"

namespace internal {
using namespace std;

static void setValuesForMZD(istream &fin, const string &myunit, float &mySumCounts,
                            int32_t &mybinsize, int32_t &myBlockBinCount,
                            int32_t &myBlockColumnCount, bool &found) {
    string unit;
    getline(fin, unit, '\0');                  // unit
    readInt32FromFile(fin);                    // Old "zoom" index -- not used
    float sumCounts = readFloatFromFile(fin);  // sumCounts
    readFloatFromFile(fin);                    // occupiedCellCount
    readFloatFromFile(fin);                    // stdDev
    readFloatFromFile(fin);                    // percent95
    int32_t binSize = readInt32FromFile(fin);
    int32_t blockBinCount = readInt32FromFile(fin);
    int32_t blockColumnCount = readInt32FromFile(fin);
    found = false;
    if (myunit == unit && mybinsize == binSize) {
        mySumCounts = sumCounts;
        myBlockBinCount = blockBinCount;
        myBlockColumnCount = blockColumnCount;
        found = true;
    }
}

static indexEntry readIndexEntry(istream &fin) {
    int64_t filePosition = readInt64FromFile(fin);
    int32_t blockSizeInBytes = readInt32FromFile(fin);
    indexEntry entry = indexEntry();
    entry.size = (int64_t)blockSizeInBytes;
    entry.position = filePosition;
    return entry;
}

static void populateBlockMap(istream &fin, int32_t nBlocks, map<int32_t, indexEntry> &blockMap) {
    for (int b = 0; b < nBlocks; b++) {
        int32_t blockNumber = readInt32FromFile(fin);
        blockMap[blockNumber] = readIndexEntry(fin);
    }
}

// reads the raw binned contact matrix at specified resolution, setting the block bin count and
// block column count
static map<int32_t, indexEntry> readMatrixZoomData(istream &fin, const string &myunit,
                                                   int32_t mybinsize, float &mySumCounts,
                                                   int32_t &myBlockBinCount,
                                                   int32_t &myBlockColumnCount, bool &found) {
    map<int32_t, indexEntry> blockMap;
    setValuesForMZD(fin, myunit, mySumCounts, mybinsize, myBlockBinCount, myBlockColumnCount,
                    found);

    int32_t nBlocks = readInt32FromFile(fin);
    if (found) {
        populateBlockMap(fin, nBlocks, blockMap);
    } else {
        fin.seekg(nBlocks * (sizeof(int32_t) + sizeof(int64_t) + sizeof(int32_t)), ios_base::cur);
    }
    return blockMap;
}

// reads the raw binned contact matrix at specified resolution, setting the block bin count and
// block column count
static map<int32_t, indexEntry> readMatrixZoomDataHttp(CURL_ptr &curl, int64_t &myFilePosition,
                                                       const string &myunit, int32_t mybinsize,
                                                       float &mySumCounts, int32_t &myBlockBinCount,
                                                       int32_t &myBlockColumnCount, bool &found) {
    map<int32_t, indexEntry> blockMap;
    int32_t header_size = 5 * sizeof(int32_t) + 4 * sizeof(float);
    auto buffer = getData(curl, myFilePosition, 1);
    if (buffer.front() == 'B') {
        header_size += 3;
    } else if (buffer.front() == 'F') {
        header_size += 5;
    } else {
        throw std::runtime_error("Unit not understood");
    }
    buffer = getData(curl, myFilePosition, header_size);
    memstream fin(buffer);
    setValuesForMZD(fin, myunit, mySumCounts, mybinsize, myBlockBinCount, myBlockColumnCount,
                    found);
    int32_t nBlocks = readInt32FromFile(fin);

    if (found) {
        int32_t chunkSize = nBlocks * (sizeof(int32_t) + sizeof(int64_t) + sizeof(int32_t));
        buffer = getData(curl, myFilePosition + header_size, chunkSize);
        fin = memstream(buffer);
        populateBlockMap(fin, nBlocks, blockMap);
    } else {
        myFilePosition = myFilePosition + header_size +
                         (nBlocks * (sizeof(int32_t) + sizeof(int64_t) + sizeof(int32_t)));
    }
    return blockMap;
}

// goes to the specified file pointer in http and finds the raw contact matrixType at specified
// resolution, calling readMatrixZoomData. sets blockbincount and blockcolumncount
static map<int32_t, indexEntry> readMatrixHttp(CURL_ptr &curl, int64_t myFilePosition,
                                               const string &unit, int32_t resolution,
                                               float &mySumCounts, int32_t &myBlockBinCount,
                                               int32_t &myBlockColumnCount) {
    int32_t size = sizeof(int32_t) * 3;
    auto buffer = getData(curl, myFilePosition, size);
    memstream bufin(buffer);

    int32_t c1 = readInt32FromFile(bufin);
    int32_t c2 = readInt32FromFile(bufin);
    int32_t nRes = readInt32FromFile(bufin);
    int32_t i = 0;
    bool found = false;
    myFilePosition = myFilePosition + size;
    map<int32_t, indexEntry> blockMap;

    while (i < nRes && !found) {
        // myFilePosition gets updated within call
        blockMap = readMatrixZoomDataHttp(curl, myFilePosition, unit, resolution, mySumCounts,
                                          myBlockBinCount, myBlockColumnCount, found);
        i++;
    }
    if (!found) {
        throw std::runtime_error("Error finding block data");
    }
    return blockMap;
}

// goes to the specified file pointer and finds the raw contact matrixType at specified resolution,
// calling readMatrixZoomData. sets blockbincount and blockcolumncount
static map<int32_t, indexEntry> readMatrix(istream &fin, int64_t myFilePosition, const string &unit,
                                           int32_t resolution, float &mySumCounts,
                                           int32_t &myBlockBinCount, int32_t &myBlockColumnCount) {
    map<int32_t, indexEntry> blockMap;

    fin.seekg(myFilePosition, ios::beg);
    int32_t c1 = readInt32FromFile(fin);
    int32_t c2 = readInt32FromFile(fin);
    int32_t nRes = readInt32FromFile(fin);
    int32_t i = 0;
    bool found = false;
    while (i < nRes && !found) {
        blockMap = readMatrixZoomData(fin, unit, resolution, mySumCounts, myBlockBinCount,
                                      myBlockColumnCount, found);
        i++;
    }
    if (!found) {
        throw std::runtime_error("Error finding block data");
    }
    return blockMap;
}

// reads the normalization vector from the file at the specified location
static vector<double> readNormalizationVector(istream &bufferin, int32_t version) {
    int64_t nValues;
    if (version > 8) {
        nValues = readInt64FromFile(bufferin);
    } else {
        nValues = (int64_t)readInt32FromFile(bufferin);
    }

    uint64_t numValues;
    numValues = static_cast<uint64_t>(nValues);
    vector<double> values(numValues);

    if (version > 8) {
        for (int i = 0; i < nValues; i++) {
            values[i] = (double)readFloatFromFile(bufferin);
        }
    } else {
        for (int i = 0; i < nValues; i++) {
            values[i] = readDoubleFromFile(bufferin);
        }
    }

    return values;
}

// gets the blocks that need to be read for this slice of the data.  needs blockbincount,
// blockcolumncount, and whether or not this is intrachromosomal.
static set<int32_t> getBlockNumbersForRegionFromBinPosition(const int64_t *regionIndices,
                                                            int32_t blockBinCount,
                                                            int32_t blockColumnCount, bool intra) {
    int32_t col1, col2, row1, row2;
    col1 = static_cast<int32_t>(regionIndices[0] / blockBinCount);
    col2 = static_cast<int32_t>((regionIndices[1] + 1) / blockBinCount);
    row1 = static_cast<int32_t>(regionIndices[2] / blockBinCount);
    row2 = static_cast<int32_t>((regionIndices[3] + 1) / blockBinCount);

    set<int32_t> blocksSet;
    // first check the upper triangular matrixType
    for (int r = row1; r <= row2; r++) {
        for (int c = col1; c <= col2; c++) {
            int32_t blockNumber = r * blockColumnCount + c;
            blocksSet.insert(blockNumber);
        }
    }
    // check region part that overlaps with lower left triangle but only if intrachromosomal
    if (intra) {
        for (int r = col1; r <= col2; r++) {
            for (int c = row1; c <= row2; c++) {
                int32_t blockNumber = r * blockColumnCount + c;
                blocksSet.insert(blockNumber);
            }
        }
    }
    return blocksSet;
}

static set<int32_t> getBlockNumbersForRegionFromBinPositionV9Intra(int64_t *regionIndices,
                                                                   int32_t blockBinCount,
                                                                   int32_t blockColumnCount) {
    // regionIndices is binX1 binX2 binY1 binY2
    set<int32_t> blocksSet;
    int32_t translatedLowerPAD, translatedHigherPAD, translatedNearerDepth, translatedFurtherDepth;
    translatedLowerPAD =
        static_cast<int32_t>((regionIndices[0] + regionIndices[2]) / 2 / blockBinCount);
    translatedHigherPAD =
        static_cast<int32_t>((regionIndices[1] + regionIndices[3]) / 2 / blockBinCount + 1);
    translatedNearerDepth = static_cast<int32_t>(
        log2(1 + abs(regionIndices[0] - regionIndices[3]) / std::sqrt(2) / blockBinCount));
    translatedFurtherDepth = static_cast<int32_t>(
        log2(1 + abs(regionIndices[1] - regionIndices[2]) / sqrt(2) / blockBinCount));

    // because code above assume above diagonal; but we could be below diagonal
    int32_t nearerDepth = min(translatedNearerDepth, translatedFurtherDepth);
    if ((regionIndices[0] > regionIndices[3] && regionIndices[1] < regionIndices[2]) ||
        (regionIndices[1] > regionIndices[2] && regionIndices[0] < regionIndices[3])) {
        nearerDepth = 0;
    }
    int32_t furtherDepth =
        max(translatedNearerDepth, translatedFurtherDepth) + 1;  // +1; integer divide rounds down

    for (int depth = nearerDepth; depth <= furtherDepth; depth++) {
        for (int pad = translatedLowerPAD; pad <= translatedHigherPAD; pad++) {
            int32_t blockNumber = depth * blockColumnCount + pad;
            blocksSet.insert(blockNumber);
        }
    }

    return blocksSet;
}

static void appendRecord(vector<contactRecord> &vector, int32_t index, int32_t binX, int32_t binY,
                         float counts) {
    contactRecord record = contactRecord();
    record.binX = binX;
    record.binY = binY;
    record.counts = counts;
    vector[index] = record;
}

static std::string decompressBlock(indexEntry idx, const std::string &compressedBytes) {
    std::string uncompressedBytes(idx.size * 10, '\0');  // biggest seen so far is 3
    z_stream infstream;
    infstream.zalloc = Z_NULL;
    infstream.zfree = Z_NULL;
    infstream.opaque = Z_NULL;
    infstream.avail_in = static_cast<uInt>(idx.size);  // size of input
    infstream.next_in = reinterpret_cast<Bytef *>(
        const_cast<char *>(&compressedBytes.front()));                  // input char array
    infstream.avail_out = static_cast<uInt>(uncompressedBytes.size());  // size of output
    infstream.next_out =
        reinterpret_cast<Bytef *>(&uncompressedBytes.front());  // output char array
    // the actual decompression work.
    inflateInit(&infstream);
    inflate(&infstream, Z_NO_FLUSH);
    inflateEnd(&infstream);

    uncompressedBytes.resize(static_cast<std::size_t>(infstream.total_out));
    return uncompressedBytes;
}

// this is the meat of reading the data.  takes in the block number and returns the set of contact
// records corresponding to that block.  the block data is compressed and must be decompressed using
// the zlib library functions
static vector<contactRecord> readBlock(const string &fileName, indexEntry idx, int32_t version) {
    if (idx.size <= 0) {
        vector<contactRecord> v;
        return v;
    }

    const auto compressedBytes = readCompressedBytesFromFile(fileName, idx);
    auto uncompressedBytes = decompressBlock(idx, compressedBytes);

    // create stream from buffer for ease of use
    memstream bufferin(uncompressedBytes);
    uint64_t nRecords;
    nRecords = static_cast<uint64_t>(readInt32FromFile(bufferin));
    vector<contactRecord> v(nRecords);
    // different versions have different specific formats
    if (version < 7) {
        for (uInt i = 0; i < nRecords; i++) {
            int32_t binX = readInt32FromFile(bufferin);
            int32_t binY = readInt32FromFile(bufferin);
            float counts = readFloatFromFile(bufferin);
            appendRecord(v, i, binX, binY, counts);
        }
    } else {
        int32_t binXOffset = readInt32FromFile(bufferin);
        int32_t binYOffset = readInt32FromFile(bufferin);
        bool useShort = readCharFromFile(bufferin) == 0;  // yes this is opposite of usual

        bool useShortBinX = true;
        bool useShortBinY = true;
        if (version > 8) {
            useShortBinX = readCharFromFile(bufferin) == 0;
            useShortBinY = readCharFromFile(bufferin) == 0;
        }

        char type = readCharFromFile(bufferin);
        int32_t index = 0;
        if (type == 1) {
            if (useShortBinX && useShortBinY) {
                int16_t rowCount = readInt16FromFile(bufferin);
                for (int i = 0; i < rowCount; i++) {
                    int32_t binY = binYOffset + readInt16FromFile(bufferin);
                    int16_t colCount = readInt16FromFile(bufferin);
                    for (int j = 0; j < colCount; j++) {
                        int32_t binX = binXOffset + readInt16FromFile(bufferin);
                        float counts;
                        if (useShort) {
                            counts = readInt16FromFile(bufferin);
                        } else {
                            counts = readFloatFromFile(bufferin);
                        }
                        appendRecord(v, index++, binX, binY, counts);
                    }
                }
            } else if (useShortBinX && !useShortBinY) {
                int32_t rowCount = readInt32FromFile(bufferin);
                for (int i = 0; i < rowCount; i++) {
                    int32_t binY = binYOffset + readInt32FromFile(bufferin);
                    int16_t colCount = readInt16FromFile(bufferin);
                    for (int j = 0; j < colCount; j++) {
                        int32_t binX = binXOffset + readInt16FromFile(bufferin);
                        float counts;
                        if (useShort) {
                            counts = readInt16FromFile(bufferin);
                        } else {
                            counts = readFloatFromFile(bufferin);
                        }
                        appendRecord(v, index++, binX, binY, counts);
                    }
                }
            } else if (!useShortBinX && useShortBinY) {
                int16_t rowCount = readInt16FromFile(bufferin);
                for (int i = 0; i < rowCount; i++) {
                    int32_t binY = binYOffset + readInt16FromFile(bufferin);
                    int32_t colCount = readInt32FromFile(bufferin);
                    for (int j = 0; j < colCount; j++) {
                        int32_t binX = binXOffset + readInt32FromFile(bufferin);
                        float counts;
                        if (useShort) {
                            counts = readInt16FromFile(bufferin);
                        } else {
                            counts = readFloatFromFile(bufferin);
                        }
                        appendRecord(v, index++, binX, binY, counts);
                    }
                }
            } else {
                int32_t rowCount = readInt32FromFile(bufferin);
                for (int i = 0; i < rowCount; i++) {
                    int32_t binY = binYOffset + readInt32FromFile(bufferin);
                    int32_t colCount = readInt32FromFile(bufferin);
                    for (int j = 0; j < colCount; j++) {
                        int32_t binX = binXOffset + readInt32FromFile(bufferin);
                        float counts;
                        if (useShort) {
                            counts = readInt16FromFile(bufferin);
                        } else {
                            counts = readFloatFromFile(bufferin);
                        }
                        appendRecord(v, index++, binX, binY, counts);
                    }
                }
            }
        } else if (type == 2) {
            int32_t nPts = readInt32FromFile(bufferin);
            int16_t w = readInt16FromFile(bufferin);

            for (int i = 0; i < nPts; i++) {
                // int32_t idx = (p.y - binOffset2) * w + (p.x - binOffset1);
                int32_t row = i / w;
                int32_t col = i - row * w;
                int32_t bin1 = binXOffset + col;
                int32_t bin2 = binYOffset + row;

                float counts;
                if (useShort) {
                    int16_t c = readInt16FromFile(bufferin);
                    if (c != -32768) {
                        appendRecord(v, index++, bin1, bin2, c);
                    }
                } else {
                    counts = readFloatFromFile(bufferin);
                    if (!isnan(counts)) {
                        appendRecord(v, index++, bin1, bin2, counts);
                    }
                }
            }
        }
    }
    return v;
}

static int32_t getNumRecordsInBlock(const string &fileName, indexEntry idx, int32_t version) {
    if (idx.size <= 0) {
        return 0;
    }
    const auto compressedBytes = readCompressedBytesFromFile(fileName, idx);
    auto uncompressedBytes = decompressBlock(idx, compressedBytes);

    // create stream from buffer for ease of use
    memstream bufferin(uncompressedBytes);
    uint64_t nRecords;
    nRecords = static_cast<uint64_t>(readInt32FromFile(bufferin));
    return nRecords;
}

MatrixZoomData::MatrixZoomData(const chromosome &chrom1, const chromosome &chrom2,
                               const string &matrixType, const string &norm, const string &unit,
                               int32_t resolution, int32_t &version, int64_t &master,
                               int64_t &totalFileSize, const string &fileName) {
    this->version = version;
    this->fileName = fileName;
    int32_t c01 = chrom1.index;
    int32_t c02 = chrom2.index;
    if (c01 <= c02) {  // default is ok
        this->c1 = c01;
        this->c2 = c02;
        this->numBins1 = static_cast<int32_t>(chrom1.length / resolution);
        this->numBins2 = static_cast<int32_t>(chrom2.length / resolution);
    } else {  // flip
        this->c1 = c02;
        this->c2 = c01;
        this->numBins1 = static_cast<int32_t>(chrom2.length / resolution);
        this->numBins2 = static_cast<int32_t>(chrom1.length / resolution);
    }
    isIntra = c1 == c2;

    this->matrixType = matrixType;
    this->norm = norm;
    this->resolution = resolution;

    HiCFileStream stream(fileName);
    indexEntry c1NormEntry{}, c2NormEntry{};

    if (stream.isRemote()) {
        readFooterURL(stream.curl(), master, version, c1, c2, matrixType, norm, unit, resolution,
                      myFilePos, c1NormEntry, c2NormEntry, expectedValues);
    } else {
        stream.fin().seekg(master, ios::beg);
        readFooter(stream.fin(), master, version, c1, c2, matrixType, norm, unit, resolution,
                   myFilePos, c1NormEntry, c2NormEntry, expectedValues);
    }

    if (norm != "NONE") {
        c1Norm = readNormalizationVectorFromFooter(c1NormEntry, version, fileName);
        if (isIntra) {
            c2Norm = c1Norm;
        } else {
            c2Norm = readNormalizationVectorFromFooter(c2NormEntry, version, fileName);
        }
    }

    stream = HiCFileStream((fileName));
    if (stream.isRemote()) {
        // readMatrix will assign blockBinCount and blockColumnCount
        blockMap = readMatrixHttp(stream.curl(), myFilePos, unit, resolution, sumCounts,
                                  blockBinCount, blockColumnCount);
    } else {
        // readMatrix will assign blockBinCount and blockColumnCount
        blockMap = readMatrix(stream.fin(), myFilePos, unit, resolution, sumCounts, blockBinCount,
                              blockColumnCount);
    }

    if (!isIntra) {
        avgCount = (sumCounts / numBins1) / numBins2;  // <= trying to avoid overflows
    }
}

vector<double> MatrixZoomData::readNormalizationVectorFromFooter(indexEntry cNormEntry,
                                                                 int32_t &version,
                                                                 const string &fileName) {
    auto buffer = readCompressedBytesFromFile(fileName, cNormEntry);
    memstream bufferin(buffer);
    vector<double> cNorm = readNormalizationVector(bufferin, version);
    return cNorm;
}

bool MatrixZoomData::isInRange(int32_t r, int32_t c, int32_t numRows, int32_t numCols) {
    return 0 <= r && r < numRows && 0 <= c && c < numCols;
}

set<int32_t> MatrixZoomData::getBlockNumbers(int64_t *regionIndices) const {
    if (version > 8 && isIntra) {
        return getBlockNumbersForRegionFromBinPositionV9Intra(regionIndices, blockBinCount,
                                                              blockColumnCount);
    } else {
        return getBlockNumbersForRegionFromBinPosition(regionIndices, blockBinCount,
                                                       blockColumnCount, isIntra);
    }
}

vector<double> MatrixZoomData::getNormVector(int32_t index) {
    if (index == c1) {
        return c1Norm;
    } else if (index == c2) {
        return c2Norm;
    }
    throw std::runtime_error("Invalid index provided: " + std::to_string(index) +
                             "\nShould be either " + std::to_string(c1) + " or " +
                             std::to_string(c2));
}

vector<double> MatrixZoomData::getExpectedValues() { return expectedValues; }

vector<contactRecord> MatrixZoomData::getRecords(int64_t gx0, int64_t gx1, int64_t gy0,
                                                 int64_t gy1) {
    int64_t origRegionIndices[] = {gx0, gx1, gy0, gy1};
    int64_t regionIndices[4];
    convertGenomeToBinPos(origRegionIndices, regionIndices, resolution);

    set<int32_t> blockNumbers = getBlockNumbers(regionIndices);
    vector<contactRecord> records;
    for (int32_t blockNumber : blockNumbers) {
        // get contacts in this block
        // cout << *it << " -- " << blockMap.size() << endl;
        // cout << blockMap[*it].size << " " <<  blockMap[*it].position << endl;
        vector<contactRecord> tmp_records = readBlock(fileName, blockMap[blockNumber], version);
        for (contactRecord rec : tmp_records) {
            int64_t x = rec.binX * resolution;
            int64_t y = rec.binY * resolution;

            if ((x >= origRegionIndices[0] && x <= origRegionIndices[1] &&
                 y >= origRegionIndices[2] && y <= origRegionIndices[3]) ||
                // or check regions that overlap with lower left
                (isIntra && y >= origRegionIndices[0] && y <= origRegionIndices[1] &&
                 x >= origRegionIndices[2] && x <= origRegionIndices[3])) {
                float c = rec.counts;
                if (norm != "NONE") {
                    c = static_cast<float>(c / (c1Norm[rec.binX] * c2Norm[rec.binY]));
                }
                if (matrixType == "oe") {
                    if (isIntra) {
                        c = static_cast<float>(
                            c / expectedValues[min(expectedValues.size() - 1,
                                                   (size_t)floor(abs(y - x) / resolution))]);
                    } else {
                        c = static_cast<float>(c / avgCount);
                    }
                } else if (matrixType == "expected") {
                    if (isIntra) {
                        c = static_cast<float>(expectedValues[min(
                            expectedValues.size() - 1, (size_t)floor(abs(y - x) / resolution))]);
                    } else {
                        c = static_cast<float>(avgCount);
                    }
                }

                if (!isnan(c) && !isinf(c)) {
                    contactRecord record = contactRecord();
                    record.binX = static_cast<int32_t>(x);
                    record.binY = static_cast<int32_t>(y);
                    record.counts = c;
                    records.push_back(record);
                }
            }
        }
    }
    return records;
}

vector<vector<float> > MatrixZoomData::getRecordsAsMatrix(int64_t gx0, int64_t gx1, int64_t gy0,
                                                          int64_t gy1) {
    vector<contactRecord> records = this->getRecords(gx0, gx1, gy0, gy1);
    if (records.empty()) {
        vector<vector<float> > res = vector<vector<float> >(1, vector<float>(1, 0));
        return res;
    }

    int64_t origRegionIndices[] = {gx0, gx1, gy0, gy1};
    int64_t regionIndices[4];
    convertGenomeToBinPos(origRegionIndices, regionIndices, resolution);

    int64_t originR = regionIndices[0];
    int64_t endR = regionIndices[1];
    int64_t originC = regionIndices[2];
    int64_t endC = regionIndices[3];
    int32_t numRows = endR - originR + 1;
    int32_t numCols = endC - originC + 1;
    float matrix[numRows][numCols];
    for (int32_t r = 0; r < numRows; r++) {
        for (int32_t c = 0; c < numCols; c++) {
            matrix[r][c] = 0;
        }
    }

    for (contactRecord cr : records) {
        if (isnan(cr.counts) || isinf(cr.counts)) continue;
        int32_t r = cr.binX / resolution - originR;
        int32_t c = cr.binY / resolution - originC;
        if (isInRange(r, c, numRows, numCols)) {
            matrix[r][c] = cr.counts;
        }
        if (isIntra) {
            r = cr.binY / resolution - originR;
            c = cr.binX / resolution - originC;
            if (isInRange(r, c, numRows, numCols)) {
                matrix[r][c] = cr.counts;
            }
        }
    }

    vector<vector<float> > finalMatrix;
    for (int32_t i = 0; i < numRows; i++) {
        vector<float> row;
        row.reserve(numCols);
        for (int32_t j = 0; j < numCols; j++) {
            row.push_back(matrix[i][j]);
        }
        finalMatrix.push_back(row);
    }
    return finalMatrix;
}

int64_t MatrixZoomData::getNumberOfTotalRecords() {
    int64_t regionIndices[4] = {0, numBins1, 0, numBins2};
    set<int32_t> blockNumbers = getBlockNumbers(regionIndices);
    int64_t total = 0;
    for (int32_t blockNumber : blockNumbers) {
        total += getNumRecordsInBlock(fileName, blockMap[blockNumber], version);
    }
    return total;
}

static void rollingMedian(vector<double> &initialValues, vector<double> &finalResult,
                          int32_t window) {
    // window is actually a ~wing-span
    if (window < 1) {
        finalResult = initialValues;
        return;
    }

    /*
    finalResult.push_back(initialValues[0]);
    int64_t length = initialValues.size();
    for (int64_t index = 1; index < length; index++) {
        int64_t initialIndex;
        int64_t finalIndex;
        if (index < window) {
            initialIndex = 0;
            finalIndex = 2 * index;
        } else {
            initialIndex = index - window;
            finalIndex = index + window;
        }

        if (finalIndex > length - 1) {
            finalIndex = length - 1;
        }

        vector<double> subVector = sliceVector(initialValues, initialIndex, finalIndex);
        finalResult.push_back(getMedian(subVector));
    }
    */
    finalResult = initialValues;
}

int64_t readThroughExpectedVectorURL(internal::CURL_ptr &curl, int64_t currentPointer,
                                     int32_t version, vector<double> &expectedValues,
                                     int64_t nValues, bool store, int32_t resolution) {
    if (store) {
        int32_t bufferSize = nValues * sizeof(double) + 10000;
        if (version > 8) {
            bufferSize = nValues * sizeof(float) + 10000;
        }
        auto buffer = internal::getData(curl, currentPointer, bufferSize);
        internal::memstream fin(buffer);

        vector<double> initialExpectedValues;
        if (version > 8) {
            populateVectorWithNumbers<float>(fin, initialExpectedValues, nValues);
        } else {
            populateVectorWithNumbers<double>(fin, initialExpectedValues, nValues);
        }

        // This seems to be copying initialValues into finalResult at the moment
        // rollingMedian(initialExpectedValues, expectedValues, window);
        // int32_t window = 5000000 / resolution;
        expectedValues = initialExpectedValues;
    }

    if (version > 8) {
        return nValues * sizeof(float);
    } else {
        return nValues * sizeof(double);
    }
}

void readThroughExpectedVector(int32_t version, istream &fin, vector<double> &expectedValues,
                               int64_t nValues, bool store, int32_t resolution) {
    if (store) {
        vector<double> initialExpectedValues;
        if (version > 8) {
            populateVectorWithNumbers<float>(fin, initialExpectedValues, nValues);
        } else {
            populateVectorWithNumbers<double>(fin, initialExpectedValues, nValues);
        }

        // This seems to be copying initialValues into finalResult at the moment
        // int32_t window = 5000000 / resolution;
        // rollingMedian(initialExpectedValues, expectedValues, window);
        expectedValues = initialExpectedValues;
    } else if (nValues > 0) {
        if (version > 8) {
            fin.seekg(nValues * sizeof(float), ios_base::cur);
        } else {
            fin.seekg(nValues * sizeof(double), ios_base::cur);
        }
    }
}

int64_t readThroughNormalizationFactorsURL(internal::CURL_ptr &curl, int64_t currentPointer,
                                           int32_t version, bool store,
                                           vector<double> &expectedValues, int32_t c1,
                                           int32_t nNormalizationFactors) {
    if (store) {
        int32_t bufferSize = nNormalizationFactors * (sizeof(int32_t) + sizeof(double)) + 10000;
        if (version > 8) {
            bufferSize = nNormalizationFactors * (sizeof(int32_t) + sizeof(float)) + 10000;
        }
        auto buffer = internal::getData(curl, currentPointer, bufferSize);
        internal::memstream fin(buffer);

        for (int j = 0; j < nNormalizationFactors; j++) {
            int32_t chrIdx = readInt32FromFile(fin);
            double v;
            if (version > 8) {
                v = readFloatFromFile(fin);
            } else {
                v = readDoubleFromFile(fin);
            }
            if (chrIdx == c1) {
                for (double &expectedValue : expectedValues) {
                    expectedValue = expectedValue / v;
                }
            }
        }
    }

    if (version > 8) {
        return nNormalizationFactors * (sizeof(int32_t) + sizeof(float));
    } else {
        return nNormalizationFactors * (sizeof(int32_t) + sizeof(double));
    }
}

void readThroughNormalizationFactors(istream &fin, int32_t version, bool store,
                                     vector<double> &expectedValues, int32_t c1) {
    int32_t nNormalizationFactors = readInt32FromFile(fin);
    if (store) {
        for (int j = 0; j < nNormalizationFactors; j++) {
            int32_t chrIdx = readInt32FromFile(fin);
            double v;
            if (version > 8) {
                v = readFloatFromFile(fin);
            } else {
                v = readDoubleFromFile(fin);
            }
            if (chrIdx == c1) {
                for (double &expectedValue : expectedValues) {
                    expectedValue = expectedValue / v;
                }
            }
        }
    } else if (nNormalizationFactors > 0) {
        if (version > 8) {
            fin.seekg(nNormalizationFactors * (sizeof(int32_t) + sizeof(float)), ios_base::cur);
        } else {
            fin.seekg(nNormalizationFactors * (sizeof(int32_t) + sizeof(double)), ios_base::cur);
        }
    }
}

// reads the footer from the master pointer location. takes in the chromosomes,
// norm, unit (BP or FRAG) and resolution or binsize, and sets the file
// position of the matrix and the normalization vectors for those chromosomes
// at the given normalization and resolution
void readFooterURL(CURL_ptr &curl, int64_t master, int32_t version, int32_t c1, int32_t c2,
                   const string &matrixType, const string &norm, const string &unit,
                   int32_t resolution, int64_t &myFilePos, indexEntry &c1NormEntry,
                   indexEntry &c2NormEntry, vector<double> &expectedValues) {
    int64_t currentPointer = master;

    auto buffer = getData(curl, currentPointer, 100);
    memstream fin(buffer);

    if (version > 8) {
        int64_t nBytes = readInt64FromFile(fin);
        currentPointer += 8;
    } else {
        int32_t nBytes = readInt32FromFile(fin);
        currentPointer += 4;
    }

    stringstream ss;
    ss << c1 << "_" << c2;
    string key = ss.str();

    int32_t nEntries = readInt32FromFile(fin);
    currentPointer += 4;

    int32_t bufferSize0 = nEntries * 50;
    buffer = getData(curl, currentPointer, bufferSize0);
    fin = memstream(buffer);

    bool found = false;
    string keyStr;
    for (int i = 0; i < nEntries; i++) {
        currentPointer += readFromFile(fin, keyStr);
        int64_t fpos = readInt64FromFile(fin);
        int32_t sizeinbytes = readInt32FromFile(fin);
        currentPointer += 12;
        if (keyStr == key) {
            myFilePos = fpos;
            found = true;
        }
    }
    if (!found) {
        throw std::runtime_error("Remote file doesn't have the given chr_chr map " + key);
    }

    if ((matrixType == "observed" && norm == "NONE") ||
        ((matrixType == "oe" || matrixType == "expected") && norm == "NONE" && c1 != c2))
        return;  // no need to read norm vector index

    // read in and ignore expected value maps; don't store; reading these to
    // get to norm vector index
    buffer = getData(curl, currentPointer, 100);
    fin = memstream(buffer);

    int32_t nExpectedValues = readInt32FromFile(fin);
    currentPointer += 4;
    string unit0;
    for (int i = 0; i < nExpectedValues; i++) {
        buffer = getData(curl, currentPointer, 1000);
        fin = memstream(buffer);

        currentPointer += readFromFile(fin, unit0);

        int32_t binSize = readInt32FromFile(fin);
        currentPointer += 4;

        int64_t nValues;
        if (version > 8) {
            nValues = readInt64FromFile(fin);
            currentPointer += 8;
        } else {
            nValues = (int64_t)readInt32FromFile(fin);
            currentPointer += 4;
        }

        bool store = c1 == c2 && (matrixType == "oe" || matrixType == "expected") &&
                     norm == "NONE" && unit0 == unit && binSize == resolution;

        currentPointer += readThroughExpectedVectorURL(curl, currentPointer, version,
                                                       expectedValues, nValues, store, resolution);

        buffer = getData(curl, currentPointer, 100);
        fin = memstream(buffer);
        int32_t nNormalizationFactors = readInt32FromFile(fin);
        currentPointer += 4;

        currentPointer += readThroughNormalizationFactorsURL(
            curl, currentPointer, version, store, expectedValues, c1, nNormalizationFactors);
    }

    if (c1 == c2 && (matrixType == "oe" || matrixType == "expected") && norm == "NONE") {
        if (expectedValues.empty()) {
            throw std::runtime_error("Remote file did not contain expected values vectors at " +
                                     std::to_string(resolution) + " " + unit);
        }
        return;
    }

    buffer = getData(curl, currentPointer, 100);
    fin = memstream(buffer);
    nExpectedValues = readInt32FromFile(fin);
    currentPointer += 4;
    string nType;
    for (int i = 0; i < nExpectedValues; i++) {
        buffer = getData(curl, currentPointer, 1000);
        fin = memstream(buffer);

        currentPointer += readFromFile(fin, nType);
        currentPointer += readFromFile(fin, unit0);

        int32_t binSize = readInt32FromFile(fin);
        currentPointer += 4;

        int64_t nValues;
        if (version > 8) {
            nValues = readInt64FromFile(fin);
            currentPointer += 8;
        } else {
            nValues = (int64_t)readInt32FromFile(fin);
            currentPointer += 4;
        }
        bool store = c1 == c2 && (matrixType == "oe" || matrixType == "expected") &&
                     nType == norm && unit0 == unit && binSize == resolution;

        currentPointer += readThroughExpectedVectorURL(curl, currentPointer, version,
                                                       expectedValues, nValues, store, resolution);

        buffer = getData(curl, currentPointer, 100);
        fin = memstream(buffer);
        int32_t nNormalizationFactors = readInt32FromFile(fin);
        currentPointer += 4;

        currentPointer += readThroughNormalizationFactorsURL(
            curl, currentPointer, version, store, expectedValues, c1, nNormalizationFactors);
    }

    if (c1 == c2 && (matrixType == "oe" || matrixType == "expected") && norm != "NONE") {
        if (expectedValues.empty()) {
            throw std::runtime_error(
                "Remote file did not contain normalized expected values vectors at " +
                std::to_string(resolution) + " " + unit);
        }
    }

    buffer = getData(curl, currentPointer, 100);
    fin = memstream(buffer);
    nEntries = readInt32FromFile(fin);
    currentPointer += 4;

    bool found1 = false;
    bool found2 = false;
    int32_t bufferSize2 = nEntries * 60;
    buffer = getData(curl, currentPointer, bufferSize2);
    fin = memstream(buffer);

    string normtype;
    for (int i = 0; i < nEntries; i++) {
        currentPointer += readFromFile(fin, normtype);

        int32_t chrIdx = readInt32FromFile(fin);
        currentPointer += 4;
        string unit1;
        currentPointer += readFromFile(fin, unit1);

        int32_t resolution1 = readInt32FromFile(fin);
        int64_t filePosition = readInt64FromFile(fin);
        currentPointer += 12;

        int64_t sizeInBytes;
        if (version > 8) {
            sizeInBytes = readInt64FromFile(fin);
            currentPointer += 8;
        } else {
            sizeInBytes = (int64_t)readInt32FromFile(fin);
            currentPointer += 4;
        }

        if (chrIdx == c1 && normtype == norm && unit1 == unit && resolution1 == resolution) {
            c1NormEntry.position = filePosition;
            c1NormEntry.size = sizeInBytes;
            found1 = true;
        }
        if (chrIdx == c2 && normtype == norm && unit1 == unit && resolution1 == resolution) {
            c2NormEntry.position = filePosition;
            c2NormEntry.size = sizeInBytes;
            found2 = true;
        }
    }
    if (!found1 || !found2) {
        throw std::runtime_error("Remote file did not contain " + norm +
                                 " normalization vectors for one or both chromosomes at " +
                                 std::to_string(resolution) + " " + unit);
    }
}

void readFooter(istream &fin, int64_t master, int32_t version, int32_t c1, int32_t c2,
                const string &matrixType, const string &norm, const string &unit,
                int32_t resolution, int64_t &myFilePos, indexEntry &c1NormEntry,
                indexEntry &c2NormEntry, vector<double> &expectedValues) {
    if (version > 8) {
        int64_t nBytes = readInt64FromFile(fin);
    } else {
        int32_t nBytes = readInt32FromFile(fin);
    }

    stringstream ss;
    ss << c1 << "_" << c2;
    string key = ss.str();

    int32_t nEntries = readInt32FromFile(fin);
    bool found = false;
    for (int i = 0; i < nEntries; i++) {
        string keyStr;
        getline(fin, keyStr, '\0');
        int64_t fpos = readInt64FromFile(fin);
        int32_t sizeinbytes = readInt32FromFile(fin);
        if (keyStr == key) {
            myFilePos = fpos;
            found = true;
        }
    }
    if (!found) {
        throw std::runtime_error("File doesn't have the given chr_chr map " + key);
    }

    if ((matrixType == "observed" && norm == "NONE") ||
        ((matrixType == "oe" || matrixType == "expected") && norm == "NONE" && c1 != c2))
        return;  // no need to read norm vector index

    // read in and ignore expected value maps; don't store; reading these to
    // get to norm vector index
    int32_t nExpectedValues = readInt32FromFile(fin);
    for (int i = 0; i < nExpectedValues; i++) {
        string unit0;
        getline(fin, unit0, '\0');  // unit
        int32_t binSize = readInt32FromFile(fin);

        int64_t nValues;
        if (version > 8) {
            nValues = readInt64FromFile(fin);
        } else {
            nValues = (int64_t)readInt32FromFile(fin);
        }

        bool store = c1 == c2 && (matrixType == "oe" || matrixType == "expected") &&
                     norm == "NONE" && unit0 == unit && binSize == resolution;
        readThroughExpectedVector(version, fin, expectedValues, nValues, store, resolution);
        readThroughNormalizationFactors(fin, version, store, expectedValues, c1);
    }

    if (c1 == c2 && (matrixType == "oe" || matrixType == "expected") && norm == "NONE") {
        if (expectedValues.empty()) {
            throw std::runtime_error("File did not contain expected values vectors at " +
                                     std::to_string(resolution) + " " + unit);
        }
        return;
    }

    nExpectedValues = readInt32FromFile(fin);
    for (int i = 0; i < nExpectedValues; i++) {
        string nType, unit0;
        getline(fin, nType, '\0');  // typeString
        getline(fin, unit0, '\0');  // unit
        int32_t binSize = readInt32FromFile(fin);

        int64_t nValues;
        if (version > 8) {
            nValues = readInt64FromFile(fin);
        } else {
            nValues = (int64_t)readInt32FromFile(fin);
        }
        bool store = c1 == c2 && (matrixType == "oe" || matrixType == "expected") &&
                     nType == norm && unit0 == unit && binSize == resolution;
        readThroughExpectedVector(version, fin, expectedValues, nValues, store, resolution);
        readThroughNormalizationFactors(fin, version, store, expectedValues, c1);
    }

    if (c1 == c2 && (matrixType == "oe" || matrixType == "expected") && norm != "NONE") {
        if (expectedValues.empty()) {
            throw std::runtime_error("File did not contain normalized expected values vectors at " +
                                     std::to_string(resolution) + " " + unit);
        }
    }

    // Index of normalization vectors
    nEntries = readInt32FromFile(fin);
    bool found1 = false;
    bool found2 = false;
    for (int i = 0; i < nEntries; i++) {
        string normtype;
        getline(fin, normtype, '\0');  // normalization type
        int32_t chrIdx = readInt32FromFile(fin);
        string unit1;
        getline(fin, unit1, '\0');  // unit
        int32_t resolution1 = readInt32FromFile(fin);
        int64_t filePosition = readInt64FromFile(fin);
        int64_t sizeInBytes;
        if (version > 8) {
            sizeInBytes = readInt64FromFile(fin);
        } else {
            sizeInBytes = (int64_t)readInt32FromFile(fin);
        }

        if (chrIdx == c1 && normtype == norm && unit1 == unit && resolution1 == resolution) {
            c1NormEntry.position = filePosition;
            c1NormEntry.size = sizeInBytes;
            found1 = true;
        }
        if (chrIdx == c2 && normtype == norm && unit1 == unit && resolution1 == resolution) {
            c2NormEntry.position = filePosition;
            c2NormEntry.size = sizeInBytes;
            found2 = true;
        }
    }
    if (!found1 || !found2) {
        throw std::runtime_error("File did not contain " + norm +
                                 " normalization vectors for one or both chromosomes at " +
                                 std::to_string(resolution) + " " + unit);
    }
}
}  // namespace internal
