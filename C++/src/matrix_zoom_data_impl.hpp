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

#ifndef MATRIX_ZOOM_DATA_IMPL_H
#define MATRIX_ZOOM_DATA_IMPL_H

#include <zlib.h>

#include <algorithm>
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

namespace internal {

inline void setValuesForMZD(std::istream &fin, const std::string &myunit, float &mySumCounts,
                            std::int32_t &mybinsize, std::int32_t &myBlockBinCount,
                            std::int32_t &myBlockColumnCount, bool &found) {
    std::string unit;
    getline(fin, unit, '\0');                  // unit
    readInt32FromFile(fin);                    // Old "zoom" index -- not used
    float sumCounts = readFloatFromFile(fin);  // sumCounts
    readFloatFromFile(fin);                    // occupiedCellCount
    readFloatFromFile(fin);                    // stdDev
    readFloatFromFile(fin);                    // percent95
    std::int32_t binSize = readInt32FromFile(fin);
    std::int32_t blockBinCount = readInt32FromFile(fin);
    std::int32_t blockColumnCount = readInt32FromFile(fin);
    found = false;
    if (myunit == unit && mybinsize == binSize) {
        mySumCounts = sumCounts;
        myBlockBinCount = blockBinCount;
        myBlockColumnCount = blockColumnCount;
        found = true;
    }
}

inline indexEntry readIndexEntry(std::istream &fin) {
    std::int64_t filePosition = readInt64FromFile(fin);
    std::int32_t blockSizeInBytes = readInt32FromFile(fin);
    indexEntry entry = indexEntry();
    entry.size = (std::int64_t)blockSizeInBytes;
    entry.position = filePosition;
    return entry;
}

inline void populateBlockMap(std::istream &fin, std::int32_t nBlocks,
                             std::map<std::int32_t, indexEntry> &blockMap) {
    for (int b = 0; b < nBlocks; b++) {
        std::int32_t blockNumber = readInt32FromFile(fin);
        blockMap[blockNumber] = readIndexEntry(fin);
    }
}

// reads the raw binned contact matrix at specified resolution, setting the block bin count and
// block column count
inline std::map<std::int32_t, indexEntry> readMatrixZoomData(
    std::istream &fin, const std::string &myunit, std::int32_t mybinsize, float &mySumCounts,
    std::int32_t &myBlockBinCount, std::int32_t &myBlockColumnCount, bool &found) {
    std::map<std::int32_t, indexEntry> blockMap;
    setValuesForMZD(fin, myunit, mySumCounts, mybinsize, myBlockBinCount, myBlockColumnCount,
                    found);

    std::int32_t nBlocks = readInt32FromFile(fin);
    if (found) {
        populateBlockMap(fin, nBlocks, blockMap);
    } else {
        fin.seekg(nBlocks * (sizeof(std::int32_t) + sizeof(std::int64_t) + sizeof(std::int32_t)),
                  std::ios::cur);
    }
    return blockMap;
}

// reads the raw binned contact matrix at specified resolution, setting the block bin count and
// block column count
inline std::map<std::int32_t, indexEntry> readMatrixZoomDataHttp(
    CURL_ptr &curl, std::int64_t &myFilePosition, const std::string &myunit, std::int32_t mybinsize,
    float &mySumCounts, std::int32_t &myBlockBinCount, std::int32_t &myBlockColumnCount,
    bool &found) {
    std::map<std::int32_t, indexEntry> blockMap;
    std::int32_t header_size = 5 * sizeof(std::int32_t) + 4 * sizeof(float);
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
    std::int32_t nBlocks = readInt32FromFile(fin);

    if (found) {
        std::int32_t chunkSize =
            nBlocks * (sizeof(std::int32_t) + sizeof(std::int64_t) + sizeof(std::int32_t));
        buffer = getData(curl, myFilePosition + header_size, chunkSize);
        fin = memstream(buffer);
        populateBlockMap(fin, nBlocks, blockMap);
    } else {
        myFilePosition =
            myFilePosition + header_size +
            (nBlocks * (sizeof(std::int32_t) + sizeof(std::int64_t) + sizeof(std::int32_t)));
    }
    return blockMap;
}

// goes to the specified file pointer in http and finds the raw contact matrixType at specified
// resolution, calling readMatrixZoomData. sets blockbincount and blockcolumncount
inline std::map<std::int32_t, indexEntry> readMatrixHttp(
    CURL_ptr &curl, std::int64_t myFilePosition, const std::string &unit, std::int32_t resolution,
    float &mySumCounts, std::int32_t &myBlockBinCount, std::int32_t &myBlockColumnCount) {
    std::int32_t size = sizeof(std::int32_t) * 3;
    auto buffer = getData(curl, myFilePosition, size);
    memstream bufin(buffer);

    std::int32_t c1 = readInt32FromFile(bufin);
    std::int32_t c2 = readInt32FromFile(bufin);
    std::int32_t nRes = readInt32FromFile(bufin);
    std::int32_t i = 0;
    bool found = false;
    myFilePosition = myFilePosition + size;
    std::map<std::int32_t, indexEntry> blockMap;

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
inline std::map<std::int32_t, indexEntry> readMatrix(std::istream &fin, std::int64_t myFilePosition,
                                                     const std::string &unit,
                                                     std::int32_t resolution, float &mySumCounts,
                                                     std::int32_t &myBlockBinCount,
                                                     std::int32_t &myBlockColumnCount) {
    std::map<std::int32_t, indexEntry> blockMap;

    fin.seekg(myFilePosition, std::ios::beg);
    std::int32_t c1 = readInt32FromFile(fin);
    std::int32_t c2 = readInt32FromFile(fin);
    std::int32_t nRes = readInt32FromFile(fin);
    std::int32_t i = 0;
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
inline std::vector<double> readNormalizationVector(std::istream &bufferin, std::int32_t version) {
    std::int64_t nValues;
    if (version > 8) {
        nValues = readInt64FromFile(bufferin);
    } else {
        nValues = (std::int64_t)readInt32FromFile(bufferin);
    }

    std::uint64_t numValues;
    numValues = static_cast<std::uint64_t>(nValues);
    std::vector<double> values(numValues);

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
inline std::set<std::int32_t> getBlockNumbersForRegionFromBinPosition(
    const std::int64_t *regionIndices, std::int32_t blockBinCount, std::int32_t blockColumnCount,
    bool intra) {
    std::int32_t col1, col2, row1, row2;
    col1 = static_cast<std::int32_t>(regionIndices[0] / blockBinCount);
    col2 = static_cast<std::int32_t>((regionIndices[1] + 1) / blockBinCount);
    row1 = static_cast<std::int32_t>(regionIndices[2] / blockBinCount);
    row2 = static_cast<std::int32_t>((regionIndices[3] + 1) / blockBinCount);

    std::set<std::int32_t> blocksSet;
    // first check the upper triangular matrixType
    for (int r = row1; r <= row2; r++) {
        for (int c = col1; c <= col2; c++) {
            std::int32_t blockNumber = r * blockColumnCount + c;
            blocksSet.insert(blockNumber);
        }
    }
    // check region part that overlaps with lower left triangle but only if intrachromosomal
    if (intra) {
        for (int r = col1; r <= col2; r++) {
            for (int c = row1; c <= row2; c++) {
                std::int32_t blockNumber = r * blockColumnCount + c;
                blocksSet.insert(blockNumber);
            }
        }
    }
    return blocksSet;
}

inline std::set<std::int32_t> getBlockNumbersForRegionFromBinPositionV9Intra(
    std::int64_t *regionIndices, std::int32_t blockBinCount, std::int32_t blockColumnCount) {
    // regionIndices is binX1 binX2 binY1 binY2
    std::set<std::int32_t> blocksSet;
    std::int32_t translatedLowerPAD, translatedHigherPAD, translatedNearerDepth,
        translatedFurtherDepth;
    translatedLowerPAD =
        static_cast<std::int32_t>((regionIndices[0] + regionIndices[2]) / 2 / blockBinCount);
    translatedHigherPAD =
        static_cast<std::int32_t>((regionIndices[1] + regionIndices[3]) / 2 / blockBinCount + 1);
    translatedNearerDepth = static_cast<std::int32_t>(std::log2(
        1 + std::abs(regionIndices[0] - regionIndices[3]) / std::sqrt(2) / blockBinCount));
    translatedFurtherDepth = static_cast<std::int32_t>(std::log2(
        1 + std::abs(regionIndices[1] - regionIndices[2]) / std::sqrt(2) / blockBinCount));

    // because code above assume above diagonal; but we could be below diagonal
    std::int32_t nearerDepth = std::min(translatedNearerDepth, translatedFurtherDepth);
    if ((regionIndices[0] > regionIndices[3] && regionIndices[1] < regionIndices[2]) ||
        (regionIndices[1] > regionIndices[2] && regionIndices[0] < regionIndices[3])) {
        nearerDepth = 0;
    }
    std::int32_t furtherDepth = std::max(translatedNearerDepth, translatedFurtherDepth) +
                                1;  // +1; integer divide rounds down

    for (int depth = nearerDepth; depth <= furtherDepth; depth++) {
        for (int pad = translatedLowerPAD; pad <= translatedHigherPAD; pad++) {
            std::int32_t blockNumber = depth * blockColumnCount + pad;
            blocksSet.insert(blockNumber);
        }
    }

    return blocksSet;
}

inline void appendRecord(std::vector<contactRecord> &vector, std::int32_t index, std::int32_t binX,
                         std::int32_t binY, float counts) {
    contactRecord record = contactRecord();
    record.binX = binX;
    record.binY = binY;
    record.counts = counts;
    vector[index] = record;
}

inline std::string decompressBlock(indexEntry idx, const std::string &compressedBytes) {
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
inline std::vector<contactRecord> readBlock(const std::string &fileName, indexEntry idx,
                                            std::int32_t version) {
    if (idx.size <= 0) {
        std::vector<contactRecord> v;
        return v;
    }

    const auto compressedBytes = readCompressedBytesFromFile(fileName, idx);
    auto uncompressedBytes = decompressBlock(idx, compressedBytes);

    // create stream from buffer for ease of use
    memstream bufferin(uncompressedBytes);
    std::uint64_t nRecords;
    nRecords = static_cast<std::uint64_t>(readInt32FromFile(bufferin));
    std::vector<contactRecord> v(nRecords);
    // different versions have different specific formats
    if (version < 7) {
        for (uInt i = 0; i < nRecords; i++) {
            std::int32_t binX = readInt32FromFile(bufferin);
            std::int32_t binY = readInt32FromFile(bufferin);
            float counts = readFloatFromFile(bufferin);
            appendRecord(v, i, binX, binY, counts);
        }
    } else {
        std::int32_t binXOffset = readInt32FromFile(bufferin);
        std::int32_t binYOffset = readInt32FromFile(bufferin);
        bool useShort = readCharFromFile(bufferin) == 0;  // yes this is opposite of usual

        bool useShortBinX = true;
        bool useShortBinY = true;
        if (version > 8) {
            useShortBinX = readCharFromFile(bufferin) == 0;
            useShortBinY = readCharFromFile(bufferin) == 0;
        }

        char type = readCharFromFile(bufferin);
        std::int32_t index = 0;
        if (type == 1) {
            if (useShortBinX && useShortBinY) {
                std::int16_t rowCount = readInt16FromFile(bufferin);
                for (int i = 0; i < rowCount; i++) {
                    std::int32_t binY = binYOffset + readInt16FromFile(bufferin);
                    std::int16_t colCount = readInt16FromFile(bufferin);
                    for (int j = 0; j < colCount; j++) {
                        std::int32_t binX = binXOffset + readInt16FromFile(bufferin);
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
                std::int32_t rowCount = readInt32FromFile(bufferin);
                for (int i = 0; i < rowCount; i++) {
                    std::int32_t binY = binYOffset + readInt32FromFile(bufferin);
                    std::int16_t colCount = readInt16FromFile(bufferin);
                    for (int j = 0; j < colCount; j++) {
                        std::int32_t binX = binXOffset + readInt16FromFile(bufferin);
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
                std::int16_t rowCount = readInt16FromFile(bufferin);
                for (int i = 0; i < rowCount; i++) {
                    std::int32_t binY = binYOffset + readInt16FromFile(bufferin);
                    std::int32_t colCount = readInt32FromFile(bufferin);
                    for (int j = 0; j < colCount; j++) {
                        std::int32_t binX = binXOffset + readInt32FromFile(bufferin);
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
                std::int32_t rowCount = readInt32FromFile(bufferin);
                for (int i = 0; i < rowCount; i++) {
                    std::int32_t binY = binYOffset + readInt32FromFile(bufferin);
                    std::int32_t colCount = readInt32FromFile(bufferin);
                    for (int j = 0; j < colCount; j++) {
                        std::int32_t binX = binXOffset + readInt32FromFile(bufferin);
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
            std::int32_t nPts = readInt32FromFile(bufferin);
            std::int16_t w = readInt16FromFile(bufferin);

            for (int i = 0; i < nPts; i++) {
                // std::int32_t idx = (p.y - binOffset2) * w + (p.x - binOffset1);
                std::int32_t row = i / w;
                std::int32_t col = i - row * w;
                std::int32_t bin1 = binXOffset + col;
                std::int32_t bin2 = binYOffset + row;

                float counts;
                if (useShort) {
                    std::int16_t c = readInt16FromFile(bufferin);
                    if (c != -32768) {
                        appendRecord(v, index++, bin1, bin2, c);
                    }
                } else {
                    counts = readFloatFromFile(bufferin);
                    if (!std::isnan(counts)) {
                        appendRecord(v, index++, bin1, bin2, counts);
                    }
                }
            }
        }
    }
    return v;
}

inline std::int32_t getNumRecordsInBlock(const std::string &fileName, indexEntry idx,
                                         std::int32_t version) {
    if (idx.size <= 0) {
        return 0;
    }
    const auto compressedBytes = readCompressedBytesFromFile(fileName, idx);
    auto uncompressedBytes = decompressBlock(idx, compressedBytes);

    // create stream from buffer for ease of use
    memstream bufferin(uncompressedBytes);
    std::uint64_t nRecords;
    nRecords = static_cast<std::uint64_t>(readInt32FromFile(bufferin));
    return nRecords;
}

inline MatrixZoomData::MatrixZoomData(const chromosome &chrom1, const chromosome &chrom2,
                                      const std::string &matrixType, const std::string &norm,
                                      const std::string &unit, std::int32_t resolution,
                                      std::int32_t &version, std::int64_t &master,
                                      std::int64_t &totalFileSize, const std::string &fileName) {
    this->version = version;
    this->fileName = fileName;
    std::int32_t c01 = chrom1.index;
    std::int32_t c02 = chrom2.index;
    if (c01 <= c02) {  // default is ok
        this->c1 = c01;
        this->c2 = c02;
        this->numBins1 = static_cast<std::int32_t>(chrom1.length / resolution);
        this->numBins2 = static_cast<std::int32_t>(chrom2.length / resolution);
    } else {  // flip
        this->c1 = c02;
        this->c2 = c01;
        this->numBins1 = static_cast<std::int32_t>(chrom2.length / resolution);
        this->numBins2 = static_cast<std::int32_t>(chrom1.length / resolution);
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
        stream.fin().seekg(master, std::ios::beg);
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

inline std::vector<double> MatrixZoomData::readNormalizationVectorFromFooter(
    indexEntry cNormEntry, std::int32_t &version, const std::string &fileName) {
    auto buffer = readCompressedBytesFromFile(fileName, cNormEntry);
    memstream bufferin(buffer);
    std::vector<double> cNorm = readNormalizationVector(bufferin, version);
    return cNorm;
}

inline bool MatrixZoomData::isInRange(std::int32_t r, std::int32_t c, std::int32_t numRows,
                                      std::int32_t numCols) {
    return 0 <= r && r < numRows && 0 <= c && c < numCols;
}

inline std::set<std::int32_t> MatrixZoomData::getBlockNumbers(std::int64_t *regionIndices) const {
    if (version > 8 && isIntra) {
        return getBlockNumbersForRegionFromBinPositionV9Intra(regionIndices, blockBinCount,
                                                              blockColumnCount);
    } else {
        return getBlockNumbersForRegionFromBinPosition(regionIndices, blockBinCount,
                                                       blockColumnCount, isIntra);
    }
}

inline std::vector<double> MatrixZoomData::getNormVector(std::int32_t index) {
    if (index == c1) {
        return c1Norm;
    } else if (index == c2) {
        return c2Norm;
    }
    throw std::runtime_error("Invalid index provided: " + std::to_string(index) +
                             "\nShould be either " + std::to_string(c1) + " or " +
                             std::to_string(c2));
}

inline std::vector<double> MatrixZoomData::getExpectedValues() { return expectedValues; }

inline std::vector<contactRecord> MatrixZoomData::getRecords(std::int64_t gx0, std::int64_t gx1,
                                                             std::int64_t gy0, std::int64_t gy1) {
    std::int64_t origRegionIndices[] = {gx0, gx1, gy0, gy1};
    std::int64_t regionIndices[4];
    convertGenomeToBinPos(origRegionIndices, regionIndices, resolution);

    std::set<std::int32_t> blockNumbers = getBlockNumbers(regionIndices);
    std::vector<contactRecord> records;
    for (std::int32_t blockNumber : blockNumbers) {
        // get contacts in this block
        // cout << *it << " -- " << blockMap.size() << endl;
        // cout << blockMap[*it].size << " " <<  blockMap[*it].position << endl;
        std::vector<contactRecord> tmp_records =
            readBlock(fileName, blockMap[blockNumber], version);
        for (contactRecord rec : tmp_records) {
            std::int64_t x = rec.binX * resolution;
            std::int64_t y = rec.binY * resolution;

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
                            c / expectedValues[std::min(
                                    expectedValues.size() - 1,
                                    (std::size_t)std::floor(std::abs(y - x) / resolution))]);
                    } else {
                        c = static_cast<float>(c / avgCount);
                    }
                } else if (matrixType == "expected") {
                    if (isIntra) {
                        c = static_cast<float>(expectedValues[std::min(
                            expectedValues.size() - 1,
                            (std::size_t)std::floor(std::abs(y - x) / resolution))]);
                    } else {
                        c = static_cast<float>(avgCount);
                    }
                }

                if (!std::isnan(c) && !std::isinf(c)) {
                    contactRecord record = contactRecord();
                    record.binX = static_cast<std::int32_t>(x);
                    record.binY = static_cast<std::int32_t>(y);
                    record.counts = c;
                    records.push_back(record);
                }
            }
        }
    }
    return records;
}

inline std::vector<std::vector<float> > MatrixZoomData::getRecordsAsMatrix(std::int64_t gx0,
                                                                           std::int64_t gx1,
                                                                           std::int64_t gy0,
                                                                           std::int64_t gy1) {
    std::vector<contactRecord> records = this->getRecords(gx0, gx1, gy0, gy1);
    if (records.empty()) {
        std::vector<std::vector<float> > res =
            std::vector<std::vector<float> >(1, std::vector<float>(1, 0));
        return res;
    }

    std::int64_t origRegionIndices[] = {gx0, gx1, gy0, gy1};
    std::int64_t regionIndices[4];
    convertGenomeToBinPos(origRegionIndices, regionIndices, resolution);

    std::int64_t originR = regionIndices[0];
    std::int64_t endR = regionIndices[1];
    std::int64_t originC = regionIndices[2];
    std::int64_t endC = regionIndices[3];
    std::int32_t numRows = endR - originR + 1;
    std::int32_t numCols = endC - originC + 1;
    float matrix[numRows][numCols];
    for (std::int32_t r = 0; r < numRows; r++) {
        for (std::int32_t c = 0; c < numCols; c++) {
            matrix[r][c] = 0;
        }
    }

    for (contactRecord cr : records) {
        if (std::isnan(cr.counts) || std::isinf(cr.counts)) continue;
        std::int32_t r = cr.binX / resolution - originR;
        std::int32_t c = cr.binY / resolution - originC;
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

    std::vector<std::vector<float> > finalMatrix;
    for (std::int32_t i = 0; i < numRows; i++) {
        std::vector<float> row;
        row.reserve(numCols);
        for (std::int32_t j = 0; j < numCols; j++) {
            row.push_back(matrix[i][j]);
        }
        finalMatrix.push_back(row);
    }
    return finalMatrix;
}

inline std::int64_t MatrixZoomData::getNumberOfTotalRecords() {
    std::int64_t regionIndices[4] = {0, numBins1, 0, numBins2};
    std::set<std::int32_t> blockNumbers = getBlockNumbers(regionIndices);
    std::int64_t total = 0;
    for (std::int32_t blockNumber : blockNumbers) {
        total += getNumRecordsInBlock(fileName, blockMap[blockNumber], version);
    }
    return total;
}

inline void rollingMedian(std::vector<double> &initialValues, std::vector<double> &finalResult,
                          std::int32_t window) {
    // window is actually a ~wing-span
    if (window < 1) {
        finalResult = initialValues;
        return;
    }

    /*
    finalResult.push_back(initialValues[0]);
    std::int64_t length = initialValues.size();
    for (std::int64_t index = 1; index < length; index++) {
        std::int64_t initialIndex;
        std::int64_t finalIndex;
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

        std::vector<double> subVector = sliceVector(initialValues, initialIndex, finalIndex);
        finalResult.push_back(getMedian(subVector));
    }
    */
    finalResult = initialValues;
}

inline std::int64_t readThroughExpectedVectorURL(internal::CURL_ptr &curl,
                                                 std::int64_t currentPointer, std::int32_t version,
                                                 std::vector<double> &expectedValues,
                                                 std::int64_t nValues, bool store,
                                                 std::int32_t resolution) {
    if (store) {
        std::int32_t bufferSize = nValues * sizeof(double) + 10000;
        if (version > 8) {
            bufferSize = nValues * sizeof(float) + 10000;
        }
        auto buffer = internal::getData(curl, currentPointer, bufferSize);
        internal::memstream fin(buffer);

        std::vector<double> initialExpectedValues;
        if (version > 8) {
            populateVectorWithNumbers<float>(fin, initialExpectedValues, nValues);
        } else {
            populateVectorWithNumbers<double>(fin, initialExpectedValues, nValues);
        }

        // This seems to be copying initialValues into finalResult at the moment
        // rollingMedian(initialExpectedValues, expectedValues, window);
        // std::int32_t window = 5000000 / resolution;
        expectedValues = initialExpectedValues;
    }

    if (version > 8) {
        return nValues * sizeof(float);
    } else {
        return nValues * sizeof(double);
    }
}

inline void readThroughExpectedVector(std::int32_t version, std::istream &fin,
                                      std::vector<double> &expectedValues, std::int64_t nValues,
                                      bool store, std::int32_t resolution) {
    if (store) {
        std::vector<double> initialExpectedValues;
        if (version > 8) {
            populateVectorWithNumbers<float>(fin, initialExpectedValues, nValues);
        } else {
            populateVectorWithNumbers<double>(fin, initialExpectedValues, nValues);
        }

        // This seems to be copying initialValues into finalResult at the moment
        // std::int32_t window = 5000000 / resolution;
        // rollingMedian(initialExpectedValues, expectedValues, window);
        expectedValues = initialExpectedValues;
    } else if (nValues > 0) {
        if (version > 8) {
            fin.seekg(nValues * sizeof(float), std::ios::cur);
        } else {
            fin.seekg(nValues * sizeof(double), std::ios::cur);
        }
    }
}

inline std::int64_t readThroughNormalizationFactorsURL(
    internal::CURL_ptr &curl, std::int64_t currentPointer, std::int32_t version, bool store,
    std::vector<double> &expectedValues, std::int32_t c1, std::int32_t nNormalizationFactors) {
    if (store) {
        std::int32_t bufferSize =
            nNormalizationFactors * (sizeof(std::int32_t) + sizeof(double)) + 10000;
        if (version > 8) {
            bufferSize = nNormalizationFactors * (sizeof(std::int32_t) + sizeof(float)) + 10000;
        }
        auto buffer = internal::getData(curl, currentPointer, bufferSize);
        internal::memstream fin(buffer);

        for (int j = 0; j < nNormalizationFactors; j++) {
            std::int32_t chrIdx = readInt32FromFile(fin);
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
        return nNormalizationFactors * (sizeof(std::int32_t) + sizeof(float));
    } else {
        return nNormalizationFactors * (sizeof(std::int32_t) + sizeof(double));
    }
}

inline void readThroughNormalizationFactors(std::istream &fin, std::int32_t version, bool store,
                                            std::vector<double> &expectedValues, std::int32_t c1) {
    std::int32_t nNormalizationFactors = readInt32FromFile(fin);
    if (store) {
        for (int j = 0; j < nNormalizationFactors; j++) {
            std::int32_t chrIdx = readInt32FromFile(fin);
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
            fin.seekg(nNormalizationFactors * (sizeof(std::int32_t) + sizeof(float)),
                      std::ios::cur);
        } else {
            fin.seekg(nNormalizationFactors * (sizeof(std::int32_t) + sizeof(double)),
                      std::ios::cur);
        }
    }
}

// reads the footer from the master pointer location. takes in the chromosomes,
// norm, unit (BP or FRAG) and resolution or binsize, and sets the file
// position of the matrix and the normalization vectors for those chromosomes
// at the given normalization and resolution
inline void readFooterURL(CURL_ptr &curl, std::int64_t master, std::int32_t version,
                          std::int32_t c1, std::int32_t c2, const std::string &matrixType,
                          const std::string &norm, const std::string &unit, std::int32_t resolution,
                          std::int64_t &myFilePos, indexEntry &c1NormEntry, indexEntry &c2NormEntry,
                          std::vector<double> &expectedValues) {
    std::int64_t currentPointer = master;

    auto buffer = getData(curl, currentPointer, 100);
    memstream fin(buffer);

    if (version > 8) {
        std::int64_t nBytes = readInt64FromFile(fin);
        currentPointer += 8;
    } else {
        std::int32_t nBytes = readInt32FromFile(fin);
        currentPointer += 4;
    }

    std::stringstream ss;
    ss << c1 << "_" << c2;
    std::string key = ss.str();

    std::int32_t nEntries = readInt32FromFile(fin);
    currentPointer += 4;

    std::int32_t bufferSize0 = nEntries * 50;
    buffer = getData(curl, currentPointer, bufferSize0);
    fin = memstream(buffer);

    bool found = false;
    std::string keyStr;
    for (int i = 0; i < nEntries; i++) {
        currentPointer += readFromFile(fin, keyStr);
        std::int64_t fpos = readInt64FromFile(fin);
        std::int32_t sizeinbytes = readInt32FromFile(fin);
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

    std::int32_t nExpectedValues = readInt32FromFile(fin);
    currentPointer += 4;
    std::string unit0;
    for (int i = 0; i < nExpectedValues; i++) {
        buffer = getData(curl, currentPointer, 1000);
        fin = memstream(buffer);

        currentPointer += readFromFile(fin, unit0);

        std::int32_t binSize = readInt32FromFile(fin);
        currentPointer += 4;

        std::int64_t nValues;
        if (version > 8) {
            nValues = readInt64FromFile(fin);
            currentPointer += 8;
        } else {
            nValues = (std::int64_t)readInt32FromFile(fin);
            currentPointer += 4;
        }

        bool store = c1 == c2 && (matrixType == "oe" || matrixType == "expected") &&
                     norm == "NONE" && unit0 == unit && binSize == resolution;

        currentPointer += readThroughExpectedVectorURL(curl, currentPointer, version,
                                                       expectedValues, nValues, store, resolution);

        buffer = getData(curl, currentPointer, 100);
        fin = memstream(buffer);
        std::int32_t nNormalizationFactors = readInt32FromFile(fin);
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
    std::string nType;
    for (int i = 0; i < nExpectedValues; i++) {
        buffer = getData(curl, currentPointer, 1000);
        fin = memstream(buffer);

        currentPointer += readFromFile(fin, nType);
        currentPointer += readFromFile(fin, unit0);

        std::int32_t binSize = readInt32FromFile(fin);
        currentPointer += 4;

        std::int64_t nValues;
        if (version > 8) {
            nValues = readInt64FromFile(fin);
            currentPointer += 8;
        } else {
            nValues = (std::int64_t)readInt32FromFile(fin);
            currentPointer += 4;
        }
        bool store = c1 == c2 && (matrixType == "oe" || matrixType == "expected") &&
                     nType == norm && unit0 == unit && binSize == resolution;

        currentPointer += readThroughExpectedVectorURL(curl, currentPointer, version,
                                                       expectedValues, nValues, store, resolution);

        buffer = getData(curl, currentPointer, 100);
        fin = memstream(buffer);
        std::int32_t nNormalizationFactors = readInt32FromFile(fin);
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
    std::int32_t bufferSize2 = nEntries * 60;
    buffer = getData(curl, currentPointer, bufferSize2);
    fin = memstream(buffer);

    std::string normtype;
    for (int i = 0; i < nEntries; i++) {
        currentPointer += readFromFile(fin, normtype);

        std::int32_t chrIdx = readInt32FromFile(fin);
        currentPointer += 4;
        std::string unit1;
        currentPointer += readFromFile(fin, unit1);

        std::int32_t resolution1 = readInt32FromFile(fin);
        std::int64_t filePosition = readInt64FromFile(fin);
        currentPointer += 12;

        std::int64_t sizeInBytes;
        if (version > 8) {
            sizeInBytes = readInt64FromFile(fin);
            currentPointer += 8;
        } else {
            sizeInBytes = (std::int64_t)readInt32FromFile(fin);
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

inline void readFooter(std::istream &fin, std::int64_t master, std::int32_t version,
                       std::int32_t c1, std::int32_t c2, const std::string &matrixType,
                       const std::string &norm, const std::string &unit, std::int32_t resolution,
                       std::int64_t &myFilePos, indexEntry &c1NormEntry, indexEntry &c2NormEntry,
                       std::vector<double> &expectedValues) {
    if (version > 8) {
        std::int64_t nBytes = readInt64FromFile(fin);
    } else {
        std::int32_t nBytes = readInt32FromFile(fin);
    }

    std::stringstream ss;
    ss << c1 << "_" << c2;
    std::string key = ss.str();

    std::int32_t nEntries = readInt32FromFile(fin);
    bool found = false;
    for (int i = 0; i < nEntries; i++) {
        std::string keyStr;
        getline(fin, keyStr, '\0');
        std::int64_t fpos = readInt64FromFile(fin);
        std::int32_t sizeinbytes = readInt32FromFile(fin);
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
    std::int32_t nExpectedValues = readInt32FromFile(fin);
    for (int i = 0; i < nExpectedValues; i++) {
        std::string unit0;
        getline(fin, unit0, '\0');  // unit
        std::int32_t binSize = readInt32FromFile(fin);

        std::int64_t nValues;
        if (version > 8) {
            nValues = readInt64FromFile(fin);
        } else {
            nValues = (std::int64_t)readInt32FromFile(fin);
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
        std::string nType, unit0;
        getline(fin, nType, '\0');  // typeString
        getline(fin, unit0, '\0');  // unit
        std::int32_t binSize = readInt32FromFile(fin);

        std::int64_t nValues;
        if (version > 8) {
            nValues = readInt64FromFile(fin);
        } else {
            nValues = (std::int64_t)readInt32FromFile(fin);
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
        std::string normtype;
        getline(fin, normtype, '\0');  // normalization type
        std::int32_t chrIdx = readInt32FromFile(fin);
        std::string unit1;
        getline(fin, unit1, '\0');  // unit
        std::int32_t resolution1 = readInt32FromFile(fin);
        std::int64_t filePosition = readInt64FromFile(fin);
        std::int64_t sizeInBytes;
        if (version > 8) {
            sizeInBytes = readInt64FromFile(fin);
        } else {
            sizeInBytes = (std::int64_t)readInt32FromFile(fin);
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

#endif
