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

#ifndef HIC_FILE_IMPL_H
#define HIC_FILE_IMPL_H

#include <curl/curl.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <ios>
#include <istream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "straw/internal/common.h"

inline std::int64_t HiCFile::readTotalFileSize(const std::string &url) {
    if (internal::StartsWith(url, "http")) {
        auto discardData = +[](void *buffer, std::size_t size, std::size_t nmemb,
                               void *userp) -> std::size_t { return size * nmemb; };

        auto curl = internal::CURL_ptr(curl_easy_init(), &curl_easy_cleanup);
        if (!curl.get()) {
            throw std::runtime_error("Unable to initialize curl");
        }
        curl_easy_setopt(curl.get(), CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl.get(), CURLOPT_HEADER, 1L);
        curl_easy_setopt(curl.get(), CURLOPT_NOBODY, 1L);
        curl_easy_setopt(curl.get(), CURLOPT_WRITEFUNCTION, discardData);
        curl_easy_setopt(curl.get(), CURLOPT_FOLLOWLOCATION, 1L);
        curl_easy_setopt(curl.get(), CURLOPT_USERAGENT, "straw");

        auto res = curl_easy_perform(curl.get());
        if (res != CURLE_OK) {
            throw std::runtime_error("Unable to fetch metadata for " + url + ": " +
                                     std::string(curl_easy_strerror(res)));
        }
        curl_off_t cl{};
        res = curl_easy_getinfo(curl.get(), CURLINFO_CONTENT_LENGTH_DOWNLOAD_T, &cl);
        if (res != CURLE_OK || cl == -1) {
            const auto reason =
                res == CURLE_OK ? std::string{} : ": " + std::string(curl_easy_strerror(res));
            throw std::runtime_error("Unable to fetch content length for " + url + reason);
        }
        if (cl == -1) {
            throw std::runtime_error("Unable to fetch content length for " + url);
        }
        return static_cast<std::int64_t>(cl);
    }
    return std::ifstream(url, std::ios::binary | std::ios::ate).tellg();
}

inline std::vector<std::int32_t> readResolutionsFromHeader(std::istream &fin) {
    int numBpResolutions = internal::readInt32FromFile(fin);
    std::vector<std::int32_t> resolutions;
    for (int i = 0; i < numBpResolutions; i++) {
        std::int32_t res = internal::readInt32FromFile(fin);
        resolutions.push_back(res);
    }
    return resolutions;
}

inline HiCFile::HiCFile(std::string fileName_)
    : fileName(std::move(fileName_)), totalFileSize(readTotalFileSize(fileName)) {
    // read header into buffer; 100K should be sufficient
    if (internal::StartsWith(fileName, "http")) {
        auto curl = internal::initCURL(fileName);
        auto buffer = internal::getData(curl, 0, 100000);
        internal::memstream bufin(buffer);
        chromosomes =
            readHeader(bufin, master, genomeID, numChromosomes, version, nviPosition, nviLength);
        resolutions = readResolutionsFromHeader(bufin);
    } else {
        std::ifstream fin(fileName, std::ios::in | std::ios::binary);
        if (!fin) {
            throw std::runtime_error("File " + fileName + " cannot be opened for reading");
        }
        chromosomes =
            readHeader(fin, master, genomeID, numChromosomes, version, nviPosition, nviLength);
        resolutions = readResolutionsFromHeader(fin);
    }
    assert(totalFileSize != 0);
}

inline const std::string &HiCFile::getGenomeID() const noexcept { return genomeID; }

inline const std::vector<std::int32_t> &HiCFile::getResolutions() const noexcept {
    return resolutions;
}

inline std::vector<chromosome> HiCFile::getChromosomes() const {
    std::vector<chromosome> flat_chroms(chromosomes.size());
    for (const auto &node : chromosomes) {
        flat_chroms[node.second.index] = node.second;
    }

    return flat_chroms;
}

inline auto HiCFile::getChromosomeMap() const noexcept -> const ChromosomeMap & {
    return chromosomes;
}

inline internal::MatrixZoomData HiCFile::getMatrixZoomData(
    const std::string &chr1, const std::string &chr2, const std::string &matrixType,
    const std::string &norm, const std::string &unit, std::int32_t resolution) {
    chromosome chrom1 = chromosomes[chr1];
    chromosome chrom2 = chromosomes[chr2];
    return {chrom1,     chrom2,  (matrixType), (norm),        (unit),
            resolution, version, master,       totalFileSize, fileName};
}

// reads the header, storing the positions of the normalization vectors and returning the
// masterIndexPosition pointer
inline auto HiCFile::readHeader(std::istream &fin, std::int64_t &masterIndexPosition,
                                std::string &genomeID, std::int32_t &numChromosomes,
                                std::int32_t &version, std::int64_t &nviPosition,
                                std::int64_t &nviLength) -> ChromosomeMap {
    std::map<std::string, chromosome> chromosomeMap;
    if (!internal::checkMagicString(fin)) {
        throw std::runtime_error(
            "Hi-C magic std::string is missing, does not appear to be a hic file");
    }

    version = internal::readInt32FromFile(fin);
    if (version < 6) {
        throw std::runtime_error("Version " + std::to_string(version) + " no longer supported");
    }
    masterIndexPosition = internal::readInt64FromFile(fin);
    getline(fin, genomeID, '\0');

    if (version > 8) {
        nviPosition = internal::readInt64FromFile(fin);
        nviLength = internal::readInt64FromFile(fin);
    }

    std::int32_t nattributes = internal::readInt32FromFile(fin);

    // reading and ignoring attribute-value dictionary
    for (int i = 0; i < nattributes; i++) {
        std::string key, value;
        getline(fin, key, '\0');
        getline(fin, value, '\0');
    }

    numChromosomes = internal::readInt32FromFile(fin);
    // chromosome map for finding matrixType
    for (int i = 0; i < numChromosomes; i++) {
        std::string name;
        std::int64_t length;
        getline(fin, name, '\0');
        if (version > 8) {
            length = internal::readInt64FromFile(fin);
        } else {
            length = (std::int64_t)internal::readInt32FromFile(fin);
        }

        chromosome chr;
        chr.index = i;
        chr.name = name;
        chr.length = length;
        chromosomeMap[name] = chr;
    }
    return chromosomeMap;
}

#endif
