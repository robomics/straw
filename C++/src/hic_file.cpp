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

#include <algorithm>

#include "straw/straw.h"

using namespace std;

size_t HiCFile::hdf(char *buffer, size_t size, size_t nitems, void *userdata) {
    assert(buffer);
    assert(userdata);
    size_t numbytes = size * nitems;
    std::string s(buffer, numbytes);

    auto &totalFileSize = *reinterpret_cast<std::int64_t *>(userdata);  // NOLINT
    totalFileSize = 0;

    int32_t found;
    found = static_cast<int32_t>(s.find("content-range"));
    if ((size_t)found == string::npos) {
        found = static_cast<int32_t>(s.find("Content-Range"));
    }
    if ((size_t)found != string::npos) {
        int32_t found2;
        found2 = static_cast<int32_t>(s.find('/'));
        // content-range: bytes 0-100000/891471462
        if ((size_t)found2 != string::npos) {
            string total = s.substr(found2 + 1);
            totalFileSize = stol(total);
        }
    }

    return numbytes;
}

internal::CURL_ptr HiCFile::oneTimeInitCURL(const std::string &url, std::int64_t &totalFileSize) {
    auto curl = internal::initCURL(url);
    curl_easy_setopt(curl.get(), CURLOPT_HEADERFUNCTION, hdf);
    curl_easy_setopt(curl.get(), CURLOPT_HEADERDATA,
                     reinterpret_cast<void *>(&totalFileSize));  // NOLINT
    return curl;
}

static vector<int32_t> readResolutionsFromHeader(istream &fin) {
    int numBpResolutions = internal::readInt32FromFile(fin);
    vector<int32_t> resolutions;
    for (int i = 0; i < numBpResolutions; i++) {
        int32_t res = internal::readInt32FromFile(fin);
        resolutions.push_back(res);
    }
    return resolutions;
}

HiCFile::HiCFile(const string &fileName) {
    this->fileName = fileName;

    // read header into buffer; 100K should be sufficient
    if (internal::StartsWith(fileName, "http")) {
        auto curl = oneTimeInitCURL(fileName, totalFileSize);
        auto buffer = internal::getData(curl, 0, 100000);
        internal::memstream bufin(buffer);
        chromosomes =
            readHeader(bufin, master, genomeID, numChromosomes, version, nviPosition, nviLength);
        resolutions = readResolutionsFromHeader(bufin);
    } else {
        ifstream fin;
        fin.open(fileName, fstream::in | fstream::binary);
        if (!fin) {
            throw std::runtime_error("File " + fileName + " cannot be opened for reading");
        }
        chromosomes =
            readHeader(fin, master, genomeID, numChromosomes, version, nviPosition, nviLength);
        resolutions = readResolutionsFromHeader(fin);
        fin.seekg(std::ios::ate);
        totalFileSize = fin.tellg();
    }
    assert(totalFileSize != 0);
}

const string &HiCFile::getGenomeID() const noexcept { return genomeID; }

const vector<int32_t> &HiCFile::getResolutions() const noexcept { return resolutions; }

vector<chromosome> HiCFile::getChromosomes() const {
    std::vector<chromosome> flat_chroms(chromosomes.size());
    for (const auto &node : chromosomes) {
        flat_chroms[node.second.index] = node.second;
    }

    return flat_chroms;
}

auto HiCFile::getChromosomeMap() const noexcept -> const ChromosomeMap & { return chromosomes; }

internal::MatrixZoomData HiCFile::getMatrixZoomData(const string &chr1, const string &chr2,
                                                    const string &matrixType, const string &norm,
                                                    const string &unit, int32_t resolution) {
    chromosome chrom1 = chromosomes[chr1];
    chromosome chrom2 = chromosomes[chr2];
    return {chrom1,     chrom2,  (matrixType), (norm),        (unit),
            resolution, version, master,       totalFileSize, fileName};
}

// reads the header, storing the positions of the normalization vectors and returning the
// masterIndexPosition pointer
map<string, chromosome> HiCFile::readHeader(istream &fin, int64_t &masterIndexPosition,
                                            string &genomeID, int32_t &numChromosomes,
                                            int32_t &version, int64_t &nviPosition,
                                            int64_t &nviLength) {
    map<string, chromosome> chromosomeMap;
    if (!internal::checkMagicString(fin)) {
        throw std::runtime_error("Hi-C magic string is missing, does not appear to be a hic file");
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

    int32_t nattributes = internal::readInt32FromFile(fin);

    // reading and ignoring attribute-value dictionary
    for (int i = 0; i < nattributes; i++) {
        string key, value;
        getline(fin, key, '\0');
        getline(fin, value, '\0');
    }

    numChromosomes = internal::readInt32FromFile(fin);
    // chromosome map for finding matrixType
    for (int i = 0; i < numChromosomes; i++) {
        string name;
        int64_t length;
        getline(fin, name, '\0');
        if (version > 8) {
            length = internal::readInt64FromFile(fin);
        } else {
            length = (int64_t)internal::readInt32FromFile(fin);
        }

        chromosome chr;
        chr.index = i;
        chr.name = name;
        chr.length = length;
        chromosomeMap[name] = chr;
    }
    return chromosomeMap;
}
