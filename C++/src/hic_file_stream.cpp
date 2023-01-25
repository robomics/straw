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

#include <curl/curl.h>

#include <cassert>
#include <fstream>
#include <ios>
#include <stdexcept>
#include <string>

#include "straw/internal/common.h"
#include "straw/straw.h"

namespace internal {

HiCFileStream::HiCFileStream(const std::string &fileName)
    : fin_(initRegularFile(fileName)), curl_(initRemoteFile(fileName)) {}

bool HiCFileStream::isLocal() const noexcept { return !curl_; }

bool HiCFileStream::isRemote() const noexcept { return !isLocal(); }

CURL_ptr &HiCFileStream::curl() noexcept {
    assert(curl_);
    return curl_;
}

const CURL_ptr &HiCFileStream::curl() const noexcept {
    assert(curl_);
    return curl_;
}

std::ifstream &HiCFileStream::fin() noexcept { return fin_; }
const std::ifstream &HiCFileStream::fin() const noexcept { return fin_; }

void HiCFileStream::readCompressedBytes(indexEntry idx, std::string &buffer) {
    if (isRemote()) {
        getData(curl_, idx.position, idx.size, buffer);
    } else {
        fin_.seekg(idx.position, std::ios::beg);
        buffer.resize(idx.size);
        fin_.read(&buffer.front(), idx.size);
    }
}

std::string HiCFileStream::readCompressedBytes(indexEntry idx) {
    std::string buffer{};
    readCompressedBytes(idx, buffer);
    return buffer;
}

CURL_ptr HiCFileStream::initRemoteFile(const std::string &url) {
    try {
        if (StartsWith(url, "http")) {
            return initCURL(url);
        }
        return {nullptr, &curl_easy_cleanup};
    } catch (const std::exception &e) {
        throw std::runtime_error("URL " + url + " cannot be opened for reading: " + e.what());
    }
}

std::ifstream HiCFileStream::initRegularFile(const std::string &path) {
    try {
        if (StartsWith(path, "http")) {
            return {};
        }
        std::ifstream ifs(path, std::ios::in | std::ios::binary);
        ifs.exceptions(ifs.exceptions() | std::ios::badbit | std::ios::failbit);
        return ifs;
    } catch (const std::exception &e) {
        throw std::runtime_error("URL " + path + " cannot be opened for reading: " + e.what());
    }
}
}  // namespace internal
