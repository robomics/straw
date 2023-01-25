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

#ifndef STRAWC_COMMON_HPP
#define STRAWC_COMMON_HPP

#include <curl/curl.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <ios>
#include <istream>
#include <memory>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <type_traits>
#include <vector>

namespace internal {
using CURL_ptr = std::unique_ptr<CURL, decltype(&curl_easy_cleanup)>;

inline bool StartsWith(const std::string &s, const std::string &prefix) {
    return s.find(prefix) == 0;
}

// callback for libcurl. data written to this buffer
inline size_t WriteMemoryCallback(void *contents, size_t size, size_t nmemb, void *userp) {
    assert(userp);
    assert(contents);
    size_t realsize = size * nmemb;
    auto buffer = *reinterpret_cast<std::string *>(userp);  // NOLINT
    try {
        buffer.reserve(realsize + 1);
        buffer.assign(static_cast<const char *>(contents), realsize);
        buffer.push_back('\0');

    } catch (const std::bad_alloc &e) {
        fprintf(stderr, "%s: out of memory!\n", e.what());
        return 0;
    }
    return buffer.size();
}

// get a buffer that can be used as an input stream from the URL
inline void getData(CURL_ptr &curl, int64_t position, int64_t chunksize, std::string &buffer) {
    assert(curl);
    const auto oss = std::to_string(position) + "-" + std::to_string(position + chunksize);
    curl_easy_setopt(curl.get(), CURLOPT_WRITEDATA, reinterpret_cast<void *>(&buffer));
    curl_easy_setopt(curl.get(), CURLOPT_RANGE, oss.c_str());
    CURLcode res = curl_easy_perform(curl.get());
    if (res != CURLE_OK) {
        throw std::runtime_error("curl_easy_perform() failed: " +
                                 std::string(curl_easy_strerror(res)));
    }
}

inline std::string getData(CURL_ptr &curl, int64_t position, int64_t chunksize) {
    std::string buffer{};
    getData(curl, position, chunksize, buffer);
    return buffer;
}

inline CURL_ptr initCURL(const std::string &url) {
    auto curl = CURL_ptr(curl_easy_init(), &curl_easy_cleanup);
    if (curl.get()) {
        curl_easy_setopt(curl.get(), CURLOPT_WRITEFUNCTION, WriteMemoryCallback);
        curl_easy_setopt(curl.get(), CURLOPT_URL, &url.front());
        curl_easy_setopt(curl.get(), CURLOPT_FOLLOWLOCATION, 1L);
        curl_easy_setopt(curl.get(), CURLOPT_USERAGENT, "straw");
    } else {
        throw std::runtime_error("Unable to initialize curl");
    }
    return curl;
}

// this is for creating a stream from a byte array for ease of use
// see
// https://stackoverflow.com/questions/41141175/how-to-implement-seekg-seekpos-on-an-in-memory-buffer
struct membuf : std::streambuf {
    inline membuf(char *first, std::size_t size) { setg(first, first, first + size); }
};

struct memstream : virtual membuf, std::istream {
    inline memstream(char *first, std::size_t size)
        : membuf(first, size), std::istream(static_cast<std::streambuf *>(this)) {}
    inline explicit memstream(std::string &s) : memstream(&s.front(), s.size()) {}

    inline std::istream::pos_type seekpos(std::istream::pos_type sp,
                                          std::ios_base::openmode which) override {
        return seekoff(sp - std::istream::pos_type(std::istream::off_type(0)), std::ios_base::beg,
                       which);
    }

    inline std::istream::pos_type seekoff(
        std::istream::off_type off, std::ios_base::seekdir dir,
        std::ios_base::openmode which = std::ios_base::in) override {
        if (dir == std::ios_base::cur)
            gbump(off);
        else if (dir == std::ios_base::end)
            setg(eback(), egptr() + off, egptr());
        else if (dir == std::ios_base::beg)
            setg(eback(), eback() + off, egptr());
        return gptr() - eback();
    }
};

inline bool checkMagicString(std::istream &fin) {
    std::string str;
    getline(fin, str, '\0');
    return str == "HIC";
}

template <typename T, typename std::enable_if<std::is_arithmetic<T>::value>::type * = nullptr>
inline T readFromFile(std::istream &fin) {
    T buffer{};
    fin.read(reinterpret_cast<char *>(&buffer), sizeof(T));  // NOLINT
    return buffer;
}

inline std::int64_t readFromFile(std::istream &fin, std::string &buffer, char delim = '\0') {
    getline(fin, buffer, delim);
    return buffer.size() + 1;
}

inline std::string readFromFile(std::istream &fin, char delim = '\0') {
    std::string buffer;
    readFromFile(fin, buffer, delim);
    return buffer;
}

inline char readCharFromFile(std::istream &fin) { return readFromFile<char>(fin); }

inline int16_t readInt16FromFile(std::istream &fin) { return readFromFile<int16_t>(fin); }

inline int32_t readInt32FromFile(std::istream &fin) { return readFromFile<int32_t>(fin); }

inline int64_t readInt64FromFile(std::istream &fin) { return readFromFile<int64_t>(fin); }

inline float readFloatFromFile(std::istream &fin) { return readFromFile<float>(fin); }

inline double readDoubleFromFile(std::istream &fin) { return readFromFile<double>(fin); }

template <typename N, typename std::enable_if<std::is_arithmetic<N>::value>::type * = nullptr>
inline void populateVectorWithNumbers(std::istream &fin, std::vector<double> &buffer,
                                      int64_t nValues) {
    buffer.resize(nValues);
    std::generate(buffer.begin(), buffer.end(),
                  [&]() { return static_cast<double>(readFromFile<N>(fin)); });
}

template <typename N, typename std::enable_if<std::is_arithmetic<N>::value>::type * = nullptr>
inline std::vector<double> populateVectorWithNumbers(std::istream &fin, int64_t nValues) {
    std::vector<double> buffer(nValues);
    populateVectorWithNumbers<N>(fin, buffer, nValues);
    return buffer;
}

}  // namespace internal

inline void convertGenomeToBinPos(const int64_t origRegionIndices[4], int64_t regionIndices[4],
                                  int32_t resolution) {
    for (uint16_t q = 0; q < 4; q++) {
        // used to find the blocks we need to access
        regionIndices[q] = origRegionIndices[q] / resolution;
    }
}

#endif  // STRAWC_COMMON_HPP
