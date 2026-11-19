/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "CodeLocation.h"

#include <deque>
#include <functional>
#include <mutex>
#include <string_view>
#include <unordered_set>

namespace atlas {

namespace {

class InternedStrings {
public:
    static InternedStrings& instance() {
        // Leaked: interned strings are referenced by timings reported during static destruction.
        static InternedStrings* strings = new InternedStrings();
        return *strings;
    }

    const char* intern(std::string_view str) {
        std::lock_guard<std::mutex> lock(mutex_);
        auto it = index_.find(str);
        if (it == index_.end()) {
            it = index_.emplace(storage_.emplace_back(str)).first;
        }
        return it->data();
    }

private:
    std::mutex mutex_;
    // deque keeps element addresses stable, so views into the strings remain valid
    std::deque<std::string> storage_;
    std::unordered_set<std::string_view> index_;
};

size_t hash_string(std::string_view s, int line) {
    // http://www.cse.yorku.ca/~oz/hash.html
    size_t h = 5381;
    for (const char c : s) {
        h = ((h << 5) + h) + c;
    }
    return h ^ (std::hash<int>{}(line) << 1);
}

}  // namespace

CodeLocation CodeLocation::temporary(const char* file, int line, const char* function) {
    CodeLocation loc(file, line, function);
    loc.temporary_ = true;
    return loc;
}

CodeLocation CodeLocation::temporary(std::string_view file, int line, const char* function) {
    CodeLocation loc(file.empty() ? "" : file.data(), line, function);
    loc.temporary_        = true;
    loc.length_delimited_ = true;
    loc.file_length_      = file.size();
    return loc;
}

CodeLocation CodeLocation::interned(const char* file, int line, const char* function) {
    auto& strings = InternedStrings::instance();
    return CodeLocation(file ? strings.intern(file) : nullptr, line, function ? strings.intern(function) : nullptr);
}

CodeLocation CodeLocation::interned(std::string_view file, int line, const char* function) {
    auto& strings = InternedStrings::instance();
    return CodeLocation(strings.intern(file), line, function ? strings.intern(function) : nullptr);
}

CodeLocation CodeLocation::persistent() const {
    if (not temporary_) {
        return *this;
    }
    CodeLocation loc = length_delimited_ ? interned(std::string_view(loc_.file(), file_length_), line(), func())
                                         : interned(file(), line(), func());
    loc.hash_ = hash_;
    return loc;
}

size_t CodeLocation::hash() const {
    if (hash_ == 0) {
        const std::string_view filename = length_delimited_ ? std::string_view(file(), file_length_)
                                                            : std::string_view(file() ? file() : "");
        hash_ = hash_string(filename, line());
    }
    return hash_;
}

std::ostream& operator<<(std::ostream& s, const CodeLocation& loc) {
    s << loc.loc_;
    return s;
}

}  // namespace atlas
