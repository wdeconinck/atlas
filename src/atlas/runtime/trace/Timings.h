/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

//-----------------------------------------------------------------------------------------------------------

namespace eckit {
class Configuration;
}
namespace atlas {
class CodeLocation;
}

namespace atlas {
namespace runtime {
namespace trace {

class CallStack;

// An unused slot is distinct from an explicitly empty label.
class Label {
public:
    Label() = default;
    Label(const char* value): Label(std::string_view(value)) {}
    Label(const std::string& value): Label(std::string_view(value)) {}
    Label(std::string_view value) {
        if (value.size() > inline_capacity) {
            auto* buffer = new char[sizeof(size_t) + value.size()];
            const size_t length = value.size();
            std::memcpy(buffer, &length, sizeof(length));
            storage_.long_value = {long_tag, buffer};
            std::memcpy(buffer + sizeof(size_t), value.data(), value.size());
        }
        else {
            storage_.short_value.size = static_cast<std::uint16_t>(value.size());
            if (!value.empty()) {
                std::memcpy(storage_.short_value.data, value.data(), value.size());
            }
        }
    }
    Label(const Label& other) {
        if (other.present()) {
            *this = Label(other.view());
        }
    }
    Label(Label&& other) noexcept { swap(other); }
    Label& operator=(Label other) noexcept {
        swap(other);
        return *this;
    }
    ~Label() {
        if (storage_.short_value.size == long_tag) {
            delete[] storage_.long_value.data;
        }
    }

    bool present() const { return storage_.short_value.size != unused; }
    std::string_view view() const {
        if (!present()) {
            return {};
        }
        if (storage_.short_value.size != long_tag) {
            return {storage_.short_value.data, storage_.short_value.size};
        }
        size_t length;
        std::memcpy(&length, storage_.long_value.data, sizeof(length));
        return {storage_.long_value.data + sizeof(size_t), length};
    }
    operator std::string_view() const { return view(); }

private:
    static constexpr size_t inline_capacity = 62;
    static constexpr std::uint16_t long_tag = inline_capacity + 1;
    static constexpr std::uint16_t unused = std::numeric_limits<std::uint16_t>::max();
    struct ShortValue {
        std::uint16_t size;
        char data[inline_capacity];
    };
    struct LongValue {
        std::uint16_t size;
        char* data;
    };
    union Storage {
        ShortValue short_value{unused, {}};
        LongValue long_value;
    } storage_;

    void swap(Label& other) noexcept {
        std::swap(storage_, other.storage_);
    }
};

static_assert(sizeof(Label) == 64, "Trace labels must occupy 64 bytes");

class Timings {
public:
    using Configuration = eckit::Configuration;
    using CodeLocation  = atlas::CodeLocation;
    using Identifier    = size_t;
    using Labels        = std::array<Label, 5>;

public:  // static methods
    static Identifier add(const CodeLocation&, const CallStack&, std::string_view title, const Labels&);

    /// Persistent copy of the location registered for id
    static const CodeLocation& location(const Identifier& id);

    static void update(const Identifier& id, double seconds);

    static std::string report();

    static std::string report(const Configuration&);
};

}  // namespace trace
}  // namespace runtime
}  // namespace atlas
