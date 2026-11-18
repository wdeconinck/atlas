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

#include <cstddef>
#include <iosfwd>
#include <string>
#include "eckit/log/CodeLocation.h"

namespace atlas {

/// Cheap, copyable source location.
///
/// The file and function strings are referenced, never owned. How long they stay valid is chosen at construction:
///
///   CodeLocation(Here()), CodeLocation(file, line, function)
///       The strings outlive the location, e.g. string literals or __FILE__.
///   CodeLocation::temporary(file, line, function)
///       The strings are only valid during the current call. Use persistent() before keeping the location.
///   CodeLocation::interned(file, line, function)
///       The strings are copied into storage that lives for the rest of the program
///       (one copy per distinct string; thread-safe).
class CodeLocation {
public:
    CodeLocation() = default;
    CodeLocation(const eckit::CodeLocation& loc): loc_(loc) {}
    CodeLocation(const char* file, int line, const char* function): loc_(file, line, function) {}

    static CodeLocation temporary(const char* file, int line, const char* function);
    static CodeLocation interned(const char* file, int line, const char* function);

    bool isTemporary() const { return temporary_; }

    /// This location, or an interned copy if it is temporary.
    CodeLocation persistent() const;

    /// Hash of file and line; computed on first use and cached.
    size_t hash() const;

    operator const eckit::CodeLocation&() const { return loc_; }

    /// conversion to bool for checking if location was set
    operator bool() const { return loc_; }

    std::string asString() const { return loc_.asString(); }
    /// accessor to line
    int line() const { return loc_.line(); }
    /// accessor to file
    const char* file() const { return loc_.file(); }
    /// accessor to function
    const char* func() const { return loc_.func(); }
    friend std::ostream& operator<<(std::ostream& s, const CodeLocation& loc);

private:
    eckit::CodeLocation loc_;
    mutable size_t hash_{0};
    bool temporary_{false};
};

}  // namespace atlas
