/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/runtime/trace/Nesting.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Config.h"

extern "C" size_t atlas_test_trace_stack_size() {
    return atlas::runtime::trace::CurrentCallStack::instance().current().size();
}

extern "C" size_t atlas_test_trace_nested_cpp() {
    atlas::Trace trace(Here(), "nested C++ trace");
    return atlas::runtime::trace::CurrentCallStack::instance().current().size();
}

extern "C" int atlas_test_trace_report_contains(const char* text) {
    return atlas::Trace::report().find(text) != std::string::npos;
}

extern "C" int atlas_test_trace_report_excludes(const char* title, const char* label) {
    auto config = atlas::util::Config("exclude", std::vector<std::string>{label});
    return atlas::Trace::report(config).find(title) == std::string::npos;
}

extern "C" int atlas_test_trace_report_contains_view(const char* text, size_t length) {
    return atlas::Trace::report().find(std::string(text, length)) != std::string::npos;
}
