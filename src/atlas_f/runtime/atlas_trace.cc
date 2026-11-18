/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas_trace.h"

#include <atomic>
#include <string>
#include <vector>

#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/runtime/trace/CodeLocation.h"

namespace atlas {

//static std::vector<std::string> Labels( int size, const char* labels[] ) {
//    std::vector<std::string> _labels;
//    _labels.reserve( size );
//    for ( int i = 0; i < size; ++i ) {
//        _labels.emplace_back( labels[i] );
//    }
//    return _labels;
//}

namespace {

class OwnedTrace final : public Trace {
public:
    using Trace::Trace;

    void retain() { count_.fetch_add(1, std::memory_order_relaxed); }

    void returnValue() {
        // Fortran constructor results are acquired by defined assignment, not finalized.
        const auto previous = count_.exchange(0, std::memory_order_relaxed);
        ATLAS_ASSERT(previous == 1, "Cannot return a shared or unowned Atlas trace");
    }

    void release() {
        const auto previous = count_.fetch_sub(1, std::memory_order_acq_rel);
        ATLAS_ASSERT(previous > 0, "Cannot release an unowned Atlas trace");
        if (previous == 1) {
            delete this;
        }
    }

private:
    std::atomic<std::size_t> count_{1};
};

}  // namespace

extern "C" {

Trace* new_atlas_Trace(const char* file, int line, const char* title) {
    return new OwnedTrace(CodeLocation(file, line, nullptr, true), std::string(title));
}

Trace* new_atlas_Trace_labels_1(const char* file, int line, const char* title, const char* label1) {
    std::vector<std::string> labels{label1};
    return new OwnedTrace(CodeLocation(file, line, nullptr, true), std::string(title), labels);
}
Trace* new_atlas_Trace_labels_2(const char* file, int line, const char* title, const char* label1, const char* label2) {
    std::vector<std::string> labels{label1, label2};
    return new OwnedTrace(CodeLocation(file, line, nullptr, true), std::string(title), labels);
}
Trace* new_atlas_Trace_labels_3(const char* file, int line, const char* title, const char* label1, const char* label2,
                                const char* label3) {
    std::vector<std::string> labels{label1, label2, label3};
    return new OwnedTrace(CodeLocation(file, line, nullptr, true), std::string(title), labels);
}
Trace* new_atlas_Trace_labels_4(const char* file, int line, const char* title, const char* label1, const char* label2,
                                const char* label3, const char* label4) {
    std::vector<std::string> labels{label1, label2, label3, label4};
    return new OwnedTrace(CodeLocation(file, line, nullptr, true), std::string(title), labels);
}
Trace* new_atlas_Trace_labels_5(const char* file, int line, const char* title, const char* label1, const char* label2,
                                const char* label3, const char* label4, const char* label5) {
    std::vector<std::string> labels{label1, label2, label3, label4, label5};
    return new OwnedTrace(CodeLocation(file, line, nullptr, true), std::string(title), labels);
}

void delete_atlas_Trace(Trace* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot delete uninitialised atlas_Trace");
    static_cast<OwnedTrace*>(This)->release();
}

void retain_atlas_Trace(Trace* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot retain uninitialised atlas_Trace");
    static_cast<OwnedTrace*>(This)->retain();
}

void return_atlas_Trace(Trace* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot return uninitialised atlas_Trace");
    static_cast<OwnedTrace*>(This)->returnValue();
}

void release_atlas_Trace(Trace* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot release uninitialised atlas_Trace");
    static_cast<OwnedTrace*>(This)->release();
}

void atlas_Trace__start(Trace* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot start uninitialised atlas_Trace");
    This->start();
}

void atlas_Trace__stop(Trace* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot stop uninitialised atlas_Trace");
    This->stop();
}

void atlas_Trace__pause(Trace* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot pause uninitialised atlas_Trace");
    This->pause();
}

void atlas_Trace__resume(Trace* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot resume uninitialised atlas_Trace");
    This->resume();
}

int atlas_Trace__running(Trace* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot check 'running' status of uninitialised atlas_Trace");
    return This->running();
}

double atlas_Trace__elapsed(Trace* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot check elapsed time of uninitialised atlas_Trace");
    return This->elapsed();
}

}  // extern C

}  // namespace atlas
