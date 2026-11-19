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

#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"

namespace atlas {

namespace {

class OwnedTrace final : public Trace {
public:
    OwnedTrace(std::string_view file, int line, std::string_view title, const Labels& labels = {}):
        Trace(file, line, title, labels) {}

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

std::string_view view(const char* data, size_t length) {
    return length ? std::string_view(data, length) : std::string_view("");
}

template <typename... LabelViews>
Trace* new_labelled_trace(std::string_view file, int line, std::string_view title, LabelViews... labels) {
    Trace::Labels owned_labels{};
    if (runtime::trace::Control::enabled()) {
        size_t index = 0;
        ((owned_labels[index++] = runtime::trace::Label(labels)), ...);
    }
    return new OwnedTrace(file, line, title, owned_labels);
}

}  // namespace

extern "C" {

Trace* new_atlas_Trace_view(const char* file, size_t file_length, int line, const char* title, size_t title_length) {
    return new OwnedTrace(view(file, file_length), line, view(title, title_length));
}

Trace* new_atlas_Trace_labels_view_1(const char* file, size_t file_length, int line, const char* title,
                                     size_t title_length, const char* label1, size_t length1) {
    return new_labelled_trace(view(file, file_length), line, view(title, title_length), view(label1, length1));
}

Trace* new_atlas_Trace_labels_view_2(const char* file, size_t file_length, int line, const char* title,
                                     size_t title_length, const char* label1, size_t length1,
                                     const char* label2, size_t length2) {
    return new_labelled_trace(view(file, file_length), line, view(title, title_length), view(label1, length1),
                              view(label2, length2));
}

Trace* new_atlas_Trace_labels_view_3(const char* file, size_t file_length, int line, const char* title,
                                     size_t title_length, const char* label1, size_t length1,
                                     const char* label2, size_t length2, const char* label3, size_t length3) {
    return new_labelled_trace(view(file, file_length), line, view(title, title_length), view(label1, length1),
                              view(label2, length2), view(label3, length3));
}

Trace* new_atlas_Trace_labels_view_4(const char* file, size_t file_length, int line, const char* title,
                                     size_t title_length, const char* label1, size_t length1,
                                     const char* label2, size_t length2, const char* label3, size_t length3,
                                     const char* label4, size_t length4) {
    return new_labelled_trace(view(file, file_length), line, view(title, title_length), view(label1, length1),
                              view(label2, length2), view(label3, length3), view(label4, length4));
}

Trace* new_atlas_Trace_labels_view_5(const char* file, size_t file_length, int line, const char* title,
                                     size_t title_length, const char* label1, size_t length1,
                                     const char* label2, size_t length2, const char* label3, size_t length3,
                                     const char* label4, size_t length4, const char* label5, size_t length5) {
    return new_labelled_trace(view(file, file_length), line, view(title, title_length), view(label1, length1),
                              view(label2, length2), view(label3, length3), view(label4, length4),
                              view(label5, length5));
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
