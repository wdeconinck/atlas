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

#include <charconv>
#include <cstring>
#include <iosfwd>
#include <limits>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/trace/CallStack.h"
#include "atlas/runtime/trace/CodeLocation.h"
#include "atlas/runtime/trace/Nesting.h"
#include "atlas/runtime/trace/StopWatch.h"
#include "atlas/runtime/trace/Timings.h"

//-----------------------------------------------------------------------------------------------------------

namespace eckit {
class Configuration;
}  // namespace eckit

//-----------------------------------------------------------------------------------------------------------

namespace atlas {
namespace runtime {
namespace trace {

//-----------------------------------------------------------------------------------------------------------

template <typename TraceTraits>
class TraceT {
public:
    using Traits   = TraceTraits;
    using Barriers = typename Traits::Barriers;
    using Tracing  = typename Traits::Tracing;
    using Labels   = std::vector<std::string>;

public:  // static methods
    static std::string report();
    static std::string report(const eckit::Configuration& config);

public:
    TraceT(const CodeLocation&);
    TraceT(const CodeLocation&, const std::string& title);
    TraceT(const CodeLocation&, const std::string& title, const Labels&);

    ~TraceT();

    bool running() const;

    void start();

    void stop();

    void pause();

    void resume();

    double elapsed() const;

protected:
    // file and title may be temporaries (e.g. from Fortran); they are only read during construction.
    TraceT(const char* file, int line, const char* title, const Labels& = {});

private:  // types
    using Identifier = Timings::Identifier;

    struct Title {
        Title() = default;
        Title(const Title& other);
        Title& operator=(const Title& other);
        ~Title();

        void format(std::string_view);
        std::string_view view() const;

    private:
        char* resize(size_t);

        static constexpr size_t inline_capacity = 63;
        union Storage {
            char short_title[inline_capacity]{};
            char* long_title;
        } storage_;
        // Sizes above inline_capacity select the owned heap buffer.
        size_t size_{0};
    };

private:  // member functions
    void barrier() const;

    void updateTimings() const;

    void registerTimer();

    void initialise(std::string_view title, const Labels&);

private:  // member data
    bool running_{false};
    bool paused_{false};
    StopWatch stopwatch_;
    CodeLocation loc_;
    Title title_;
    Identifier id_;
    Labels labels_;
};

//-----------------------------------------------------------------------------------------------------------
// Definitions

template <typename TraceTraits>
inline TraceT<TraceTraits>::Title::Title(const Title& other) {
    const auto text = other.view();
    std::memcpy(resize(text.size()), text.data(), text.size());
}

template <typename TraceTraits>
inline auto TraceT<TraceTraits>::Title::operator=(const Title& other) -> Title& {
    Title copy(other);
    std::swap(storage_, copy.storage_);
    std::swap(size_, copy.size_);
    return *this;
}

template <typename TraceTraits>
inline TraceT<TraceTraits>::Title::~Title() {
    if (size_ > inline_capacity) {
        delete[] storage_.long_title;
    }
}

template <typename TraceTraits>
inline char* TraceT<TraceTraits>::Title::resize(size_t size) {
    Storage next;
    if (size > inline_capacity) {
        next.long_title = new char[size];
    }
    if (size_ > inline_capacity) {
        delete[] storage_.long_title;
    }
    storage_ = next;
    size_ = size;
    return size_ <= inline_capacity ? storage_.short_title : storage_.long_title;
}

template <typename TraceTraits>
inline void TraceT<TraceTraits>::Title::format(std::string_view base) {
    const std::string_view barriers = Barriers::state() ? " [b]" : "";
    std::string_view thread_prefix;
    std::string_view thread_number;
    std::string_view thread_end;
    char number[std::numeric_limits<int>::digits10 + 2];
    if (atlas_omp_get_num_threads() > 1) {
        const auto result = std::to_chars(number, number + sizeof(number), atlas_omp_get_thread_num());
        thread_prefix = " @thread[";
        thread_number = {number, static_cast<std::size_t>(result.ptr - number)};
        thread_end = "]";
    }
    const std::string_view parts[] = {base, barriers, thread_prefix, thread_number, thread_end};
    const auto size = base.size() + barriers.size() + thread_prefix.size() + thread_number.size() + thread_end.size();
    auto* out = resize(size);
    for (auto part : parts) {
        if (!part.empty()) {
            std::memcpy(out, part.data(), part.size());
            out += part.size();
        }
    }
}

template <typename TraceTraits>
inline std::string_view TraceT<TraceTraits>::Title::view() const {
    if (size_ <= inline_capacity) {
        return {storage_.short_title, size_};
    }
    return {storage_.long_title, size_};
}

template <typename TraceTraits>
inline TraceT<TraceTraits>::TraceT(const CodeLocation& loc, const std::string& title): TraceT(loc, title, Labels{}) {}

template <typename TraceTraits>
inline TraceT<TraceTraits>::TraceT(const CodeLocation& loc): TraceT(loc, loc ? loc.func() : "", Labels{}) {}

template <typename TraceTraits>
inline TraceT<TraceTraits>::TraceT(const CodeLocation& loc, const std::string& title, const Labels& labels): loc_(loc) {
    initialise(title, labels);
}

template <typename TraceTraits>
inline TraceT<TraceTraits>::TraceT(const char* file, int line, const char* title, const Labels& labels):
    loc_(CodeLocation::temporary(file, line, nullptr)) {
    initialise(title, labels);
}

template <typename TraceTraits>
inline void TraceT<TraceTraits>::initialise(std::string_view title, const Labels& labels) {
    if (Control::enabled()) {
        title_.format(title);
        labels_ = labels;
        start();
    }
    else if (loc_.isTemporary()) {
        // This trace is not recorded on this thread; drop the caller's strings rather than keep dangling pointers.
        loc_ = CodeLocation("", loc_.line(), nullptr);
    }
}

template <typename TraceTraits>
inline TraceT<TraceTraits>::~TraceT() {
    stop();
}

template <typename TraceTraits>
inline void TraceT<TraceTraits>::barrier() const {
    Barriers::execute();
}

template <typename TraceTraits>
inline void TraceT<TraceTraits>::registerTimer() {
    id_ = Timings::add(loc_, CurrentCallStack::instance().current(), title_.view(), labels_);
    if (loc_.isTemporary()) {
        // The registry keeps a persistent copy of the location
        loc_ = Timings::location(id_);
    }
}

template <typename TraceTraits>
inline void TraceT<TraceTraits>::updateTimings() const {
    Timings::update(id_, stopwatch_.elapsed());
}

template <typename TraceTraits>
inline bool TraceT<TraceTraits>::running() const {
    return running_;
}

template <typename TraceTraits>
inline void TraceT<TraceTraits>::start() {
    if (Control::enabled()) {
        ATLAS_ASSERT(not running_, "Cannot start an already running trace: " + std::string(title_.view()));
        running_ = true;
        CurrentCallStack::instance().push(loc_, title_.view());
        registerTimer();
        Tracing::start(title_.view());
        barrier();
        stopwatch_.start();
    }
}

template <typename TraceTraits>
inline void TraceT<TraceTraits>::stop() {
    if (running_) {
        barrier();
        stopwatch_.stop();
        if (not paused_) {
            CurrentCallStack::instance().pop();
        }
        updateTimings();
        Tracing::stop(title_.view(), stopwatch_.elapsed());
        running_ = false;
        paused_  = false;
    }
}

template <typename TraceTraits>
inline void TraceT<TraceTraits>::pause() {
    if (not Control::enabled()) {
        return;
    }
    ATLAS_ASSERT(running_, "Cannot pause a stopped trace: " + std::string(title_.view()));
    ATLAS_ASSERT(not paused_, "Cannot pause an already paused trace: " + std::string(title_.view()));
    barrier();
    stopwatch_.stop();
    CurrentCallStack::instance().pop();
    paused_ = true;
}

template <typename TraceTraits>
inline void TraceT<TraceTraits>::resume() {
    if (not Control::enabled()) {
        return;
    }
    ATLAS_ASSERT(running_, "Cannot resume a stopped trace: " + std::string(title_.view()));
    ATLAS_ASSERT(paused_, "Cannot resume a trace that is not paused: " + std::string(title_.view()));
    barrier();
    CurrentCallStack::instance().push(loc_, title_.view());
    stopwatch_.start();
    paused_ = false;
}

template <typename TraceTraits>
inline double TraceT<TraceTraits>::elapsed() const {
    return stopwatch_.elapsed();
}

template <typename TraceTraits>
inline std::string TraceT<TraceTraits>::report() {
    return Timings::report() + Barriers::report();
}

template <typename TraceTraits>
inline std::string TraceT<TraceTraits>::report(const eckit::Configuration& config) {
    return Timings::report(config) + Barriers::report();
}

//-----------------------------------------------------------------------------------------------------------

}  // namespace trace
}  // namespace runtime
}  // namespace atlas
