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

private:  // member functions
    void barrier() const;

    void updateTimings() const;

    void registerTimer();

    void initialise(std::string_view title, const Labels&);

    void formatTitle(std::string_view);

    std::string_view title() const;

private:  // member data
    bool running_{false};
    bool paused_{false};
    StopWatch stopwatch_;
    CodeLocation loc_;
    char short_title_[63]{};
    // A size above the inline capacity selects long_title_; neither representation borrows storage.
    unsigned char short_title_size_{0};
    std::string long_title_;
    Identifier id_;
    Labels labels_;
};

//-----------------------------------------------------------------------------------------------------------
// Definitions

template <typename TraceTraits>
inline void TraceT<TraceTraits>::formatTitle(std::string_view base) {
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
    if (size <= sizeof(short_title_)) {
        auto* out = short_title_;
        for (auto part : parts) {
            if (!part.empty()) {
                std::memcpy(out, part.data(), part.size());
                out += part.size();
            }
        }
        short_title_size_ = static_cast<unsigned char>(size);
    }
    else {
        short_title_size_ = sizeof(short_title_) + 1;
        long_title_.reserve(size);
        for (auto part : parts) {
            long_title_.append(part);
        }
    }
}

template <typename TraceTraits>
inline std::string_view TraceT<TraceTraits>::title() const {
    if (short_title_size_ <= sizeof(short_title_)) {
        return {short_title_, short_title_size_};
    }
    return long_title_;
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
        formatTitle(title);
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
    id_ = Timings::add(loc_, CurrentCallStack::instance().current(), title(), labels_);
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
        ATLAS_ASSERT(not running_, "Cannot start an already running trace: " + std::string(title()));
        running_ = true;
        CurrentCallStack::instance().push(loc_, title());
        registerTimer();
        Tracing::start(title());
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
        Tracing::stop(title(), stopwatch_.elapsed());
        running_ = false;
        paused_  = false;
    }
}

template <typename TraceTraits>
inline void TraceT<TraceTraits>::pause() {
    if (not Control::enabled()) {
        return;
    }
    ATLAS_ASSERT(running_, "Cannot pause a stopped trace: " + std::string(title()));
    ATLAS_ASSERT(not paused_, "Cannot pause an already paused trace: " + std::string(title()));
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
    ATLAS_ASSERT(running_, "Cannot resume a stopped trace: " + std::string(title()));
    ATLAS_ASSERT(paused_, "Cannot resume a trace that is not paused: " + std::string(title()));
    barrier();
    CurrentCallStack::instance().push(loc_, title());
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
