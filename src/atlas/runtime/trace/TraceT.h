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

#include <iosfwd>
#include <string>
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

private:  // types
    using Identifier = Timings::Identifier;

private:  // member functions
    void barrier() const;

    void updateTimings() const;

    void registerTimer();

    static std::string formatTitle(const std::string&);

private:  // member data
    bool running_{false};
    bool paused_{false};
    StopWatch stopwatch_;
    CodeLocation loc_;
    std::string title_;
    Identifier id_;
    Labels labels_;
};

//-----------------------------------------------------------------------------------------------------------
// Definitions

template <typename TraceTraits>
inline std::string TraceT<TraceTraits>::formatTitle(const std::string& _title) {
    const bool barriers = Barriers::state();
    const bool threaded = atlas_omp_get_num_threads() > 1;
    std::string title;
    title.reserve(_title.size() + (barriers ? 4 : 0) + (threaded ? 16 : 0));
    title.append(_title);
    if (barriers) {
        title.append(" [b]");
    }
    if (threaded) {
        title.append(" @thread[");
        title.append(std::to_string(atlas_omp_get_thread_num()));
        title.push_back(']');
    }
    return title;
}

template <typename TraceTraits>
inline TraceT<TraceTraits>::TraceT(const CodeLocation& loc, const std::string& title): TraceT(loc, title, Labels{}) {}

template <typename TraceTraits>
inline TraceT<TraceTraits>::TraceT(const CodeLocation& loc): TraceT(loc, loc ? loc.func() : "", Labels{}) {}

template <typename TraceTraits>
inline TraceT<TraceTraits>::TraceT(const CodeLocation& loc, const std::string& title, const Labels& labels): loc_(loc) {
    if (Control::enabled()) {
        title_  = formatTitle(title);
        labels_ = labels;
        start();
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
    id_ = Timings::add(loc_, CurrentCallStack::instance().current(), title_, labels_);
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
        ATLAS_ASSERT(not running_, "Cannot start an already running trace: " + title_);
        running_ = true;
        CurrentCallStack::instance().push(loc_, title_);
        registerTimer();
        Tracing::start(title_);
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
        Tracing::stop(title_, stopwatch_.elapsed());
        running_ = false;
        paused_  = false;
    }
}

template <typename TraceTraits>
inline void TraceT<TraceTraits>::pause() {
    if (not Control::enabled()) {
        return;
    }
    ATLAS_ASSERT(running_, "Cannot pause a stopped trace: " + title_);
    ATLAS_ASSERT(not paused_, "Cannot pause an already paused trace: " + title_);
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
    ATLAS_ASSERT(running_, "Cannot resume a stopped trace: " + title_);
    ATLAS_ASSERT(paused_, "Cannot resume a trace that is not paused: " + title_);
    barrier();
    CurrentCallStack::instance().push(loc_, title_);
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
