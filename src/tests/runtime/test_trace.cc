/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <chrono>
#include <cmath>
#include <memory>
#include <sstream>
#include <thread>

#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Config.h"
#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {

static void work() {
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
}


CASE("test elapsed") {
    auto trace = Trace(Here());

    EXPECT(trace.running());

    EXPECT(trace.elapsed() == 0.);

    work();

    trace.pause();

    EXPECT(trace.running());

    double elapsed = trace.elapsed();
    EXPECT(elapsed != 0.);
    EXPECT(trace.elapsed() == elapsed);

    trace.resume();

    work();

    EXPECT(trace.running());

    trace.stop();

    EXPECT(trace.elapsed() != elapsed);
}

CASE("test trace OpenMP") {
    atlas_omp_parallel_for(int i = 0; i < 10; ++i) {
        auto trace = Trace(Here(), "loop");
        if (ATLAS_HAVE_OMP) {
            work();

            trace.stop();
            if (atlas_omp_get_thread_num() > 0) {
                EXPECT(trace.elapsed() == 0.);
            }
            else {
                EXPECT(trace.elapsed() != 0.);
            }
        }
    }
}

CASE("test paused trace finalization and restart") {
    runtime::trace::Barriers barriers(false);
    auto& stack             = runtime::trace::CurrentCallStack::instance();
    const auto initial_size = stack.current().size();
    const auto initial_hash = stack.current().hash();
    {
        Trace trace(Here(), "paused trace");
        for (int i = 0; i < 32; ++i) {
            trace.pause();
            EXPECT(stack.current().size() == initial_size);
            EXPECT_THROWS_AS(trace.pause(), eckit::AssertionFailed);
            trace.resume();
            EXPECT(stack.current().size() == initial_size + 1);
            EXPECT_THROWS_AS(trace.resume(), eckit::AssertionFailed);
        }
        trace.pause();
        trace.stop();
        EXPECT(not trace.running());
        EXPECT_THROWS_AS(trace.resume(), eckit::AssertionFailed);
        {
            Trace parent(Here(), "restart parent");
            trace.start();
            trace.pause();
            trace.stop();
        }
        trace.start();
        trace.pause();
    }
    EXPECT(stack.current().size() == initial_size);
    EXPECT(stack.current().hash() == initial_hash);
    EXPECT_NO_THROW(Trace::report());
}

CASE("test call stack hash restoration and copying") {
    runtime::trace::CallStack stack;
    stack.push(Here(), "parent");
    const auto parent_hash = stack.hash();
    stack.push(Here(), "child");
    const auto child_hash = stack.hash();
    auto copied           = stack;
    runtime::trace::CallStack assigned;
    assigned = stack;
    stack.pop();
    EXPECT(stack.hash() == parent_hash);
    EXPECT(copied.hash() == child_hash);
    EXPECT(assigned.hash() == child_hash);
    copied.pop();
    assigned.pop();
    EXPECT(copied.hash() == parent_hash);
    EXPECT(assigned.hash() == parent_hash);
    stack.pop();
    EXPECT(stack.hash() == 0);
    EXPECT(stack.size() == 0);
}

CASE("test timing sample variance") {
    using runtime::trace::CallStack;
    using runtime::trace::Timings;
    auto loc = Here();
    CallStack stack;
    stack.push(loc, "variance regression");
    auto id = Timings::add(loc, stack, "variance regression", {});
    Timings::update(id, 1.);
    Timings::update(id, 3.);
    Timings::update(id, 2.);
    const auto report = Timings::report(util::Config("decimals", 8)("header", false));
    const auto row    = report.find("variance regression");
    EXPECT(row != std::string::npos);
    const auto column = report.find("std: ", row);
    EXPECT(column != std::string::npos);
    std::istringstream value(report.substr(column + 5));
    double deviation = 0.;
    value >> deviation;
    EXPECT(std::abs(deviation - 1.) < 1.e-8);
}

CASE("test trace title suffixes") {
    runtime::trace::Barriers barriers(true);
    {
        Trace plain(Here(), "plain suffix");
        Trace labelled(Here(), "labelled suffix", {"suffix"});
    }
    const auto report = Trace::report();
    const std::string suffix = Trace::Barriers::state() ? " [b]" : "";
    EXPECT(report.find("plain suffix" + suffix + " [b]") == std::string::npos);
    EXPECT(report.find("plain suffix" + suffix) != std::string::npos);
    EXPECT(report.find("labelled suffix" + suffix) != std::string::npos);
}

CASE("test inline and long trace titles survive restart") {
    runtime::trace::Barriers barriers(false);
    auto& stack = runtime::trace::CurrentCallStack::instance();
    const std::string inline_title = "inline title";
    const std::string full_inline_title(31, 'a');
    const std::string first_long_title(32, 'b');
    const std::string long_title = "a considerably longer trace title that needs heap storage";
    for (const auto& title : {std::string{}, inline_title, full_inline_title, first_long_title, long_title}) {
        std::string input = title;
        Trace trace(Here(), input);
        const auto hash = stack.current().hash();
        input.assign(input.size(), 'x');
        trace.pause();
        trace.resume();
        EXPECT(stack.current().hash() == hash);
        trace.stop();
        trace.start();
        EXPECT(stack.current().hash() == hash);
    }
    const auto report = Trace::report();
    EXPECT(report.find(inline_title) != std::string::npos);
    EXPECT(report.find(full_inline_title) != std::string::npos);
    EXPECT(report.find(first_long_title) != std::string::npos);
    EXPECT(report.find(long_title) != std::string::npos);
}

CASE("test trace titles at the 63 byte inline capacity") {
    runtime::trace::Barriers barriers(false);
    auto& stack = runtime::trace::CurrentCallStack::instance();
    const std::string full_inline_title(63, 'c');
    const std::string first_long_title(64, 'd');
    const std::string long_title(200, 'e');
    for (const auto& title : {full_inline_title, first_long_title, long_title}) {
        std::string input = title;
        Trace trace(Here(), input);
        const auto hash = stack.current().hash();
        input.assign(input.size(), 'x');
        trace.pause();
        trace.resume();
        EXPECT(stack.current().hash() == hash);
        trace.stop();
        trace.start();
        EXPECT(stack.current().hash() == hash);
    }
    const auto report = Trace::report();
    EXPECT(report.find(full_inline_title) != std::string::npos);
    EXPECT(report.find(first_long_title) != std::string::npos);
    EXPECT(report.find(long_title) != std::string::npos);
}

CASE("test trace does not keep temporary location strings") {
    runtime::trace::Barriers barriers(false);
    auto& stack = runtime::trace::CurrentCallStack::instance();
    for (bool temporary : {false, true}) {
        const std::string name = temporary ? "temporary_location_file.cc" : "interned_location_file.cc";
        std::string file       = "/a/temporary/path/" + name;
        std::string func       = "location_function";
        std::unique_ptr<Trace> trace;
        size_t hash;
        {
            const auto loc = temporary ? CodeLocation::temporary(file.c_str(), 4242, func.c_str())
                                       : CodeLocation::interned(file.c_str(), 4242, func.c_str());
            trace          = std::make_unique<Trace>(loc, name);
            hash           = stack.current().hash();
        }
        file.assign(file.size(), 'x');
        func.assign(func.size(), 'x');
        trace->pause();
        trace->resume();
        EXPECT(stack.current().hash() == hash);
        trace->stop();
        trace->start();
        EXPECT(stack.current().hash() == hash);
        trace.reset();
        EXPECT(Trace::report().find(name + " +4242") != std::string::npos);
    }
}

CASE("test length-delimited temporary locations") {
    std::string file = "/a/long/temporary/location.F90 trailing bytes";
    const std::string expected = "/a/long/temporary/location.F90";
    auto temporary = CodeLocation::temporary(std::string_view(file.data(), expected.size()), 731, nullptr);
    auto reference = CodeLocation::temporary(expected.c_str(), 731, nullptr);
    EXPECT(temporary.hash() == reference.hash());
    const auto cached_hash = temporary.hash();
    file[0] = 'x';
    EXPECT(temporary.hash() == cached_hash);
    file[0] = '/';
    auto persistent = temporary.persistent();
    file.assign(file.size(), 'x');
    EXPECT(std::string(persistent.file()) == expected);
    EXPECT(persistent.hash() == reference.hash());
    auto empty = CodeLocation::temporary(std::string_view{}, 731, nullptr);
    EXPECT(empty.hash() == CodeLocation::temporary("", 731, nullptr).hash());
    EXPECT(std::string(empty.persistent().file()).empty());

    std::string zero_file = "a";
    auto zero = CodeLocation::temporary(std::string_view(zero_file), 88835, nullptr);
    EXPECT(zero.hash() == 0);
    EXPECT(zero.persistent().hash() == 0);
}

CASE("test trace labels own inline and long strings") {
    using runtime::trace::Label;
    for (const auto& text : {std::string{}, std::string(62, 's'), std::string(63, 'l'),
                             std::string(65536, 'x')}) {
        Label label(std::string_view(text.data(), text.size()));
        Label copy = label;
        label = Label("replaced");
        EXPECT(copy.present());
        EXPECT(copy.view() == text);
        Label moved = std::move(copy);
        EXPECT(moved.view() == text);
        copy = moved;
        EXPECT(copy.view() == text);
    }
    EXPECT(!Label{}.present());
    EXPECT(Label("").present());
}

CASE("test copied trace titles own their storage") {
    runtime::trace::Barriers barriers(false);
    auto& stack = runtime::trace::CurrentCallStack::instance();
    for (const auto& title : {std::string{}, std::string(63, 'f'), std::string(64, 'g'), std::string(200, 'h')}) {
        std::unique_ptr<Trace> copy;
        size_t hash;
        {
            Trace source(Here(), title);
            hash = stack.current().hash();
            source.stop();
            copy = std::make_unique<Trace>(source);
        }
        copy->start();
        EXPECT(stack.current().hash() == hash);
        copy->stop();
        for (const auto& initial : {std::string("short"), std::string(200, 'i')}) {
            Trace assigned(Here(), initial);
            assigned.stop();
            assigned = *copy;
            assigned.start();
            EXPECT(stack.current().hash() == hash);
            assigned.stop();
        }
    }
}

CASE("test barrier") {
    EXPECT(runtime::trace::Barriers::state() == Library::instance().traceBarriers());
    {
        runtime::trace::Barriers set_barriers(true);
        EXPECT(runtime::trace::Barriers::state() == true);
        {
            runtime::trace::Barriers set_barriers(false);
            EXPECT(runtime::trace::Barriers::state() == false);
        }
        EXPECT(runtime::trace::Barriers::state() == true);
    }
    EXPECT(runtime::trace::Barriers::state() == Library::instance().traceBarriers());
}

// --------------------------------------------------------------------------


void haloexchange() {
    ATLAS_TRACE();
}
void wind_const() {
    ATLAS_TRACE();
}
void wind_next() {
    ATLAS_TRACE();
    haloexchange();
}
static int count = 0;
void dp_meth() {
    ATLAS_TRACE();
    if (count == 0) {
        wind_const();
    }
    else {
        wind_next();
    }
    ++count;
}
void tracer_interpolation() {
    ATLAS_TRACE();
    haloexchange();
}
void execute_sladv() {
    ATLAS_TRACE();
    dp_meth();
    tracer_interpolation();
}

CASE("test report") {
    SECTION("1") {
        for (int i = 0; i < 3; ++i) {
            execute_sladv();
        }
    }
    SECTION("2") {
        count = 0;
        for (int i = 0; i < 3; ++i) {
            execute_sladv();
        }
    }
    Log::info() << atlas::Trace::report() << std::endl;
}

// --------------------------------------------------------------------------


}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
