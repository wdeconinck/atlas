/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include "eckit/testing/Setup.h"

#include "atlas/library/Library.h"
#include "eckit/runtime/Main.h"
#include "eckit/mpi/Comm.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace test {


struct AtlasFixture : public eckit::testing::Setup {

    AtlasFixture() {
        eckit::Main::instance().taskID( eckit::mpi::comm("world").rank() );
        if( eckit::Main::instance().taskID() != 0 ) Log::reset();
        atlas::Library::instance().initialise();
    }

    ~AtlasFixture() {
        atlas::Library::instance().finalise();
    }
};

//----------------------------------------------------------------------------------------------------------------------

} // end namespace test
} // end namespace atlas