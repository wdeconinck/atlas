/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Runtime.h"

#include "pluto/pluto_config.h"
#include "hic/hic.h"

namespace pluto {

std::size_t devices() {
    int num_devices = 0;
    if constexpr(PLUTO_HAVE_HIC) {
        auto err = hicGetDeviceCount(&num_devices);
        if (err) {
            num_devices = 0;
        }
    }
    return static_cast<std::size_t>(num_devices);
}

}
