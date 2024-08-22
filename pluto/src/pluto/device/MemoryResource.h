/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */
#pragma once

#include <string_view>
#include <memory_resource>

namespace pluto::device {

// --------------------------------------------------------------------------------------------------------

using memory_resource = std::pmr::memory_resource;

// --------------------------------------------------------------------------------------------------------

std::pmr::memory_resource* get_default_resource();
void set_default_resource(std::pmr::memory_resource*);
void set_default_resource(std::string_view name);

// --------------------------------------------------------------------------------------------------------

class DefaultResource {
public:

    DefaultResource(std::string_view name);

    DefaultResource(std::pmr::memory_resource* mr) :
        saved_(get_default_resource()) {
        device::set_default_resource(mr);
    }

    ~DefaultResource() {
        device::set_default_resource(saved_);
    }
private:
    std::pmr::memory_resource* saved_;
};

// --------------------------------------------------------------------------------------------------------

}
