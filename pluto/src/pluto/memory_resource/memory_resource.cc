/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "memory_resource.h"

#include "pluto/util/Registry.h"

namespace pluto {

// --------------------------------------------------------------------------------------------------------

using MemoryResourceRegistry = Registry<std::pmr::memory_resource>;

std::pmr::memory_resource* register_resource(std::string_view name, std::pmr::memory_resource* mr) {
    return &MemoryResourceRegistry::instance().enregister(name, *mr);
}

std::pmr::memory_resource* register_resource(std::string_view name, std::unique_ptr<std::pmr::memory_resource>&& mr) {
    return &MemoryResourceRegistry::instance().enregister(name, std::move(mr));
}

void unregister_resource(std::string_view name) {
    MemoryResourceRegistry::instance().unregister(name);
}

std::pmr::memory_resource* get_registered_resource(std::string_view name) {
    return &MemoryResourceRegistry::instance().get(name);
}

bool has_registered_resource(std::string_view name) {
    return MemoryResourceRegistry::instance().has(name);
}

void unregister_resources() {
    MemoryResourceRegistry::instance().clear();
}

// --------------------------------------------------------------------------------------------------------

}
