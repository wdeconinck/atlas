/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <string_view>
#include <memory>

namespace atlas {
    class Field;
    class FieldSet;
    class FunctionSpace;
}

namespace atlas {

class Remap {
public:
    Remap();
    ~Remap();
    bool has_interpolation(std::string_view method, const FunctionSpace& source, const FunctionSpace& target);
    void setup_interpolation(std::string_view method, const FunctionSpace& source, const FunctionSpace& target);
    void interpolate(std::string_view method, const Field& source, Field& target, bool on_device);
    void interpolate(std::string_view method, const FieldSet& source, FieldSet& target, bool on_device);
    void clear();

private:
    struct InterpolationsMap;
    std::unique_ptr<InterpolationsMap> interpolations_;
};

} // atlas
