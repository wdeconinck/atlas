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

#include <iosfwd>

#include "eckit/memory/SharedPtr.h"
#include "atlas/mesh/detail/MeshImpl.h"

//----------------------------------------------------------------------------------------------------------------------
// Forward declarations

namespace atlas {
namespace grid {
    class Projection;
} }

namespace atlas {
namespace util {
    class Metadata;
} }

namespace atlas {
namespace mesh {
    class Nodes;
    class HybridElements;
    typedef HybridElements Edges;
    typedef HybridElements Cells;
} }

namespace atlas {
namespace meshgenerator {
    class MeshGeneratorImpl;
} }

//----------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace mesh {

//----------------------------------------------------------------------------------------------------------------------

class Mesh {

public:

    using Implementation = detail::MeshImpl;

public:

    Mesh();
    Mesh( const Mesh& );
    Mesh( const Implementation* );

    /// @brief Construct a mesh from a Stream (serialization)
    explicit Mesh(eckit::Stream&);

    /// @brief Serialization to Stream
    void encode(eckit::Stream& s) const { return impl_->encode(s); }

    void print(std::ostream& out) const { impl_->print(out); }

    const util::Metadata& metadata() const { return impl_->metadata(); }
          util::Metadata& metadata()       { return impl_->metadata(); }

    const Nodes& nodes() const { return impl_->nodes(); }
          Nodes& nodes()       { return impl_->nodes(); }

    const Cells& cells() const { return impl_->cells(); }
          Cells& cells()       { return impl_->cells();; }

    const Edges& edges() const { return impl_->edges(); }
          Edges& edges()       { return impl_->edges(); }

    const HybridElements& facets() const { return impl_->facets(); }
          HybridElements& facets()       { return impl_->facets(); }

    const HybridElements& ridges() const { return impl_->ridges(); }
          HybridElements& ridges()       { return impl_->ridges(); }

    const HybridElements& peaks() const { return impl_->peaks(); }
          HybridElements& peaks()       { return impl_->peaks(); }

    bool generated() const { return impl_->generated(); }

    /// @brief Return the memory footprint of the mesh
    size_t footprint() const { return impl_->footprint(); }

    size_t nb_partitions() const { return impl_->nb_partitions(); }

    void cloneToDevice() const { impl_->cloneToDevice(); }

    void cloneFromDevice() const { impl_->cloneFromDevice(); }

    void syncHostDevice() const { impl_->syncHostDevice(); }

    const grid::Projection& projection() const { return impl_->projection(); }

    const Implementation* get() const { return impl_.get(); }
          Implementation* get()       { return impl_.get(); }

private:  // methods

    friend std::ostream& operator<<(std::ostream& s, const Mesh& p) {
        p.print(s);
        return s;
    }

    friend class meshgenerator::MeshGeneratorImpl;
    void setProjection(const grid::Projection& p) { impl_->setProjection(p); }

private:

    eckit::SharedPtr<Implementation> impl_;

};

//----------------------------------------------------------------------------------------------------------------------

} // namespace mesh
} // namespace atlas
