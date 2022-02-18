/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#pragma once

#include "atlas/interpolation/method/knn/KNearestNeighboursBase.h"

#include <forward_list>

#include "atlas/functionspace.h"
#include "atlas/util/ConvexSphericalPolygon.h"


namespace atlas {
namespace interpolation {
namespace method {


class ConservativeMethod : public Method {

private:
    struct InterpolationParameters {
        std::vector<idx_t> tcell_id;
        std::vector<PointXYZ> centroids;
        std::vector<double> weights;
        std::vector<double> sweights;
    };

public:
    struct RemapStat {
        bool setup_computed = false;
        bool remap_computed = false;
        enum Counts { 
            SRC_PLG = 0,    // index, number of source polygons
            TGT_PLG,    // index, number of target polygons
            INT_PLG,    // index, number of intersection polygons
            UNCVR_SRC   // index, number of uncovered source polygons
        };
        std::array<int, 4> counts;
        enum Errors {
            SRC_PLG_L1 = 0,   // index, over/undershoot in source subpolygon creation
            SRC_PLG_LINF, 
            TGT_PLG_L1,   // index, over/untershoot in target subpolygon creation
            TGT_PLG_LINF, 
            GEO_L1,       // index, cumulative area mismatch in polygon intersections
            GEO_LINF,     // index, like GEO_L1 but in L_infinity norm
            GEO_DIFF,     // index, difference in earth area coverages
            REMAP_CONS,   // index, error in mass conservation
            REMAP_L2,     // index, error accuracy for given analytical function
            REMAP_LINF    // index, like REMAP_L2 but in L_infinity norm
        };
        std::array<double, 10> errors;
    };

public:
    typedef util::ConvexSphericalPolygon CSPolygon;
    typedef std::vector<std::pair<CSPolygon, int>> PolygonArray;
    typedef std::vector<std::tuple<CSPolygon, int>> CSPolygonArray;
    typedef array::ArrayView<double, 1> FieldArray;

    ConservativeMethod(const Config& = util::NoConfig());

    using Method::do_setup;
    void do_setup(const FunctionSpace& src_fs, const FunctionSpace& tgt_fs);
    void do_setup(const Grid& src_grid, const Grid& tgt_grid);
    void do_setup(const Grid& src_grid, const Grid& tgt_grid, const Cache&) { ATLAS_NOTIMPLEMENTED; }
    void do_execute(const Field& src_field, Field& tgt_field);

    void set_order(int order);
    void setup_stat() const;
    void remap_stat(const FieldArray& src_field, const FieldArray& tgt_field,
                            FieldArray* diff_field, double func(const PointLonLat&)) const;
    void print(std::ostream& out) const { out << "ConservativeMethod[]"; }

    RemapStat& remap_stat() const;
    bool src_cell_data() const { return src_cell_data_; }
    bool tgt_cell_data() const { return tgt_cell_data_; }
    const FunctionSpace& source() const { return src_fs_; }
    const FunctionSpace& target() const { return tgt_fs_; }
    Mesh src_mesh() const { return src_mesh_; }
    Mesh tgt_mesh() const { return tgt_mesh_; }
    int normalise_intersections() const { return normalise_intersections_; }
    int order() const { return order_; }
    int matrix_free() const { return matrix_free_; }
    inline const std::vector<InterpolationParameters>& iparam() const { return iparam_; }
    inline const PointXYZ& src_points(size_t id) const { return src_points_[id]; }
    inline const PointXYZ& tgt_points(size_t id) const { return tgt_points_[id]; }

protected:
    void intersect_polygons(const CSPolygonArray& src_csp, const CSPolygonArray& tgt_scp);
    void setup_1st_order_matrix();
    void setup_2nd_order_matrix();
    template <class TargetCellsIDs>
    void dump_intersection(const util::ConvexSphericalPolygon& s_csp, const CSPolygonArray& tgt_csp,
                           const TargetCellsIDs& tgt_cells) const;
    PointXYZ get_point(idx_t node, const Mesh& mesh) const; 
    PointXYZ get_point(idx_t node, const Mesh& mesh, PointLonLat& pll) const; 
    std::vector<idx_t> sort_cell_edges(Mesh& mesh, idx_t cell_id) const;
    std::vector<idx_t> sort_node_edges(Mesh& mesh, idx_t cell_id) const;
    std::vector<idx_t> get_cell_neighbours(Mesh& mesh, idx_t jcell) const;
    std::vector<idx_t> get_node_neighbours(Mesh& mesh, idx_t jcell) const;
    CSPolygonArray get_polygons_celldata(Mesh& mesh) const;
    CSPolygonArray get_polygons_nodedata(Mesh& mesh, std::vector<idx_t>& csp2node,
                                         std::vector<std::vector<idx_t>>& node2csp,
                                         std::array<double,2>& errors) const;

private:
    int next_index(int current_index, int size, int offset = 1) const; 
    int prev_index(int current_index, int size, int offset = 1) const; 

protected:
    bool src_cell_data_;
    bool tgt_cell_data_;
    FunctionSpace src_fs_;
    FunctionSpace tgt_fs_;
    Mesh src_mesh_;
    Mesh tgt_mesh_;
    int normalise_intersections_;
    int order_;
    idx_t n_spoints_;
    idx_t n_tpoints_;
    bool matrix_free_;
    std::vector<PointXYZ> src_points_;
    std::vector<PointXYZ> tgt_points_;
    Field src_areas_;
    Field tgt_areas_;
    mutable RemapStat remap_stat_;
    std::vector<InterpolationParameters> iparam_;   // TODO: remove after setup
    std::vector<idx_t> src_csp2node_;               // TODO: remove
    std::vector<idx_t> tgt_csp2node_;               // TODO: remove
    std::vector<std::vector<idx_t>> src_node2csp_;  // TODO: remove
    std::vector<std::vector<idx_t>> tgt_node2csp_;  // TODO: remove
};


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
