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
public:
    struct InterpolationParameters {      // one polygon intersection
        std::vector<idx_t> cell_idx;      // target cells used for intersection
        std::vector<PointXYZ> centroids;  // intersection cell centroids
        std::vector<double> src_weights;  // intersection cell areas
                                          // TODO: tgt_weights can be computed on the fly
        std::vector<double> tgt_weights;  // (intersection cell areas) / (target cell area)
    };

    class CachableData : public InterpolationCacheEntry {
    public:
        ~CachableData() override = default;
        size_t footprint() const override;
        static std::string static_type() { return "ConservativeMethod"; }
        std::string type() const override { return static_type(); }

        void print(std::ostream& out) const;

    private:
        friend class ConservativeMethod;

        // position and effective area of data points
        std::vector<PointXYZ> src_points_;
        std::vector<PointXYZ> tgt_points_;
        std::vector<double> src_areas_;
        std::vector<double> tgt_areas_;

        // indexing of subpolygons
        std::vector<idx_t> src_csp2node_;
        std::vector<idx_t> tgt_csp2node_;
        std::vector<std::vector<idx_t>> src_node2csp_;
        std::vector<std::vector<idx_t>> tgt_node2csp_;

        std::vector<InterpolationParameters> src_iparam_;  // TODO: remove after setup


        // Reconstructible if need be
        FunctionSpace src_fs_;
        FunctionSpace tgt_fs_;
    };

public:
    class Cache final : public interpolation::Cache {
    public:
        Cache() = default;
        Cache(const interpolation::Cache& c);
        //        Cache(const Interpolation&);

        operator bool() const { return entry_; }
        size_t footprint() const;
        const CachableData* get() const { return entry_; }

    private:
        friend class ConservativeMethod;
        Cache(std::shared_ptr<InterpolationCacheEntry> entry);
        const CachableData* entry_{nullptr};
    };

    struct RemapStat {
        bool setup_computed = false;
        bool remap_computed = false;
        enum Counts
        {
            SRC_PLG = 0,  // index, number of source polygons
            TGT_PLG,      // index, number of target polygons
            INT_PLG,      // index, number of intersection polygons
            UNCVR_SRC     // index, number of uncovered source polygons
        };
        std::array<int, 4> counts;
        enum Errors
        {
            SRC_PLG_L1 = 0,  // index, over/undershoot in source subpolygon creation
            SRC_PLG_LINF,
            TGT_PLG_L1,  // index, over/untershoot in target subpolygon creation
            TGT_PLG_LINF,
            GEO_L1,      // index, cumulative area mismatch in polygon intersections
            GEO_LINF,    // index, like GEO_L1 but in L_infinity norm
            GEO_DIFF,    // index, difference in earth area coverages
            REMAP_CONS,  // index, error in mass conservation
            REMAP_L2,    // index, error accuracy for given analytical function
            REMAP_LINF   // index, like REMAP_L2 but in L_infinity norm
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
    void do_setup(const FunctionSpace& src_fs, const FunctionSpace& tgt_fs) override;
    void do_setup_impl(const Grid& src_grid, const Grid& tgt_grid);
    void do_setup(const Grid& src_grid, const Grid& tgt_grid, const interpolation::Cache&) override;
    void do_execute(const Field& src_field, Field& tgt_field) const override;

    void set_order(int order);
    void setup_stat() const;
    void remap_stat(const FieldArray& src_field, const FieldArray& tgt_field, FieldArray* diff_field,
                    double func(const PointLonLat&)) const;
    void print(std::ostream& out) const override { out << "ConservativeMethod[]"; }

    RemapStat& remap_stat() const;
    bool src_cell_data() const { return src_cell_data_; }
    bool tgt_cell_data() const { return tgt_cell_data_; }
    const FunctionSpace& source() const override { return cachable_data_->src_fs_; }
    const FunctionSpace& target() const override { return cachable_data_->tgt_fs_; }
    Mesh src_mesh() const { return src_mesh_; }
    Mesh tgt_mesh() const { return tgt_mesh_; }
    int normalise_intersections() const { return normalise_intersections_; }
    int order() const { return order_; }
    int matrix_free() const { return matrix_free_; }
    inline const std::vector<InterpolationParameters>& iparam() const { return cachable_data_->src_iparam_; }
    inline const PointXYZ& src_points(size_t id) const { return cachable_data_->src_points_[id]; }
    inline const PointXYZ& tgt_points(size_t id) const { return cachable_data_->tgt_points_[id]; }

    interpolation::Cache createCache() const override {
        interpolation::Cache cache;
        if (not matrix_free_) {
            cache.add(matrix_cache_);
        }
        cache.add(cache_);
        return cache;
    }

protected:
    void intersect_polygons(const CSPolygonArray& src_csp, const CSPolygonArray& tgt_scp);
    Matrix compute_1st_order_matrix();
    Matrix compute_2nd_order_matrix();
    void dump_intersection(const CSPolygon& plg_1, const CSPolygonArray& plg_2_array,
                           const std::vector<idx_t>& plg_2_idx_array) const;
    template <class TargetCellsIDs>
    void dump_intersection(const CSPolygon& plg_1, const CSPolygonArray& plg_2_array,
                           const TargetCellsIDs& plg_2_idx_array) const;
    std::vector<idx_t> sort_cell_edges(Mesh& mesh, idx_t cell_id) const;
    std::vector<idx_t> sort_node_edges(Mesh& mesh, idx_t cell_id) const;
    std::vector<idx_t> get_cell_neighbours(Mesh& mesh, idx_t jcell) const;
    std::vector<idx_t> get_node_neighbours(Mesh& mesh, idx_t jcell) const;
    CSPolygonArray get_polygons_celldata(Mesh& mesh) const;
    CSPolygonArray get_polygons_nodedata(Mesh& mesh, std::vector<idx_t>& csp2node,
                                         std::vector<std::vector<idx_t>>& node2csp,
                                         std::array<double, 2>& errors) const;

private:
    int next_index(int current_index, int size, int offset = 1) const;
    int prev_index(int current_index, int size, int offset = 1) const;

protected:
    bool src_cell_data_;
    bool tgt_cell_data_;
    FunctionSpace src_fs_;
    FunctionSpace tgt_fs_;
    mutable Mesh src_mesh_;
    mutable Mesh tgt_mesh_;
    int normalise_intersections_;
    int order_;
    bool matrix_free_;
    mutable RemapStat remap_stat_;

    std::shared_ptr<CachableData> cachable_data_shared_;
    const CachableData* cachable_data_;
    Cache cache_;

    // position and effective area of data points
    idx_t n_spoints_;
    idx_t n_tpoints_;
};


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
