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
#include "atlas/mesh/actions/BuildEdges.h"

#include "atlas/util/ConvexSphericalPolygon.h"


namespace atlas {
namespace interpolation {
namespace method {


class ConservativeMethod : public Method {
public:
    typedef util::ConvexSphericalPolygon CSPolygon;
    typedef std::vector<std::pair<CSPolygon, int>> PolygonArray;
    typedef std::vector<std::tuple<CSPolygon, int>> CSPolygonArray;
    typedef array::ArrayView<double, 1> FieldArray;

    struct InterpolationParameters {
        std::vector<idx_t> tcell_id;
        std::vector<PointXYZ> centroids;
        std::vector<double> weights;
        std::vector<double> sweights;
    };

    ConservativeMethod(const Config& = util::NoConfig());

    using Method::do_setup;
    void do_setup(const FunctionSpace& src_fs, const FunctionSpace& tgt_fs);
    void do_setup(const Grid& src_grid, const Grid& tgt_grid);
    void do_setup(const Grid& src_grid, const Grid& tgt_grid, const Cache&) { ATLAS_NOTIMPLEMENTED; }
    void do_execute(const Field& src_field, Field& tgt_field);

    void set_order(int order) {
        order_ = order;
        if (order == 2) {
            mesh::actions::build_edges(src_mesh_, util::Config("pole_edges", false));
            setup_2nd_order_matrix();
        }
        if (order == 1) {
            setup_1st_order_matrix();
        }
    }
    void setup_stat(double& geo_create_err) const;
    void remap_stat(const FieldArray& src_field, const FieldArray& tgt_field, FieldArray& diff_field,
                    double& global_cons_err, double func(const PointLonLat&), double& remap_error_l2,
                    double& remap_error_linf) const;
    void print(std::ostream& out) const { out << "ConservativeMethod[]"; }

    bool src_cell_data() const { return src_cell_data_; }
    bool tgt_cell_data() const { return tgt_cell_data_; }
    const FunctionSpace& source() const { return src_fs_; }
    const FunctionSpace& target() const { return tgt_fs_; }
    int order() const { return order_; }
    Mesh src_mesh() const { return src_mesh_; }
    Mesh tgt_mesh() const { return tgt_mesh_; }
    double remap_err_l1() const { return remap_err_l1_; }
    double remap_err_linf() const { return remap_err_linf_; }
    double geo_err_intsc_l1() const { return geo_err_intsc_l1_; }
    double geo_err_intsc_linf() const { return geo_err_intsc_linf_; }
    inline const std::vector<InterpolationParameters>& iparam() const { return iparam_; }
    inline const PointXYZ& src_points(size_t id) const { return src_points_[id]; }
    inline const PointXYZ& tgt_points(size_t id) const { return tgt_points_[id]; }

protected:
    void intersect_polygons(const CSPolygonArray& src_csp, const CSPolygonArray& tgt_scp);
    void setup_1st_order_matrix();
    void setup_2nd_order_matrix();

private:
    template <class TargetCellsIDs>
    void dump_intersection(const util::ConvexSphericalPolygon& s_csp, const CSPolygonArray& tgt_csp,
                           const TargetCellsIDs& tgt_cells) const;

    std::vector<idx_t> sort_cell_edges(Mesh& mesh, idx_t cell_id) const;
    std::vector<idx_t> sort_node_edges(Mesh& mesh, idx_t cell_id) const;
    std::vector<idx_t> get_cell_neighbours(Mesh& mesh, idx_t jcell) const;
    std::vector<idx_t> get_node_neighbours(Mesh& mesh, idx_t jcell) const;
    CSPolygonArray get_polygons_celldata(Mesh& mesh) const;
    CSPolygonArray get_polygons_nodedata(Mesh& mesh, std::vector<idx_t>& csp2node,
                                         std::vector<std::vector<idx_t>>& node2csp) const;

protected:
    bool src_cell_data_;
    bool tgt_cell_data_;
    FunctionSpace src_fs_;
    FunctionSpace tgt_fs_;
    Mesh src_mesh_;
    Mesh tgt_mesh_;
    int normalise_intersections_;
    int order_;
    int fvtype_;
    idx_t n_spoints_;
    idx_t n_tpoints_;
    bool matrix_free_;
    std::vector<PointXYZ> src_points_;
    std::vector<PointXYZ> tgt_points_;
    Field src_areas_;
    Field tgt_areas_;
    double geo_err_intsc_l1_;                       // error in polygon intersections
    double geo_err_intsc_linf_;                     // error in polygon intersections
    double remap_err_l1_;                           // error in remapping
    double remap_err_linf_;                         // error in remapping
    std::vector<InterpolationParameters> iparam_;   // TODO: remove
    std::vector<idx_t> src_csp2node_;               // TODO: remove
    std::vector<idx_t> tgt_csp2node_;               // TODO: remove
    std::vector<std::vector<idx_t>> src_node2csp_;  // TODO: remove
    std::vector<std::vector<idx_t>> tgt_node2csp_;  // TODO: remove
};


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
