#include "atlas/interpolation/Remap.h"

#include <map>
#include <string>

#include "atlas/functionspace.h"
#include "atlas/interpolation.h"

#include "atlas/util/relayout.h"

namespace atlas {

namespace {
Interpolation make_interpolation(const std::string& interpolation_method, FunctionSpace src_fs, FunctionSpace tgt_fs) {
    auto src_grid = src_fs.grid();
    auto tgt_grid = tgt_fs.grid();
    if (interpolation_method == "bilinear") {
        auto src_nonblocked_fs = functionspace::StructuredColumns{src_grid, grid::MatchingPartitioner(src_fs), option::halo(1)};
        auto tgt_nonblocked_fs = functionspace::PointCloud{tgt_grid, grid::MatchingPartitioner(src_fs) };
        return Interpolation(option::type("structured-bilinear"),  src_nonblocked_fs, tgt_nonblocked_fs);
    }
    else if (interpolation_method == "bicubic") {
        auto src_nonblocked_fs = functionspace::StructuredColumns{src_grid, grid::MatchingPartitioner(src_fs), option::halo(2)};
        auto tgt_nonblocked_fs = functionspace::PointCloud{tgt_grid, grid::MatchingPartitioner(src_fs) };
        return Interpolation(option::type("structured-bicubic"),  src_nonblocked_fs, tgt_nonblocked_fs);
    }
    else if (interpolation_method == "cons") {
        auto src_mesh = Mesh(src_grid, grid::MatchingPartitioner(src_fs));
        auto tgt_mesh = Mesh(tgt_grid, grid::MatchingPartitioner(src_mesh));
        auto src_nonblocked_fs = functionspace::NodeColumns{src_mesh, option::halo(10)};
        auto tgt_nonblocked_fs = functionspace::NodeColumns{tgt_mesh, option::halo(1)};
        return Interpolation(option::type("conservative-spherical-polygon"),  src_nonblocked_fs, tgt_nonblocked_fs);
    }
    else if (interpolation_method == "cons2") {
        auto src_mesh = Mesh(src_grid, grid::MatchingPartitioner(src_fs));
        auto tgt_mesh = Mesh(tgt_grid, grid::MatchingPartitioner(src_mesh));
        auto src_nonblocked_fs = functionspace::NodeColumns{src_mesh, option::halo(10)};
        auto tgt_nonblocked_fs = functionspace::NodeColumns{tgt_mesh, option::halo(1)};
        return Interpolation(option::type("conservative-spherical-polygon") | util::Config("order",2),  src_nonblocked_fs, tgt_nonblocked_fs);
    }
    else {
        throw_NotImplemented("Interpolation method "+interpolation_method+" not implemented",Here());
    }
}

void execute_interpolation(Interpolation interpolation, const FieldSet& src_blocked_fields, FieldSet& tgt_blocked_fields, bool on_device) {
    FieldSet src_nonblocked_fields;
    FieldSet tgt_nonblocked_fields;
    auto src_nonblocked_fs = interpolation.source();
    auto tgt_nonblocked_fs = interpolation.target();
    ATLAS_TRACE_SCOPE("create nonblocked fields "+std::string(on_device?"[device]":"[host]")){
        for(int i=0; i<src_blocked_fields.size(); ++i) {
            src_nonblocked_fields.add(src_nonblocked_fs.createField(src_blocked_fields[i]));
            tgt_nonblocked_fields.add(tgt_nonblocked_fs.createField(src_blocked_fields[i]));
        }
        if (on_device) {
            src_nonblocked_fields.allocateDevice();
            tgt_nonblocked_fields.allocateDevice();
        }
    }
    ATLAS_TRACE_SCOPE("copy_blocked_to_nonblocked "+std::string(on_device?"[device]":"[host]")){
        on_device ? src_blocked_fields.syncDevice() : src_blocked_fields.syncHost();
        copy_blocked_to_nonblocked(src_blocked_fields, src_nonblocked_fields, on_device);
        src_nonblocked_fields.haloExchange(on_device);
    }
    ATLAS_TRACE_SCOPE("nonblocked interpolation "+std::string(on_device?"[device]":"[host]")) {
        // TODO: add on_device argument to interpolation.execute
        interpolation.execute(src_nonblocked_fields, tgt_nonblocked_fields /*, on_device*/);
    }
    ATLAS_TRACE_SCOPE("copy_nonblocked_to_blocked "+std::string(on_device?"[device]":"[host]")){
        copy_nonblocked_to_blocked(tgt_nonblocked_fields, tgt_blocked_fields, on_device);
        tgt_blocked_fields.set_dirty(true);
        tgt_blocked_fields.setDeviceNeedsUpdate(on_device);
        tgt_blocked_fields.setHostNeedsUpdate(!on_device);
    }
}

std::string interpolation_key(std::string_view method, const FunctionSpace& source, const FunctionSpace& target) {
    std::string sgrid_name = source.grid().name();
    std::string tgrid_name = target.grid().name();
    std::string uid = "<uid>"; // could be based on source and target functionspace properties that are relevant for interpolation, e.g. halo size, blocking size, etc.
    return "{method="+std::string(method)+",source={type="+source.type()+",grid="+sgrid_name+"},target={type="+target.type()+",grid="+tgrid_name+"},uid="+uid+"}";
}

} // namespace

struct Remap::InterpolationsMap : public std::map<std::string, Interpolation> {
};

Remap::Remap() : 
    interpolations_(new InterpolationsMap()) {}

Remap::~Remap() = default;

bool Remap::has_interpolation(std::string_view method, const FunctionSpace& source, const FunctionSpace& target) {
    std::string key = interpolation_key(method, source, target);
    return interpolations_->count(key) > 0;
}

void Remap::setup_interpolation(std::string_view method, const FunctionSpace& source, const FunctionSpace& target) {
    std::string key = interpolation_key(method, source, target);
    if (interpolations_->count(key) == 0) {
        ATLAS_TRACE("setup_interpolation("+std::string(method)+")");
        interpolations_->emplace(key, make_interpolation(std::string(method), source, target));
    }
}

void Remap::interpolate(std::string_view method, const Field& source, Field& target, bool on_device) {
    FieldSet source_fields;
    FieldSet target_fields;
    source_fields.add(source);
    target_fields.add(target);
    interpolate(method, source_fields, target_fields, on_device);
}

void Remap::interpolate(std::string_view method, const FieldSet& source, FieldSet& target, bool on_device) {
    std::string key = interpolation_key(method, source[0].functionspace(), target[0].functionspace());
    if (interpolations_->count(key) == 0) {
        setup_interpolation(method, source[0].functionspace(), target[0].functionspace());
    }
    ATLAS_TRACE("execute_interpolation("+std::string(method)+")");
    execute_interpolation(interpolations_->at(key), source, target, on_device);
}

void Remap::clear() {
    interpolations_->clear();
}

extern "C" {
atlas::Remap* atlas__Remap__new() {
    return new atlas::Remap();
}

void atlas__Remap__delete(atlas::Remap* ptr) {
    delete ptr;
}

void atlas__Remap__setup_interpolation(atlas::Remap* remap, const char* method, size_t method_len, const atlas::functionspace::FunctionSpaceImpl* source, const atlas::functionspace::FunctionSpaceImpl* target) {
    remap->setup_interpolation(std::string_view(method, method_len), FunctionSpace(source), FunctionSpace(target));
}

void atlas__Remap__interpolate_field(atlas::Remap* remap, const char* method,
                                                  size_t method_len, const field::FieldImpl* source,
                                                  field::FieldImpl* target, int on_device) {
    const Field source_field(source);
    Field target_field(target);
    remap->interpolate(std::string_view(method, method_len), source_field, target_field, on_device);
}

void atlas__Remap__interpolate_fieldset(atlas::Remap* remap,
                                      const char* method, size_t method_len,
                                      const field::FieldSetImpl* source,
                                      field::FieldSetImpl* target, int on_device) {
    const FieldSet source_fields(source);
    FieldSet target_fields(target);
    remap->interpolate(std::string_view(method, method_len), source_fields, target_fields, on_device);
}
}



}