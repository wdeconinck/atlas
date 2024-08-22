/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "MemoryPoolResource.h"

#include "GatorMemoryResource.h"
#include "TraceMemoryResource.h"
#include "pluto/util/Trace.h"

namespace pluto{

std::pmr::memory_resource* MemoryPoolResource::resource(std::size_t bytes) {
    constexpr std::size_t MB = 1024*1024;
    constexpr std::size_t GB = 1024*MB;
    if (options_.largest_required_pool_block > 0 && bytes > options_.largest_required_pool_block) {
        return upstream_;
    }
    if (pools_.empty()) {
        if (options_.largest_required_pool_block != 0) {
            pool_block_sizes_ = {options_.largest_required_pool_block};
        }
        else {
            pool_block_sizes_ = {32*MB, 64*MB, 128*MB, 256*MB, 512*MB, 1*GB, 2*GB, 4*GB, 8*GB, 16*GB, 32*GB, 64*GB, 128*GB, 256*GB, 512*GB, 1024*GB};
        }
        pools_.resize(pool_block_sizes_.size());
        pool_block_size_ = 0;
    }

    auto upper_or_equal = std::upper_bound(pool_block_sizes_.begin(), pool_block_sizes_.end(), bytes-1);
    if (upper_or_equal == pool_block_sizes_.end()) {
        return upstream_;
    }

    if (*upper_or_equal > pool_block_size_) {
        pool_block_size_ = *upper_or_equal;
        std::size_t pool_index = upper_or_equal - pool_block_sizes_.begin();
        std::size_t blocks_per_chunk = std::max(options_.max_blocks_per_chunk, static_cast<std::size_t>(1));
        GatorOptions options;
        options.initial_size = blocks_per_chunk*pool_block_size_;
        options.grow_size    = blocks_per_chunk*pool_block_size_;
        if (TraceOptions::instance().enabled) {
            std::string size_str;
            if( pool_block_size_ >= GB ) {
                size_str = std::to_string(pool_block_size_/GB) +"GB"; 
            }
            else {
                size_str = std::to_string(pool_block_size_/MB) +"MB"; 
            }
            pools_[pool_index] =
                std::make_unique<TraceMemoryResource>("gator["+size_str+"]",
                std::make_unique<GatorMemoryResource>(options, upstream_) );
        }
        else {
            pools_[pool_index] = std::make_unique<GatorMemoryResource>(options, upstream_);
        }
        pool_ = pools_[pool_index].get();
    }
    return pool_;
}

void MemoryPoolResource::cleanup_unused_gators() {
    for (int i=0; i<pool_block_sizes_.size(); ++i) {
        if (pools_[i] && pool_block_size_ > pool_block_sizes_[i]) {
            pools_[i].reset();
        }
    }
}


}
