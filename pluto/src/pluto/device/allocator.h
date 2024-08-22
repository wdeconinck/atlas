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

#include <memory_resource>

#include "hic/hic.h"
#include "pluto/offload/wait.h"

namespace pluto::device {

// --------------------------------------------------------------------------------------------------------


    template <class U, class... Args>
    HIC_GLOBAL void new_on_device(U* p, Args... args) {
        // printf("new_on_device %f\n",args...);
        new (p) U(args...);
    }


    template <class U>
    HIC_GLOBAL void delete_on_device(U* p) {
        p->~U();
    }

template<typename T>
class allocator {
public:
    using value_type = T;
    allocator() :
        memory_resource_(get_default_resource()) {}
    
    allocator(const allocator& other) :
        memory_resource_(other.memory_resource_) {}

    allocator(memory_resource* mr) :
        memory_resource_(mr) {}

    value_type* allocate(std::size_t size) const {
        return static_cast<value_type*>(memory_resource_->allocate(size * sizeof(value_type), 256));
    }

    void deallocate(value_type* ptr, std::size_t size) const {
        memory_resource_->deallocate(ptr, size * sizeof(value_type), 256);
    }

    template <class U, class... Args>
    void construct(U* p, Args&&... args) {
#if HIC_COMPILER
        new_on_device<<<1, 1>>>(p, std::forward<Args>(args)...);
        pluto::wait();
#else
        new_on_device(p, args...);
#endif
    }

    template <class U>
    void destroy(U* p) {
#if HIC_COMPILER
        delete_on_device<<<1, 1>>>(p);
        pluto::wait();
#else
        delete_on_device(p);
#endif
    }
private:
    std::pmr::memory_resource* memory_resource_{nullptr};
};

// --------------------------------------------------------------------------------------------------------

}

