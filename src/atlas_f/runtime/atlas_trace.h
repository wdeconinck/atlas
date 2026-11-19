/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <cstddef>

namespace atlas {
class Trace;
}

extern "C" {

atlas::Trace* new_atlas_Trace_view(const char* file, size_t file_length, int line, const char* title,
                                   size_t title_length);
atlas::Trace* new_atlas_Trace_labels_view_1(const char* file, size_t file_length, int line, const char* title,
                                            size_t title_length, const char* label1, size_t length1);
atlas::Trace* new_atlas_Trace_labels_view_2(const char* file, size_t file_length, int line, const char* title,
                                            size_t title_length, const char* label1, size_t length1,
                                            const char* label2, size_t length2);
atlas::Trace* new_atlas_Trace_labels_view_3(const char* file, size_t file_length, int line, const char* title,
                                            size_t title_length, const char* label1, size_t length1,
                                            const char* label2, size_t length2, const char* label3, size_t length3);
atlas::Trace* new_atlas_Trace_labels_view_4(const char* file, size_t file_length, int line, const char* title,
                                            size_t title_length, const char* label1, size_t length1,
                                            const char* label2, size_t length2, const char* label3, size_t length3,
                                            const char* label4, size_t length4);
atlas::Trace* new_atlas_Trace_labels_view_5(const char* file, size_t file_length, int line, const char* title,
                                            size_t title_length, const char* label1, size_t length1,
                                            const char* label2, size_t length2, const char* label3, size_t length3,
                                            const char* label4, size_t length4, const char* label5, size_t length5);
// Ownership functions require a live pointer returned by new_atlas_Trace*.
void retain_atlas_Trace(atlas::Trace* This);
void return_atlas_Trace(atlas::Trace* This);
void release_atlas_Trace(atlas::Trace* This);
void atlas_Trace__start(atlas::Trace* This);
void atlas_Trace__stop(atlas::Trace* This);
void atlas_Trace__pause(atlas::Trace* This);
void atlas_Trace__resume(atlas::Trace* This);
int atlas_Trace__running(atlas::Trace* This);
double atlas_Trace__elapsed(atlas::Trace* This);

}  // extern C
