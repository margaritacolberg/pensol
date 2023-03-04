// Copyright (c) 2018-2023 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PEN_SNAPSHOT_H
#define PEN_SNAPSHOT_H

#include "hardspheres.h"
#include <string>

void write_snapshot(const std::string snapshot_name,
                    const std::vector<Vec3> &pos, Random &mt,
                    UpdateConfig &update_config);

void read_snapshot(const std::string snapshot_name, std::vector<Vec3> &pos,
                   Random &mt, UpdateConfig &update_config);

void read_input(const std::string input_name, std::vector<Vec3> &pos);

#endif
