// Copyright (c) 2018-2023 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PEN_SOLVENT_H
#define PEN_SOLVENT_H

#include "hardspheres.h"

#if NDEBUG
#define LOG_DEBUG(x)
#else
#define LOG_DEBUG(x)                                                           \
  std::cout << __FILE__ << ":" << __LINE__ << ": " << x << std::endl
#endif

void collide_sol(System &sys, Random &mt, const Param &p);

#endif
