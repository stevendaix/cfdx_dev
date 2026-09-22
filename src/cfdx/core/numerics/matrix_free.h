#pragma once

// M0.8-T10 — Public matrix-free FVM operator API.
// The numerical implementation lives in the linear-algebra layer so that
// physics code can consume the same LinearOperatorBase abstraction as
// assembled operators. This header is the numerics-facing entry point.

#include "cfdx/core/linalg/matrix_free_fv_operator.h"
#include "cfdx/core/geometry/geometry_cache.h"

namespace cfdx::core {

using MatrixFreeFvDiffusionOperator = FvDiffusionOperator;

} // namespace cfdx::core
