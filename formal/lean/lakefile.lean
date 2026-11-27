import Lake
open Lake DSL

package «cfdx-formal» where
  -- Keep the formal project independent from the C++/CMake build.

require mathlib from git
  "https://github.com/leanprover-community/mathlib4" @ "v4.19.0"

lean_lib CFDX
