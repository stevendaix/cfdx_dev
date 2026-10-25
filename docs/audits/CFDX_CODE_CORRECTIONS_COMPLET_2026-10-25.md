# CFDX — Code source intégral du candidat de corrections et portes de validation

**Référence vérifiée de l'audit :** `ca01b3f4dd9ae0621f05bb0463ac3eb03f938726` (25 octobre 2026). **Commit local du candidat :** `c8694bdd`. **IMPORTANT :** ce document contient le code réel archivé du candidat, pas du code fictif. Il ne représente pas un correctif final et validé des six verrous : plusieurs restent ouverts. Aucun test n'a été exécuté dans cette nouvelle session ; les résultats indiqués sont ceux des journaux de l'audit antérieur. Le master GitHub peut avoir évolué depuis ce SHA.

## Mode d'emploi

Le **patch Git directement applicable** figure en partie I. La **version complète, non tronquée, de chacun des 28 fichiers impactés** figure en partie II, avec son chemin. Pour appliquer le candidat sur un clone placé exactement au commit de référence, extraire le bloc `diff` de la partie I vers un fichier `.patch`, puis utiliser `git apply --check` et `git apply`. Le patch a été vérifié sur ce SHA lors de l'audit précédent, pas dans la session présente.

## Portes de qualification — état réel et travaux encore nécessaires

### P0 — S9 et COUPLED/Rhie–Chow : **NON FERMÉ**

**Code dans le document :** restauration de l'exécution du test S9 et ses assertions (`tests/validation/test_n9_s9_rhie_chow_consistency.cpp`). **Reste à coder et prouver :** contrat unique de vitesse physique `U` versus `H/A`, cohérence des deux API de flux et des corrections retardées, tests de pression affine avec `dA` variable et conservation locale. Sur le candidat : 13 tests S9 réussis, 2 assertions échouent. Ne pas les masquer pour obtenir un vert artificiel.

### P0 — Reconstruction de pression de production : **NON FERMÉ**

**Code dans le document :** modifications archivées du solveur stationnaire et contrôles géométriques. **Reste à coder et prouver :** reconstruction WLS/Gauss sélectionnable, BC explicites, reconstruction de face cohérente avec l'opérateur implicite SIMPLE/COUPLED, résidus et forces. L'overlay WLS exploratoire qui améliore les 82 premières itérations n'est PAS intégré au candidat ; l'essai n'a pas convergé à l'arrêt. Ne pas affirmer qu'il s'agit d'un correctif NACA final.

### P0 — Solveurs linéaires et AMG : **PARTIEL**

**Code dans le document :** vrai résidu de diagnostic à l'échec CG, gardes de finitude, régressions (`src/cfdx/core/linalg/cg_solver.h`, tests). **Reste à coder :** certification systématique du résidu vrai avec seuil absolu/relatif robuste, normes extrêmes, gestion FP32, invariance MPI, contrôle spectral/positivité du préconditionneur et politique de reconstruction AMG. La PR #1092 possédait des correctifs complémentaires, mais n'est pas contenue en entier ici. Le refus du CG sur NACA capturé reste distinct de la divergence de l'écoulement.

### P0 — NACA sur 3 grilles : **NON VALIDÉ**

**Code dans le document :** réparations import NASA/TMR, frontière sillage, fixtures/tests Python, métriques de paroi et journal de diagnostic. **Reste à exécuter après les corrections numériques :** n64/n128/n256 convergés sans assouplir les tolérances, vérification Cl/Cd/Cp/Cf, vérification GCI et indépendance au maillage, référence adaptée au régime physique. L'audit reproduit n64 convergé et n128 échoué à l'itération 20 ; n256 non atteint.

### P1 — Distances aux parois : **PARTIEL**

**Code dans le document :** gardes sur résultats invalides, correction Godunov/Eikonal, résidu advection–diffusion correspondant à l'opérateur assemblé, tests et benchmarks. **Reste :** obstacles non convexes, graphes mesh-wave traversant les parois, biais mathématique de Poisson (qui ne peut pas être effacé en itérant), tests géométriques sur raffinement ; score d'erreur indépendant du résidu de l'algorithme. La méthode BVH géométrique reste la référence exacte sur la surface discrète.

### P1 — RANS, LES, thermique, multiphasique : **NON VALIDÉS GLOBALEMENT**

**Code dans le document :** corrections temporelles, tests restaurés et protections numériques transverses. **Reste à intégrer/tester :** sélection de turbulence depuis les cas utilisateur, RANS fermé jusqu'aux forces et à `y+`, LES avec fermeture et convergence statistique, CHT interface conforme, P1/S2S/DOM avec bilan indépendant, VOF multiphasique complet et cas compressibles si revendiqués. Les composants testés ne prouvent pas la validation physique bout en bout.

## Résultats de référence, sans extrapolation

- Sur le candidat ancien `c8694bdd` : **25/26 cibles CTest ciblées** (S9 rouge, volontairement conservée).
- Python sur master `ca01b3f4` : **663 réussis, 49 ignorés**.
- Ghia et canal laminaire : campagnes longues réussies, BFS toujours diagnostique.
- Aucun résultat d'essai complet de ces six correctifs ne peut être honnêtement appelé « tous les modèles validés ».

## Partie I — Patch Git réel et intégral

**Attention :** la partie II répète les fichiers complets pour consultation, mais seul le patch de cette partie est destiné à être appliqué.

```diff
diff --git a/CMakeLists.txt b/CMakeLists.txt
index fcb70ed1..58582f20 100644
--- a/CMakeLists.txt
+++ b/CMakeLists.txt
@@ -419,6 +419,46 @@ if(CFDX_BUILD_TESTS)
         add_test(NAME ${test_name} COMMAND ${test_name})
     endmacro()
 
+    # Recovered standalone regressions: these sources previously had no CTest
+    # target, so API drift in them was invisible to normal builds.
+    add_cfdx_test(test_temporal tests/numerical/test_temporal.cpp)
+    add_cfdx_test(test_discretization_verification tests/validation/test_discretization_verification.cpp)
+    add_cfdx_test(test_matrix_diagnostics tests/unit/test_matrix_diagnostics.cpp)
+    add_cfdx_test(test_application_layer tests/unit/test_application_layer.cpp)
+    add_cfdx_test(test_memory_budget tests/unit/test_memory_budget.cpp)
+    add_cfdx_test(test_memory_planner tests/unit/test_memory_planner.cpp)
+    set_tests_properties(test_temporal test_discretization_verification PROPERTIES
+        LABELS "numerical;verification;recovered")
+    set_tests_properties(test_matrix_diagnostics test_application_layer
+        test_memory_budget test_memory_planner PROPERTIES LABELS "unit;recovered")
+
+    if(CFDX_ENABLE_MPI AND MPI_FOUND)
+        set(CFDX_RECOVERED_MPIEXEC_NUMPROCS_FLAG "${MPIEXEC_NUMPROCS_FLAG}")
+        if(CFDX_RECOVERED_MPIEXEC_NUMPROCS_FLAG STREQUAL "")
+            set(CFDX_RECOVERED_MPIEXEC_NUMPROCS_FLAG "-n")
+        endif()
+        add_executable(test_mpi_deterministic_fused tests/unit/test_mpi_deterministic_fused.cpp)
+        target_link_libraries(test_mpi_deterministic_fused PRIVATE cfdx_core MPI::MPI_CXX)
+        add_test(NAME test_mpi_deterministic_fused
+            COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_RECOVERED_MPIEXEC_NUMPROCS_FLAG} 2
+                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_deterministic_fused> ${MPIEXEC_POSTFLAGS})
+        set_tests_properties(test_mpi_deterministic_fused PROPERTIES
+            LABELS "unit;mpi;recovered" PROCESSORS 2)
+
+        # This source calls MPI-IO entry points directly and requires two ranks.
+        # Do not substitute a serial HDF5 build or a one-rank smoke test.
+        if(CFDX_ENABLE_PARALLEL_HDF5 AND HDF5_FOUND AND HDF5_IS_PARALLEL)
+            add_executable(test_parallel_hdf5 tests/unit/test_parallel_hdf5.cpp)
+            target_link_libraries(test_parallel_hdf5 PRIVATE cfdx_core MPI::MPI_CXX)
+            target_include_directories(test_parallel_hdf5 PRIVATE ${CMAKE_SOURCE_DIR}/tests)
+            add_test(NAME test_parallel_hdf5
+                COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_RECOVERED_MPIEXEC_NUMPROCS_FLAG} 2
+                    ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_parallel_hdf5> ${MPIEXEC_POSTFLAGS})
+            set_tests_properties(test_parallel_hdf5 PROPERTIES
+                LABELS "unit;mpi;hdf5;recovered" PROCESSORS 2)
+        endif()
+    endif()
+
     # N8 physical/production qualification reuses the existing validation
     # executables and records a machine-readable campaign report. It is
     # intentionally a separate CTest entry so the individual physical or
@@ -661,6 +701,9 @@ set_tests_properties(test_phase12_radiation_vv PROPERTIES LABELS "phase12;valida
     add_cfdx_test(test_n11_bound_reporting_hardening tests/unit/test_n11_bound_reporting_hardening.cpp)
     add_cfdx_test(test_phase10_model_registry tests/validation/test_phase10_model_registry.cpp)
     add_cfdx_test(test_wall_distance_methods tests/unit/test_wall_distance_methods.cpp)
+    add_cfdx_test(test_wall_distance_invalid_results tests/unit/test_wall_distance_invalid_results.cpp)
+    add_cfdx_test(test_wall_distance_godunov_ridges tests/unit/test_wall_distance_godunov_ridges.cpp)
+    add_cfdx_test(test_wall_distance_advection_residual tests/unit/test_wall_distance_advection_residual.cpp)
     add_test(NAME benchmark_wall_distance_naca0012
         COMMAND cfdx_naca0012_wall_distance_comparison
                 ${CMAKE_SOURCE_DIR}/tests/data/su2/mesh_NACA0012_inv.su2
@@ -999,6 +1042,18 @@ add_cfdx_test(test_n11_fw2_positivity_clamp_counter tests/unit/test_n11_fw2_posi
         if(CFDX_MPIEXEC_NUMPROCS_FLAG STREQUAL "")
             set(CFDX_MPIEXEC_NUMPROCS_FLAG "-n")
         endif()
+        # Only Open MPI accepts --oversubscribe. MPICH/Hydra rejects it
+        # before starting the restart reader. Preserve three-rank readers on
+        # two-core Open MPI CI workers without sending a vendor flag to all MPI.
+        set(CFDX_MPIEXEC_OVERSUBSCRIBE_FLAG "")
+        execute_process(COMMAND "${MPIEXEC_EXECUTABLE}" --version
+            OUTPUT_VARIABLE CFDX_MPIEXEC_VERSION_OUTPUT
+            ERROR_VARIABLE CFDX_MPIEXEC_VERSION_ERROR
+            RESULT_VARIABLE CFDX_MPIEXEC_VERSION_RESULT)
+        if(CFDX_MPIEXEC_VERSION_RESULT EQUAL 0 AND
+           "${CFDX_MPIEXEC_VERSION_OUTPUT}${CFDX_MPIEXEC_VERSION_ERROR}" MATCHES "Open MPI|OpenRTE")
+            set(CFDX_MPIEXEC_OVERSUBSCRIBE_FLAG "--oversubscribe")
+        endif()
         add_test(NAME test_mpi_deterministic
                  COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                          ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_deterministic> ${MPIEXEC_POSTFLAGS})
@@ -1038,7 +1093,7 @@ add_cfdx_test(test_n11_fw2_positivity_clamp_counter tests/unit/test_n11_fw2_posi
                              --reference ${CMAKE_BINARY_DIR}/production_poisson_reference.h5
                              ${MPIEXEC_POSTFLAGS})
             add_test(NAME test_production_mpi_poisson_reader
-                     COMMAND ${MPIEXEC_EXECUTABLE} --oversubscribe ${CFDX_MPIEXEC_NUMPROCS_FLAG} 3
+                     COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_OVERSUBSCRIBE_FLAG} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 3
                              ${MPIEXEC_PREFLAGS} $<TARGET_FILE:cfdx_mpi_poisson_restart>
                              --mode reader
                              --checkpoint ${CMAKE_BINARY_DIR}/production_poisson_checkpoint.h5
@@ -1050,7 +1105,7 @@ add_cfdx_test(test_n11_fw2_positivity_clamp_counter tests/unit/test_n11_fw2_posi
                      COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                              ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_phase5_restart_writer> ${MPIEXEC_POSTFLAGS})
             add_test(NAME test_phase5_restart_reader
-                     COMMAND ${MPIEXEC_EXECUTABLE} --oversubscribe ${CFDX_MPIEXEC_NUMPROCS_FLAG} 3
+                     COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_OVERSUBSCRIBE_FLAG} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 3
                              ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_phase5_restart_reader> ${MPIEXEC_POSTFLAGS})
             set_tests_properties(test_phase5_restart_reader PROPERTIES DEPENDS test_phase5_restart_writer)
         endif()
@@ -1117,6 +1172,13 @@ add_cfdx_test(test_n11_fw2_positivity_clamp_counter tests/unit/test_n11_fw2_posi
         add_cfdx_test(test_cross_import_equivalence tests/unit/test_cross_import_equivalence.cpp)
     endif()
     if(Python3_Interpreter_FOUND)
+        add_test(NAME test_cpp_harness_execution_contract
+            COMMAND ${Python3_EXECUTABLE}
+                    ${CMAKE_SOURCE_DIR}/scripts/check_cpp_test_execution.py
+                    --root ${CMAKE_SOURCE_DIR} --self-test)
+        set_tests_properties(test_cpp_harness_execution_contract PROPERTIES
+            LABELS "audit;test-infrastructure"
+            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
         add_test(NAME test_n9_s8_naca0012_contract
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_SOURCE_DIR}/tests/python/test_n9_s8_naca0012_contract.py)
diff --git a/apps/verification/naca0012_wall_distance_comparison.cpp b/apps/verification/naca0012_wall_distance_comparison.cpp
index c47bf1af..c424d060 100644
--- a/apps/verification/naca0012_wall_distance_comparison.cpp
+++ b/apps/verification/naca0012_wall_distance_comparison.cpp
@@ -286,6 +286,9 @@ int main(int argc,char** argv) {
                       << row.elapsed_ms << ',' << (row.converged ? "true" : "false") << ','
                       << row.stopping_reason << '\n';
         std::cout << "csv=" << output_path << '\n';
+        for(const auto& row:rows)
+            if(row.invalid_samples!=0)
+                throw std::runtime_error("invalid fluid distance field for "+row.method);
         return 0;
     } catch(const std::exception& error) {
         std::cerr << "NACA0012 wall-distance comparison failed: " << error.what() << '\n';
diff --git a/apps/verification/wall_distance_complex_benchmark.cpp b/apps/verification/wall_distance_complex_benchmark.cpp
index f851a450..0393efe1 100644
--- a/apps/verification/wall_distance_complex_benchmark.cpp
+++ b/apps/verification/wall_distance_complex_benchmark.cpp
@@ -204,5 +204,8 @@ int main(int argc,char** argv) {
     // gate for the benchmark infrastructure, not a physical validation claim.
     if(rows.front().l2>1e-14 || rows.front().near_l2>1e-14)
         throw std::runtime_error("exact wall-distance reference is not self-consistent");
+    for(const auto& row:rows)
+        if(row.invalid!=0)
+            throw std::runtime_error("invalid fluid distance field for "+row.method);
     return 0;
 }
diff --git a/docs/validation/NACA0012_IMPORT_REPAIR_2026-10-25.md b/docs/validation/NACA0012_IMPORT_REPAIR_2026-10-25.md
new file mode 100644
index 00000000..faea418f
--- /dev/null
+++ b/docs/validation/NACA0012_IMPORT_REPAIR_2026-10-25.md
@@ -0,0 +1,266 @@
+# NACA0012: importer repair and verification boundary, 2026-10-25
+
+Baseline: `71858a20fbfc67006824936d0646b3dd1cb79394`. This record distinguishes
+verified input/topology repairs from physical flow qualification. No NACA force
+or turbulence status is promoted by this change.
+
+## The different NACA cases
+
+| Path | Declared problem | Evidence supplied by that path |
+|---|---|---|
+| `tests/data/su2/inv_NACA0012_basic.cfg` | Euler, Mach 0.8, incidence 1.25 degrees | SU2 format-reader fixture; no CFDX compressible production solver |
+| `tests/data/su2/incomp_NACA0012.cfg` | Incompressible RANS, SA, incidence 10 degrees; source values imply Re approximately 6 million | Imported case metadata; production executable rejects RANS |
+| `test_naca0012_qualification` | Incompressible laminar Re=1000, incidence zero, generated O-grid | SIMPLE + upwind solve in the long test; the quick test does not solve the flow |
+| N9-S8 | Laminar Re=1000, incidence zero on imported NASA/TMR meshes | Contract and planned six-algorithm/three-mesh campaign; not a turbulent experimental comparison |
+| TUR-NACA0012 | Turbulent NACA validation against independent data | Still requires a production RANS path and retained physical evidence |
+
+The NASA turbulent validation problem uses Re=6 million and a modified sharp
+trailing edge. Its grid may be reused in a declared laminar study, but the
+experimental turbulent drag is not the oracle for a Re=1000 laminar solve.
+The generated O-grid uses the finite-thickness original trailing edge; those
+geometry definitions must also be recorded before interpreting a drag difference.
+
+Primary references:
+
+- https://tmbwg.github.io/turbmodels/naca0012_grids.html
+- https://tmbwg.github.io/turbmodels/naca0012_val.html
+- https://tmbwg.github.io/turbmodels/NACA0012_grids/n0012_225-65.nmf
+
+## Reproduced defects in the original NASA/TMR importer
+
+The official `n0012_225-65.p2dfmt` starts with three integers: `1`, `225`, `65`.
+The first integer is the number of blocks. The old parser instead interpreted
+`1` and `225` as the grid dimensions and rejected the actual file with:
+
+```text
+ValueError: expected 452 numeric tokens for 1x225, got 29253
+```
+
+The advertised `.gz` files were also read as plain ASCII.
+
+After bypassing only that header defect for diagnosis, the original `convert()`
+assigned every j=0 edge to `airfoil`. On the official 225x65 grid this produced
+224 wall faces instead of 128. Ninety-six of these faces are the two wake cuts;
+their centres extend to x=453.135486 chords. The neutral map explicitly assigns
+wall points 49..177 and a reversed one-to-one connection 1..49 to 225..177.
+These are **48 pairs of internal wake faces**, not solid walls. The original
+converter kept their 49 coincident vertex pairs topologically disconnected.
+
+The i=0 and i=ni-1 sides were named `wake_lower` and `wake_upper`, although the
+neutral map identifies them as downstream farfield. Extrusion patches were
+named `_front_cells` and `_back_cells`, receiving generic patch type 7 from the
+writer instead of EMPTY type 6.
+
+The previous `test_n9_s8_naca0012_contract.py` checked filename dimensions and
+words in a chapter. It did not open a PLOT3D file or inspect its converted
+connectivity, so none of these defects were exercised.
+
+## Repair
+
+`import_tmr_naca0012.py` now:
+
+1. reads the single-block formatted PLOT3D header, validates token count and
+   finite coordinates, and supports gzip;
+2. identifies the reversed coincident j=0 wake endpoint interval rather than
+   inventing a wall-index ratio from the filename;
+3. validates that the detected interval terminates at the unit-chord trailing
+   edge, so a partly corrupted wake cannot silently become extra solid wall;
+4. joins the paired vertices in both extrusion planes, creating ordinary
+   internal owner/neighbour faces along the wake;
+5. exports only the airfoil as WALL, combines all outer sides as farfield, and
+   exports `front` and `back` as EMPTY;
+6. retains the source SHA256 and topology counts in the HDF5 provenance. The
+   status remains `imported-provenance-only`.
+
+The implementation is deliberately scoped to the advertised single-block,
+unit-chord NASA NACA C-grid family. Other topology is rejected explicitly.
+It does not claim to be a generic PLOT3D multiblock connector.
+
+## Verification actually executed
+
+Nine added importer regression cases fail **8/9 on the baseline** and pass
+**9/9 after repair**. They cover the real header form, compressed input,
+malformed input, nonfinite coordinates, an independent miniature C-grid with
+exact face-count oracles, missing wake matching, partial wake corruption and
+provenance. They are collected by the existing full `tests/python` CI job.
+
+The official 225x65 grid was then imported and read through the actual native
+CFDX HDF5 reader, `validate_mesh` and `build_fv_geometry`, linked against the
+baseline production core. The following are measured values, not expected
+values copied into a numerical solution:
+
+| Native check | Measured result |
+|---|---:|
+| Cells | 14,336 |
+| Points after wake joining | 29,152 |
+| Faces | 57,584 |
+| Internal faces | 28,432 |
+| Internal wake faces | 48 |
+| Airfoil wall faces | 128 |
+| Farfield faces | 352 |
+| Front / back EMPTY faces | 14,336 each |
+| Unassigned physical boundary faces | 0 |
+| Minimum cell volume | 3.3070392178011831e-9 |
+| Minimum oriented face dot product Sf·d | 1.6535196089005912e-9 |
+| Geometry errors | 0 |
+| Geometry warnings | 1, large aspect ratio |
+| Maximum nonorthogonality | 12.06419507085366 degrees |
+| Maximum skewness | 0.15920642368411458 |
+| Maximum absolute surface closure | 3.5527136788005009e-15 |
+
+Every cell volume and every face-normal orientation was checked for finiteness
+and strict positivity. All boundary faces belonged to exactly one patch.
+
+Official source hashes:
+
+- `n0012_225-65.p2dfmt`: `66817d6bbc9353cb0ca9867f958da1bcce6485f7d68ec7038362a341bdb0fa8c`
+- `n0012_225-65.nmf`: `c6c32668f19bbd1a89510dfeb3a43c924137ae90afa7772eeb83875f9d9b34b0`
+
+Reproduction, with assets downloaded deliberately before the command:
+
+```bash
+python scripts/import_tmr_naca0012.py \
+  --input /path/to/n0012_225-65.p2dfmt \
+  --output /path/to/n0012_225-65.h5
+PYTHONPATH=src/cfdx/python python -m pytest -q \
+  tests/python/test_tmr_naca0012_import.py \
+  tests/python/test_su2_naca0012_tutorial.py \
+  tests/python/test_su2_naca0012_validation.py
+```
+
+The three selected Python files contain **12 passing tests** after the complete
+change. A passed mesh import is not a passed flow simulation.
+
+## Honest SU2 tutorial scope
+
+The original Python command stopped at canonical conversion: JST has no faithful
+registered CFDX convection counterpart. Its error printer also classified
+`unsupported_nonblocking` entries as blocking because it searched for the
+substring `block`. Its unused "theory" comparison treated a freestream static
+pressure as stagnation pressure and did not provide an aerodynamic oracle.
+The shell wrapper used unsupported options `--case`, `--output`, `--max-iter`.
+The production executable accepts `--mesh`, `--output-dir`, `--iterations` and
+explicitly rejects the compressible energy model.
+
+The Python command now retains distinct reader, boundary/initial-condition,
+canonical-conversion and physical-run verdicts in `reader_report.json`; the
+shell command delegates to it. Its default exit code remains nonzero (2) when
+full aerodynamic validation is unavailable. An explicit `--reader-only` permits
+success for the restricted input-reader scope:
+
+```bash
+python validation/run_su2_naca0012.py --reader-only --output-dir /path/to/evidence
+```
+
+No JST substitution is introduced and no solver result is manufactured.
+
+The current SU2 adapter still has these independently observed metadata defects,
+which this importer-focused patch does not repair:
+
+- `MARKER_EULER=(airfoil)` becomes `wall_no_slip` through a name-based fallback;
+- `farfield` becomes `outlet`/`fixed` with zero velocity;
+- initial scalar speed is 272 m/s, but the vector contains the unit incidence
+  direction `[0.9997620270799091, 0.02181488503456112, 0]`.
+
+The updated report exposes all three mismatches. Existing format-reader tests
+check patch names and cardinalities, not faithful boundary/initial-condition
+semantics. These defects must be fixed before that importer can be described as
+physically equivalent to the original SU2 case.
+
+## Remaining numerical limitations relevant to the laminar NACA solve
+
+The existing generated O-grid n64 is topologically valid, but native geometry
+inspection gives nonorthogonality 80.550987 degrees, maximum skewness 0.986432
+and 82 warnings. Its wall-adjacent centre distance has median 0.0268116 chord,
+minimum 0.00511106 and maximum 0.0326396. It is a coarse, highly skewed starting
+mesh. The indexing of its circular outer contour also breaks the reflection
+pairing of the airfoil surface: the maximum corresponding reflection difference
+in the first layer is 0.00319299 chord. These are measured quality limitations;
+their contribution to Cd/Cl must be quantified by controlled flow comparisons.
+
+On both actual meshes, the exact manufactured affine field phi=x+2y was supplied
+to the production `gauss_gradient_with_boundary`. Only cells with no physical
+boundary face were included, so the result is independent of boundary-value
+approximations. WLS received explicit exact Dirichlet values at physical faces
+and the correct zero spanwise derivative.
+
+| Mesh | Interior cells | Production Gauss gradient max absolute error | WLS max absolute error |
+|---|---:|---:|---:|
+| Generated O-grid n64 | 2,838 | 0.10987995423605813 | 1.1082784313214814e-14 |
+| Repaired TMR 225x65 | 13,858 | 0.17283073646457187 | 5.4393680889484996e-9 |
+| Repaired TMR 449x129 | 56,386 | 0.10357401848014633 | 1.3019123605904252e-8 |
+| Repaired TMR 897x257 | 227,458 | 0.07351694706579924 | 2.9871163766996353e-8 |
+
+This demonstrates a reconstruction error before any nonlinear flow solve. The
+production helper still uses Euclidean interpolation distances, separately from
+`core/numerics/gradient.h`, and has no face-skew correction. Merely correcting a
+standalone gradient kernel cannot certify this production path. These numbers
+are **not** drag errors and do not by themselves establish the whole cause of
+an aerodynamic discrepancy.
+
+Further requirements remain unchanged: fix and independently verify the
+pressure/velocity face-flux contract; validate the actual production gradient
+and nonorthogonal diffusion; freeze consistent farfield conditions; compare
+pressure and viscous force parts and surface Cp/Cf; execute the declared three
+meshes and six algorithms; retain true residuals, independent balances and an
+appropriate physical reference. A unit test, quick contract, converged linear
+solve or this mesh repair cannot replace those results.
+
+## Refined official grids and boundary-membership cost
+
+The same native reader/geometry checks were also completed on the two finer
+official files. Their gzip payloads were retrieved from the NASA archive,
+`https://www.nasa.gov/wp-content/uploads/2026/02/naca0012-grids.zip`, and checked
+against its ZIP member CRC and length before decompression.
+
+| Check | 449x129 | 897x257 |
+|---|---:|---:|
+| Cells | 57,344 | 229,376 |
+| Points after joining | 115,648 | 460,672 |
+| Faces | 229,856 | 918,464 |
+| Internal faces | 114,208 | 457,792 |
+| Internal wake faces | 96 | 192 |
+| Airfoil wall faces | 256 | 512 |
+| Farfield faces | 704 | 1,408 |
+| Front/back EMPTY faces each | 57,344 | 229,376 |
+| Unassigned boundary faces | 0 | 0 |
+| Minimum cell volume | 5.1820719983179846e-10 | 9.579601570460766e-11 |
+| Minimum oriented face dot product | 2.5910359991589918e-10 | 4.7897995986945982e-11 |
+| Geometry errors | 0 | 0 |
+| Geometry warnings | 1, aspect ratio | 1, aspect ratio |
+| Maximum nonorthogonality, degrees | 9.2209192083970688 | 6.1731507752808854 |
+| Maximum skewness | 0.10395424486342719 | 0.055277414197076311 |
+| Maximum absolute surface closure | 1.7763568394002505e-15 | 4.4408920985006262e-16 |
+
+Source SHA256:
+
+- `n0012_449-129.p2dfmt.gz`: `6d60dfa6ef9c7544ee156276fcc5487d7eacddadbdf6e81cc6f10e78fb03eff8`
+- `n0012_897-257.p2dfmt.gz`: `c603ea77017ae8b97af581237b81a8251fd2d3af1c49c8c6f2cd2687aafcb28e`
+
+During these checks, `BoundaryPatches::is_consistent_with_mesh` was found to
+rescan every patch face for each physical boundary face. With the two large
+EMPTY patches, this is quadratic in the number of cells. The native HDF5
+reader invokes this validation, and subsequent `validate_mesh` invokes it again.
+The repair retains the duplicate/range and internal-face checks, records listed
+faces in a bitmap, and then scans the boundary once. It is linear in the number
+of faces plus patch entries.
+
+A benchmark on the three actual imported grids called the exact original
+function body and its replacement on the same loaded mesh. Both returned true.
+These are single-process wall-clock observations for that function alone,
+not an end-to-end CFD speedup guarantee:
+
+| Grid | Boundary faces | Original, seconds | Replacement, seconds |
+|---|---:|---:|---:|
+| 225x65 | 29,152 | 0.121096512 | 0.000140853 |
+| 449x129 | 115,648 | 2.34260247 | 0.000537623 |
+| 897x257 | 460,672 | 31.2930124 | 0.002588804 |
+
+Five membership regression cases explicitly check a complete patch cover,
+an absent face, a repeated face, an internal face in a patch, and an out-of-range
+face. The expanded `test_mesh` passes 16/16 both with the original and repaired
+function, demonstrating preserved acceptance/rejection decisions. The existing
+`test_boundary` passes 11/11 with the repair. The complete native geometry and
+manufactured-gradient diagnostics on all three actual grids are byte-for-byte
+identical before and after this membership-only change. No timing threshold is
+used as a unit-test oracle.
diff --git a/docs/validation/NACA0012_RUNTIME_AUDIT_2026-10-25.md b/docs/validation/NACA0012_RUNTIME_AUDIT_2026-10-25.md
new file mode 100644
index 00000000..e7a29f64
--- /dev/null
+++ b/docs/validation/NACA0012_RUNTIME_AUDIT_2026-10-25.md
@@ -0,0 +1,263 @@
+# NACA0012 runtime audit, 2026-10-25
+
+Baseline: `71858a20fbfc67006824936d0646b3dd1cb79394`.
+This is an execution record for the existing laminar qualification, separate
+from the NASA mesh importer repair. It does not promote the NACA case to
+physically validated status.
+
+## Actual baseline outcome
+
+The complete baseline CTest campaign ran `test_naca0012_qualification` on the
+three meshes generated by the repository. It aborted after 20.01 seconds.
+`test_naca0012_wall_shear_audit` subsequently failed by the same mechanism after
+19.78 seconds. Their quick counterparts passed because they do not execute the
+full three-level flow campaign.
+
+The controls in `tests/validation/test_naca0012_geometry.cpp` were unmodified:
+laminar incompressible Re=1000, zero incidence, SIMPLE, upwind convection,
+linear tolerance 1e-10 and limit 4000, nonlinear tolerance 1e-8 and limit 8000,
+velocity relaxation 0.6 and pressure relaxation 0.25. The source prescribes
+U=(1,0,0) and p=0 at the farfield, no slip on the airfoil, and spanwise EMPTY
+faces. This is not the turbulent NASA Re=6-million validation case.
+
+| Level | Actual outcome |
+|---|---|
+| n64, 3,096 cells | Solver converged in 435 iterations |
+| n128, 12,336 cells | Pressure solve returned NOT_APPLICABLE at nonlinear iteration 20; three lower-relaxation retries also failed |
+| n256, 49,248 cells | Not reached by either full test |
+
+The converged n64 output was Cd=0.0903395, Cl=1.55115e-5,
+Cm=-0.00094537, Cd_pressure=0.0345473, Cd_viscous=0.0557922 and
+continuity_linf=1.50029e-10. A converged coarse solution does not establish
+mesh convergence or agreement with a physical reference. The final test itself
+labels its prospective successful outcome DIAGNOSTIC, with the physical
+reference gate still pending.
+
+The wall-shear audit applied three postprocessors to the same n64 solution:
+
+| Reconstruction | Viscous Cd |
+|---|---:|
+| Production Gauss | 0.0557922 |
+| Weighted least squares | 0.0581681 |
+| One-dimensional wall-normal estimate | 0.106364 |
+
+Its reported relative spread was 0.688604. This is a measured sensitivity of
+postprocessing on one coarse field, not three independent flow solutions and
+not proof that the largest value is the correct one.
+
+## The default failure message hides the actual residual
+
+The baseline exception prints:
+
+```text
+pressure-correction solve did not converge
+(status=3, iterations=3, residual=0.000000, relative=0.000000)
+```
+
+Status 3 is `SolverStatus::NOT_APPLICABLE`, not exhaustion of the 4000-iteration
+linear budget. The automatic pressure plan resolves to CG with
+`NativeBoomerAMGPreconditioner`, using `ReusableCgContext`.
+
+A diagnostic-only local include overlay captured the real matrix, RHS and CG
+iterate on that failure. It changed neither the flow formulas nor the controls.
+A second diagnostic run printed every nonlinear iteration and captured the
+first pressure matrix as well. These runs reproduced the original failure;
+they are follow-up observations of the same defect, not additional independent
+qualification cases.
+
+At the first failed pressure solve:
+
+- RHS norm: 5.5156653682823174e12;
+- true residual norm, recomputed as ||b-Ax||: 1.7332238406933466e11;
+- true relative residual: 0.031423658343383246;
+- requested relative tolerance: 1e-10;
+- stored `SolverResult.residual` and `residual_relative`: both zero;
+- initial preconditioned product r^T M^-1 r: -4.3863563470325308e37.
+
+The early NOT_APPLICABLE exits in `cg_solver.h` do not consistently populate
+the residual fields. An independently instrumented replay identified the
+specific exit here: the nonpositive `replaced_rsnew` check after residual
+replacement, at baseline `cg_solver.h:361`. The initial r^T M^-1 r had already
+been negative; CG only checked its finiteness and took a negative first step
+coefficient. Both the missing true-residual reporting and the missing initial
+positivity check deserve explicit regression coverage.
+
+## Nonlinear divergence precedes the linear refusal
+
+Printing every iteration instead of only 1 and 200 exposed the earlier growth:
+
+| Nonlinear iteration | Momentum RHS-relative residual | Momentum raw residual | Pressure linear relative residual |
+|---|---:|---:|---:|
+| 1 | 0.02170991661 | 0.01063918177 | 5.796087420e-11 |
+| 5 | 0.2776201489 | 0.1360274329 | 2.989268378e-11 |
+| 10 | 1.377742198 | 11.22542016 | 1.022865973e-11 |
+| 13 | 51.84263613 | 6062.070165 | 2.017762218e-11 |
+| 16 | 8479.466128 | 2.491988208e9 | 2.603139310e-11 |
+| 19 | 2.760647353e9 | 5.450210718e24 | 5.888560796e-11 |
+
+At iteration 19 the flux/velocity mismatch was 1.8018721310623572e14.
+The linear pressure solves still met their requested tolerance up to that
+iteration. Therefore the final AMG/CG refusal is a later failure in a nonlinear
+calculation that was already diverging. Replacing the last pressure solver
+alone cannot certify a physically meaningful flow.
+
+The earliest largest momentum errors are adjacent to the airfoil. This narrows
+where to investigate, but does not isolate which production discretization
+term initiates the instability. The existing generated O-grids have the
+following independently measured geometry/reconstruction limitations:
+
+| O-grid level | Max nonorthogonality | Max skewness | Geometry warnings | Affine production-gradient maximum error |
+|---|---:|---:|---:|---:|
+| n64 | 80.550987 degrees | 0.986432 | 82 | 0.1098799542 |
+| n128 | 81.280277 degrees | 0.988442 | 163 | 0.1671984660 |
+| n256 | 81.647616 degrees | 0.989393 | 351 | 0.1594223136 |
+
+All three meshes have positive finite volumes, correctly oriented faces and
+zero topological/geometry errors. The affine test uses phi=x+2y, includes only
+cells without a physical boundary face, and precedes any flow solve. Its
+corresponding WLS errors are 1.11e-14, 6.14e-14 and 4.79e-13. These are
+reconstruction errors, not measured drag errors. A controlled production
+discretization comparison is still needed to assign the primary nonlinear
+instability to a particular term.
+
+### One controlled relaxation diagnostic
+
+A separate n128 A/B experiment changed only the initial velocity relaxation
+from 0.6 to 0.3. Mesh, boundary conditions, schemes, alpha_p=0.25, residual
+tolerances and the configured 8000/4000 nonlinear/linear iteration budgets
+were preserved. The wrapper declared an external 300-second time budget and
+an early-stop threshold of momentum RHS-relative residual greater than 1e6.
+Printing each iteration and collecting the initial pressure snapshot were
+observational changes only.
+
+This run was stopped at iteration 48, after 8.019055 seconds, for manifest
+nonlinear divergence. The last measured momentum RHS-relative residual was
+3.2768351750380591e7, raw momentum residual 1.1839691475182645e19,
+normalized continuity 1.9507736565740657, and flux/velocity mismatch
+2.2336081815654251e10. The pressure linear solve still converged to
+2.3469556868328088e-11 at that iteration, below the unchanged 1e-10 tolerance.
+
+Thus halving this relaxation delays but does not remove the nonlinear
+instability. This single comparison does not rule out other stable controls;
+no parameter sweep was used and no force qualification is inferred. The
+independent pressure-solver failure remains a later issue, rather than the
+only reason the original flow calculation fails.
+
+## The captured pressure matrix remains grounded and SPD
+
+The matrix has 12,336 rows and 61,166 entries. It is exactly symmetric, has a
+positive diagonal, strictly negative off-diagonal entries, one connected
+component and 257 strictly dominant rows attached to prescribed farfield
+pressure. There is no missing pressure gauge in this snapshot.
+
+The diagonal ranges from 1.6798584831606658e-16 to
+0.011604898446567188, a ratio of about 6.91e13. The independent normalized
+analysis of B=D^-1/2 A D^-1/2 found:
+
+- lambda_min = 1.21312018e-9;
+- lambda_max approximately 1.999999999;
+- normalized condition number approximately 1.649e9;
+- positive LDL^T pivots, smallest 7.707e-8;
+- factor symmetry defect 1.67e-16.
+
+The fine matrix's SPD character is strongly supported, while its scaling and
+conditioning have become extreme. A separate SciPy sparse direct solve gave
+relative residual 2.41565e-10 and a maximum pressure correction 2.80923e24;
+this is another observation of the damaged nonlinear state, not an accepted
+flow solution.
+
+## Reusing the AMG hierarchy breaks its CG-compatible action
+
+At the failure, the reusable context had one hierarchy setup and 19 numeric
+updates. An isolated native replay built AMG from the first pressure matrix
+and then updated it directly to the captured iteration-20 matrix. It reproduced
+the original negative r^T M^-1 r and the complete failed CG result bit for bit.
+Thus intermediate flow states are not needed to reproduce this particular
+AMG reuse defect.
+
+All following results use the same captured A and b, zero initial correction,
+linear tolerance 1e-10 and budget 4000:
+
+| Linear path | Initial r^T M^-1 r | Outcome | Iterations | Recomputed relative residual |
+|---|---:|---|---:|---:|
+| Original first-step AMG hierarchy + numeric update, CG | -4.386356347e37 | NOT_APPLICABLE | 3 | 0.03142365834 |
+| Newly built AMG, CG | 2.113141873e33 | MAX_ITER_REACHED / exhausted useful recurrence | 126 | 5.482475639e-10 |
+| Newly built AMG, FGMRES | 2.113141873e33 | CONVERGED | 517 | 9.873355974e-11 |
+| Jacobi, CG | 1.530768460e33 | MAX_ITER_REACHED | 4000 | 6.824463987e-9 |
+
+FGMRES convergence is a real linear snapshot result. No force, conservation or
+nonlinear-convergence claim is inferred from it.
+
+### Measured mechanism inside the reused hierarchy
+
+The native AMG implementation uses six pre- and six post-Jacobi sweeps with
+fixed omega=0.7. `update_values` retains the old interpolation P and rebuilds
+the Galerkin operators. It does not requalify the stability of that fixed
+smoother on the refreshed levels.
+
+A read-only diagnostic getter exported the actual stored level matrices.
+Its direct `MatrixFreeVcyclePreconditioner` replay first reproduced the native
+wrapper's r^T M^-1 r exactly, so the level analysis concerns the same hierarchy.
+For a normalized symmetric level operator, Jacobi's eigenmode factor is
+1-omega*lambda. Stability requires lambda_max < 2/omega = 2.857142857.
+
+| Reused hierarchy level | Unknowns | Measured lambda_max | High-mode factor 1-0.7*lambda_max |
+|---|---:|---:|---:|
+| 3 | 1297 | 2.976252949 | -1.083377064 |
+| 4 | 566 | 3.518576285 | -1.463003399 |
+| 5 | 288 | 3.987227820 | -1.791059474 |
+| 6 | 163 | 4.455448843 | -2.118814190 |
+| 7 | 116 | 4.173164525 | -1.921215168 |
+| 8 | 102 | 3.337720948 | -1.336404663 |
+
+The level-6 eigenpair residual is 2.20e-15. These modes are amplified on each
+smoothing sweep. In contrast, every coarse level of newly built AMG is below
+the stability bound; its largest observed coarse lambda_max is 1.988241902.
+The issue is not an indefinite fine pressure matrix. It is an old hierarchy
+whose refreshed levels no longer support the configured fixed Jacobi damping.
+
+A repair must preserve the explicit numeric-only contract of `update_values`:
+either requalify/adapt smoothing at refresh, or report lost applicability so
+the caller can deliberately rebuild a hierarchy. It must also retain true
+residual checks. Simply increasing linear budgets, weakening tolerances, or
+declaring the n64 force validated would not address these observations.
+
+## Retained evidence and reproduction
+
+Workspace evidence is in `audit_naca/flow/` next to the baseline checkout:
+
+- `naca128_diagnostic.cpp` and `include/`: diagnostic-only wrapper and overlay;
+- `naca128_diagnostic.log`: original failure and four true residual snapshots;
+- `pressure_failure_1.mtx`: failing iteration-20 matrix in Matrix Market format;
+- `pressure_failure_1.vectors`: columns RHS and failed CG correction, 17 digits;
+- `trace/naca128_trace.log`: every accepted nonlinear iteration;
+- `trace/pressure_failure_1.mtx`: initial iteration-1 matrix;
+- `relaxed/`: source fixture with the single alpha_u change, bounded runner,
+  and complete `naca128_alpha03.log`;
+- `replay_pressure_snapshot.cpp` and `.log`: identical-system linear comparison;
+- `export_amg_levels.cpp`, `amg_include/`, `amg_levels/`: actual refreshed/fresh
+  level matrices and spectral analysis script, JSON and log.
+
+The separate solver audit retains `naca_pressure_matrix_analysis.*` and the
+CG exit trace. The original baseline CTest log is
+`audit_output/ctest_release_all.log`; the full NACA baseline was not rerun as a
+second qualification campaign.
+
+To replay the isolated linear failure, from the retained evidence directory:
+
+```bash
+./replay_pressure_snapshot \
+  pressure_failure_1.mtx pressure_failure_1.vectors \
+  trace/pressure_failure_1.mtx
+
+cd amg_levels
+../export_amg_levels ../trace/pressure_failure_1.mtx \
+  ../pressure_failure_1.mtx ../pressure_failure_1.vectors
+python analyze_levels.py
+```
+
+The executable sources compile with C++17 against the unchanged baseline
+headers and `libcfdx_core.a`. The flow diagnostic additionally needs the native
+HDF5 build. Neither diagnostic overlay is a production change. The importer
+and boundary-membership fixes are recorded separately in
+`NACA0012_IMPORT_REPAIR_2026-10-25.md`.
diff --git a/docs/validation/SU2_NACA0012_READER.md b/docs/validation/SU2_NACA0012_READER.md
index ec30cb75..fe4dd027 100644
--- a/docs/validation/SU2_NACA0012_READER.md
+++ b/docs/validation/SU2_NACA0012_READER.md
@@ -83,3 +83,22 @@ The reader qualification passes only when:
 
 No CFDX numerical result is compared here. Numerical validation remains a
 separate solver-level campaign.
+
+## Reproducible command and remaining fidelity limits (2026-10-25)
+
+Run the explicitly limited reader scope with:
+
+```bash
+python validation/run_su2_naca0012.py --reader-only --output-dir reader-evidence
+```
+
+The JSON report separates input parsing from boundary/initial-condition
+fidelity, canonical numerical conversion and physical simulation. Without
+`--reader-only`, the legacy command exits nonzero when the full aerodynamic
+validation cannot run. It never substitutes another scheme for unsupported JST
+and never presents an isentropic static-pressure expression as a NACA Cp oracle.
+
+The current adapter's Euler wall, farfield velocity and velocity-vector scaling
+still require correction; a reader PASS only covers the explicitly listed
+format checks. See [the importer repair and audit record](NACA0012_IMPORT_REPAIR_2026-10-25.md)
+for measured counterexamples and the separate NASA/TMR topology repair.
diff --git a/scripts/check_cpp_test_execution.py b/scripts/check_cpp_test_execution.py
new file mode 100644
index 00000000..beceb75c
--- /dev/null
+++ b/scripts/check_cpp_test_execution.py
@@ -0,0 +1,120 @@
+#!/usr/bin/env python3
+"""Reject deferred CFDX harness registrations without a run_all() call.
+
+This narrow static guard recognizes literal-name/lambda run_case registrations
+from common/test_harness.h. Local numerical helpers named run_case, comments,
+and string literals are not test registrations. It complements execution; it
+does not prove that an arbitrary call is reachable or that an oracle is sound.
+"""
+
+from __future__ import annotations
+
+import argparse
+from pathlib import Path
+import re
+
+
+def cpp_tokens(source: str) -> list[tuple[str, str]]:
+    tokens: list[tuple[str, str]] = []
+    i = 0
+    while i < len(source):
+        if source[i].isspace():
+            i += 1
+        elif source.startswith("//", i):
+            end = source.find("\n", i + 2)
+            i = len(source) if end < 0 else end + 1
+        elif source.startswith("/*", i):
+            end = source.find("*/", i + 2)
+            i = len(source) if end < 0 else end + 2
+        elif source.startswith('R"', i):
+            opening = source.find("(", i + 2, i + 19)
+            if opening < 0:
+                i += 1
+                continue
+            closing = ")" + source[i + 2 : opening] + '"'
+            end = source.find(closing, opening + 1)
+            end = len(source) if end < 0 else end + len(closing)
+            tokens.append(("string", source[i:end]))
+            i = end
+        elif source[i] in "\"'":
+            quote, start = source[i], i
+            i += 1
+            while i < len(source):
+                if source[i] == "\\":
+                    i += 2
+                elif source[i] == quote:
+                    i += 1
+                    break
+                else:
+                    i += 1
+            tokens.append(("string", source[start:i]))
+        elif source[i].isalpha() or source[i] == "_":
+            match = re.match(r"[A-Za-z_][A-Za-z_0-9]*", source[i:])
+            assert match is not None
+            value = match.group(0)
+            tokens.append(("identifier", value))
+            i += len(value)
+        else:
+            tokens.append(("punctuation", source[i]))
+            i += 1
+    return tokens
+
+
+def inspect_source(source: str) -> tuple[int, bool]:
+    tokens = cpp_tokens(source)
+    values = [token[1] for token in tokens]
+    uses_harness = any(
+        values[i : i + 3] == ["#", "include", '"common/test_harness.h"']
+        for i in range(len(tokens) - 2)
+    )
+    if not uses_harness:
+        return 0, False
+    registrations = sum(
+        values[i] == "run_case"
+        and values[i + 1] == "("
+        and tokens[i + 2][0] == "string"
+        and values[i + 3 : i + 5] == [",", "["]
+        for i in range(len(tokens) - 4)
+    )
+    executes = any(
+        values[i : i + 3] == ["run_all", "(", ")"]
+        for i in range(len(tokens) - 2)
+    )
+    return registrations, registrations > 0 and not executes
+
+
+def self_test() -> None:
+    prefix = '#include "common/test_harness.h"\n'
+    registration = 'int main() { run_case("example", [] { fail(); }); return 0; }'
+    assert inspect_source(prefix + registration) == (1, True)
+    assert inspect_source(prefix + registration.replace("return 0", "return run_all()")) == (1, False)
+    assert inspect_source(prefix + '// run_case("ignored", [] {});\n') == (0, False)
+    assert inspect_source(prefix + 'const char *x = R"tag(run_case("ignored", [] {});)tag";') == (0, False)
+    assert inspect_source(prefix + 'double run_case(int n) { return n; } int main() { run_case(16); }') == (0, False)
+    assert inspect_source(prefix + registration + '// run_all()\nconst char *x="run_all()";') == (1, True)
+    print("Harness guard self-tests: 6 passed")
+
+
+def main() -> int:
+    parser = argparse.ArgumentParser(description=__doc__)
+    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
+    parser.add_argument("--self-test", action="store_true")
+    args = parser.parse_args()
+    if args.self_test:
+        self_test()
+    files, cases, missing = 0, 0, []
+    for path in sorted((args.root / "tests").rglob("*.cpp")):
+        count, absent = inspect_source(path.read_text(encoding="utf-8"))
+        if count:
+            files += 1
+            cases += count
+        if absent:
+            missing.append((path.relative_to(args.root), count))
+    for path, count in missing:
+        print(f"FAIL: {path}: {count} deferred case(s) but no run_all() call")
+    print(f"Harness execution guard: {files} files, {cases} registered cases, {len(missing)} missing executions")
+    return 1 if missing else 0
+
+
+if __name__ == "__main__":
+    raise SystemExit(main())
diff --git a/scripts/import_tmr_naca0012.py b/scripts/import_tmr_naca0012.py
index 9fe07702..d6667b27 100644
--- a/scripts/import_tmr_naca0012.py
+++ b/scripts/import_tmr_naca0012.py
@@ -1,68 +1,173 @@
 #!/usr/bin/env python3
-"""Import a NASA/TMR NACA0012 PLOT3D 2-D C-grid into CFDX-HDF5."""
+"""Import the single-block NASA/TMR NACA0012 PLOT3D C-grid.
+
+The j=0 boundary contains the airfoil AND two coincident wake cuts.  The
+wake cuts are joined as internal faces; only the closed airfoil is a wall.
+This is an input/topology conversion, not a physical qualification.
+"""
 from __future__ import annotations
+
 import argparse
+import gzip
+import hashlib
 from pathlib import Path
 import re
+
 import h5py
 import numpy as np
+
 from generate_2d_validation_meshes import write_mesh
 
-TMR_URL="https://turbmodels.larc.nasa.gov/naca0012_grids.html"
-TMR_ARCHIVE_URL="https://www.nasa.gov/wp-content/uploads/2026/02/naca0012-grids.zip"
+TMR_URL = "https://tmbwg.github.io/turbmodels/naca0012_grids.html"
+TMR_ARCHIVE_URL = "https://www.nasa.gov/wp-content/uploads/2026/02/naca0012-grids.zip"
+
 
 def read_plot3d_2d(path: Path):
-    tokens=path.read_text(encoding="ascii").replace(","," ").split()
-    if len(tokens)<2: raise ValueError(f"PLOT3D file is empty: {path}")
-    ni,nj=int(tokens[0]),int(tokens[1])
-    expected=2+2*ni*nj
-    if len(tokens)!=expected: raise ValueError(f"expected {expected} numeric tokens for {ni}x{nj}, got {len(tokens)}")
-    values=np.asarray([float(v) for v in tokens[2:]],dtype=np.float64)
-    return values[:ni*nj].reshape((nj,ni)),values[ni*nj:].reshape((nj,ni))
-
-def convert(x,y):
-    nj,ni=x.shape; n=ni*nj
-    points=[(float(x[j,i]),float(y[j,i]),0.0) for j in range(nj) for i in range(ni)]
-    points += [(float(x[j,i]),float(y[j,i]),1.0) for j in range(nj) for i in range(ni)]
-    def p(j,i,layer=0): return layer*n+j*ni+i
-    cells=[]
-    for j in range(nj-1):
-        for i in range(ni-1):
-            a,b,c,d=p(j,i),p(j,i+1),p(j+1,i+1),p(j+1,i)
-            cells.append((a,b,c,d,a+n,b+n,c+n,d+n))
-    airfoil=[(p(0,i+1),p(0,i),p(0,i,1),p(0,i+1,1)) for i in range(ni-1)]
-    farfield=[(p(nj-1,i),p(nj-1,i+1),p(nj-1,i+1,1),p(nj-1,i,1)) for i in range(ni-1)]
-    wake_lower=[(p(j,0),p(j+1,0),p(j+1,0,1),p(j,0,1)) for j in range(nj-1)]
-    wake_upper=[(p(j+1,ni-1),p(j,ni-1),p(j,ni-1,1),p(j+1,ni-1,1)) for j in range(nj-1)]
-    front=[tuple(c[:4]) for c in cells]
-    back=[tuple(reversed(tuple(v+n for v in c[:4]))) for c in cells]
-    return points,cells,{"airfoil":airfoil,"farfield":farfield,"wake_lower":wake_lower,"wake_upper":wake_upper,"_front_cells":front,"_back_cells":back}
+    """Read NASA's formatted multigrid PLOT3D, with one 2-D block.
+
+    The first record is nbl, not ni.  Compressed files in the advertised
+    grid family use the same format and are read without manual unpacking.
+    """
+    opener = gzip.open if path.suffix == ".gz" else open
+    with opener(path, "rt", encoding="ascii") as stream:
+        tokens = stream.read().replace(",", " ").split()
+    if len(tokens) < 3:
+        raise ValueError(f"PLOT3D block header is missing: {path}")
+    nblocks = int(tokens[0])
+    if nblocks != 1:
+        raise ValueError(f"only single-block NASA/TMR PLOT3D is supported, got {nblocks}")
+    ni, nj = int(tokens[1]), int(tokens[2])
+    if ni < 5 or nj < 2:
+        raise ValueError(f"invalid NASA/TMR C-grid dimensions: {ni}x{nj}")
+    expected = 3 + 2 * ni * nj
+    if len(tokens) != expected:
+        raise ValueError(f"expected {expected} numeric tokens for {ni}x{nj}, got {len(tokens)}")
+    values = np.asarray(
+        [float(value.replace("D", "E").replace("d", "e")) for value in tokens[3:]],
+        dtype=np.float64,
+    )
+    if not np.isfinite(values).all():
+        raise ValueError("NASA/TMR grid coordinates must be finite")
+    return values[:ni * nj].reshape((nj, ni)), values[ni * nj:].reshape((nj, ni))
+
+
+def convert(x, y):
+    """Join the coincident j=0 wake and extrude with EMPTY front/back.
+
+    NASA's neutral map declares a reversed one-to-one connection between
+    the leading and trailing index intervals on j=0.  Determine that
+    interval from coincident endpoint pairs, including the trailing-edge
+    pair.  No wall-index ratio is guessed from the filename.  Reject a
+    different topology instead of manufacturing a wall in its wake.
+    """
+    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
+    if x.ndim != 2 or x.shape != y.shape or min(x.shape) < 2:
+        raise ValueError("C-grid x/y arrays must have the same two-dimensional shape")
+    if not np.isfinite(x).all() or not np.isfinite(y).all():
+        raise ValueError("C-grid coordinates must be finite")
+    nj, ni = x.shape
+    # Absolute tolerance is in the unit-chord NASA grid coordinates.  It
+    # is deliberately independent of the distant farfield radius.
+    paired_points = 0
+    for i in range(ni // 2):
+        opposite = ni - 1 - i
+        if max(abs(x[0, i] - x[0, opposite]),
+               abs(y[0, i] - y[0, opposite])) > 1e-12:
+            break
+        paired_points += 1
+    first_wall = paired_points - 1
+    last_wall = ni - paired_points
+    if paired_points < 2 or last_wall - first_wall < 3:
+        raise ValueError("expected a reversed coincident j=0 wake and a nondegenerate airfoil")
+
+    # The supplied NASA family is normalized to unit chord.  Check the
+    # independently known geometry before accepting the detected seam: a
+    # partly corrupted wake must not be reclassified as extra solid wall.
+    wall_x = x[0, first_wall:last_wall + 1]
+    wall_y = y[0, first_wall:last_wall + 1]
+    geometry_tolerance = 1e-6  # admits the documented original TMR TE typo.
+    if (abs(wall_x[0] - 1.0) > geometry_tolerance or
+            abs(wall_y[0]) > geometry_tolerance or
+            abs(wall_x.min()) > geometry_tolerance or
+            wall_x.max() > 1.0 + geometry_tolerance or
+            (np.diff(x[0, :paired_points]) >= 0.0).any()):
+        raise ValueError("matched wake must terminate at the unit-chord NACA trailing edge")
+
+    # Merge only the explicitly matched wake vertices.  This gives each
+    # wake face two owners through write_mesh's ordinary connectivity
+    # matching; it is not a periodic/no-slip physical boundary.
+    representatives = np.arange(ni * nj)
+    for i in range(paired_points):
+        representatives[ni - 1 - i] = i
+    unique, inverse = np.unique(representatives, return_inverse=True)
+    nplane = len(unique)
+    flat_x, flat_y = x.ravel(), y.ravel()
+    points = [(float(flat_x[i]), float(flat_y[i]), z)
+              for z in (0.0, 1.0) for i in unique]
+
+    def p(j, i, layer=0):
+        return int(inverse[j * ni + i]) + layer * nplane
+
+    cells = []
+    for j in range(nj - 1):
+        for i in range(ni - 1):
+            lower = (p(j, i), p(j, i + 1), p(j + 1, i + 1), p(j + 1, i))
+            cell = lower + tuple(v + nplane for v in lower)
+            if len(set(cell)) != 8:
+                raise ValueError("wake stitching created a degenerate hexahedron")
+            cells.append(cell)
+    airfoil = [(p(0, i + 1), p(0, i), p(0, i, 1), p(0, i + 1, 1))
+               for i in range(first_wall, last_wall)]
+    farfield = [(p(nj - 1, i), p(nj - 1, i + 1),
+                 p(nj - 1, i + 1, 1), p(nj - 1, i, 1))
+                for i in range(ni - 1)]
+    # i=0 and i=ni-1 are downstream farfield boundaries in the neutral
+    # map.  The actual wake connection is along j=0, handled above.
+    farfield += [(p(j, 0), p(j + 1, 0), p(j + 1, 0, 1), p(j, 0, 1))
+                 for j in range(nj - 1)]
+    farfield += [(p(j + 1, ni - 1), p(j, ni - 1),
+                  p(j, ni - 1, 1), p(j + 1, ni - 1, 1))
+                 for j in range(nj - 1)]
+    front = [tuple(cell[:4]) for cell in cells]
+    back = [tuple(reversed(cell[4:])) for cell in cells]
+    return points, cells, {"airfoil": airfoil, "farfield": farfield,
+                           "front": front, "back": back}
+
 
 def infer_dimensions(path: Path):
-    m=re.search(r"_(\d+)-(\d+)\.p2dfmt(?:\.gz)?$",path.name)
-    if not m: raise ValueError(f"cannot infer TMR grid dimensions from {path.name}")
-    return int(m.group(1)),int(m.group(2))
+    match = re.search(r"_(\d+)-(\d+)\.p2dfmt(?:\.gz)?$", path.name)
+    if not match:
+        raise ValueError(f"cannot infer TMR grid dimensions from {path.name}")
+    return int(match.group(1)), int(match.group(2))
+
 
 def main():
-    ap=argparse.ArgumentParser()
-    ap.add_argument("--input",type=Path,required=True)
-    ap.add_argument("--output",type=Path,required=True)
-    args=ap.parse_args()
-    ni_name,nj_name=infer_dimensions(args.input)
-    x,y=read_plot3d_2d(args.input)
-    nj,ni=x.shape
-    if (ni,nj)!=(ni_name,nj_name): raise ValueError(f"filename declares {ni_name}x{nj_name}, file contains {ni}x{nj}")
-    points,cells,boundaries=convert(x,y)
-    write_mesh(args.output,points,cells,boundaries)
-    with h5py.File(args.output,"r+") as h5:
-        h5.attrs["validation_source"]="NASA TMR NACA0012"
-        h5.attrs["validation_source_url"]=TMR_URL
-        h5.attrs["validation_archive_url"]=TMR_ARCHIVE_URL
-        h5.attrs["validation_grid"]=f"{ni}x{nj}"
-        h5.attrs["validation_topology"]="C-grid"
-        h5.attrs["validation_physics"]="N9-S8 laminar incompressible campaign"
-        h5.attrs["validation_status"]="imported-provenance-only"
-    print(f"Imported {args.input} -> {args.output} ({ni}x{nj})")
+    parser = argparse.ArgumentParser(description=__doc__)
+    parser.add_argument("--input", type=Path, required=True)
+    parser.add_argument("--output", type=Path, required=True)
+    args = parser.parse_args()
+    ni_name, nj_name = infer_dimensions(args.input)
+    x, y = read_plot3d_2d(args.input)
+    nj, ni = x.shape
+    if (ni, nj) != (ni_name, nj_name):
+        raise ValueError(f"filename declares {ni_name}x{nj_name}, file contains {ni}x{nj}")
+    points, cells, boundaries = convert(x, y)
+    write_mesh(args.output, points, cells, boundaries)
+    with h5py.File(args.output, "r+") as h5:
+        h5.attrs["validation_source"] = "NASA TMR NACA0012"
+        h5.attrs["validation_source_url"] = TMR_URL
+        h5.attrs["validation_archive_url"] = TMR_ARCHIVE_URL
+        h5.attrs["validation_source_sha256"] = hashlib.sha256(args.input.read_bytes()).hexdigest()
+        h5.attrs["validation_grid"] = f"{ni}x{nj}"
+        h5.attrs["validation_topology"] = "C-grid; j=0 reversed wake stitched"
+        h5.attrs["validation_wall_faces"] = len(boundaries["airfoil"])
+        h5.attrs["validation_wake_pairs"] = (ni - 1 - len(boundaries["airfoil"])) // 2
+        h5.attrs["validation_physics"] = "N9-S8 laminar incompressible campaign"
+        h5.attrs["validation_status"] = "imported-provenance-only"
+    print(f"Imported {args.input} -> {args.output} ({ni}x{nj}); "
+          f"airfoil_faces={len(boundaries['airfoil'])}; physical qualification NOT RUN")
     return 0
 
-if __name__=="__main__": raise SystemExit(main())
+
+if __name__ == "__main__":
+    raise SystemExit(main())
diff --git a/src/cfdx/core/linalg/cg_solver.h b/src/cfdx/core/linalg/cg_solver.h
index c21ecaa3..16daccf1 100644
--- a/src/cfdx/core/linalg/cg_solver.h
+++ b/src/cfdx/core/linalg/cg_solver.h
@@ -91,6 +91,26 @@ inline SolverResult solve_cg_impl(
             return result;
         }
     if (null_space) null_space->remove(x);
+    // A valid system can fail before the usual convergence/final-residual
+    // reporting path (preconditioner failure, non-positive curvature, or
+    // rejected residual replacement). Preserve its current b-A*x diagnostic
+    // instead of returning SolverResult's default zero residuals. This helper
+    // is used only on failure; it does not change an iterate or a stopping gate.
+    const auto not_applicable = [&](std::size_t iterations) {
+        result.status = SolverStatus::NOT_APPLICABLE;
+        result.iterations = iterations;
+        Vector failure_residual(n);
+        mixed_precision_true_residual(A, b, x, failure_residual);
+        result.residual =
+            krylov_norm2(failure_residual, SolverPrecision::FP64, reduction);
+        const double rhs_norm =
+            krylov_norm2(b, SolverPrecision::FP64, reduction);
+        result.residual_relative = rhs_norm > 0.0
+            ? result.residual / rhs_norm
+            : (result.residual == 0.0
+                ? 0.0 : std::numeric_limits<double>::infinity());
+        return result;
+    };
     const double* Av = A.values_data();
     const bool mp32 = precision.enabled && precision.operator_precision == SolverPrecision::FP32;
     const SolverPrecision redp = precision.enabled ? precision.reduction_precision : SolverPrecision::FP64;
@@ -102,8 +122,7 @@ inline SolverResult solve_cg_impl(
     std::vector<double> M(n, 0.0);
     if (preconditioner) {
         if (setup_preconditioner && !preconditioner->setup(A)) {
-            result.status = SolverStatus::NOT_APPLICABLE;
-            return result;
+            return not_applicable(0);
         }
     } else {
         for (std::size_t i = 0; i < n; ++i) {
@@ -114,8 +133,7 @@ inline SolverResult solve_cg_impl(
                 }
             }
             if (!(M[i] > 0.0) || !std::isfinite(M[i])) {
-                result.status = SolverStatus::NOT_APPLICABLE;
-                return result;
+                return not_applicable(0);
             }
         }
     }
@@ -144,8 +162,7 @@ inline SolverResult solve_cg_impl(
         return true;
     };
     if (!apply_preconditioner()) {
-        result.status = SolverStatus::NOT_APPLICABLE;
-        return result;
+        return not_applicable(0);
     }
 
     // p = z
@@ -210,9 +227,7 @@ inline SolverResult solve_cg_impl(
         }
         const double pAp = krylov_dot(p_vector, Ap_vector, redp, reduction);
         if (!(pAp > 0.0) || !std::isfinite(pAp)) {
-            result.status = SolverStatus::NOT_APPLICABLE;
-            result.iterations = iter;
-            return result;
+            return not_applicable(iter);
         }
         const double alpha = rsold / pAp;
         if (!std::isfinite(alpha)) { result.status=SolverStatus::DIVERGED; result.iterations=iter; return result; }
@@ -233,9 +248,7 @@ inline SolverResult solve_cg_impl(
 
         // z = M^{-1} r
         if (!apply_preconditioner()) {
-            result.status = SolverStatus::NOT_APPLICABLE;
-            result.iterations = iter;
-            return result;
+            return not_applicable(iter);
         }
 
         // rsnew = r · z
@@ -347,9 +360,7 @@ inline SolverResult solve_cg_impl(
                 for (std::size_t i = 0; i < n; ++i) r[i] = rv(i);
                 if (null_space) null_space->remove(r);
                 if (!apply_preconditioner()) {
-                    result.status = SolverStatus::NOT_APPLICABLE;
-                    result.iterations = iter;
-                    return result;
+                    return not_applicable(iter);
                 }
                 for (std::size_t i = 0; i < n; ++i) {
                     r_vector(i) = r[i];
@@ -358,9 +369,7 @@ inline SolverResult solve_cg_impl(
                 const double replaced_rsnew =
                     krylov_dot(r_vector, z_vector, redp, reduction);
                 if (!std::isfinite(replaced_rsnew) || !(replaced_rsnew > 0.0)) {
-                    result.status = SolverStatus::NOT_APPLICABLE;
-                    result.iterations = iter;
-                    return result;
+                    return not_applicable(iter);
                 }
                 // The replacement recomputes b-Ax after the normal convergence
                 // check. Re-test the true residual here so a converged iterate
diff --git a/src/cfdx/core/mesh/boundary.cpp b/src/cfdx/core/mesh/boundary.cpp
index 945885ca..e0aedb00 100644
--- a/src/cfdx/core/mesh/boundary.cpp
+++ b/src/cfdx/core/mesh/boundary.cpp
@@ -10,36 +10,25 @@ bool BoundaryPatches::is_consistent_with_mesh(const Mesh& m) const {
     if (!is_consistent(m.n_faces())) return false;
 
     const FaceOwnership& own = m.ownership();
-    std::vector<bool> is_boundary_face(m.n_faces(), false);
+    std::vector<bool> covered(m.n_faces(), false);
 
-    // Marquer les faces de frontière selon owner/neighbour
-    for (std::size_t f = 0; f < m.n_faces(); ++f) {
-        is_boundary_face[f] = (own.neighbour(f) == FaceOwnership::BOUNDARY);
-    }
-
-    // Vérifier que chaque face dans un patch est bien une face de frontière
+    // is_consistent() has already rejected duplicate/out-of-range ids.
+    // Mark each listed boundary face once, retaining the internal-face gate.
     for (const auto& p : patches_) {
         for (const auto fid : p.face_ids) {
-            if (!is_boundary_face[fid]) {
-                return false;  // Face interne dans un patch de frontière
+            if (own.neighbour(fid) != FaceOwnership::BOUNDARY) {
+                return false;
             }
+            covered[fid] = true;
         }
     }
 
-    // Vérifier que toutes les faces de frontière sont dans un patch
+    // Membership must be tested in O(1), not by rescanning every patch for
+    // every boundary face. In a one-cell-thick 2-D extrusion, the two EMPTY
+    // patches contain O(n_cells) faces and the former search was quadratic.
     for (std::size_t f = 0; f < m.n_faces(); ++f) {
-        if (is_boundary_face[f]) {
-            bool found = false;
-            for (const auto& p : patches_) {
-                for (const auto fid : p.face_ids) {
-                    if (fid == static_cast<FaceIndex>(f)) {
-                        found = true;
-                        break;
-                    }
-                }
-                if (found) break;
-            }
-            if (!found) return false;  // Face de frontière orpheline
+        if (own.neighbour(f) == FaceOwnership::BOUNDARY && !covered[f]) {
+            return false;
         }
     }
 
diff --git a/src/cfdx/core/numerics/temporal.h b/src/cfdx/core/numerics/temporal.h
index f2a70743..75cce998 100644
--- a/src/cfdx/core/numerics/temporal.h
+++ b/src/cfdx/core/numerics/temporal.h
@@ -135,6 +135,17 @@ inline Field<double, Location::CELL> advance_time(
     if (ctx && (ctx->n_cells != n_cells || ctx->field_dim != dim)) {
         throw std::invalid_argument("advance_time: time-integration context does not match field shape");
     }
+    const auto require_finite = [&](const Field<double, Location::CELL>& field,
+                                    const char* message) {
+        for (std::size_t d = 0; d < dim; ++d) {
+            for (std::size_t cell = 0; cell < n_cells; ++cell) {
+                if (!std::isfinite(field.component_data(d)[cell])) {
+                    throw std::runtime_error(message);
+                }
+            }
+        }
+    };
+    require_finite(phi, "advance_time: input field contains a non-finite value");
     Field<double, Location::CELL> phi_new(
         n_cells, "phi_new", phi.metadata().unit, dim);
 
@@ -160,6 +171,7 @@ inline Field<double, Location::CELL> advance_time(
         predictor();
         for (int iteration = 0; iteration < max_iterations; ++iteration) {
             rhs_func(phi_new, rhs_iter);
+            require_finite(rhs_iter, "advance_time: implicit RHS contains a non-finite value");
             double max_delta = 0.0;
             for (std::size_t d = 0; d < dim; ++d) {
                 const double* old_data = phi_new.component_data(d);
@@ -182,6 +194,9 @@ inline Field<double, Location::CELL> advance_time(
                         candidate = (a0 * phi_data[cell] - a1 * prev[cell] +
                              b * rhs_data[cell]) / denom;
                     }
+                    if (!std::isfinite(candidate)) {
+                        throw std::runtime_error("advance_time: implicit candidate contains a non-finite value");
+                    }
                     max_delta = std::max(max_delta, std::abs(candidate - new_data[cell]));
                     new_data[cell] = candidate;
                 }
@@ -242,6 +257,7 @@ inline Field<double, Location::CELL> advance_time(
                 }
                 for (int iteration = 0; iteration < max_iterations; ++iteration) {
                     rhs_func(phi_new, rhs_iter);
+                    require_finite(rhs_iter, "advance_time: BDF2 bootstrap RHS contains a non-finite value");
                     double max_delta = 0.0;
                     for (std::size_t d = 0; d < dim; ++d) {
                         const double* phi_data = phi.component_data(d);
@@ -249,6 +265,9 @@ inline Field<double, Location::CELL> advance_time(
                         double* new_data = phi_new.component_data(d);
                         for (std::size_t cell = 0; cell < n_cells; ++cell) {
                             const double candidate = phi_data[cell] + dt * rhs_data[cell];
+                            if (!std::isfinite(candidate)) {
+                                throw std::runtime_error("advance_time: BDF2 bootstrap candidate contains a non-finite value");
+                            }
                             max_delta = std::max(max_delta, std::abs(candidate - new_data[cell]));
                             new_data[cell] = candidate;
                         }
@@ -264,6 +283,8 @@ inline Field<double, Location::CELL> advance_time(
             break;
     }
 
+    // Reject arithmetic overflow and invalid iterates before committing history.
+    require_finite(phi_new, "advance_time: output field contains a non-finite value");
     if (ctx) {
         ctx->dt_prev = dt;
         ctx->shift(phi_new);
@@ -450,4 +471,4 @@ inline void apply_local_time_stepping(
 }
 
 }  // namespace core
-}  // namespace cfdx
\ No newline at end of file
+}  // namespace cfdx
diff --git a/src/cfdx/physics/steady_incompressible_solver.h b/src/cfdx/physics/steady_incompressible_solver.h
index 2bfa725d..32594af0 100644
--- a/src/cfdx/physics/steady_incompressible_solver.h
+++ b/src/cfdx/physics/steady_incompressible_solver.h
@@ -3392,13 +3392,17 @@ inline IncompressibleSolveResult solve_steady_incompressible(
             }
             pressure_residual = rp.residual_relative;
             pressure_iterations = rp.iterations;
-            if (rp.status != cfdx::core::SolverStatus::CONVERGED)
-                throw NonlinearRetryableFailure(
-                    "solve_steady_incompressible: pressure-correction solve did not converge "
-                    "(status=" + std::to_string(static_cast<int>(rp.status)) +
-                    ", iterations=" + std::to_string(rp.iterations) +
-                    ", residual=" + std::to_string(rp.residual) +
-                    ", relative=" + std::to_string(rp.residual_relative) + ")");
+            if (rp.status != cfdx::core::SolverStatus::CONVERGED) {
+                std::ostringstream message;
+                message << "solve_steady_incompressible: pressure-correction solve did not converge "
+                        << "(status=" << static_cast<int>(rp.status)
+                        << ", iterations=" << rp.iterations
+                        << ", residual=" << std::scientific
+                        << std::setprecision(std::numeric_limits<double>::max_digits10)
+                        << rp.residual
+                        << ", relative=" << rp.residual_relative << ")";
+                throw NonlinearRetryableFailure(message.str());
+            }
 
             // Relax the physical pressure exactly once. The reference is a
             // gauge: enforce it by a uniform shift, never by overwriting one
diff --git a/src/cfdx/physics/wall_distance.h b/src/cfdx/physics/wall_distance.h
index a88dc057..d96d3a72 100644
--- a/src/cfdx/physics/wall_distance.h
+++ b/src/cfdx/physics/wall_distance.h
@@ -41,7 +41,8 @@ enum class WallDistanceConvergenceStatus {
     MAX_ITER,
     RESIDUAL_TOO_HIGH,
     RESIDUAL_STAGNATION,
-    NOT_REPORTED
+    NOT_REPORTED,
+    INVALID_FIELD
 };
 
 struct WallDistanceResult {
@@ -951,8 +952,10 @@ inline double godunov_gradient_at(const WallDistanceBvh& bvh,
             if(k>0) face(g.index(i,j,k-1),true,g.spacing.z,4);
             if(k+1<g.nz) face(g.index(i,j,k+1),false,g.spacing.z,5);
         }
-        return std::sqrt(std::max(dm,0.0)*std::max(dm,0.0)
-                       + std::min(dp,0.0)*std::min(dp,0.0));
+        // The upwind Eikonal Hamiltonian chooses one direction per axis.
+        // Adding both squared slopes double-counts a medial ridge and
+        // disagrees with eikonal_update(), which uses the nearest neighbour.
+        return std::max({dm,-dp,0.0});
     };
     const double gx=component(0), gy=component(1), gz=component(2);
     return std::sqrt(gx*gx+gy*gy+gz*gz);
@@ -997,8 +1000,10 @@ inline double godunov_gradient_at(const std::vector<double>& d,const WallDistanc
                 else if(std::isfinite(d[p])) dp=(d[p]-d[id])/g.spacing.z;
             }
         }
-        return std::sqrt(std::max(dm,0.0)*std::max(dm,0.0)
-                       + std::min(dp,0.0)*std::min(dp,0.0));
+        // The upwind Eikonal Hamiltonian chooses one direction per axis.
+        // Adding both squared slopes double-counts a medial ridge and
+        // disagrees with eikonal_update(), which uses the nearest neighbour.
+        return std::max({dm,-dp,0.0});
     };
     const double gx=component(0), gy=component(1), gz=component(2);
     return std::sqrt(gx*gx+gy*gy+gz*gz);
@@ -2079,13 +2084,28 @@ inline double wall_distance_pde_residual_inf(WallDistanceMethod method,
                                                         const std::vector<unsigned char>* fixed,
                                                         double advection_diffusion_gamma,
                                                         const WallDistanceBvh* bvh) {
+    if(d.size()!=g.points.size() || g.solid.size()!=g.points.size() ||
+       (fixed && fixed->size()!=g.points.size()))
+        throw std::invalid_argument("wall_distance_pde_residual_inf: field/grid size mismatch");
     double rmax=0.0;
+    std::size_t fluid_samples=0;
+    const auto accumulate=[&](double residual) {
+        rmax=std::isfinite(residual)
+            ? std::max(rmax,std::abs(residual))
+            : std::numeric_limits<double>::infinity();
+    };
     for(std::size_t id=0;id<d.size();++id) {
-        if(g.solid[id] || !std::isfinite(d[id]) || (fixed && (*fixed)[id])) continue;
+        if(g.solid[id]) continue;
+        ++fluid_samples;
+        // A missing or corrupt value is not a zero PDE residual. Check fixed
+        // seeds too: otherwise an invalid boundary field can certify a solve.
+        if(!std::isfinite(d[id]) || d[id]<0.0)
+            return std::numeric_limits<double>::infinity();
+        if(fixed && (*fixed)[id]) continue;
         if(method==WallDistanceMethod::EIKONAL) {
             const double grad=bvh ? godunov_gradient_at(*bvh,d,g,id)
                                   : godunov_gradient_at(d,g,id);
-            rmax=std::max(rmax,std::abs(grad-1.0));
+            accumulate(grad-1.0);
             continue;
         }
         if(method==WallDistanceMethod::HAMILTON_JACOBI ||
@@ -2095,7 +2115,7 @@ inline double wall_distance_pde_residual_inf(WallDistanceMethod method,
                                   : godunov_gradient_at(d,g,id);
             const double lap=bvh ? laplacian_at(*bvh,d,g,id)
                                  : laplacian_at(d,g,id);
-            rmax=std::max(rmax,std::abs(grad-1.0-gamma*lap));
+            accumulate(grad-1.0-gamma*lap);
             continue;
         }
         if(method==WallDistanceMethod::ADVECTION_DIFFUSION) {
@@ -2117,30 +2137,34 @@ inline double wall_distance_pde_residual_inf(WallDistanceMethod method,
             };
             const double gx=one(0), gy=one(1), gz=one(2);
             const double gn=std::sqrt(gx*gx+gy*gy+gz*gz);
-            if(!(gn>1e-14)) continue;
-            const double ux=gx/gn, uy=gy/gn, uz=gz/gn;
-            auto upwind=[&](double u,int axis)->double {
-                if(axis==0) {
-                    if(u>0.0 && i>0) return u*(d[id]-value(g.index(i-1,j,k)))/g.spacing.x;
-                    if(u<0.0 && i+1<g.nx) return u*(value(g.index(i+1,j,k))-d[id])/g.spacing.x;
-                } else if(axis==1) {
-                    if(u>0.0 && j>0) return u*(d[id]-value(g.index(i,j-1,k)))/g.spacing.y;
-                    if(u<0.0 && j+1<g.ny) return u*(value(g.index(i,j+1,k))-d[id])/g.spacing.y;
-                } else {
-                    if(u>0.0 && k>0) return u*(d[id]-value(g.index(i,j,k-1)))/g.spacing.z;
-                    if(u<0.0 && k+1<g.nz) return u*(value(g.index(i,j,k+1))-d[id])/g.spacing.z;
-                }
-                return 0.0;
+            const double ux=gn>1e-14 ? gx/gn : 0.0;
+            const double uy=gn>1e-14 ? gy/gn : 0.0;
+            const double uz=gn>1e-14 ? gz/gn : 0.0;
+            // Certify the same discrete equation as the in-place update in
+            // advection_diffusion_distance(): normalized signed upwind U,
+            // first-order face transport, and gamma/(h*delta) diffusion.
+            // Godunov |grad d| plus a Shortley-Weller Laplacian is a different
+            // operator and can stay O(1) after this solve reaches roundoff.
+            double lhs=0.0;
+            const auto face=[&](std::size_t q,bool exists,bool minus_side,double h,double u) {
+                if(!exists) return;
+                const double delta=(bvh && g.solid[q])
+                    ? poisson_wall_offset(*bvh,g.points[id],g.points[q],h) : h;
+                const double difference=d[id]-value(q);
+                lhs+=advection_diffusion_gamma*difference/(h*delta);
+                if((minus_side && u>0.0) || (!minus_side && u<0.0))
+                    lhs+=std::abs(u)*difference/delta;
             };
-            const double transport=bvh ? godunov_gradient_at(*bvh,d,g,id)
-                                        : upwind(ux,0)+upwind(uy,1)+upwind(uz,2);
-            const double lap=bvh ? laplacian_at(*bvh,d,g,id)
-                                 : laplacian_at(d,g,id);
-            rmax=std::max(rmax,std::abs(
-                transport-1.0-advection_diffusion_gamma*lap));
+            face(i>0?g.index(i-1,j,k):0,i>0,true,g.spacing.x,ux);
+            face(i+1<g.nx?g.index(i+1,j,k):0,i+1<g.nx,false,g.spacing.x,ux);
+            face(j>0?g.index(i,j-1,k):0,j>0,true,g.spacing.y,uy);
+            face(j+1<g.ny?g.index(i,j+1,k):0,j+1<g.ny,false,g.spacing.y,uy);
+            face(k>0?g.index(i,j,k-1):0,k>0,true,g.spacing.z,uz);
+            face(k+1<g.nz?g.index(i,j,k+1):0,k+1<g.nz,false,g.spacing.z,uz);
+            accumulate(lhs-1.0);
         }
     }
-    return rmax;
+    return fluid_samples>0 ? rmax : std::numeric_limits<double>::infinity();
 }
 
 // Default arguments deliberately not applied: callers must state the method
@@ -2251,8 +2275,21 @@ inline WallDistanceResult compute_wall_distance(WallDistanceMethod method,const
             break;
         }
     }
-    for(std::size_t i=0;i<r.distance.size();++i)
-        if(g.solid[i] || !std::isfinite(r.distance[i]) || r.distance[i]<0.0) r.valid[i]=0;
+    std::size_t fluid_samples=0, invalid_fluid_samples=0;
+    for(std::size_t i=0;i<r.distance.size();++i) {
+        if(g.solid[i]) { r.valid[i]=0; continue; }
+        ++fluid_samples;
+        if(!std::isfinite(r.distance[i]) || r.distance[i]<0.0) {
+            r.valid[i]=0;
+            ++invalid_fluid_samples;
+        }
+    }
+    if(fluid_samples==0 || invalid_fluid_samples!=0) {
+        r.converged=false;
+        r.residual_inf=std::numeric_limits<double>::infinity();
+        r.stopping_reason=fluid_samples==0 ? "no_fluid_samples" : "invalid_distance_field";
+        r.convergence_status=WallDistanceConvergenceStatus::INVALID_FIELD;
+    }
     return r;
 }
 
@@ -2272,24 +2309,41 @@ struct WallDistanceBenchmarkMetrics {
     double near_wall_l2_relative{0.0};
     double monotonicity_violations{0.0};
     std::size_t samples{0};
+    std::size_t invalid_samples{0};
+    std::size_t near_wall_samples{0};
 };
 
 inline WallDistanceBenchmarkMetrics compare_wall_distance(const WallDistanceGrid& g,
                                                           const std::vector<double>& ref,
                                                           const std::vector<double>& value,
                                                           double near_wall) {
+    if(ref.size()!=g.points.size() || value.size()!=g.points.size() ||
+       g.solid.size()!=g.points.size() || g.neighbours.size()!=g.points.size())
+        throw std::invalid_argument("compare_wall_distance: field/grid size mismatch");
+    if(!std::isfinite(near_wall) || near_wall<0.0)
+        throw std::invalid_argument("compare_wall_distance: invalid near-wall width");
     double sum2=0.0, ref2=0.0, maxe=0.0, max_ref=0.0, near_sum2=0.0, near_ref2=0.0;
-    std::size_t n=0, near_n=0, violations=0;
+    std::size_t n=0, near_n=0, violations=0, invalid=0;
     for(std::size_t i=0;i<ref.size();++i) {
-        if(g.solid[i] || !std::isfinite(ref[i]) || !std::isfinite(value[i])) continue;
+        if(g.solid[i]) continue;
+        if(!std::isfinite(ref[i]) || !std::isfinite(value[i]) || ref[i]<0.0 || value[i]<0.0) {
+            ++invalid;
+            continue;
+        }
         const double e=std::abs(value[i]-ref[i]); sum2+=e*e; ref2+=ref[i]*ref[i]; maxe=std::max(maxe,e); max_ref=std::max(max_ref,ref[i]); ++n;
         if(ref[i]<=near_wall) { near_sum2+=e*e; near_ref2+=ref[i]*ref[i]; ++near_n; }
         for(auto nb:g.neighbours[i]) if(!g.solid[nb] && ref[nb]>ref[i] && value[nb]+1e-12<value[i]) ++violations;
     }
+    const double unavailable=std::numeric_limits<double>::infinity();
+    // Never turn an empty set or a partially invalid field into an accuracy
+    // success. Counts are retained so the caller can diagnose the failure.
+    if(n==0 || invalid!=0 || !std::isfinite(sum2) || !std::isfinite(ref2) ||
+       !std::isfinite(near_sum2) || !std::isfinite(near_ref2))
+        return {unavailable,unavailable,unavailable,static_cast<double>(violations),n,invalid,near_n};
     return {std::sqrt(sum2/std::max(1e-30,ref2)),
             maxe/std::max(1e-30,max_ref),
-            std::sqrt(near_sum2/std::max(1e-30,near_ref2)),
-            static_cast<double>(violations),n};
+            near_n>0 ? std::sqrt(near_sum2/std::max(1e-30,near_ref2)) : unavailable,
+            static_cast<double>(violations),n,invalid,near_n};
 }
 
 } // namespace cfdx::physics
diff --git a/tests/numerical/test_temporal.cpp b/tests/numerical/test_temporal.cpp
index 31421da4..6fab26a2 100644
--- a/tests/numerical/test_temporal.cpp
+++ b/tests/numerical/test_temporal.cpp
@@ -12,50 +12,54 @@
 using namespace cfdx::core;
 using namespace cfdx::testing;
 
-// Helper: build a simple 1D mesh with 4 cells (thin 3D slab)
+// Four closed hexahedra: dx=1/4, dy=dz=1/10, V=1/400 each.
+// The former coplanar fixture had zero volume and only exercised the 1e30
+// fallback in the CFL/local-step routines.
 Mesh make_test_mesh() {
-    Mesh m;
-
-    m.points().resize(10);
-    m.points().set(0, 0.0, 0.0, 0.0);
-    m.points().set(1, 0.25, 0.0, 0.0);
-    m.points().set(2, 0.5, 0.0, 0.0);
-    m.points().set(3, 0.75, 0.0, 0.0);
-    m.points().set(4, 1.0, 0.0, 0.0);
-    m.points().set(5, 0.0, 0.0, 0.1);
-    m.points().set(6, 0.25, 0.0, 0.1);
-    m.points().set(7, 0.5, 0.0, 0.1);
-    m.points().set(8, 0.75, 0.0, 0.1);
-    m.points().set(9, 1.0, 0.0, 0.1);
-
-    m.faces().push_face({0, 5, 6, 1});
-    m.faces().push_face({1, 6, 7, 2});
-    m.faces().push_face({2, 7, 8, 3});
-    m.faces().push_face({3, 8, 9, 4});
-    m.faces().push_face({0, 1, 6, 5});
-    m.faces().push_face({4, 9, 8, 3});
-
-    m.ownership().resize(6);
-    m.ownership().set_owner(0, 0); m.ownership().set_neighbour(0, 1);
-    m.ownership().set_owner(1, 1); m.ownership().set_neighbour(1, 2);
-    m.ownership().set_owner(2, 2); m.ownership().set_neighbour(2, 3);
-    m.ownership().set_owner(3, 3); m.ownership().set_neighbour(3, FaceOwnership::BOUNDARY);
-    m.ownership().set_owner(4, 0); m.ownership().set_neighbour(4, FaceOwnership::BOUNDARY);
-    m.ownership().set_owner(5, 3); m.ownership().set_neighbour(5, FaceOwnership::BOUNDARY);
-
-    m.cells().push_cell({0, 4});
-    m.cells().push_cell({0, 1});
-    m.cells().push_cell({1, 2});
-    m.cells().push_cell({2, 3, 5});
-
-    BoundaryPatches bp;
-    bp.add_patch("inlet", PatchType::INLET);
-    bp.patch(0).face_ids = {4};
-    bp.add_patch("outlet", PatchType::OUTLET);
-    bp.patch(1).face_ids = {5};
-    m.set_boundary(bp);
-
-    return m;
+    constexpr std::size_t cells=4;
+    Mesh mesh;
+    mesh.points().resize(4*(cells+1));
+    for(std::size_t i=0;i<=cells;++i) {
+        const auto b=4*i;
+        const double x=.25*static_cast<double>(i);
+        mesh.points().set(b,  x,0.0,0.0);
+        mesh.points().set(b+1,x,.1,0.0);
+        mesh.points().set(b+2,x,.1,.1);
+        mesh.points().set(b+3,x,0.0,.1);
+        if(i==0) mesh.faces().push_face({b,b+3,b+2,b+1});
+        else mesh.faces().push_face({b,b+1,b+2,b+3});
+    }
+    for(std::size_t i=0;i<cells;++i) {
+        const auto b=4*i;
+        const auto side=mesh.n_faces();
+        mesh.faces().push_face({b,b+4,b+7,b+3});
+        mesh.faces().push_face({b+1,b+2,b+6,b+5});
+        mesh.faces().push_face({b,b+1,b+5,b+4});
+        mesh.faces().push_face({b+3,b+7,b+6,b+2});
+        mesh.cells().push_cell({i,i+1,side,side+1,side+2,side+3});
+    }
+    mesh.ownership().resize(mesh.n_faces());
+    for(std::size_t i=0;i<=cells;++i) {
+        mesh.ownership().set_owner(i,i==0 ? 0 : i-1);
+        mesh.ownership().set_neighbour(i,i>0 && i<cells
+            ? static_cast<CellIndex>(i) : FaceOwnership::BOUNDARY);
+    }
+    for(std::size_t i=0;i<cells;++i)
+        for(std::size_t side=0;side<4;++side) {
+            const auto face=cells+1+4*i+side;
+            mesh.ownership().set_owner(face,i);
+            mesh.ownership().set_neighbour(face,FaceOwnership::BOUNDARY);
+        }
+    BoundaryPatches patches;
+    patches.add_patch("inlet",PatchType::INLET);
+    patches.patch(0).face_ids={0};
+    patches.add_patch("outlet",PatchType::OUTLET);
+    patches.patch(1).face_ids={cells};
+    patches.add_patch("walls",PatchType::WALL);
+    for(std::size_t face=cells+1;face<mesh.n_faces();++face)
+        patches.patch(2).face_ids.push_back(face);
+    mesh.set_boundary(patches);
+    return mesh;
 }
 
 // Simple RHS for testing: dφ/dt = -φ (decay)
@@ -71,37 +75,47 @@ void decay_rhs(const Field<double, Location::CELL>& phi, Field<double, Location:
 
 int main() {
     run_case("euler_explicit_decay", [&]() {
-        Mesh mesh = make_test_mesh();
         Field<double, Location::CELL> phi(4, "phi", "1", 1);
         phi.fill(1.0);
 
         double dt = 0.1;
-        auto phi_new = advance_time(TimeScheme::EULER_EXPLICIT, phi, dt, decay_rhs);
+        auto phi_new = advance_time(phi, dt, decay_rhs, TimeScheme::EULER_EXPLICIT);
 
         // Exact: φ(t) = φ0 * exp(-t)
         // Euler: φ^{n+1} = φ^n * (1 - dt)
         EXPECT_NEAR(phi_new(0), 1.0 * (1.0 - dt), 1e-12);
     });
 
+    run_case("euler_implicit_decay", []() {
+        Field<double,Location::CELL> phi(4,"phi","1",2);
+        for(std::size_t cell=0;cell<phi.size();++cell) {
+            phi(cell,0)=1.0+.1*static_cast<double>(cell);
+            phi(cell,1)=.5+.2*static_cast<double>(cell);
+        }
+        constexpr double dt=.1;
+        const auto result=advance_time(phi,dt,decay_rhs,TimeScheme::EULER_IMPLICIT);
+        for(std::size_t cell=0;cell<phi.size();++cell)
+            for(std::size_t component=0;component<phi.dimension();++component)
+                EXPECT_NEAR(result(cell,component),phi(cell,component)/(1.0+dt),1e-12);
+    });
+
     run_case("crank_nicolson_decay", [&]() {
-        Mesh mesh = make_test_mesh();
         Field<double, Location::CELL> phi(4, "phi", "1", 1);
         phi.fill(1.0);
 
         double dt = 0.1;
-        auto phi_new = advance_time(TimeScheme::CRANK_NICOLSON, phi, dt, decay_rhs);
-
-        // CN predictor-corrector:
-        // Predictor: φ* = φ^n + Δt * RHS(φ^n) = 1.0 + 0.1 * (-1.0) = 0.9
-        // Corrector: φ^{n+1} = φ^n + 0.5*Δt * (RHS(φ^n) + RHS(φ*))
-        // = 1.0 + 0.5*0.1 * (-1.0 + -0.9) = 1.0 - 0.05 * 1.9 = 1.0 - 0.095 = 0.905
-        double phi_star = 1.0 + dt * (-1.0);
-        double expected = 1.0 + 0.5 * dt * (-1.0 + -phi_star);
-        EXPECT_NEAR(phi_new(0), expected, 1e-12);
+        auto phi_new = advance_time(phi, dt, decay_rhs, TimeScheme::CRANK_NICOLSON);
+
+        // Solve (1+dt/2) phi^{n+1} = (1-dt/2) phi^n.
+        // The former oracle was Heun's explicit predictor-corrector, not CN.
+        const double expected=(1.0-.5*dt)/(1.0+.5*dt);
+        for(std::size_t cell=0;cell<phi.size();++cell) {
+            EXPECT_NEAR(phi_new(cell),expected,1e-12);
+            EXPECT_NEAR(phi_new(cell)-phi(cell)+.5*dt*(phi(cell)+phi_new(cell)),0.0,1e-12);
+        }
     });
 
     run_case("bdf2_requires_history", [&]() {
-        Mesh mesh = make_test_mesh();
         Field<double, Location::CELL> phi(4, "phi", "1", 1);
         phi.fill(1.0);
 
@@ -109,16 +123,40 @@ int main() {
         ctx.initialize(phi);
 
         double dt = 0.1;
-        // First step (falls back to Euler)
-        auto phi1 = advance_time(TimeScheme::BDF2, phi, dt, decay_rhs, &ctx);
-        EXPECT_NEAR(phi1(0), 1.0 * (1.0 - dt), 1e-12);
-
-        // Second step (uses BDF2 explicit):
-        // φ2 = (4φ1 - φ0 + 2Δt * RHS(φ1)) / 3
-        // RHS(φ1) = -φ1
-        auto phi2 = advance_time(TimeScheme::BDF2, phi1, dt, decay_rhs, &ctx);
-        double expected = (4.0 * phi1(0) - 1.0 + 2.0 * dt * (-phi1(0))) / 3.0;
-        EXPECT_NEAR(phi2(0), expected, 1e-12);
+        // Bootstrap solves (1+dt) phi1 = phi0 by implicit Euler.
+        EXPECT_FALSE(ctx.has_prev);
+        const auto phi1=advance_time(phi,dt,decay_rhs,TimeScheme::BDF2,&ctx);
+        EXPECT_TRUE(ctx.has_prev);
+        EXPECT_NEAR(ctx.dt_prev,dt,1e-15);
+        for(std::size_t cell=0;cell<phi.size();++cell) {
+            EXPECT_NEAR(phi1(cell),1.0/(1.0+dt),1e-12);
+            EXPECT_NEAR(ctx.phi_prev(cell),phi(cell),1e-15);
+            EXPECT_NEAR(ctx.phi_curr(cell),phi1(cell),1e-15);
+        }
+
+        // BDF2: (3+2dt) phi2 = 4 phi1 - phi0.
+        const auto phi2=advance_time(phi1,dt,decay_rhs,TimeScheme::BDF2,&ctx);
+        for(std::size_t cell=0;cell<phi.size();++cell) {
+            const double expected=(4.0*phi1(cell)-phi(cell))/(3.0+2.0*dt);
+            EXPECT_NEAR(phi2(cell),expected,1e-12);
+            EXPECT_NEAR((3.0+2.0*dt)*phi2(cell)-4.0*phi1(cell)+phi(cell),0.0,1e-12);
+            EXPECT_NEAR(ctx.phi_prev(cell),phi1(cell),1e-15);
+            EXPECT_NEAR(ctx.phi_curr(cell),phi2(cell),1e-15);
+        }
+    });
+
+    run_case("bdf2_variable_step_implicit_decay", []() {
+        Field<double,Location::CELL> phi(4,"phi","1",1);
+        phi.fill(1.0);
+        TimeIntegrationContext ctx(4,1,"phi");
+        ctx.initialize(phi);
+        const auto first=advance_time(phi,.1,decay_rhs,TimeScheme::BDF2,&ctx);
+        // dt_new / dt_old = 2: 5*phi2 - 9*phi1 + 4*phi0 = .6*RHS(phi2).
+        // This closed-form oracle does not call bdf2_coefficients().
+        const auto second=advance_time(first,.2,decay_rhs,TimeScheme::BDF2,&ctx);
+        for(std::size_t cell=0;cell<phi.size();++cell)
+            EXPECT_NEAR(second(cell),(9.0*first(cell)-4.0*phi(cell))/5.6,1e-12);
+        EXPECT_NEAR(ctx.dt_prev,.2,1e-15);
     });
 
     run_case("compute_cfl_time_step", [&]() {
@@ -133,9 +171,11 @@ int main() {
 
         double dt = compute_cfl_time_step(mesh, velocity, 0.5, 0.0);
 
-        // dt should be positive and finite (mesh dependent)
-        EXPECT_TRUE(dt > 0.0);
-        EXPECT_TRUE(std::isfinite(dt));
+        // The existing API uses h=V^(1/3), with V=.25*.1*.1.
+        const double expected=.5*std::cbrt(.0025)/10.0;
+        EXPECT_NEAR(dt,expected,1e-12);
+        EXPECT_NEAR(compute_cfl_time_step(mesh,velocity,.5,2.0),
+                    .5*std::cbrt(.0025)/12.0,1e-12);
     });
 
     run_case("adaptive_time_stepper", [&]() {
@@ -171,11 +211,10 @@ int main() {
 
         auto dt_local = compute_local_time_steps(mesh, velocity, 0.5);
 
-        // All dt should be positive and finite
-        for (std::size_t c = 0; c < 4; ++c) {
-            EXPECT_TRUE(dt_local[c] > 0.0);
-            EXPECT_TRUE(std::isfinite(dt_local[c]));
-        }
+        // Equal positive volumes make the expected local dt independent of
+        // the geometry routine under test. A 1e30 fallback must fail.
+        for(std::size_t c=0;c<4;++c)
+            EXPECT_NEAR(dt_local[c],.5*std::cbrt(.0025)/velocity(c,0),1e-12);
 
         // Apply local time stepping
         Field<double, Location::CELL> phi(4, "phi", "1", 1);
diff --git a/tests/python/test_su2_naca0012_tutorial.py b/tests/python/test_su2_naca0012_tutorial.py
new file mode 100644
index 00000000..24f7bf9c
--- /dev/null
+++ b/tests/python/test_su2_naca0012_tutorial.py
@@ -0,0 +1,28 @@
+"""The public NACA reader command cannot claim a CFDX flow validation."""
+import importlib.util
+from pathlib import Path
+
+ROOT = Path(__file__).resolve().parents[2]
+spec = importlib.util.spec_from_file_location("naca_reader_cli", ROOT / "validation/run_su2_naca0012.py")
+module = importlib.util.module_from_spec(spec)
+spec.loader.exec_module(module)
+
+
+def test_real_fixture_reports_reader_and_production_separately(tmp_path):
+    report = module.qualify_reader(tmp_path)
+    assert report["reader_passed"]
+    assert all(report["reader_checks"].values())
+    assert not report["canonical_conversion_passed"]
+    blocking = [item for item in report["canonical_conversion_findings"]["findings"]
+                if item["severity"] == "unsupported_blocking"]
+    assert any(item["feature"] == "canonical_selection" and "JST" in item["detail"]
+               for item in blocking)
+    assert report["boundary_initial_condition_equivalence_passed"] == all(report["physical_metadata_checks"].values())
+    assert not report["production_supported"]
+    assert report["physical_validation"] == "NOT_RUN"
+    assert not list(tmp_path.glob("*.cfdx.h5"))
+
+
+def test_default_command_cannot_exit_success_for_an_unrun_flow(tmp_path):
+    assert module.main(["--output-dir", str(tmp_path)]) == 2
+    assert module.main(["--output-dir", str(tmp_path), "--reader-only"]) == 0
diff --git a/tests/python/test_tmr_naca0012_import.py b/tests/python/test_tmr_naca0012_import.py
new file mode 100644
index 00000000..0d1f935c
--- /dev/null
+++ b/tests/python/test_tmr_naca0012_import.py
@@ -0,0 +1,118 @@
+"""Executable topology regressions for the NASA/TMR NACA C-grid importer."""
+from pathlib import Path
+import gzip
+import sys
+
+import h5py
+import numpy as np
+import pytest
+
+ROOT = Path(__file__).resolve().parents[2]
+sys.path.insert(0, str(ROOT / "scripts"))
+from import_tmr_naca0012 import convert, read_plot3d_2d, main
+from generate_2d_validation_meshes import write_mesh
+
+
+def miniature_c_grid():
+    # Two paired wake edges on each side, a four-edge airfoil, one radial
+    # layer. This is independent of NASA's ni/nj ratios and is deliberately
+    # small enough for topology counts to be an exact independent oracle.
+    x = np.array([[3, 2, 1, .5, 0, .5, 1, 2, 3],
+                  [3, 2, 1, 0, -1, 0, 1, 2, 3]], dtype=float)
+    y = np.array([[0, 0, 0, -.1, 0, .1, 0, 0, 0],
+                  [-1, -1, -1, -1, 0, 1, 1, 1, 1]], dtype=float)
+    return x, y
+
+
+def plot3d_text(x, y):
+    nj, ni = x.shape
+    return f"1\n{ni} {nj}\n" + "\n".join(str(v) for v in np.r_[x.ravel(), y.ravel()])
+
+
+@pytest.mark.parametrize("compressed", [False, True])
+def test_reads_official_single_block_header_and_gzip(tmp_path, compressed):
+    x, y = miniature_c_grid()
+    path = tmp_path / ("n0012_9-2.p2dfmt.gz" if compressed else "n0012_9-2.p2dfmt")
+    text = plot3d_text(x, y)
+    if compressed:
+        with gzip.open(path, "wt", encoding="ascii") as f:
+            f.write(text)
+    else:
+        path.write_text(text, encoding="ascii")
+    actual_x, actual_y = read_plot3d_2d(path)
+    np.testing.assert_array_equal(actual_x, x)
+    np.testing.assert_array_equal(actual_y, y)
+
+
+@pytest.mark.parametrize("mutation,diagnostic", [
+    (lambda s: s.replace("1\n", "2\n", 1), "single-block"),
+    (lambda s: s + "\n123", "numeric tokens"),
+    (lambda s: s.replace("3.0", "nan", 1), "finite"),
+])
+def test_rejects_unsupported_or_malformed_plot3d(tmp_path, mutation, diagnostic):
+    x, y = miniature_c_grid()
+    path = tmp_path / "n0012_9-2.p2dfmt"
+    path.write_text(mutation(plot3d_text(x, y)), encoding="ascii")
+    with pytest.raises(ValueError, match=diagnostic):
+        read_plot3d_2d(path)
+
+
+def test_wake_is_internal_and_extrusion_is_empty(tmp_path):
+    points, cells, boundaries = convert(*miniature_c_grid())
+    assert len(points) == 30  # 18 plane points minus 3 duplicate pairs, twice.
+    assert len(cells) == 8
+    assert {name: len(faces) for name, faces in boundaries.items()} == {
+        "airfoil": 4, "farfield": 10, "front": 8, "back": 8,
+    }
+    path = tmp_path / "miniature.h5"
+    write_mesh(path, points, cells, boundaries)
+    with h5py.File(path, "r") as f:
+        neighbours = f["neighbour"][:]
+        assert np.count_nonzero(neighbours >= 0) == 9  # 7 adjacent + 2 joined wake.
+        assert np.count_nonzero(neighbours == -1) == 30
+        patch_ids = f["patch_face_ids"][:]
+        assert len(np.unique(patch_ids)) == len(patch_ids)
+        assert set(patch_ids) == set(np.flatnonzero(neighbours == -1))
+        patch_types = {v.split(":")[0]: int(v.split(":")[-1])
+                       for v in f.attrs["boundary_patches"].split(";")}
+        assert patch_types == {"airfoil": 0, "farfield": 7, "front": 6, "back": 6}
+        points = f["points"][:]
+        face_vertices, face_offsets = f["face_vertices"][:], f["face_offsets"][:]
+        centres = np.array([points[face_vertices[face_offsets[i]:face_offsets[i+1]]].mean(axis=0)
+                            for i in range(len(neighbours))])
+        on_wake = (centres[:, 0] > 1) & (np.abs(centres[:, 1]) < 1e-14) & (centres[:, 2] == .5)
+        assert np.count_nonzero(on_wake) == 2
+        assert (neighbours[on_wake] >= 0).all()
+        # Signed polygon area need not share a winding convention; its
+        # nonzero magnitude is independently checked for every input cell.
+        for cell in cells:
+            q = points[np.asarray(cell[:4]), :2]
+            area = .5 * np.sum(q[:, 0] * np.roll(q[:, 1], -1) - q[:, 1] * np.roll(q[:, 0], -1))
+            assert abs(area) > 0
+
+
+def test_rejects_nonmatching_wake():
+    x, y = miniature_c_grid()
+    y[0, -1] = .01
+    with pytest.raises(ValueError, match="coincident"):
+        convert(x, y)
+
+
+def test_command_records_provenance_without_physical_promotion(tmp_path, monkeypatch):
+    source = tmp_path / "n0012_9-2.p2dfmt"
+    output = tmp_path / "grid.h5"
+    source.write_text(plot3d_text(*miniature_c_grid()), encoding="ascii")
+    monkeypatch.setattr(sys, "argv", ["import_tmr_naca0012", "--input", str(source), "--output", str(output)])
+    assert main() == 0
+    with h5py.File(output, "r") as f:
+        assert f.attrs["validation_wake_pairs"] == 2
+        assert f.attrs["validation_wall_faces"] == 4
+        assert f.attrs["validation_status"] == "imported-provenance-only"
+        assert len(f.attrs["validation_source_sha256"]) == 64
+
+
+def test_partial_wake_corruption_cannot_become_solid_wall():
+    x, y = miniature_c_grid()
+    y[0, -3] = .01  # break only the trailing-edge pair, after two good pairs.
+    with pytest.raises(ValueError, match="trailing edge"):
+        convert(x, y)
diff --git a/tests/unit/test_cg_solver.cpp b/tests/unit/test_cg_solver.cpp
index dcd2c90a..d34faf23 100644
--- a/tests/unit/test_cg_solver.cpp
+++ b/tests/unit/test_cg_solver.cpp
@@ -22,7 +22,106 @@ SparseMatrix make_poisson_1d(std::size_t n) {
     return A;
 }
 
+namespace {
+enum class FailurePoint { Setup, InitialApply, IterationApply, ReplacementDot };
+
+class FailingPreconditioner final : public Preconditioner {
+public:
+    explicit FailingPreconditioner(FailurePoint point) : point_(point) {}
+    bool setup(const SparseMatrix&) override { return point_ != FailurePoint::Setup; }
+    bool apply(const Vector& r, Vector& z) const override {
+        ++applications;
+        if ((point_ == FailurePoint::InitialApply && applications == 1) ||
+            (point_ == FailurePoint::IterationApply && applications == 2))
+            return false;
+        const double sign =
+            point_ == FailurePoint::ReplacementDot && applications >= 2 ? -1.0 : 1.0;
+        for (std::size_t i = 0; i < r.size(); ++i) z(i) = sign * r(i);
+        return true;
+    }
+    const char* name() const override { return "Failure injection"; }
+    mutable std::size_t applications = 0;
+private:
+    FailurePoint point_;
+};
+
+void check_preconditioner_failure_residual(FailurePoint point) {
+    SparseMatrix A(2, 2);
+    A.push_back(0, 0, 1.0);
+    A.push_back(1, 1, 2.0);
+    A.finalize();
+    Vector b(2, 1.0), x(2, 0.0);
+    const bool before_iteration =
+        point == FailurePoint::Setup || point == FailurePoint::InitialApply;
+    if (before_iteration) {
+        x(0) = 0.25;
+        x(1) = -0.5;
+    }
+    FailingPreconditioner pc(point);
+    const auto result = solve_cg(A, b, x, pc, 10, 1e-12);
+    EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE);
+    EXPECT_TRUE(result.iterations == (before_iteration ? 0u : 1u));
+    if (point == FailurePoint::ReplacementDot) {
+        // The third application follows residual replacement and rejects
+        // r^T M^-1 r < 0, as in the captured NACA pressure failure.
+        EXPECT_TRUE(pc.applications == 3);
+        EXPECT_TRUE(result.residual_replacements == 1);
+    }
+    if (!before_iteration) {
+        EXPECT_NEAR(x(0), 2.0 / 3.0, 1e-15);
+        EXPECT_NEAR(x(1), 2.0 / 3.0, 1e-15);
+    }
+    // Independent two-equation oracle, evaluated at the returned iterate.
+    const double expected = std::hypot(1.0 - x(0), 1.0 - 2.0 * x(1));
+    EXPECT_TRUE(expected > 0.0);
+    EXPECT_NEAR(result.residual, expected, 1e-14);
+    EXPECT_NEAR(result.residual_relative, expected / std::sqrt(2.0), 1e-14);
+}
+} // namespace
+
 int main() {
+    run_case("cg_failed_setup_reports_current_true_residual", [] {
+        check_preconditioner_failure_residual(FailurePoint::Setup);
+    });
+    run_case("cg_failed_initial_apply_reports_current_true_residual", [] {
+        check_preconditioner_failure_residual(FailurePoint::InitialApply);
+    });
+    run_case("cg_failed_iteration_apply_reports_updated_true_residual", [] {
+        check_preconditioner_failure_residual(FailurePoint::IterationApply);
+    });
+    run_case("cg_failed_replacement_reports_updated_true_residual", [] {
+        check_preconditioner_failure_residual(FailurePoint::ReplacementDot);
+    });
+    run_case("cg_failed_zero_rhs_reports_nonzero_true_residual", [] {
+        SparseMatrix A(1, 1);
+        A.push_back(0, 0, 2.0);
+        A.finalize();
+        Vector b(1, 0.0), x(1, 3.0);
+        FailingPreconditioner pc(FailurePoint::Setup);
+        const auto result = solve_cg(A, b, x, pc, 10, 1e-12);
+        EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE);
+        EXPECT_TRUE(result.iterations == 0);
+        EXPECT_NEAR(x(0), 3.0, 0.0);
+        EXPECT_NEAR(result.residual, 6.0, 1e-14);
+        EXPECT_TRUE(std::isinf(result.residual_relative));
+        EXPECT_TRUE(result.residual_relative > 0.0);
+    });
+    run_case("cg_nonpositive_curvature_reports_true_residual", [] {
+        SparseMatrix A(2, 2);
+        A.push_back(0, 0, 1.0); A.push_back(0, 1, 2.0);
+        A.push_back(1, 0, 2.0); A.push_back(1, 1, 1.0);
+        A.finalize();
+        Vector b(2), x(2, 0.0);
+        b(0) = 1.0; b(1) = -1.0;
+        const auto result = solve_cg(A, b, x, 10, 1e-12);
+        EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE);
+        EXPECT_TRUE(result.iterations == 1);
+        EXPECT_NEAR(x(0), 0.0, 0.0);
+        EXPECT_NEAR(x(1), 0.0, 0.0);
+        EXPECT_NEAR(result.residual, std::sqrt(2.0), 1e-14);
+        EXPECT_NEAR(result.residual_relative, 1.0, 1e-14);
+    });
+
     run_case("cg_solve_diagonal", []() {
         // A = 2I, b = [2, 2, 2], solution attendue x = [1, 1, 1].
         SparseMatrix A(3, 3);
@@ -177,4 +276,4 @@ int main() {
     });
 
     return run_all();
-}
\ No newline at end of file
+}
diff --git a/tests/unit/test_core_temporal.cpp b/tests/unit/test_core_temporal.cpp
index 928e9809..40a007a2 100644
--- a/tests/unit/test_core_temporal.cpp
+++ b/tests/unit/test_core_temporal.cpp
@@ -75,5 +75,59 @@ int main() {
         EXPECT_TRUE(second(0) > 0.0);
     });
 
+    for (const auto scheme : {TimeScheme::EULER_IMPLICIT,
+                              TimeScheme::CRANK_NICOLSON, TimeScheme::BDF2}) {
+        run_case(std::string("intermediate_nonfinite_rhs_rejected_") + to_string(scheme), [scheme] {
+            ScalarCellField phi(1, "phi", "1", 1);
+            phi(0) = 1.0;
+            TimeIntegrationContext ctx(1, 1, "phi");
+            ctx.initialize(phi);
+            const RhsFunction rhs = [](const ScalarCellField& in, ScalarCellField& out) {
+                out(0) = -std::sqrt(in(0));
+            };
+            const double dt = scheme == TimeScheme::CRANK_NICOLSON ? 4.0 : 2.0;
+            EXPECT_THROW(advance_time(phi, dt, rhs, scheme, &ctx), std::runtime_error);
+            EXPECT_NEAR(ctx.phi_curr(0), 1.0, 0.0);
+            EXPECT_FALSE(ctx.has_prev);
+        });
+    }
+
+    run_case("nonfinite_initial_field_rejected_with_finite_rhs", [] {
+        ScalarCellField phi(1, "phi", "1", 1);
+        phi(0) = std::numeric_limits<double>::quiet_NaN();
+        const RhsFunction rhs = [](const ScalarCellField&, ScalarCellField& out) { out(0) = 0.0; };
+        EXPECT_THROW(advance_time(phi, 0.1, rhs, TimeScheme::EULER_EXPLICIT), std::runtime_error);
+    });
+
+    run_case("explicit_overflow_rejected_before_history_commit", [] {
+        ScalarCellField phi(1, "phi", "1", 1);
+        phi(0) = std::numeric_limits<double>::max();
+        TimeIntegrationContext ctx(1, 1, "phi");
+        ctx.initialize(phi);
+        const RhsFunction rhs = [](const ScalarCellField&, ScalarCellField& out) {
+            out(0) = std::numeric_limits<double>::max();
+        };
+        EXPECT_THROW(advance_time(phi, 1.0, rhs, TimeScheme::EULER_EXPLICIT, &ctx), std::runtime_error);
+        EXPECT_NEAR(ctx.phi_curr(0), phi(0), 0.0);
+        EXPECT_FALSE(ctx.has_prev);
+    });
+
+    run_case("bdf2_history_path_rejects_nonfinite_rhs_without_committing", [] {
+        ScalarCellField phi(1, "phi", "1", 1);
+        phi(0) = 1.0;
+        TimeIntegrationContext ctx(1, 1, "phi");
+        ctx.initialize(phi);
+        ctx.shift(phi);
+        ctx.dt_prev = 2.0;
+        const RhsFunction rhs = [](const ScalarCellField& in, ScalarCellField& out) {
+            out(0) = -std::sqrt(in(0));
+        };
+        EXPECT_THROW(advance_time(phi, 2.0, rhs, TimeScheme::BDF2, &ctx), std::runtime_error);
+        EXPECT_NEAR(ctx.phi_curr(0), 1.0, 0.0);
+        EXPECT_NEAR(ctx.phi_prev(0), 1.0, 0.0);
+        EXPECT_NEAR(ctx.dt_prev, 2.0, 0.0);
+        EXPECT_TRUE(ctx.has_prev);
+    });
+
     return cfdx::testing::run_all();
 }
diff --git a/tests/unit/test_matrix_characteristics.cpp b/tests/unit/test_matrix_characteristics.cpp
index d5f64c25..f9847af4 100644
--- a/tests/unit/test_matrix_characteristics.cpp
+++ b/tests/unit/test_matrix_characteristics.cpp
@@ -50,5 +50,5 @@ int main() {
         EXPECT_TRUE(characteristics.coefficient_range > 1.0e10);
     });
 
-    return 0;
+    return run_all();
 }
diff --git a/tests/unit/test_mesh.cpp b/tests/unit/test_mesh.cpp
index cc0e8597..2c1c2104 100644
--- a/tests/unit/test_mesh.cpp
+++ b/tests/unit/test_mesh.cpp
@@ -117,6 +117,37 @@ int main() {
         EXPECT_FALSE(m.boundary().is_consistent(m.n_faces()));
     });
 
+    run_case("boundary_mesh_membership_complete", []() {
+        const Mesh m = make_unit_cube();
+        EXPECT_TRUE(m.boundary().is_consistent_with_mesh(m));
+    });
+
+    run_case("boundary_mesh_membership_orphan", []() {
+        Mesh m = make_unit_cube();
+        m.boundary().remove_patch(0);
+        EXPECT_FALSE(m.boundary().is_consistent_with_mesh(m));
+    });
+
+    run_case("boundary_mesh_membership_duplicate", []() {
+        Mesh m = make_unit_cube();
+        m.boundary().patch(0).face_ids.push_back(0);
+        EXPECT_FALSE(m.boundary().is_consistent_with_mesh(m));
+    });
+
+    run_case("boundary_mesh_membership_internal_face", []() {
+        Mesh m = make_unit_cube();
+        // Only the boundary-membership contract is under test here: a
+        // nonnegative neighbour makes this a non-boundary face.
+        m.ownership().set_neighbour(0, 1);
+        EXPECT_FALSE(m.boundary().is_consistent_with_mesh(m));
+    });
+
+    run_case("boundary_mesh_membership_out_of_range", []() {
+        Mesh m = make_unit_cube();
+        m.boundary().patch(0).face_ids[0] = static_cast<FaceIndex>(m.n_faces());
+        EXPECT_FALSE(m.boundary().is_consistent_with_mesh(m));
+    });
+
     run_case("validate_mesh_full_ok", []() {
         Mesh m = make_unit_cube();
         auto report = validate_mesh(m);
diff --git a/tests/unit/test_wall_distance_advection_residual.cpp b/tests/unit/test_wall_distance_advection_residual.cpp
new file mode 100644
index 00000000..7a796605
--- /dev/null
+++ b/tests/unit/test_wall_distance_advection_residual.cpp
@@ -0,0 +1,55 @@
+#include "cfdx/physics/wall_distance.h"
+#include <cmath>
+#include <iostream>
+using namespace cfdx::physics;
+int main() {
+    int failed=0;
+    const auto expect=[&](bool ok,const char* message){if(!ok){++failed;std::cerr<<message<<'\n';}};
+    // Manufactured stationary ridge for the actual discrete transport rule.
+    // Opposing derivatives cancel in U. With h=1/4, gamma=1/20 and
+    // d_E=d_W=1/4, gamma*(2*d_P-d_E-d_W)/h^2=1 gives d_P=7/8.
+    auto g=make_wall_distance_grid(5,3,3,{0,-1,-1},{.25,1,1},
+        [](const WallDistanceVec3& p){return p.x<=0 || p.x>=1;});
+    std::vector<double> d(g.points.size(),0.0);
+    std::vector<unsigned char> fixed(g.points.size(),1);
+    for(std::size_t id=0;id<d.size();++id)
+        if(!g.solid[id]) d[id]=g.points[id].x==.5 ? .875 : .25;
+    const auto centre=g.index(2,1,1);fixed[centre]=0;
+    const double stationary=wall_distance_pde_residual_inf(
+        WallDistanceMethod::ADVECTION_DIFFUSION,d,g,&fixed,.05);
+    expect(stationary<1e-13,"stationary transport row was measured with another operator");
+    // Zero field must retain the unit source even when U is undefined/zero.
+    std::fill(d.begin(),d.end(),0.0);
+    expect(std::abs(wall_distance_pde_residual_inf(
+        WallDistanceMethod::ADVECTION_DIFFUSION,d,g,&fixed,.05)-1.0)<1e-13,
+        "zero gradient erased the unit source");
+
+    // Fractional wall at x=.1: at P=.25, d=(x-.1)^2 gives d_P=.0225,
+    // d_E=.16, delta_W=.15. U=1, transport=.15 and the assembled diffusion
+    // Laplacian is -.0225/(.25*.15)+(.16-.0225)/.25^2=1.6.
+    // Hence |.15-1-.05*1.6|=.93, versus .95 with Shortley-Weller.
+    WallSurface wall;wall.points={{.1,-4,-4},{.1,4,-4},{.1,4,4},{.1,-4,4}};
+    wall.triangles={{{0,1,2}},{{0,2,3}}};WallDistanceBvh bvh(wall);
+    auto cut=make_wall_distance_grid(5,3,3,{0,-1,-1},{.25,1,1},
+        [](const WallDistanceVec3& p){return p.x<.1;});
+    std::vector<double> quadratic(cut.points.size(),0.0);
+    std::vector<unsigned char> mask(cut.points.size(),1);
+    for(std::size_t id=0;id<quadratic.size();++id)
+        if(!cut.solid[id])quadratic[id]=std::pow(cut.points[id].x-.1,2);
+    mask[cut.index(1,1,1)]=0;
+    const double fractional=wall_distance_pde_residual_inf(
+        WallDistanceMethod::ADVECTION_DIFFUSION,quadratic,cut,&mask,.05,&bvh);
+    expect(std::abs(fractional-.93)<1e-13,"fractional-wall residual changed the assembled diffusion stencil");
+    // A stationary fractional-cut row with nonzero U also distinguishes the
+    // operators without relying on the old zero-gradient skip: d_P=.2,
+    // d_E=.95 gives transport=4/3 and the assembled Laplacian=20/3.
+    // Thus transport-gamma*lap=1, while the S-W residual is 1/12.
+    for(std::size_t id=0;id<quadratic.size();++id)
+        if(!cut.solid[id])quadratic[id]=.2+3.0*(cut.points[id].x-.25);
+    const double stationary_cut=wall_distance_pde_residual_inf(
+        WallDistanceMethod::ADVECTION_DIFFUSION,quadratic,cut,&mask,.05,&bvh);
+    expect(stationary_cut<1e-13,"stationary cut row was certified with a different diffusion operator");
+    std::cout<<"stationary_residual="<<stationary<<" fractional_residual="<<fractional
+             <<" stationary_cut_residual="<<stationary_cut<<" failed="<<failed<<'\n';
+    return failed==0?0:1;
+}
diff --git a/tests/unit/test_wall_distance_godunov_ridges.cpp b/tests/unit/test_wall_distance_godunov_ridges.cpp
new file mode 100644
index 00000000..9f4a71cc
--- /dev/null
+++ b/tests/unit/test_wall_distance_godunov_ridges.cpp
@@ -0,0 +1,38 @@
+#include "cfdx/physics/wall_distance.h"
+#include <cmath>
+#include <iostream>
+#include <stdexcept>
+using namespace cfdx::physics;
+int main() {
+    WallSurface wall;
+    for(double x:{0.0,2.0}) {
+        const std::size_t b=wall.points.size();
+        wall.points.insert(wall.points.end(),{{x,-4,-4},{x,4,-4},{x,4,4},{x,-4,4}});
+        wall.triangles.push_back({{b,b+1,b+2}});
+        wall.triangles.push_back({{b,b+2,b+3}});
+    }
+    const auto grid=make_wall_distance_grid(9,3,3,{0,-1,-1},{.25,1,1},
+        [](const WallDistanceVec3& p){return p.x<=0.0 || p.x>=2.0;});
+    const WallDistanceBvh bvh(wall);
+    const WallDistanceCutFaceOffsetCache cache(bvh,grid);
+    const auto expected=exact_reference(wall,grid);
+    const auto result=compute_wall_distance(WallDistanceMethod::EIKONAL,wall,grid,50);
+    const std::size_t ridge=grid.index(4,1,1);
+    const auto require=[](bool ok,const char* message){if(!ok)throw std::runtime_error(message);};
+    require(std::abs(result.distance[ridge]-1.0)<1e-13,"channel ridge distance is not exact");
+    for(const std::size_t i:{2u,4u,6u}) {
+        const std::size_t id=grid.index(i,1,1);
+        require(std::abs(godunov_gradient_at(expected,grid,id)-1.0)<1e-13,
+                "Godunov gradient double-counted opposite directions");
+        require(std::abs(godunov_gradient_at(bvh,expected,grid,id)-1.0)<1e-13,
+                "BVH Godunov gradient double-counted opposite directions");
+        require(std::abs(godunov_gradient_at(bvh,expected,grid,id,&cache)-1.0)<1e-13,
+                "cached Godunov gradient is inconsistent");
+    }
+    const auto metrics=compare_wall_distance(grid,expected,result.distance,.5);
+    require(metrics.l2_relative<1e-13,"channel distance differs from independent two-wall oracle");
+    require(result.converged && result.residual_inf<1e-13,
+            "exact channel distance is incorrectly reported as unconverged");
+    std::cout << "Godunov ridge and monotone branches: PASS, L2=" << metrics.l2_relative
+              << " residual=" << result.residual_inf << '\n';
+}
diff --git a/tests/unit/test_wall_distance_invalid_results.cpp b/tests/unit/test_wall_distance_invalid_results.cpp
new file mode 100644
index 00000000..3dc07ac4
--- /dev/null
+++ b/tests/unit/test_wall_distance_invalid_results.cpp
@@ -0,0 +1,76 @@
+#include "cfdx/physics/wall_distance.h"
+
+#include <cmath>
+#include <iostream>
+#include <limits>
+#include <vector>
+
+using namespace cfdx::physics;
+
+int main() {
+    int failures=0, checks=0;
+    const auto expect=[&](bool ok,const char* label) {
+        ++checks;
+        if(!ok) { ++failures; std::cerr << "FAIL: " << label << '\n'; }
+    };
+    const auto grid=make_wall_distance_grid(4,4,4,{10,0,0},{1,1,1},
+        [](const WallDistanceVec3&){ return false; });
+    std::vector<double> reference(grid.points.size(),1.0), distance=reference;
+    auto metrics=compare_wall_distance(grid,reference,distance,2.0);
+    expect(metrics.samples==grid.points.size() && metrics.l2_relative==0.0,
+           "valid identical fields retain zero error");
+    distance[0]=std::numeric_limits<double>::quiet_NaN();
+    metrics=compare_wall_distance(grid,reference,distance,2.0);
+    expect(std::isinf(metrics.l2_relative) && std::isinf(metrics.linf_relative),
+           "one NaN cannot be excluded from the accuracy verdict");
+    distance.assign(grid.points.size(),std::numeric_limits<double>::infinity());
+    metrics=compare_wall_distance(grid,reference,distance,2.0);
+    expect(metrics.samples==0 && std::isinf(metrics.l2_relative) &&
+           std::isinf(metrics.near_wall_l2_relative),
+           "a wholly invalid field cannot have zero global or near-wall error");
+    distance=reference;
+    reference[1]=std::numeric_limits<double>::quiet_NaN();
+    metrics=compare_wall_distance(grid,reference,distance,2.0);
+    expect(std::isinf(metrics.l2_relative),"a corrupt oracle cannot certify accuracy");
+    reference=distance;
+    distance[2]=-1.0;
+    metrics=compare_wall_distance(grid,reference,distance,2.0);
+    expect(std::isinf(metrics.l2_relative),"negative distance is an invalid field");
+    distance=reference;
+    metrics=compare_wall_distance(grid,reference,distance,0.5);
+    expect(std::isinf(metrics.near_wall_l2_relative),
+           "an empty near-wall sample band cannot report zero error");
+
+    WallSurface wall;
+    wall.points={{0,-4,-4},{0,4,-4},{0,4,4},{0,-4,4}};
+    wall.triangles={{{0,1,2}},{{0,2,3}}};
+    // The wall is outside the grid and its seed band. Geometric distance is
+    // well defined, but propagation methods cannot reach any of these nodes.
+    for(const auto method:{WallDistanceMethod::MESH_WAVE,
+                           WallDistanceMethod::DIRECTIONAL_MESH_WAVE,
+                           WallDistanceMethod::EIKONAL,
+                           WallDistanceMethod::HAMILTON_JACOBI,
+                           WallDistanceMethod::ADVECTION_DIFFUSION}) {
+        const auto result=compute_wall_distance(method,wall,grid,20);
+        expect(!result.converged && std::isinf(result.residual_inf),
+               "a propagation result with no reachable sample cannot converge");
+    }
+    auto solid=grid;
+    solid.solid.assign(grid.points.size(),1);
+    const auto empty_result=compute_wall_distance(WallDistanceMethod::EXACT_GEOMETRIC,wall,solid,20);
+    expect(!empty_result.converged,"an all-solid domain is not an accuracy certificate");
+    metrics=compare_wall_distance(solid,reference,distance,2.0);
+    expect(metrics.samples==0 && std::isinf(metrics.l2_relative),
+           "an all-solid comparison has no accuracy evidence");
+
+    distance.assign(grid.points.size(),1.0);
+    distance[0]=std::numeric_limits<double>::quiet_NaN();
+    std::vector<unsigned char> fixed(grid.points.size(),1);
+    expect(std::isinf(wall_distance_pde_residual_inf(WallDistanceMethod::EIKONAL,distance,grid,&fixed)),
+           "invalid fixed seeds invalidate the PDE residual");
+    distance.assign(grid.points.size(),0.0);
+    expect(wall_distance_pde_residual_inf(WallDistanceMethod::ADVECTION_DIFFUSION,distance,grid)>0.99,
+           "zero gradient does not erase the unit source residual");
+    std::cout << "checks=" << checks << " failures=" << failures << '\n';
+    return failures==0 ? 0 : 1;
+}
diff --git a/tests/validation/test_n11_scalar_conservation.cpp b/tests/validation/test_n11_scalar_conservation.cpp
index 2a63e624..373e1551 100644
--- a/tests/validation/test_n11_scalar_conservation.cpp
+++ b/tests/validation/test_n11_scalar_conservation.cpp
@@ -129,4 +129,5 @@ int main()
         EXPECT_TRUE(bounds.above == 0);
         EXPECT_TRUE(bounds.nonfinite == 1);
     });
+    return run_all();
 }
diff --git a/tests/validation/test_n9_s10_pressure_reference_policies.cpp b/tests/validation/test_n9_s10_pressure_reference_policies.cpp
index 4bc997c3..e8ed58a6 100644
--- a/tests/validation/test_n9_s10_pressure_reference_policies.cpp
+++ b/tests/validation/test_n9_s10_pressure_reference_policies.cpp
@@ -83,13 +83,30 @@ int main()
         EXPECT_NEAR(p(1), 3.0, 1e-12);
     });
 
-    run_case("n9_s10_zero_mean_policy", [] {
+    run_case("n9_s10_zero_mean_shift_preserves_pressure_differences", [] {
         const Mesh mesh = make_two_cell_channel();
-        Field<double, Location::CELL> U(2, "U", "m/s", 3);
+        const auto geometry = build_fv_geometry(mesh);
         Field<double, Location::CELL> p(2, "p", "Pa", 1);
-        U.fill(0.0);
         p(0) = 2.0;
         p(1) = 6.0;
+
+        // Gauge removal changes only the additive constant. This is an
+        // operator contract, independent of the momentum/continuity solve.
+        apply_pressure_gauge_policy(
+            p, geometry, PressureGaugePolicy::ZERO_MEAN, 0, 3.0, false);
+        EXPECT_NEAR(0.5 * (p(0) + p(1)), 0.0, 1e-12);
+        EXPECT_NEAR(p(1) - p(0), 4.0, 1e-12);
+    });
+
+    run_case("n9_s10_zero_mean_equilibrium_solver", [] {
+        const Mesh mesh = make_two_cell_channel();
+        Field<double, Location::CELL> U(2, "U", "m/s", 3);
+        Field<double, Location::CELL> p(2, "p", "Pa", 1);
+        U.fill(0.0);
+        // A closed, unforced box at rest has constant pressure. A nonzero
+        // initial pressure gradient is not an equilibrium that a steady solve
+        // should preserve; that property belongs only to the gauge helper.
+        p.fill(7.0);
         VelocityBoundaryConditions ubc;
         ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
         ScalarBoundaryConditions pbc;
@@ -99,7 +116,10 @@ int main()
         const auto result = solve_steady_incompressible(mesh, U, p, ubc, pbc, c);
         EXPECT_TRUE(result.converged);
         EXPECT_NEAR(0.5 * (p(0) + p(1)), 0.0, 1e-12);
-        EXPECT_NEAR(p(1) - p(0), 4.0, 1e-12);
+        EXPECT_NEAR(p(1) - p(0), 0.0, 1e-12);
+        for (std::size_t cell = 0; cell < U.size(); ++cell)
+            for (std::size_t component = 0; component < 3; ++component)
+                EXPECT_NEAR(U(cell, component), 0.0, 1e-12);
     });
 
     run_case("n9_s10_none_policy_preserves_pressure_level", [] {
@@ -120,5 +140,5 @@ int main()
         EXPECT_NEAR(p(0), 7.0, 1e-12);
         EXPECT_NEAR(p(1), 7.0, 1e-12);
     });
-    return 0;
+    return run_all();
 }
diff --git a/tests/validation/test_n9_s9_rhie_chow_consistency.cpp b/tests/validation/test_n9_s9_rhie_chow_consistency.cpp
index 5e3150f9..355247d3 100644
--- a/tests/validation/test_n9_s9_rhie_chow_consistency.cpp
+++ b/tests/validation/test_n9_s9_rhie_chow_consistency.cpp
@@ -788,5 +788,5 @@ int main()
         EXPECT_TRUE(identity_stale > 1e-12);
     });
 
-    return 0;
+    return run_all();
 }
diff --git a/validation/run_su2_naca0012.py b/validation/run_su2_naca0012.py
index 2880ef90..e12e2154 100755
--- a/validation/run_su2_naca0012.py
+++ b/validation/run_su2_naca0012.py
@@ -1,94 +1,115 @@
 #!/usr/bin/env python3
-"""SU2 NACA0012 validation: convert -> solve -> compare with theory.
+"""Qualify reading the SU2 NACA0012 Euler fixture, with explicit scope.
 
-Usage: python3 validation/run_su2_naca0012.py
+This fixture is not a CFDX aerodynamic validation.  Its compressible Euler
+physics and JST discretization have no production path in the current
+incompressible executable.  The faithful converter must keep reporting
+unsupported numerical selections instead of silently replacing them.
 """
+from __future__ import annotations
+
+import argparse
+import json
+import math
 import shutil
-import subprocess
 import sys
 import tempfile
 from pathlib import Path
 
 ROOT = Path(__file__).resolve().parent.parent
 DATA_DIR = ROOT / "tests" / "data" / "su2"
-BUILD_DIR = ROOT / "build"
-SOLVER = BUILD_DIR / "cfdx_production_solver"
-
 sys.path.insert(0, str(ROOT / "src" / "cfdx" / "python"))
 
+from cfdx.io.adapters.su2 import Su2Adapter  # noqa: E402
 from cfdx.io.converter import convert  # noqa: E402
+from cfdx.io.gap_analysis import GapAnalysisReport  # noqa: E402
+from cfdx.io.interfaces import ConversionResult, SourceInfo  # noqa: E402
 
 
-def main() -> int:
-    print("=== SU2 NACA0012 Validation ===")
-    print("Case: 2-D NACA0012, Euler, Mach 0.8, AoA 1.25 deg")
-    print("Mesh: 5233 vertices, 10216 triangles")
-    print()
+def qualify_reader(directory: Path) -> dict:
+    """Retain separate adapter, canonical-conversion and physics verdicts."""
+    directory.mkdir(parents=True, exist_ok=True)
+    shutil.copy2(DATA_DIR / "mesh_NACA0012_inv.su2", directory / "mesh.su2")
+    shutil.copy2(DATA_DIR / "inv_NACA0012_basic.cfg", directory / "config.cfg")
+    result = ConversionResult(source=SourceInfo())
+    read_ok = Su2Adapter().convert(directory, result)
+    setup = result.setup
+    checks = {
+        "reader_completed": bool(read_ok),
+        "no_blocking_reader_finding": not result.gap_report.has_blocking(),
+        "points": result.mesh is not None and result.mesh["points"].shape == (5233, 3),
+        "triangles": result.mesh is not None and len(result.mesh["cells"].get("triangle", [])) == 10216,
+        "physics": setup is not None and setup.physics_model == "compressible_euler",
+        "energy": setup is not None and setup.energy_model == "compressible",
+        "boundaries": setup is not None and {bc.patch_name for bc in setup.boundary_conditions} == {"airfoil", "farfield"},
+        "pressure": setup is not None and setup.initial_condition.pressure == 101325.0,
+        "temperature": setup is not None and setup.initial_condition.temperature == 288.15,
+    }
+    # Keep the stricter converter result visible.  A reader PASS must not
+    # be interpreted as acceptance of unsupported JST or compressible flow.
+    canonical = convert(directory)
+    boundary_by_name = {bc.patch_name: bc for bc in setup.boundary_conditions} if setup else {}
+    airfoil = boundary_by_name.get("airfoil")
+    farfield = boundary_by_name.get("farfield")
+    velocity = setup.initial_condition if setup else None
+    physical_metadata_checks = {
+        "euler_wall_slip": airfoil is not None and airfoil.value_type.value == "wall_slip",
+        "initial_velocity_vector_matches_speed": velocity is not None and math.isclose(
+            math.sqrt(sum(value * value for value in velocity.velocity_vector)),
+            velocity.velocity, rel_tol=1e-12, abs_tol=1e-12),
+        "farfield_velocity_matches_initial_speed": velocity is not None and farfield is not None and math.isclose(
+            math.sqrt(sum(value * value for value in farfield.velocity_vector)) if farfield.velocity_vector
+            else farfield.velocity_magnitude,
+            velocity.velocity, rel_tol=1e-12, abs_tol=1e-12),
+    }
+    return {
+        "scope": "SU2 NACA0012 input-reader qualification",
+        "reader_passed": all(checks.values()),
+        "reader_checks": checks,
+        "physical_metadata_checks": physical_metadata_checks,
+        "boundary_initial_condition_equivalence_passed": all(physical_metadata_checks.values()),
+        "reader_findings": json.loads(GapAnalysisReport(result.gap_report).to_json()),
+        "canonical_conversion_passed": canonical.success,
+        "canonical_conversion_findings": json.loads(GapAnalysisReport(canonical.gap_analysis).to_json()),
+        "production_supported": False,
+        "production_reason": "The CFDX production executable supports an incompressible carrier; this fixture declares compressible Euler at Mach 0.8.",
+        "physical_validation": "NOT_RUN",
+    }
 
-    tmp = Path(tempfile.mkdtemp(prefix="cfdx_naca0012_"))
-    print(f"Working dir: {tmp}")
 
-    # 1. Convert SU2 -> CFDX
-    print("[1/4] Converting SU2 -> CFDX...")
-    shutil.copy2(DATA_DIR / "mesh_NACA0012_inv.su2", tmp / "mesh.su2")
-    shutil.copy2(DATA_DIR / "inv_NACA0012_basic.cfg", tmp / "config.cfg")
-
-    result = convert(tmp, output=tmp / "case.cfdx.h5")
-    if not result.success:
-        print("Conversion FAILED:")
-        for f in result.gap_analysis.findings:
-            if "block" in str(f.severity).lower():
-                print(f"  BLOCKING: {f.feature} - {f.detail}")
+def main(argv=None) -> int:
+    parser = argparse.ArgumentParser(description=__doc__)
+    parser.add_argument("--output-dir", type=Path,
+                        help="directory for input copies and the reader report")
+    parser.add_argument("--reader-only", action="store_true",
+                        help="request only format-reader checks; permit success without a flow solve")
+    args = parser.parse_args(argv)
+    directory = args.output_dir or Path(tempfile.mkdtemp(prefix="cfdx_naca0012_reader_"))
+    report = qualify_reader(directory)
+    report_path = directory / "reader_report.json"
+    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
+    print("=== SU2 NACA0012 input-reader qualification ===")
+    print("Reader:", "PASS" if report["reader_passed"] else "FAIL")
+    print("Boundary/initial-condition equivalence:",
+          "PASS" if report["boundary_initial_condition_equivalence_passed"] else "FAIL")
+    for check, passed in report["physical_metadata_checks"].items():
+        if not passed:
+            print(f"  MISMATCH: {check}")
+    print("Canonical conversion:", "PASS" if report["canonical_conversion_passed"] else "BLOCKED")
+    for finding in report["canonical_conversion_findings"]["findings"]:
+        if finding["severity"] == "unsupported_blocking":
+            print(f"  BLOCKING: {finding['feature']} - {finding['detail']}")
+    print("Production: UNSUPPORTED for this compressible Euler fixture")
+    print("Physical validation: NOT_RUN; no Cd, Cl or Cp prediction is asserted")
+    print(f"Reader evidence: {report_path}")
+    if not report["reader_passed"]:
         return 1
-    print("  OK Conversion complete")
-    print(f"  Mesh: {result.mesh['points'].shape[0]} points, "
-          f"{sum(len(v) for v in result.mesh['cells'].values())} cells")
-
-    case = result.case
-    print(f"  BCs: {[b.patch_name for b in case.boundary_conditions]}")
-    print(f"  Physics: {case.physics_model}")
-    print(f"  Initial pressure: {case.initial_condition.pressure} Pa")
-    print(f"  Initial temperature: {case.initial_condition.temperature} K")
-
-    # 2. Theory check
-    print("[2/4] Checking theory...")
-    p0 = 101325.0
-    M = 0.8
-    gamma = 1.4
-    p_theory = p0 * (1 + 0.2 * M**2)**(-3.5)
-    print(f"  Freestream pressure: {p0} Pa")
-    print(f"  Theory (isentropic): {p_theory:.2f} Pa")
-    print(f"  Initial condition matches freestream: "
-          f"{abs(case.initial_condition.pressure - p0) < 0.1}")
-
-    # 3. Mesh validation
-    print("[3/4] Validating mesh...")
-    points = result.mesh["points"]
-    cells = result.mesh["cells"]
-    print(f"  Points: {points.shape[0]} x {points.shape[1]}")
-    print(f"  Cells: {sum(len(v) for v in cells.values())} total")
-    for cell_type, conn in cells.items():
-        print(f"    {cell_type}: {len(conn)}")
-    print(f"  Boundary patches: {len(case.boundary_conditions)}")
-
-    # 4. Summary
-    print()
-    print("=== Validation Summary ===")
-    print("  Case:        SU2 NACA0012 Euler")
-    print("  Converter:   cfdx.io.converter (Python)")
-    print("  Mesh:        5233 vertices, 10216 triangles")
-    print("  BCs:         airfoil, farfield")
-    print("  Physics:     compressible_euler")
-    print("  Output:      case.cfdx.h5")
-    print()
-    print("  NOTE: SU2 solver not installed. To run SU2:")
-    print("    cd", tmp, "&& su2 config.cfg")
-    print()
-    print("  To run CFDX solver (after fixing meshio import):")
-    print(f"    {SOLVER} --mesh {tmp}/case.cfdx.h5 --output-dir {tmp} --iterations 1000")
-    print()
-    print(f"TMP_DIR={tmp}")
-    return 0
+    if args.reader_only:
+        return 0
+    # Preserve a failing default verdict for legacy callers that expected
+    # this command to validate the full aerodynamic run.
+    print("Full aerodynamic validation is unavailable. Use --reader-only for the explicitly limited reader scope.")
+    return 2
 
 
 if __name__ == "__main__":
diff --git a/validation/run_su2_naca0012.sh b/validation/run_su2_naca0012.sh
index 4a0bd336..5de82d70 100755
--- a/validation/run_su2_naca0012.sh
+++ b/validation/run_su2_naca0012.sh
@@ -1,58 +1,5 @@
 #!/bin/bash
-# SU2 NACA0012 validation: convert -> solve -> compare with theory
+# Input-reader qualification only; the Python entry point records its scope.
 set -euo pipefail
-
 ROOT="$(cd "$(dirname "$0")/.." && pwd)"
-DATA_DIR="$ROOT/tests/data/su2"
-TMP_DIR="$(mktemp -d /tmp/cfdx_naca0012_XXXXXX)"
-CLI="$ROOT/build/cfdx_convert"
-SOLVER="$ROOT/build/cfdx_production_solver"
-
-echo "=== SU2 NACA0012 Validation ==="
-echo "Case: 2-D NACA0012, Euler, Mach 0.8, AoA 1.25 deg"
-echo "Mesh: 5233 vertices, 10216 triangles"
-echo ""
-
-# 1. Convert SU2 -> CFDX
-echo "[1/4] Converting SU2 -> CFDX..."
-cp "$DATA_DIR/mesh_NACA0012_inv.su2" "$TMP_DIR/mesh.su2"
-cp "$DATA_DIR/inv_NACA0012_basic.cfg" "$TMP_DIR/config.cfg"
-
-"$CLI" convert "$TMP_DIR" "$TMP_DIR/case.cfdx.h5"
-echo "  OK Conversion complete"
-
-# 2. Run CFDX solver
-echo "[2/4] Running CFDX solver..."
-cd "$TMP_DIR"
-"$SOLVER" --case case.cfdx.h5 --output solution.cfdx.h5 --max-iter 1000 2>&1 | tail -5
-echo "  OK Solver complete"
-
-# 3. Extract and compare results
-echo "[3/4] Extracting results..."
-python3 -c "
-import h5py, numpy as np
-with h5py.File('$TMP_DIR/solution.cfdx.h5', 'r') as f:
-    p = f['fields/scalar/pressure'][:]
-    U = f['fields/vector/velocity'][:]
-    print(f'Pressure range: {p.min():.2f} - {p.max():.2f} Pa')
-    print(f'Velocity range: {np.linalg.norm(U, axis=1).min():.1f} - {np.linalg.norm(U, axis=1).max():.1f} m/s')
-    p0 = 101325.0; M = 0.8; g = 1.4
-    p_theory = p0 * (1 + 0.2 * M**2)**(-3.5)
-    print(f'Freestream pressure: {p0} Pa')
-    print(f'Theory (isentropic): {p_theory:.2f} Pa')
-    print(f'Computed min pressure: {p.min():.2f} Pa')
-    print(f'Ratio computed/theory: {p.min()/p_theory:.4f}')
-"
-
-# 4. Summary
-echo ""
-echo "=== Validation Summary ==="
-echo "  Case:        SU2 NACA0012 Euler"
-echo "  Converter:   cfdx_convert (Python bridge)"
-echo "  Solver:      cfdx_production_solver"
-echo "  Output:      $TMP_DIR/solution.cfdx.h5"
-echo ""
-echo "To compare with SU2 results, run:"
-echo "  cd $TMP_DIR && su2 config.cfg"
-echo ""
-echo "TMP_DIR=$TMP_DIR"
+exec "${PYTHON:-python3}" "$ROOT/validation/run_su2_naca0012.py" "$@"

```

## Partie II — Sources complètes de tous les fichiers modifiés

### `CMakeLists.txt`

SHA256 `b3b1cb8fe9dcfb99ac82e6de2fd313aa4c1a67bb377cee65140a8456471e7ebf`

```cmake
cmake_minimum_required(VERSION 3.16)

project(CFDX
    VERSION 0.7
    LANGUAGES CXX C
    DESCRIPTION "CFDX - High-Performance General-Purpose CFD Framework"
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---------------------------------------------------------------------------
# Build options
# ---------------------------------------------------------------------------
option(CFDX_BUILD_TESTS "Build CFDX tests" ON)
option(CFDX_BUILD_PYTHON "Build CFDX Python bindings" OFF)
option(CFDX_ENABLE_GPU "Enable CUDA backend (requires full CUDA toolkit)" OFF)
option(CFDX_ENABLE_MPI "Enable MPI support" ON)
option(CFDX_USE_HDF5_SERIAL "Use locally-built HDF5 serial" ON)
option(CFDX_ENABLE_PARALLEL_HDF5 "Enable MPI-enabled HDF5 distributed checkpoints" OFF)
option(CFDX_BUILD_RUNTIME_TESTS "Build execution/GPU-OOC foundation tests" ON)
option(CFDX_ENABLE_LONG_VALIDATION "Register full long-running validation campaigns" OFF)

if(CFDX_ENABLE_GPU)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 14)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

find_package(Python3 COMPONENTS Interpreter QUIET)

# ---------------------------------------------------------------------------
# Sanitizer build type
# ---------------------------------------------------------------------------
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -DNDEBUG")

# Sanitizer build type (DebugSanitizers)
if(CMAKE_BUILD_TYPE STREQUAL "DebugSanitizers")
    set(CMAKE_CXX_FLAGS_DEBUGSANITIZERS "-g -O1 -fsanitize=address,undefined -fno-omit-frame-pointer")
    set(CMAKE_C_FLAGS_DEBUGSANITIZERS "-g -O1 -fsanitize=address,undefined -fno-omit-frame-pointer")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUGSANITIZERS "-fsanitize=address,undefined")
    set(CMAKE_SHARED_LINKER_FLAGS_DEBUGSANITIZERS "-fsanitize=address,undefined")
    set(CMAKE_MODULE_LINKER_FLAGS_DEBUGSANITIZERS "-fsanitize=address,undefined")
endif()

# ---------------------------------------------------------------------------
# Compiler flags
# ---------------------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion)
    if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug" AND NOT CMAKE_BUILD_TYPE STREQUAL "DebugSanitizers")
        add_compile_options(-O2 -DNDEBUG)
    endif()
endif()

# ---------------------------------------------------------------------------
# HDF5 detection (before any target uses it)
# ---------------------------------------------------------------------------
if(CFDX_USE_HDF5_SERIAL)
    set(CFDX_HDF5_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party/hdf5")
endif()

find_package(HDF5 QUIET)
if(HDF5_FOUND)
    set(CFDX_HAS_HDF5 ON)
    message(STATUS "Using system HDF5 ${HDF5_VERSION}")
    if(CFDX_ENABLE_PARALLEL_HDF5 AND NOT HDF5_IS_PARALLEL)
        message(FATAL_ERROR "CFDX_ENABLE_PARALLEL_HDF5=ON but detected HDF5 is not parallel-enabled")
    endif()
    # System HDF5 provides these via find_package
    if(NOT HDF5_INCLUDE_DIRS)
        set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
    endif()
    if(NOT HDF5_C_LIBRARIES)
        set(HDF5_C_LIBRARIES ${HDF5_C_LIBRARIES})
    endif()
    if(NOT HDF5_HL_LIBRARIES)
        set(HDF5_HL_LIBRARIES ${HDF5_HL_LIBRARIES})
    endif()
elseif(EXISTS "${CFDX_HDF5_DIR}/include/H5public.h")
    set(CFDX_HAS_HDF5 ON)
    message(STATUS "Using locally-built HDF5 (static) at ${CFDX_HDF5_DIR}")
    set(HDF5_INCLUDE_DIRS "${CFDX_HDF5_DIR}/include")
    set(HDF5_C_LIBRARIES "${CFDX_HDF5_DIR}/lib/libhdf5.a")
    set(HDF5_HL_LIBRARIES "${CFDX_HDF5_DIR}/lib/libhdf5_hl.a")
else()
    set(CFDX_HAS_HDF5 OFF)
    message(WARNING "HDF5 not found. HDF5-dependent features disabled.")
endif()

if(CFDX_ENABLE_PARALLEL_HDF5 AND NOT HDF5_FOUND)
    message(FATAL_ERROR "CFDX_ENABLE_PARALLEL_HDF5=ON requires a parallel HDF5 installation")
endif()

# ---------------------------------------------------------------------------
# CFDX core library - Sources organized by module
# ---------------------------------------------------------------------------

# Mesh module
set(CFDX_MESH_SOURCES
    src/cfdx/core/mesh/point.cpp
    src/cfdx/core/mesh/face.cpp
    src/cfdx/core/mesh/ownership.cpp
    src/cfdx/core/mesh/cell.cpp
    src/cfdx/core/mesh/mesh.cpp
    src/cfdx/core/mesh/boundary.cpp
)

# Importers module
set(CFDX_LINALG_SOURCES
    src/cfdx/core/linalg/hypre_amg.cpp
)

# Memory sources
set(CFDX_MEMORY_SOURCES
    src/cfdx/core/memory/memory_planner.cpp
    src/cfdx/runtime/gpu/gpu_kernels.cpp
)

# Geometry module
set(CFDX_GEOMETRY_SOURCES
    src/cfdx/core/geometry/face_geometry.cpp
    src/cfdx/core/geometry/cell_geometry.cpp
    src/cfdx/core/geometry/mesh_quality.cpp
    src/cfdx/core/geometry/mesh_validator.cpp
)

# Field module
set(CFDX_FIELD_SOURCES
    src/cfdx/core/field/field.cpp
    src/cfdx/core/field/storage.cpp
)

# IO module
set(CFDX_IO_SOURCES
    src/cfdx/io/vtu/vtu_writer.cpp
    src/cfdx/io/openfoam/openfoam_importer.cpp
    src/cfdx/io/gmsh/gmsh_importer.cpp
    src/cfdx/io/mesh/mesh_importer.cpp
    src/cfdx/io/cfdx_io/gap_analysis.cpp
    src/cfdx/io/cfdx_io/mapping_rules.cpp
)
if(CFDX_HAS_HDF5)
    list(APPEND CFDX_IO_SOURCES
        src/cfdx/io/hdf5/hdf5_writer.cpp
        src/cfdx/io/hdf5/hdf5_reader.cpp
    )
    list(APPEND CFDX_IO_SOURCES
        src/cfdx/io/hdf5/case_hdf5_io.cpp
    )
else()
    list(APPEND CFDX_IO_SOURCES src/cfdx/io/hdf5/hdf5_unavailable.cpp)
endif()

# Physics sources
set(CFDX_PHYSICS_SOURCES
    src/cfdx/physics/equation_of_state.h
    src/cfdx/physics/multi_component.h
    src/cfdx/physics/source_term_linearization.h
    src/cfdx/physics/transport_models.h
    src/cfdx/physics/temporal.h
    src/cfdx/physics/low_storage_time_integration.h
    src/cfdx/physics/adaptive_cfl.h
    src/cfdx/physics/timestep_control.h
    src/cfdx/physics/dual_time_physical_lifecycle.h
    src/cfdx/physics/dual_time_navier_stokes.h
    src/cfdx/physics/advanced_convergence.h
    src/cfdx/physics/continuation.h
    src/cfdx/physics/forces.h
    src/cfdx/physics/local_time_stepping.h
    src/cfdx/physics/low_mach.h
    src/cfdx/physics/compressible_flux.h
    src/cfdx/physics/radiation_models.h
    src/cfdx/physics/radiation_s2s.h
    src/cfdx/physics/boussinesq.h
    src/cfdx/physics/spalart_allmaras.h
    src/cfdx/physics/wall_distance.h
    src/cfdx/physics/wall_distance_turbulence.h
    src/cfdx/physics/radiation_streaming.h
    src/cfdx/thermodynamics/thermo_state.h
    src/cfdx/thermodynamics/thermo_cache.h
)

# Thermodynamics sources (P2 reorganization)
set(CFDX_THERMODYNAMICS_SOURCES
    src/cfdx/thermodynamics/equation_of_state.h
)

# Chemistry sources (P2 reorganization)
set(CFDX_CHEMISTRY_SOURCES
    src/cfdx/chemistry/source_terms.h
)

# Transport sources (P2 reorganizization)
set(CFDX_TRANSPORT_SOURCES
    src/cfdx/transport/viscosity/viscosity.h
    src/cfdx/transport/conductivity/conductivity.h
    src/cfdx/transport/diffusivity/diffusivity.h
)

set(CFDX_CORE_SOURCES
    src/cfdx/application/simulation_controller.cpp
    ${CFDX_MESH_SOURCES}
    ${CFDX_GEOMETRY_SOURCES}
    ${CFDX_FIELD_SOURCES}
    ${CFDX_IO_SOURCES}
    ${CFDX_PHYSICS_SOURCES}
    ${CFDX_THERMODYNAMICS_SOURCES}
    ${CFDX_CHEMISTRY_SOURCES}
    ${CFDX_TRANSPORT_SOURCES}
    ${CFDX_LINALG_SOURCES}
    ${CFDX_MEMORY_SOURCES}
)

# Headers (for IDE visibility, not compiled)
set(CFDX_CORE_HEADERS
    src/cfdx/core/mesh/point.h
    src/cfdx/core/mesh/face.h
    src/cfdx/core/mesh/ownership.h
    src/cfdx/core/mesh/cell.h
    src/cfdx/core/mesh/boundary.h
    src/cfdx/core/mesh/mesh.h
    src/cfdx/core/geometry/face_geometry.h
    src/cfdx/core/geometry/cell_geometry.h
    src/cfdx/core/geometry/mesh_quality.h
    src/cfdx/core/geometry/mesh_validator.h
    src/cfdx/core/field/field.h
    src/cfdx/core/field/storage.h
    src/cfdx/core/boundary/boundary_field.h
    src/cfdx/core/boundary/boundary.h
    src/cfdx/core/boundary/boundary_role.h
    src/cfdx/core/boundary/boundary_condition.h
    src/cfdx/core/boundary/boundary_constraint.h
    src/cfdx/core/boundary/boundary_validator.h
    src/cfdx/core/boundary/mathematical_condition.h
    src/cfdx/core/boundary/value_provider.h
    src/cfdx/core/boundary/legacy_patch_field_adapter.h
    src/cfdx/core/boundary/flow_boundary_conditions.h
    src/cfdx/core/numerics/interpolation.h
    src/cfdx/core/numerics/gradient.h
    src/cfdx/core/numerics/gradient_stencil.h
    src/cfdx/core/numerics/divergence.h
    src/cfdx/core/numerics/laplacian.h
    src/cfdx/core/numerics/flux.h
    src/cfdx/core/numerics/integrate.h
    src/cfdx/core/numerics/source_term.h
    src/cfdx/core/numerics/temporal.h
    src/cfdx/core/numerics/convection.h
    src/cfdx/core/numerics/backend_equivalence.h
    src/cfdx/physics/incompressible.h
    src/cfdx/physics/pressure_velocity.h
    src/cfdx/physics/pressure_velocity_algorithms.h
    src/cfdx/physics/turbulence.h
    src/cfdx/physics/thermal.h
    src/cfdx/physics/thermophysical_models.h
    src/cfdx/physics/turbulence_models.h
    src/cfdx/physics/radiation.h
    src/cfdx/physics/m1_m4_models.h
    src/cfdx/physics/solver_control.h
    src/cfdx/physics/finite_volume_transport.h
    src/cfdx/physics/boundary_constraint_fvm.h
    src/cfdx/physics/steady_incompressible_solver.h
    src/cfdx/physics/dual_time_driver.h
    src/cfdx/physics/dual_time_adaptive.h
    src/cfdx/physics/dual_time_field_driver.h
    src/cfdx/physics/dual_time_field_integrator.h
    src/cfdx/physics/dual_time_navier_stokes.h
    src/cfdx/physics/turbulence_transport.h
    src/cfdx/physics/turbulence_solver.h
    src/cfdx/physics/sst_solver.h
    src/cfdx/physics/energy_solver.h
    src/cfdx/physics/thermal_zones.h
    src/cfdx/physics/production_energy_stage.h
    src/cfdx/physics/cht_solver.h
    src/cfdx/physics/radiation_solver.h
    src/cfdx/physics/radiation_s2s.h
    src/cfdx/physics/equation_of_state.h
    src/cfdx/physics/transport_models.h
    src/cfdx/physics/multi_component.h
    src/cfdx/core/linalg/sparse_matrix.h
    src/cfdx/core/linalg/vector.h
    src/cfdx/core/linalg/linear_system.h
    src/cfdx/core/linalg/cg_solver.h
    src/cfdx/core/linalg/bicgstab_solver.h
    src/cfdx/core/linalg/gmres_solver.h
    src/cfdx/core/linalg/linear_operator.h
    src/cfdx/core/linalg/solver_workspace.h
    src/cfdx/core/linalg/linear_solver_context.h
    src/cfdx/core/linalg/linear_solver_models.h
    src/cfdx/core/linalg/linear_solver_dispatch.h
    src/cfdx/core/linalg/null_space.h
    src/cfdx/core/linalg/matrix_free_fv_operator.h
    src/cfdx/core/linalg/amg_setup_footprint.h
    src/cfdx/core/linalg/amg_preconditioner.h
    src/cfdx/core/linalg/hypre_amg.h
    src/cfdx/core/linalg/coupled_amg_schur.h
    src/cfdx/core/linalg/native_fieldsplit.h
    src/cfdx/core/linalg/mixed_precision.h
    src/cfdx/core/linalg/chebyshev_smoother.h
    src/cfdx/core/linalg/fv_operator.h
    src/cfdx/core/linalg/krylov_controls.h
    src/cfdx/core/linalg/communication_avoiding.h
    src/cfdx/core/parallel/distributed_execution.h
    src/cfdx/core/parallel/distributed_poisson.h
    src/cfdx/runtime/gpu/cuda_poisson.h
    src/cfdx/core/memory/reuse_pool.h
    src/cfdx/core/memory/field_lifetime.h
    src/cfdx/core/linalg/mpi_overlap.h
    src/cfdx/core/linalg/block_schur.h
    src/cfdx/core/mesh/solver_mesh.h
    src/cfdx/core/mesh/sfc_ordering.h
    src/cfdx/core/parallel/restart_mapping.h
    src/cfdx/thermodynamics/equation_of_state.h
    src/cfdx/chemistry/source_terms.h
    src/cfdx/transport/viscosity/viscosity.h
    src/cfdx/transport/conductivity/conductivity.h
    src/cfdx/transport/diffusivity/diffusivity.h
    src/cfdx/io/hdf5/hdf5_writer.h
    src/cfdx/io/hdf5/hdf5_reader.h
    src/cfdx/io/hdf5/case_hdf5_io.h
    src/cfdx/io/hdf5/mini_json.h
    src/cfdx/io/vtu/vtu_writer.h
)

add_executable(cfdx_production_solver apps/production_solver/main.cpp)
target_link_libraries(cfdx_production_solver PRIVATE cfdx_core)

add_executable(cfdx_thermal_solver apps/thermal_solver/main.cpp)
target_link_libraries(cfdx_thermal_solver PRIVATE cfdx_core)

add_executable(cfdx_convert
    apps/cfdx_convert/main.cpp
    apps/cfdx_convert/PythonBridge.cpp
)
target_link_libraries(cfdx_convert PRIVATE cfdx_core)
if(CFDX_HAS_HDF5)
    target_include_directories(cfdx_convert PRIVATE ${HDF5_INCLUDE_DIRS})
    target_link_libraries(cfdx_convert PRIVATE ${HDF5_C_LIBRARIES} ${HDF5_HL_LIBRARIES} z)
    target_compile_definitions(cfdx_convert PRIVATE CFDX_HAS_HDF5=1)
else()
    target_compile_definitions(cfdx_convert PRIVATE CFDX_HAS_HDF5=0)
endif()

add_executable(cfdx_wall_distance_complex_benchmark apps/verification/wall_distance_complex_benchmark.cpp)
target_link_libraries(cfdx_wall_distance_complex_benchmark PRIVATE cfdx_core)

add_executable(cfdx_naca0012_wall_distance_comparison
    apps/verification/naca0012_wall_distance_comparison.cpp)
target_link_libraries(cfdx_naca0012_wall_distance_comparison PRIVATE cfdx_core)

add_executable(cfdx_wall_distance_controlled_benchmark apps/verification/wall_distance_controlled_benchmark.cpp)
target_link_libraries(cfdx_wall_distance_controlled_benchmark PRIVATE cfdx_core)

add_executable(cfdx_wall_distance_poisson_manufactured apps/verification/wall_distance_poisson_manufactured.cpp)
target_link_libraries(cfdx_wall_distance_poisson_manufactured PRIVATE cfdx_core)

add_executable(cfdx_wall_distance_poisson_decomposition apps/verification/wall_distance_poisson_decomposition.cpp)
target_link_libraries(cfdx_wall_distance_poisson_decomposition PRIVATE cfdx_core)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(cfdx_wall_distance_poisson_manufactured PRIVATE
        -Werror=multichar -Werror=unused-parameter)
endif()

add_library(cfdx_core STATIC ${CFDX_CORE_SOURCES})
target_include_directories(cfdx_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_compile_definitions(cfdx_core PRIVATE CFDX_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
if(Python3_Interpreter_FOUND)
    target_compile_definitions(cfdx_core PRIVATE CFDX_PYTHON_EXECUTABLE="${Python3_EXECUTABLE}")
endif()

if(CFDX_HAS_HDF5)
    target_include_directories(cfdx_core PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(cfdx_core PUBLIC ${HDF5_C_LIBRARIES} ${HDF5_HL_LIBRARIES} z)
    target_compile_definitions(cfdx_core PUBLIC CFDX_HAS_HDF5=1)
else()
    target_compile_definitions(cfdx_core PUBLIC CFDX_HAS_HDF5=0)
endif()

if(CFDX_ENABLE_GPU)
    target_sources(cfdx_core PRIVATE
        src/cfdx/runtime/gpu/cuda_kernels.cu
        src/cfdx/runtime/gpu/cuda_poisson.cu
    )
    target_compile_definitions(cfdx_core PUBLIC CFDX_ENABLE_GPU=1)
endif()

if(CFDX_ENABLE_MPI)
    find_package(MPI QUIET)
    if(MPI_FOUND)
        target_link_libraries(cfdx_core PUBLIC MPI::MPI_CXX)
        target_compile_definitions(cfdx_core PUBLIC CFDX_HAS_MPI=1)
    else()
        target_compile_definitions(cfdx_core PUBLIC CFDX_HAS_MPI=0)
    endif()
else()
    target_compile_definitions(cfdx_core PUBLIC CFDX_HAS_MPI=0)
endif()

if(CFDX_ENABLE_PARALLEL_HDF5)
    target_compile_definitions(cfdx_core PUBLIC CFDX_ENABLE_PARALLEL_HDF5=1)
endif()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(CFDX_BUILD_TESTS)
    enable_testing()

    # Helper function for tests
    macro(add_cfdx_test test_name source_file)
        add_executable(${test_name} ${source_file})
        target_link_libraries(${test_name} PRIVATE cfdx_core)
        target_include_directories(${test_name} PRIVATE ${CMAKE_SOURCE_DIR}/tests)
        add_test(NAME ${test_name} COMMAND ${test_name})
    endmacro()

    # Recovered standalone regressions: these sources previously had no CTest
    # target, so API drift in them was invisible to normal builds.
    add_cfdx_test(test_temporal tests/numerical/test_temporal.cpp)
    add_cfdx_test(test_discretization_verification tests/validation/test_discretization_verification.cpp)
    add_cfdx_test(test_matrix_diagnostics tests/unit/test_matrix_diagnostics.cpp)
    add_cfdx_test(test_application_layer tests/unit/test_application_layer.cpp)
    add_cfdx_test(test_memory_budget tests/unit/test_memory_budget.cpp)
    add_cfdx_test(test_memory_planner tests/unit/test_memory_planner.cpp)
    set_tests_properties(test_temporal test_discretization_verification PROPERTIES
        LABELS "numerical;verification;recovered")
    set_tests_properties(test_matrix_diagnostics test_application_layer
        test_memory_budget test_memory_planner PROPERTIES LABELS "unit;recovered")

    if(CFDX_ENABLE_MPI AND MPI_FOUND)
        set(CFDX_RECOVERED_MPIEXEC_NUMPROCS_FLAG "${MPIEXEC_NUMPROCS_FLAG}")
        if(CFDX_RECOVERED_MPIEXEC_NUMPROCS_FLAG STREQUAL "")
            set(CFDX_RECOVERED_MPIEXEC_NUMPROCS_FLAG "-n")
        endif()
        add_executable(test_mpi_deterministic_fused tests/unit/test_mpi_deterministic_fused.cpp)
        target_link_libraries(test_mpi_deterministic_fused PRIVATE cfdx_core MPI::MPI_CXX)
        add_test(NAME test_mpi_deterministic_fused
            COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_RECOVERED_MPIEXEC_NUMPROCS_FLAG} 2
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_deterministic_fused> ${MPIEXEC_POSTFLAGS})
        set_tests_properties(test_mpi_deterministic_fused PROPERTIES
            LABELS "unit;mpi;recovered" PROCESSORS 2)

        # This source calls MPI-IO entry points directly and requires two ranks.
        # Do not substitute a serial HDF5 build or a one-rank smoke test.
        if(CFDX_ENABLE_PARALLEL_HDF5 AND HDF5_FOUND AND HDF5_IS_PARALLEL)
            add_executable(test_parallel_hdf5 tests/unit/test_parallel_hdf5.cpp)
            target_link_libraries(test_parallel_hdf5 PRIVATE cfdx_core MPI::MPI_CXX)
            target_include_directories(test_parallel_hdf5 PRIVATE ${CMAKE_SOURCE_DIR}/tests)
            add_test(NAME test_parallel_hdf5
                COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_RECOVERED_MPIEXEC_NUMPROCS_FLAG} 2
                    ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_parallel_hdf5> ${MPIEXEC_POSTFLAGS})
            set_tests_properties(test_parallel_hdf5 PROPERTIES
                LABELS "unit;mpi;hdf5;recovered" PROCESSORS 2)
        endif()
    endif()

    # N8 physical/production qualification reuses the existing validation
    # executables and records a machine-readable campaign report. It is
    # intentionally a separate CTest entry so the individual physical or
    # algebraic tests remain independently visible in CI.
    if(Python3_Interpreter_FOUND)
        add_test(NAME test_n8_physical_qualification
                 COMMAND ${Python3_EXECUTABLE}
                         ${CMAKE_SOURCE_DIR}/scripts/n8_physical_qualification.py
                         --build-dir ${CMAKE_BINARY_DIR})
        set_tests_properties(test_n8_physical_qualification PROPERTIES
            LABELS "n8;qualification;validation;sanitizer-slow"
            RUN_SERIAL TRUE)
    endif()

    # N11 conservation/boundedness qualification reuses the existing
    # executable evidence and writes a machine-readable campaign report.
    # It is intentionally separate from the individual tests so their
    # independent evidence remains visible in CI.
    if(Python3_Interpreter_FOUND)
        add_test(NAME test_n11_conservation_boundedness
                 COMMAND ${Python3_EXECUTABLE}
                         ${CMAKE_SOURCE_DIR}/scripts/n11_conservation_boundedness.py
                         --build-dir ${CMAKE_BINARY_DIR})
        set_tests_properties(test_n11_conservation_boundedness PROPERTIES
            LABELS "n11;qualification;validation;long"
            RUN_SERIAL TRUE)

        add_test(NAME test_n11_campaign_integrity
                 COMMAND ${Python3_EXECUTABLE}
                         ${CMAKE_SOURCE_DIR}/tests/python/test_n11_campaign_integrity.py)
        set_tests_properties(test_n11_campaign_integrity PROPERTIES
            LABELS "n11;unit;validation")
    endif()

    # Unit tests - Core
    add_cfdx_test(test_point tests/unit/test_point.cpp)
    add_cfdx_test(test_face tests/unit/test_face.cpp)
    add_cfdx_test(test_ownership tests/unit/test_ownership.cpp)
    add_cfdx_test(test_cell tests/unit/test_cell.cpp)
    add_cfdx_test(test_boundary tests/unit/test_boundary.cpp)
    add_cfdx_test(test_mesh tests/unit/test_mesh.cpp)
    add_cfdx_test(test_field tests/unit/test_field.cpp)
    add_cfdx_test(test_field_storage tests/unit/test_field_storage.cpp)
    add_cfdx_test(test_field_wall_hardening tests/unit/test_field_wall_hardening.cpp)
    add_cfdx_test(test_boundary_field tests/unit/test_boundary_field.cpp)
    add_cfdx_test(test_boundary_contract tests/unit/test_boundary_contract.cpp)
    add_cfdx_test(test_boundary_fvm_constraints tests/unit/test_boundary_fvm_constraints.cpp)
    add_cfdx_test(test_boundary_validator tests/unit/test_boundary_validator.cpp)
    add_cfdx_test(test_interpolation tests/unit/test_interpolation.cpp)
    add_cfdx_test(test_gradient tests/unit/test_gradient.cpp)
    add_cfdx_test(test_pressure_velocity_system tests/unit/test_pressure_velocity_system.cpp)
    add_cfdx_test(test_oriented_numerics tests/unit/test_oriented_numerics.cpp)
    add_cfdx_test(test_divergence tests/unit/test_divergence.cpp)
    add_cfdx_test(test_flux tests/unit/test_flux.cpp)
    add_cfdx_test(test_spatial_scheme_contract tests/unit/test_spatial_scheme_contract.cpp)
    add_cfdx_test(test_conservation_boundedness tests/unit/test_conservation_boundedness.cpp)
    add_cfdx_test(test_transport_conservation tests/unit/test_transport_conservation.cpp)
    add_cfdx_test(test_conservation_assembly tests/unit/test_conservation_assembly.cpp)
    add_cfdx_test(test_integrate tests/unit/test_integrate.cpp)
    add_cfdx_test(test_forces tests/unit/test_forces.cpp)
    add_cfdx_test(test_probe_csv tests/unit/test_probe_csv.cpp)

    # Unit tests - Linear Algebra
    add_cfdx_test(test_sparse_matrix tests/unit/test_sparse_matrix.cpp)
    add_cfdx_test(test_vector tests/unit/test_vector.cpp)
    add_cfdx_test(test_linear_system tests/unit/test_linear_system.cpp)
    add_cfdx_test(test_cg_solver tests/unit/test_cg_solver.cpp)
    add_cfdx_test(test_null_space tests/unit/test_null_space.cpp)
    add_cfdx_test(test_n8_pinned_pressure_gauge tests/unit/test_n8_pinned_pressure_gauge.cpp)
    add_cfdx_test(test_bicgstab_solver tests/unit/test_bicgstab_solver.cpp)
    add_cfdx_test(test_gmres_solver tests/unit/test_gmres_solver.cpp)
    add_cfdx_test(test_fgmres_solver tests/unit/test_fgmres_solver.cpp)
    add_cfdx_test(test_krylov_preconditioning tests/unit/test_krylov_preconditioning.cpp)
    add_cfdx_test(test_krylov_diagnostics tests/unit/test_krylov_diagnostics.cpp)
    add_cfdx_test(test_stationary_solvers tests/unit/test_stationary_solvers.cpp)
    add_cfdx_test(test_performance_foundations tests/unit/test_performance_foundations.cpp)
    add_cfdx_test(test_performance_runtime tests/unit/test_performance_runtime.cpp)
    add_cfdx_test(test_phase6_memory_benchmark tests/benchmark/benchmark_phase6_memory.cpp)
    add_cfdx_test(test_amg_chebyshev tests/unit/test_amg_chebyshev.cpp)
    add_cfdx_test(test_native_external_methods tests/unit/test_native_external_methods.cpp)
    add_cfdx_test(test_matrix_characteristics tests/unit/test_matrix_characteristics.cpp)
    add_cfdx_test(test_amg_preconditioner_qualification tests/validation/test_amg_preconditioner_qualification.cpp)
    add_cfdx_test(test_matrix_free_fv_operator tests/unit/test_matrix_free_fv_operator.cpp)
    add_cfdx_test(test_n11_scaling_equivalence tests/validation/test_n11_scaling_equivalence.cpp)
    add_cfdx_test(test_n11_production_scaling tests/validation/test_n11_production_scaling.cpp)
    add_cfdx_test(test_n11_unit_invariance_robustness tests/validation/test_n11_unit_invariance_robustness.cpp)
    set_tests_properties(test_n11_unit_invariance_robustness PROPERTIES LABELS "n11;qualification;validation")
    add_cfdx_test(test_backend_equivalence tests/unit/test_backend_equivalence.cpp)
    add_cfdx_test(test_n13_restart_equivalence tests/validation/test_n13_restart_equivalence.cpp)

    # Unit tests - Numerical
    add_cfdx_test(test_laplacian tests/unit/test_laplacian.cpp)
    add_cfdx_test(test_diffusion_diagnostics tests/unit/test_diffusion_diagnostics.cpp)
    add_cfdx_test(test_phase123_regression tests/unit/test_phase123_regression.cpp)
add_cfdx_test(test_phase12_radiation_vv tests/validation/test_phase12_radiation_vv.cpp)
set_tests_properties(test_phase12_radiation_vv PROPERTIES LABELS "phase12;validation")
    add_cfdx_test(test_scalar_diffusion tests/unit/test_scalar_diffusion.cpp)
    add_cfdx_test(test_advanced_preconditioners tests/unit/test_advanced_preconditioners.cpp)
    add_cfdx_test(test_schur_preconditioner tests/unit/test_schur_preconditioner.cpp)

    # Numerical tests - Geometry
    add_cfdx_test(test_face_geometry tests/numerical/test_face_geometry.cpp)
    add_cfdx_test(test_cell_geometry tests/numerical/test_cell_geometry.cpp)
    add_cfdx_test(test_mesh_quality tests/numerical/test_mesh_quality.cpp)
    add_cfdx_test(test_geometry_invariants tests/unit/test_geometry_invariants.cpp)
    add_cfdx_test(test_dat_restart_identity tests/unit/test_dat_restart_identity.cpp)
    add_cfdx_test(test_mesh_family_volumes tests/numerical/test_mesh_family_volumes.cpp)
    add_cfdx_test(test_source_term tests/numerical/test_source_term.cpp)
    add_cfdx_test(test_temporal_physics tests/unit/test_temporal.cpp)
    add_cfdx_test(test_dual_time_stepping tests/unit/test_dual_time_stepping.cpp)
    add_cfdx_test(test_dual_time_adaptive tests/unit/test_dual_time_adaptive.cpp)
    add_cfdx_test(test_dual_time_field_driver tests/unit/test_dual_time_field_driver.cpp)
    add_cfdx_test(test_dual_time_field_integrator tests/unit/test_dual_time_field_integrator.cpp)
    add_cfdx_test(test_dual_time_physical_lifecycle tests/unit/test_dual_time_physical_lifecycle.cpp)
    add_cfdx_test(test_dual_time_navier_stokes tests/unit/test_dual_time_navier_stokes.cpp)
    add_cfdx_test(test_case_numerics_report tests/unit/test_case_numerics_report.cpp)
    add_cfdx_test(test_temporal_history tests/unit/test_temporal_history.cpp)
    add_cfdx_test(test_core_temporal tests/unit/test_core_temporal.cpp)
    add_cfdx_test(test_newton_field tests/unit/test_newton_field.cpp)
    add_cfdx_test(test_simulation_controller tests/unit/test_simulation_controller.cpp)
    add_cfdx_test(test_adaptive_cfl tests/unit/test_adaptive_cfl.cpp)
    add_cfdx_test(test_timestep_control tests/unit/test_timestep_control.cpp)
    add_cfdx_test(test_continuation tests/unit/test_continuation.cpp)
    add_cfdx_test(test_local_time_stepping tests/unit/test_local_time_stepping.cpp)
    add_cfdx_test(test_incompressible tests/unit/test_incompressible.cpp)
    add_cfdx_test(test_m1_m4_physics tests/unit/test_m1_m4_physics.cpp)
    add_cfdx_test(test_m1_m4_validation tests/validation/test_m1_m4_validation.cpp)
    add_cfdx_test(test_solver_control tests/unit/test_solver_control.cpp)
    add_cfdx_test(test_convergence_monitor tests/unit/test_convergence_monitor.cpp)
    add_cfdx_test(test_initialization tests/unit/test_initialization.cpp)
    add_executable(test_n7_convergence_qualification tests/validation/test_n7_convergence_qualification.cpp)
    target_link_libraries(test_n7_convergence_qualification PRIVATE cfdx_core)
    target_include_directories(test_n7_convergence_qualification PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    add_test(NAME test_n7_convergence_qualification COMMAND test_n7_convergence_qualification)
    set_tests_properties(test_n7_convergence_qualification PROPERTIES
        LABELS "n7;qualification;validation;sanitizer-slow")
    add_executable(test_n7_adaptive_relaxation_re400
        tests/validation/test_n7_adaptive_relaxation_re400.cpp)
    target_link_libraries(test_n7_adaptive_relaxation_re400 PRIVATE cfdx_core)
    target_include_directories(test_n7_adaptive_relaxation_re400 PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    add_test(NAME test_n7_adaptive_relaxation_re400
             COMMAND test_n7_adaptive_relaxation_re400)
    set_tests_properties(test_n7_adaptive_relaxation_re400 PROPERTIES
        LABELS "n7;qualification;validation;sanitizer-slow"
        TIMEOUT 1800)
    add_cfdx_test(test_n7_causal_overshoot_regression tests/unit/test_n7_causal_overshoot_regression.cpp)
    set_tests_properties(test_n7_causal_overshoot_regression PROPERTIES
        LABELS "n7;regression;unit")
    add_executable(test_n11_scalar_conservation tests/validation/test_n11_scalar_conservation.cpp)
    target_link_libraries(test_n11_scalar_conservation PRIVATE cfdx_core)
    target_include_directories(test_n11_scalar_conservation PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    add_test(NAME test_n11_scalar_conservation COMMAND test_n11_scalar_conservation)
    set_tests_properties(test_n11_scalar_conservation PROPERTIES
        LABELS "n11;qualification;validation")
    add_cfdx_test(test_finite_volume_transport tests/unit/test_finite_volume_transport.cpp)
    add_cfdx_test(test_vof_curvature_contract tests/unit/test_vof_curvature_contract.cpp)
    add_cfdx_test(test_steady_incompressible_solver tests/validation/test_steady_incompressible_solver.cpp)
    add_cfdx_test(test_vtu_writer_validation tests/validation/test_vtu_writer_validation.cpp)
    add_cfdx_test(test_solver_vtu_provenance tests/validation/test_solver_vtu_provenance.cpp)
    add_executable(test_ghia_cavity tests/validation/test_ghia_cavity.cpp)
    target_link_libraries(test_ghia_cavity PRIVATE cfdx_core)
    target_include_directories(test_ghia_cavity PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    add_test(NAME test_ghia_cavity_quick COMMAND test_ghia_cavity --quick)
    set_tests_properties(test_ghia_cavity_quick PROPERTIES
        LABELS "phase13;validation")
    # Production PCD is an explicit N8 qualification test. It must never be
    # pulled into the quick smoke path because its coupled Krylov/preconditioner
    # cost is materially different from the SIMPLE Ghia regression.
    if(CFDX_ENABLE_LONG_VALIDATION)
        # Run a bounded low-resolution preflight first. The dependency prevents
        # the physical PCD qualification from starting when the nested solver
        # already exhibits runaway behaviour on the smaller case.
        add_test(NAME test_ghia_cavity_n8_pcd_preflight COMMAND test_ghia_cavity --pcd-preflight)
        set_tests_properties(test_ghia_cavity_n8_pcd_preflight PROPERTIES
            LABELS "n8;qualification;validation;long;preflight"
            TIMEOUT 120)
        add_test(NAME test_ghia_cavity_n8_pcd COMMAND test_ghia_cavity --pcd)
        set_tests_properties(test_ghia_cavity_n8_pcd PROPERTIES
            LABELS "n8;qualification;validation;long"
            TIMEOUT 720
            DEPENDS test_ghia_cavity_n8_pcd_preflight)
        # PCD Schur-action A/B benchmark on the REAL assembled Ghia coupled
        # systems (evidence-only replay through the coupled_assembly_callback
        # hook; the test header documents the five action variants). This is a
        # measurement test: the physical case remains the responsibility of
        # test_ghia_cavity_n8_pcd, no production tolerance, gate, budget or
        # oracle is modified, and the baseline variant must reproduce the
        # production 1e-10 convergence on every captured stage or the test
        # fails (capture-fidelity guard). The ctest budget bounds the capture
        # run (five nonlinear iterations) plus the ten replayed coupled
        # solves, each capped by the production 5000 coupled budget.
        add_cfdx_test(test_pcd_schur_action_benchmark
            tests/validation/test_pcd_schur_action_benchmark.cpp)
        set_tests_properties(test_pcd_schur_action_benchmark PROPERTIES
            LABELS "n8;schur;benchmark;validation;long"
            TIMEOUT 1200)
    endif()
    # Merge-gated observed-order tripwire. Three resolutions at 16/32/64 give a
    # Richardson-style order in about 90s, so observed-order evidence now gates
    # every merge instead of running only in the skipped long-validation
    # workflow. The 32/64/128 campaign stays behind CFDX_ENABLE_LONG_VALIDATION.
    add_test(NAME test_ghia_cavity_order COMMAND test_ghia_cavity --gate)
    set_tests_properties(test_ghia_cavity_order PROPERTIES
        LABELS "phase13;validation;campaign"
        TIMEOUT 900)
    if(CFDX_ENABLE_LONG_VALIDATION)
        add_test(NAME test_ghia_cavity COMMAND test_ghia_cavity)
        set_tests_properties(test_ghia_cavity PROPERTIES
            LABELS "phase13;validation;long")
    endif()
    add_executable(test_phase9_acceptance tests/validation/test_phase9_acceptance.cpp)
    target_link_libraries(test_phase9_acceptance PRIVATE cfdx_core)
    target_include_directories(test_phase9_acceptance PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    # Full N8 pressure-velocity matrix: deliberately not the --quick smoke test.
    # This exercises SIMPLE, SIMPLEC, PISO, PIMPLE, fractional-step and COUPLED.
    add_test(NAME test_n8_pressure_velocity_matrix COMMAND test_phase9_acceptance)
    set_tests_properties(test_n8_pressure_velocity_matrix PROPERTIES
        LABELS "n8;qualification;validation;sanitizer-slow"
        RUN_SERIAL TRUE)
    add_test(NAME test_couette_quick COMMAND test_phase9_acceptance --quick)
    # Phase 9 diagnostics must be deterministic under CTest. The validation
    # workflow may set these variables at the shell level, but CTest test
    # properties are the authoritative per-test environment and avoid losing
    # the microscope when the test is launched through CTest.
    set_tests_properties(test_couette_quick PROPERTIES
        LABELS "phase9;validation;coupled-ns"
        ENVIRONMENT "CFDX_DEBUG_CELL=33;CFDX_FREEZE_STATE=1;CFDX_DEBUG_COUPLED=1")
    if(CFDX_ENABLE_LONG_VALIDATION)
        add_test(NAME test_phase9_acceptance COMMAND test_phase9_acceptance)
        set_tests_properties(test_phase9_acceptance PROPERTIES
            LABELS "phase9;validation;coupled-ns;long"
            ENVIRONMENT "CFDX_DEBUG_CELL=33;CFDX_FREEZE_STATE=1;CFDX_DEBUG_COUPLED=1")
    endif()
    add_cfdx_test(test_m2_m4_solver_validation tests/validation/test_m2_m4_solver_validation.cpp)
    add_cfdx_test(test_phase10_turbulence_hardening tests/validation/test_phase10_turbulence_hardening.cpp)
    add_cfdx_test(test_turbulence_qualification_equations tests/unit/test_turbulence_qualification_equations.cpp)
    add_cfdx_test(test_production_turbulence_contract tests/unit/test_production_turbulence_contract.cpp)
    add_cfdx_test(test_production_radiation_contract tests/unit/test_production_radiation_contract.cpp)
    add_cfdx_test(test_n11_sa_negative_nu_tilde_policy tests/unit/test_n11_sa_negative_nu_tilde_policy.cpp)
    add_cfdx_test(test_n11_bound_reporting_hardening tests/unit/test_n11_bound_reporting_hardening.cpp)
    add_cfdx_test(test_phase10_model_registry tests/validation/test_phase10_model_registry.cpp)
    add_cfdx_test(test_wall_distance_methods tests/unit/test_wall_distance_methods.cpp)
    add_cfdx_test(test_wall_distance_invalid_results tests/unit/test_wall_distance_invalid_results.cpp)
    add_cfdx_test(test_wall_distance_godunov_ridges tests/unit/test_wall_distance_godunov_ridges.cpp)
    add_cfdx_test(test_wall_distance_advection_residual tests/unit/test_wall_distance_advection_residual.cpp)
    add_test(NAME benchmark_wall_distance_naca0012
        COMMAND cfdx_naca0012_wall_distance_comparison
                ${CMAKE_SOURCE_DIR}/tests/data/su2/mesh_NACA0012_inv.su2
                ${CMAKE_BINARY_DIR}/wall_distance_naca0012.csv 25)
    set_tests_properties(benchmark_wall_distance_naca0012 PROPERTIES
        LABELS "phase13;validation;benchmark;wall-distance"
        TIMEOUT 180)
    add_cfdx_test(test_n11_fw1_nearest_normal_fallback tests/unit/test_n11_fw1_nearest_normal_fallback.cpp)
add_cfdx_test(test_n11_fw2_positivity_clamp_counter tests/unit/test_n11_fw2_positivity_clamp_counter.cpp)
    add_cfdx_test(test_n11_b1_b2_b3_linalg_guards tests/unit/test_n11_b1_b2_b3_linalg_guards.cpp)
    add_cfdx_test(test_n11_b4_b5_amg_footprint tests/unit/test_n11_b4_b5_amg_footprint.cpp)
    add_test(NAME benchmark_wall_distance_complex
             COMMAND cfdx_wall_distance_complex_benchmark ${CMAKE_BINARY_DIR}/wall_distance_complex_benchmark.csv)
    set_tests_properties(benchmark_wall_distance_complex PROPERTIES LABELS "phase13;validation;benchmark")
    add_cfdx_test(test_phase10_sst_limiter tests/validation/test_phase10_sst_limiter.cpp)
    add_cfdx_test(test_analytical_benchmarks tests/validation/test_analytical_benchmarks.cpp)
    add_cfdx_test(test_vmfl036_reference_case tests/validation/test_vmfl036_reference_case.cpp)
    add_cfdx_test(test_axisymmetric_vmfl036_kernel tests/validation/test_axisymmetric_vmfl036_kernel.cpp)
    add_executable(test_vmfl036_axisymmetric tests/validation/test_vmfl036_axisymmetric.cpp)
    target_link_libraries(test_vmfl036_axisymmetric PRIVATE cfdx_core)
    target_include_directories(test_vmfl036_axisymmetric PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    target_compile_definitions(test_vmfl036_axisymmetric PRIVATE CFDX_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
    if(Python3_Interpreter_FOUND)
        target_compile_definitions(test_vmfl036_axisymmetric PRIVATE CFDX_PYTHON_EXECUTABLE="${Python3_EXECUTABLE}")
    endif()
    add_test(NAME test_vmfl036_axisymmetric_quick COMMAND test_vmfl036_axisymmetric --quick)
    set_tests_properties(test_vmfl036_axisymmetric_quick PROPERTIES LABELS "phase14;validation;axisymmetric")
    if(CFDX_ENABLE_LONG_VALIDATION)
        add_test(NAME test_vmfl036_axisymmetric COMMAND test_vmfl036_axisymmetric)
        set_tests_properties(test_vmfl036_axisymmetric PROPERTIES LABELS "phase14;validation;axisymmetric;long")
        add_test(NAME test_vmfl036_alpha_sensitivity
                 COMMAND test_vmfl036_axisymmetric --alpha-sweep
                         ${CMAKE_BINARY_DIR}/validation-report/vmfl036-alpha-sweep.json)
        set_tests_properties(test_vmfl036_alpha_sensitivity PROPERTIES
            LABELS "phase14;validation;axisymmetric;long;validation-total")
    endif()
    set_tests_properties(test_axisymmetric_vmfl036_kernel PROPERTIES LABELS "phase14;validation;axisymmetric")
    target_compile_definitions(test_vmfl036_reference_case PRIVATE
        CFDX_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
    )
    if(Python3_Interpreter_FOUND)
        target_compile_definitions(test_vmfl036_reference_case PRIVATE
            CFDX_PYTHON_EXECUTABLE="${Python3_EXECUTABLE}"
        )
    endif()
    set_tests_properties(test_vmfl036_reference_case PROPERTIES LABELS "phase14;validation;reference")
    add_executable(test_bfs_qualification tests/validation/test_bfs_driver_smoke.cpp)
    target_link_libraries(test_bfs_qualification PRIVATE cfdx_core)
    add_test(NAME test_bfs_quick COMMAND test_bfs_qualification --quick ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/bfs_re200_n16.h5)
    set_tests_properties(test_bfs_quick PROPERTIES LABELS "qualification;validation")
    add_executable(test_naca0012_qualification tests/validation/test_naca0012_geometry.cpp)
    target_link_libraries(test_naca0012_qualification PRIVATE cfdx_core)
    add_test(NAME test_naca0012_quick COMMAND test_naca0012_qualification --quick ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/naca0012_laminar_n64.h5)
    set_tests_properties(test_naca0012_quick PROPERTIES LABELS "qualification;validation")
    if(Python3_Interpreter_FOUND)
        set_tests_properties(test_bfs_quick test_naca0012_quick PROPERTIES DEPENDS test_generate_2d_validation_meshes)
    endif()
    if(CFDX_ENABLE_LONG_VALIDATION)
        add_test(NAME test_bfs_qualification COMMAND test_bfs_qualification ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/bfs_re200_n16.h5 ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/bfs_re200_n32.h5 ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/bfs_re200_n64.h5)
        set_tests_properties(test_bfs_qualification PROPERTIES LABELS "qualification;validation;long")
        add_test(NAME test_naca0012_qualification COMMAND test_naca0012_qualification ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/naca0012_laminar_n64.h5 ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/naca0012_laminar_n128.h5 ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/naca0012_laminar_n256.h5)
        set_tests_properties(test_naca0012_qualification PROPERTIES LABELS "qualification;validation;long")
        if(Python3_Interpreter_FOUND)
            set_tests_properties(test_bfs_qualification test_naca0012_qualification PROPERTIES DEPENDS test_generate_2d_validation_meshes)
        endif()
    endif()
    # Wall-shear reconstruction audit (GG vs WLS vs 1D) - post-processing
    # only, no solver method change. The quick test runs the wall MMS on the
    # channel mesh without any solve; the long campaign adds the NACA0012
    # Re=1000 three-level comparison where the spread is the measurement
    # (explicitly NOT gated).
    add_executable(test_naca0012_wall_shear_audit tests/validation/test_naca0012_wall_shear_audit.cpp)
    target_link_libraries(test_naca0012_wall_shear_audit PRIVATE cfdx_core)
    add_test(NAME test_naca0012_wall_shear_audit_quick COMMAND test_naca0012_wall_shear_audit --quick ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/channel_n16.h5)
    set_tests_properties(test_naca0012_wall_shear_audit_quick PROPERTIES LABELS "qualification;validation")
    if(Python3_Interpreter_FOUND)
        set_tests_properties(test_naca0012_wall_shear_audit_quick PROPERTIES DEPENDS test_generate_2d_validation_meshes)
    endif()
    if(CFDX_ENABLE_LONG_VALIDATION)
        add_test(NAME test_naca0012_wall_shear_audit COMMAND test_naca0012_wall_shear_audit ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/channel_n16.h5 ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/naca0012_laminar_n64.h5 ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/naca0012_laminar_n128.h5 ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/naca0012_laminar_n256.h5)
        set_tests_properties(test_naca0012_wall_shear_audit PROPERTIES LABELS "qualification;validation;long")
        if(Python3_Interpreter_FOUND)
            set_tests_properties(test_naca0012_wall_shear_audit PROPERTIES DEPENDS test_generate_2d_validation_meshes)
        endif()
    endif()
    if(Python3_Interpreter_FOUND)
        add_test(NAME test_generate_2d_validation_meshes COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/generate_2d_validation_meshes.py all --output-dir ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes)
        set_tests_properties(test_generate_2d_validation_meshes PROPERTIES LABELS "qualification;mesh;validation")
    endif()
    add_executable(test_channel_qualification tests/validation/test_channel_qualification.cpp)
    target_link_libraries(test_channel_qualification PRIVATE cfdx_core)
    target_include_directories(test_channel_qualification PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    add_test(NAME test_channel_quick COMMAND test_channel_qualification --quick ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/channel_n16.h5)
    set_tests_properties(test_channel_quick PROPERTIES LABELS "qualification;validation;sanitizer-slow")
    if(Python3_Interpreter_FOUND)
        set_tests_properties(test_channel_quick PROPERTIES DEPENDS test_generate_2d_validation_meshes)
    endif()
    if(CFDX_ENABLE_LONG_VALIDATION)
        add_test(NAME test_channel_qualification COMMAND test_channel_qualification ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/channel_n16.h5 ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/channel_n32.h5 ${CMAKE_CURRENT_BINARY_DIR}/validation_meshes/channel_n64.h5)
        set_tests_properties(test_channel_qualification PROPERTIES LABELS "qualification;validation;long")
        if(Python3_Interpreter_FOUND)
            set_tests_properties(test_channel_qualification PROPERTIES DEPENDS test_generate_2d_validation_meshes)
        endif()
    endif()
    add_cfdx_test(test_numerical_model_verification tests/validation/test_numerical_model_verification.cpp)
    add_executable(test_poiseuille_diagnostics tests/validation/test_poiseuille_diagnostics.cpp)
    target_link_libraries(test_poiseuille_diagnostics PRIVATE cfdx_core)
    target_include_directories(test_poiseuille_diagnostics PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    add_test(NAME test_poiseuille_quick COMMAND test_poiseuille_diagnostics --quick)
    set_tests_properties(test_poiseuille_quick PROPERTIES
        LABELS "phase13;validation")
    if(CFDX_ENABLE_LONG_VALIDATION)
        add_test(NAME test_poiseuille_diagnostics COMMAND test_poiseuille_diagnostics)
        set_tests_properties(test_poiseuille_diagnostics PROPERTIES
            LABELS "phase13;validation;long")
    endif()
    add_cfdx_test(test_mms_scalar_diffusion tests/validation/test_mms_scalar_diffusion.cpp)
    add_executable(test_mesh_refinement_order tests/validation/test_mesh_refinement_order.cpp)
    target_link_libraries(test_mesh_refinement_order PRIVATE cfdx_core)
    target_include_directories(test_mesh_refinement_order PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    if(CFDX_ENABLE_LONG_VALIDATION)
        add_test(NAME test_mesh_refinement_order COMMAND test_mesh_refinement_order)
        set_tests_properties(test_mesh_refinement_order PROPERTIES
            LABELS "phase13;validation;long")
    endif()
    add_cfdx_test(test_nonorthogonal_skew_campaign tests/validation/test_nonorthogonal_skew_campaign.cpp)
    add_cfdx_test(test_nonorthogonal_laplacian_campaign tests/validation/test_nonorthogonal_laplacian_campaign.cpp)
    add_cfdx_test(test_transient_mms tests/validation/test_transient_mms.cpp)
    add_cfdx_test(test_transient_mms_closure tests/validation/test_transient_mms_closure.cpp)
    add_cfdx_test(test_transient_mms_3d_stiff tests/validation/test_transient_mms_3d_stiff.cpp)
    add_cfdx_test(test_n12_momentum_mms tests/validation/test_n12_momentum_mms.cpp)
    add_cfdx_test(test_n12_pressure_continuity_mms tests/validation/test_n12_pressure_continuity_mms.cpp)
    add_cfdx_test(test_n12_coupled_pressure_velocity_mms tests/validation/test_n12_coupled_pressure_velocity_mms.cpp)
    set_tests_properties(
        test_mms_scalar_diffusion
        test_transient_mms
        test_transient_mms_closure
        test_transient_mms_3d_stiff
        test_n12_momentum_mms
        test_n12_pressure_continuity_mms
        PROPERTIES LABELS "n12;mms;validation")
    set_tests_properties(
        test_n12_coupled_pressure_velocity_mms
        PROPERTIES LABELS "n12;mms;validation;long;validation-total")
    add_cfdx_test(test_polyhedral_laplacian_campaign tests/validation/test_polyhedral_laplacian_campaign.cpp)
    add_cfdx_test(test_polyhedral_poisson_convergence tests/validation/test_polyhedral_poisson_convergence.cpp)
    add_cfdx_test(test_convection_scheme_verification tests/validation/test_convection_scheme_verification.cpp)
    add_cfdx_test(test_convection_3d_verification tests/validation/test_convection_3d_verification.cpp)
    add_cfdx_test(test_convection_blending tests/validation/test_convection_blending.cpp)
    add_cfdx_test(test_n16_weno5z tests/validation/test_n16_weno5z.cpp)
    set_tests_properties(test_n16_weno5z PROPERTIES LABELS "n16;validation;experimental")
    add_cfdx_test(test_temporal_order_matrix tests/validation/test_temporal_order_matrix.cpp)
    add_cfdx_test(test_gradient_verification tests/validation/test_gradient_verification.cpp)
    add_cfdx_test(test_polyhedral_gradient_campaign tests/validation/test_polyhedral_gradient_campaign.cpp)
    add_cfdx_test(test_n10_difficult_mesh_robustness tests/validation/test_n10_difficult_mesh_robustness.cpp)
    add_executable(test_n10_external_fixture_import tests/validation/test_n10_external_fixture_import.cpp)
    target_link_libraries(test_n10_external_fixture_import PRIVATE cfdx_core)
    target_include_directories(test_n10_external_fixture_import PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    add_executable(test_n10_fixture_numerical_qualification tests/validation/test_n10_fixture_numerical_qualification.cpp)
    target_link_libraries(test_n10_fixture_numerical_qualification PRIVATE cfdx_core)
    target_include_directories(test_n10_fixture_numerical_qualification PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    add_executable(test_n10_fixture_production_case tests/validation/test_n10_fixture_production_case.cpp)
    target_link_libraries(test_n10_fixture_production_case PRIVATE cfdx_core)
    target_include_directories(test_n10_fixture_production_case PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    if(Python3_Interpreter_FOUND)
        add_test(NAME test_n10_qualification
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/n10_qualification.py
                    --build-dir ${CMAKE_BINARY_DIR})
        set_tests_properties(test_n10_qualification PROPERTIES
            LABELS "n10;qualification;validation"
            RUN_SERIAL TRUE)
        add_test(NAME test_n10_mesh_fixture_manifest
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/validate_mesh_fixture_manifest.py)
        set_tests_properties(test_n10_mesh_fixture_manifest PROPERTIES
            LABELS "n10;mesh;provenance;validation"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    endif()
    set_tests_properties(
        test_n10_difficult_mesh_robustness
        test_nonorthogonal_skew_campaign
        test_nonorthogonal_laplacian_campaign
        test_polyhedral_laplacian_campaign
        test_gradient_verification
        test_polyhedral_gradient_campaign
        test_polyhedral_poisson_convergence
        PROPERTIES LABELS "validation;qualification;n10")
    # The Poisson solution-convergence assembly is a black-box indicator
    # campaign (hundreds of production operator applications per level);
    # keep it out of the DebugSanitizers pass like the other
    # sanitizer-slow qualification tests. Release coverage stays
    # merge-gated.
    set_tests_properties(test_polyhedral_poisson_convergence PROPERTIES
        LABELS "n10;qualification;validation;sanitizer-slow"
        TIMEOUT 1800)
    add_cfdx_test(test_convection_polyhedral_campaign tests/validation/test_convection_polyhedral_campaign.cpp)
    add_cfdx_test(test_convection_solver_qualification tests/validation/test_convection_solver_qualification.cpp)
    set_tests_properties(test_convection_solver_qualification PROPERTIES LABELS "phase13;validation;qualification")
    set_tests_properties(test_convection_polyhedral_campaign PROPERTIES LABELS "phase13;validation;qualification")
    add_cfdx_test(test_fluent_vmfl_reference tests/validation/test_fluent_vmfl_reference.cpp)
    add_cfdx_test(test_benchmark_matrix tests/validation/test_benchmark_matrix.cpp)
    add_cfdx_test(test_level_b_reference_benchmarks tests/validation/test_level_b_reference_benchmarks.cpp)
    add_cfdx_test(test_level_c_coupled_verification tests/validation/test_level_c_coupled_verification.cpp)
    set_tests_properties(test_level_c_coupled_verification PROPERTIES LABELS "validation;coupled-ns")

    # Canonical engineering qualification cases. These labels make the
    # qualification subset independently selectable from ordinary unit tests.
    add_cfdx_test(test_cht_validation tests/validation/test_cht_validation.cpp)
    add_cfdx_test(test_thermal_vv tests/validation/test_thermal_vv.cpp)
    add_cfdx_test(test_thermal_zones_vv tests/validation/test_thermal_zones_vv.cpp)
    add_cfdx_test(test_production_energy_stage_vv tests/validation/test_production_energy_stage_vv.cpp)
    add_cfdx_test(test_production_s2s_energy_stage_vv tests/validation/test_production_s2s_energy_stage_vv.cpp)
    add_cfdx_test(test_production_s2s_energy_outer_loop_vv tests/validation/test_production_s2s_energy_outer_loop_vv.cpp)
    add_cfdx_test(test_production_buoyancy_outer_loop_vv tests/validation/test_production_buoyancy_outer_loop_vv.cpp)
    add_cfdx_test(test_radiation_vv tests/validation/test_radiation_vv.cpp)
    add_cfdx_test(test_s2s_radiation_vv tests/validation/test_s2s_radiation_vv.cpp)
    add_cfdx_test(test_thermal_radiation_model_matrix tests/validation/test_thermal_radiation_model_matrix.cpp)
    add_cfdx_test(test_thermal_radiation_regression tests/validation/test_thermal_radiation_regression.cpp)
    add_cfdx_test(test_thermal_challenge_vv tests/validation/test_thermal_challenge_vv.cpp)
    add_cfdx_test(test_thermophysical_models_vv tests/validation/test_thermophysical_models_vv.cpp)
    set_tests_properties(
        test_ghia_cavity_quick
        test_couette_quick
        test_poiseuille_quick
        test_mms_scalar_diffusion
        test_cht_validation
        test_thermal_vv
        test_radiation_vv
        test_s2s_radiation_vv
        PROPERTIES LABELS "qualification;validation")
    set_tests_properties(test_couette_quick test_ghia_cavity_quick
        PROPERTIES LABELS "qualification;validation;sanitizer-slow")

    add_test(NAME test_turbulence_qualification_matrix
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/validate_turbulence_qualification_matrix.py
                ${CMAKE_SOURCE_DIR}/docs/validation/TURBULENCE_QUALIFICATION_MATRIX.json)
    set_tests_properties(test_turbulence_qualification_matrix PROPERTIES
        LABELS "qualification;validation;turbulence")

    # Engineering qualification registry contract. This is intentionally a
    # metadata gate: numerical qualification remains the responsibility of
    # the executable case and its retained calculation evidence.
    if(Python3_Interpreter_FOUND)
        add_test(NAME test_cfdx_qualification_registry
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/validate_qualification_registry.py
                    ${CMAKE_SOURCE_DIR}/docs/validation/CFDX_QUALIFICATION_REGISTRY.json)
        set_tests_properties(test_cfdx_qualification_registry PROPERTIES
            LABELS "qualification;validation")
    endif()

    # Phase 13 verification labels. The campaign driver selects these tests
    # without conflating unit/regression coverage with numerical verification.
    set_tests_properties(test_backend_equivalence test_n13_restart_equivalence
        PROPERTIES LABELS "phase13;validation;n13")
    set_tests_properties(test_mms_scalar_diffusion
        test_m1_m4_validation test_m2_m4_solver_validation test_cht_validation
        test_thermal_vv test_radiation_vv test_s2s_radiation_vv test_thermal_radiation_model_matrix
        test_thermal_challenge_vv
        test_phase10_turbulence_hardening test_thermophysical_models_vv
        test_analytical_benchmarks test_benchmark_matrix
        test_level_b_reference_benchmarks test_level_c_coupled_verification
        test_steady_incompressible_solver test_matrix_free_fv_operator
        test_temporal_physics test_core_temporal
        PROPERTIES LABELS "phase13;validation")
    add_cfdx_test(test_schur_infrastructure tests/unit/test_schur_infrastructure.cpp)
    add_cfdx_test(test_exact_schur tests/unit/test_exact_schur.cpp)
    add_cfdx_test(test_n9_schur_oracle tests/unit/test_n9_schur_oracle.cpp)
    add_cfdx_test(test_simplerc_schur tests/unit/test_simplerc_schur.cpp)
    add_cfdx_test(test_lsc_bfbt_schur tests/unit/test_lsc_bfbt_schur.cpp)
    add_cfdx_test(test_pcd_schur tests/unit/test_pcd_schur.cpp)
    add_cfdx_test(test_coupled_schur_selection tests/unit/test_coupled_schur_selection.cpp)
    add_cfdx_test(test_lsc_bfbt_schur_null_space tests/unit/test_lsc_bfbt_schur_null_space.cpp)
    add_cfdx_test(test_schur_approximation_comparison tests/unit/test_schur_approximation_comparison.cpp)
    add_cfdx_test(test_schur_quantitative_qualification tests/unit/test_schur_quantitative_qualification.cpp)
    add_cfdx_test(test_mgr_preconditioner tests/unit/test_mgr_preconditioner.cpp)
    add_cfdx_test(test_wall_distance_turbulence tests/unit/test_wall_distance_turbulence.cpp)
    add_cfdx_test(test_geometry_cache tests/unit/test_geometry_cache.cpp)
    add_cfdx_test(test_coupled_block_schur_amg tests/validation/test_coupled_block_schur_amg.cpp)
    add_cfdx_test(test_n8_schur_production_benchmark tests/validation/test_n8_schur_production_benchmark.cpp)
    add_cfdx_test(test_n8_schur_representative_cfd tests/validation/test_n8_schur_representative_cfd.cpp)
    add_cfdx_test(test_n8_schur_complete_krylov tests/validation/test_n8_schur_complete_krylov.cpp)
    add_cfdx_test(test_n8_amg_anisotropy_scaling tests/validation/test_n8_amg_anisotropy_scaling.cpp)
    add_cfdx_test(test_n9_common_schur_qualification tests/validation/test_n9_common_schur_qualification.cpp)
    add_cfdx_test(test_n9_s10_pressure_reference_policies tests/validation/test_n9_s10_pressure_reference_policies.cpp)
    add_cfdx_test(test_n9_physical_matrix tests/validation/test_n9_physical_matrix.cpp)
    set_tests_properties(test_n9_physical_matrix PROPERTIES
        LABELS "n9;qualification;validation;long;validation-total"
        TIMEOUT 1800)
    set_tests_properties(test_n9_s10_pressure_reference_policies PROPERTIES
        LABELS "n9;qualification;validation;long;validation-total"
        TIMEOUT 300
        RUN_SERIAL TRUE)
    add_cfdx_test(test_n9_s9_rhie_chow_consistency tests/validation/test_n9_s9_rhie_chow_consistency.cpp)
    set_tests_properties(test_n9_s9_rhie_chow_consistency PROPERTIES
        LABELS "n9;qualification;validation;long;validation-total"
        TIMEOUT 300
        RUN_SERIAL TRUE)
    add_cfdx_test(test_n9_negative_diagnostics tests/validation/test_n9_negative_diagnostics.cpp)
    set_tests_properties(test_n9_negative_diagnostics PROPERTIES
        LABELS "n9;qualification;validation;validation-total"
        TIMEOUT 600
        RUN_SERIAL TRUE)
    if(CFDX_ENABLE_MPI AND MPI_FOUND)
        add_executable(test_mpi_partition tests/unit/test_mpi_partition.cpp)
        add_executable(test_n13_mpi_equivalence tests/validation/test_n13_mpi_equivalence.cpp)
        add_executable(test_mpi_deterministic tests/unit/test_mpi_deterministic.cpp)
        add_executable(test_mpi_fused_reduction tests/unit/test_mpi_fused_reduction.cpp)
        add_executable(test_mpi_krylov_reduction tests/unit/test_mpi_krylov_reduction.cpp)
        add_executable(test_phase8_mpi_poisson tests/validation/test_phase8_mpi_poisson.cpp)
        add_executable(test_n8_mpi_boundary tests/validation/test_n8_mpi_boundary.cpp)
        if(CFDX_ENABLE_PARALLEL_HDF5 AND HDF5_FOUND)
            add_executable(cfdx_mpi_poisson_restart apps/mpi_poisson_restart/main.cpp)
            target_link_libraries(cfdx_mpi_poisson_restart PRIVATE cfdx_core MPI::MPI_CXX)
            add_executable(test_phase5_distributed tests/unit/test_phase5_distributed.cpp)
            add_executable(test_phase5_restart_writer tests/unit/test_phase5_restart_writer.cpp)
            add_executable(test_phase5_restart_reader tests/unit/test_phase5_restart_reader.cpp)
            target_link_libraries(test_phase5_distributed PRIVATE cfdx_core MPI::MPI_CXX)
            target_link_libraries(test_phase5_restart_writer PRIVATE cfdx_core MPI::MPI_CXX)
            target_link_libraries(test_phase5_restart_reader PRIVATE cfdx_core MPI::MPI_CXX)
            target_include_directories(test_phase5_distributed PRIVATE ${CMAKE_SOURCE_DIR}/tests)
            target_include_directories(test_phase5_restart_writer PRIVATE ${CMAKE_SOURCE_DIR}/tests)
            target_include_directories(test_phase5_restart_reader PRIVATE ${CMAKE_SOURCE_DIR}/tests)
        endif()
        target_link_libraries(test_mpi_deterministic PRIVATE cfdx_core MPI::MPI_CXX)
        target_include_directories(test_mpi_deterministic PRIVATE ${CMAKE_SOURCE_DIR}/tests)
        target_link_libraries(test_mpi_partition PRIVATE cfdx_core MPI::MPI_CXX)
        target_link_libraries(test_n13_mpi_equivalence PRIVATE cfdx_core MPI::MPI_CXX)
        target_link_libraries(test_mpi_fused_reduction PRIVATE cfdx_core MPI::MPI_CXX)
        target_link_libraries(test_mpi_krylov_reduction PRIVATE cfdx_core MPI::MPI_CXX)
        target_link_libraries(test_phase8_mpi_poisson PRIVATE cfdx_core MPI::MPI_CXX)
        target_link_libraries(test_n8_mpi_boundary PRIVATE cfdx_core MPI::MPI_CXX)
        target_include_directories(test_mpi_partition PRIVATE ${CMAKE_SOURCE_DIR}/tests)
        target_include_directories(test_n13_mpi_equivalence PRIVATE ${CMAKE_SOURCE_DIR}/tests)
        target_include_directories(test_mpi_fused_reduction PRIVATE ${CMAKE_SOURCE_DIR}/tests)
        target_include_directories(test_mpi_krylov_reduction PRIVATE ${CMAKE_SOURCE_DIR}/tests)
        set(CFDX_MPIEXEC_NUMPROCS_FLAG "${MPIEXEC_NUMPROCS_FLAG}")
        if(CFDX_MPIEXEC_NUMPROCS_FLAG STREQUAL "")
            set(CFDX_MPIEXEC_NUMPROCS_FLAG "-n")
        endif()
        # Only Open MPI accepts --oversubscribe. MPICH/Hydra rejects it
        # before starting the restart reader. Preserve three-rank readers on
        # two-core Open MPI CI workers without sending a vendor flag to all MPI.
        set(CFDX_MPIEXEC_OVERSUBSCRIBE_FLAG "")
        execute_process(COMMAND "${MPIEXEC_EXECUTABLE}" --version
            OUTPUT_VARIABLE CFDX_MPIEXEC_VERSION_OUTPUT
            ERROR_VARIABLE CFDX_MPIEXEC_VERSION_ERROR
            RESULT_VARIABLE CFDX_MPIEXEC_VERSION_RESULT)
        if(CFDX_MPIEXEC_VERSION_RESULT EQUAL 0 AND
           "${CFDX_MPIEXEC_VERSION_OUTPUT}${CFDX_MPIEXEC_VERSION_ERROR}" MATCHES "Open MPI|OpenRTE")
            set(CFDX_MPIEXEC_OVERSUBSCRIBE_FLAG "--oversubscribe")
        endif()
        add_test(NAME test_mpi_deterministic
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_deterministic> ${MPIEXEC_POSTFLAGS})
        add_test(NAME test_n13_mpi_equivalence
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_n13_mpi_equivalence> ${MPIEXEC_POSTFLAGS})
        set_tests_properties(test_n13_mpi_equivalence PROPERTIES
            LABELS "phase13;validation;n13"
            ENVIRONMENT "LSAN_OPTIONS=suppressions=${CMAKE_SOURCE_DIR}/tests/sanitizers/lsan_mpi.supp")
        add_test(NAME test_mpi_partition
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_partition> ${MPIEXEC_POSTFLAGS})
        add_test(NAME test_mpi_fused_reduction
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_fused_reduction> ${MPIEXEC_POSTFLAGS})
        add_test(NAME test_mpi_krylov_reduction
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_krylov_reduction> ${MPIEXEC_POSTFLAGS})
        add_test(NAME test_phase8_mpi_poisson
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_phase8_mpi_poisson> ${MPIEXEC_POSTFLAGS})
        add_test(NAME test_n8_mpi_boundary
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_n8_mpi_boundary> ${MPIEXEC_POSTFLAGS})
        set_tests_properties(test_n8_mpi_boundary PROPERTIES
            LABELS "n8;qualification;validation;mpi")
        set_tests_properties(test_phase8_mpi_poisson PROPERTIES LABELS "phase13;validation")
        if(CFDX_ENABLE_PARALLEL_HDF5 AND HDF5_FOUND)
            add_test(NAME test_phase5_distributed
                     COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                             ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_phase5_distributed> ${MPIEXEC_POSTFLAGS})
            add_test(NAME test_production_mpi_poisson_writer
                     COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                             ${MPIEXEC_PREFLAGS} $<TARGET_FILE:cfdx_mpi_poisson_restart>
                             --mode writer
                             --checkpoint ${CMAKE_BINARY_DIR}/production_poisson_checkpoint.h5
                             --reference ${CMAKE_BINARY_DIR}/production_poisson_reference.h5
                             ${MPIEXEC_POSTFLAGS})
            add_test(NAME test_production_mpi_poisson_reader
                     COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_OVERSUBSCRIBE_FLAG} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 3
                             ${MPIEXEC_PREFLAGS} $<TARGET_FILE:cfdx_mpi_poisson_restart>
                             --mode reader
                             --checkpoint ${CMAKE_BINARY_DIR}/production_poisson_checkpoint.h5
                             --reference ${CMAKE_BINARY_DIR}/production_poisson_reference.h5
                             ${MPIEXEC_POSTFLAGS})
            set_tests_properties(test_production_mpi_poisson_reader PROPERTIES
                DEPENDS test_production_mpi_poisson_writer)
            add_test(NAME test_phase5_restart_writer
                     COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                             ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_phase5_restart_writer> ${MPIEXEC_POSTFLAGS})
            add_test(NAME test_phase5_restart_reader
                     COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_OVERSUBSCRIBE_FLAG} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 3
                             ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_phase5_restart_reader> ${MPIEXEC_POSTFLAGS})
            set_tests_properties(test_phase5_restart_reader PROPERTIES DEPENDS test_phase5_restart_writer)
        endif()
        if(CMAKE_BUILD_TYPE STREQUAL "DebugSanitizers")
            # OpenMPI/hwloc allocates process-global state during MPI_Init which
            # is reported as a leak after MPI_Finalize by LeakSanitizer. Keep
            # ASan/UBSan active for the MPI regression, but exclude those
            # external-library teardown allocations from the leak check.
            # N13 has the same external OpenMPI/hwloc teardown behavior as the
            # other MPI sanitizer tests. Its numerical equivalence and reduction
            # assertions remain active; only LeakSanitizer's process-global
            # external-library teardown check is disabled for these MPI tests.
            set_tests_properties(test_n13_mpi_equivalence test_mpi_deterministic
                test_mpi_partition test_mpi_fused_reduction test_mpi_krylov_reduction
                test_phase8_mpi_poisson test_n8_mpi_boundary PROPERTIES
                ENVIRONMENT "LSAN_OPTIONS=detect_leaks=0")
        endif()
    endif()
    add_cfdx_test(test_convergence_history tests/unit/test_convergence_history.cpp)
    add_cfdx_test(test_execution_summary tests/unit/test_execution_summary.cpp)
    add_cfdx_test(test_restart_mapping tests/unit/test_restart_mapping.cpp)
    add_cfdx_test(test_transport_models tests/unit/test_transport_models.cpp)
    add_cfdx_test(test_equation_of_state_consistency tests/unit/test_equation_of_state_consistency.cpp)
    add_cfdx_test(test_sst_blending tests/unit/test_sst_blending.cpp)
    add_cfdx_test(test_thermo_transport_contracts tests/unit/test_thermo_transport_contracts.cpp)
    add_cfdx_test(test_thermophysical_models tests/unit/test_thermophysical_models.cpp)
    add_cfdx_test(test_m5_m7_foundations tests/unit/test_m5_m7_foundations.cpp)
    if(CFDX_BUILD_RUNTIME_TESTS)
        add_cfdx_test(test_runtime_execution_ooc tests/unit/test_runtime_execution_ooc.cpp)
        add_cfdx_test(test_runtime_pipeline tests/unit/test_runtime_pipeline.cpp)
        add_cfdx_test(test_phase6_runtime_acceptance tests/validation/test_phase6_runtime_acceptance.cpp)
        add_cfdx_test(test_ooc_pipeline tests/unit/test_ooc_pipeline.cpp)
        add_cfdx_test(test_phase7_ooc_validation tests/validation/test_phase7_ooc_validation.cpp)
    endif()
    if(CFDX_ENABLE_GPU)
        add_cfdx_test(test_cuda_runtime tests/unit/test_cuda_runtime.cpp)
        add_cfdx_test(test_phase8_cuda_poisson tests/validation/test_phase8_cuda_poisson.cpp)
        add_executable(test_n13_cuda_equivalence tests/validation/test_n13_cuda_equivalence.cpp)
        target_link_libraries(test_n13_cuda_equivalence PRIVATE cfdx_core)
        add_test(NAME test_n13_cuda_equivalence COMMAND test_n13_cuda_equivalence)
        set_tests_properties(test_n13_cuda_equivalence PROPERTIES LABELS "phase13;validation;n13;hardware")
    endif()

    # HDF5 tests require the native reader/writer implementation.
    if(CFDX_HAS_HDF5)
        add_cfdx_test(test_hdf5_writer tests/unit/test_hdf5_writer.cpp)
        add_cfdx_test(test_hdf5_roundtrip tests/unit/test_hdf5_roundtrip.cpp)
        add_cfdx_test(test_field_roundtrip tests/unit/test_field_roundtrip.cpp)
        add_cfdx_test(test_case_hdf5_io tests/unit/test_case_hdf5_io.cpp)
    else()
        add_cfdx_test(test_hdf5_unavailable tests/unit/test_hdf5_unavailable.cpp)
    endif()

    # Memory Planner tests
    add_cfdx_test(test_memory_planner_validation tests/unit/test_memory_planner_validation.cpp)

    # Plan completion tests
    add_cfdx_test(test_plan_complet tests/unit/test_plan_complet.cpp)

    # Module validation tests (M0.8-T06 to T10, M0.9-T05, M0.10-T01)
    add_cfdx_test(test_validation_simple tests/unit/test_validation_simple.cpp)
    add_cfdx_test(test_openfoam_importer tests/unit/test_openfoam_importer.cpp)
    if(CFDX_HAS_HDF5)
        add_cfdx_test(test_cross_import_equivalence tests/unit/test_cross_import_equivalence.cpp)
    endif()
    if(Python3_Interpreter_FOUND)
        add_test(NAME test_cpp_harness_execution_contract
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_SOURCE_DIR}/scripts/check_cpp_test_execution.py
                    --root ${CMAKE_SOURCE_DIR} --self-test)
        set_tests_properties(test_cpp_harness_execution_contract PROPERTIES
            LABELS "audit;test-infrastructure"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
        add_test(NAME test_n9_s8_naca0012_contract
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_SOURCE_DIR}/tests/python/test_n9_s8_naca0012_contract.py)
        set_tests_properties(test_n9_s8_naca0012_contract PROPERTIES
            LABELS "n9;qualification;validation;total-validation"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    endif()
    if(Python3_Interpreter_FOUND AND CFDX_HAS_HDF5)
        add_test(NAME test_n9_coupled_schur_case_selection
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_SOURCE_DIR}/tests/python/test_n9_coupled_schur_case_selection.py)
        set_tests_properties(test_n9_coupled_schur_case_selection PROPERTIES
            LABELS "n9;schur;integration;validation;validation-total"
            ENVIRONMENT "CFDX_PRODUCTION_SOLVER=\$<TARGET_FILE:cfdx_production_solver>;PYTHONPATH=${CMAKE_SOURCE_DIR}/src/cfdx/python"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
        add_test(NAME test_n9_poiseuille_body_force_production
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_SOURCE_DIR}/tests/python/test_n9_poiseuille_body_force_production.py)
        set_tests_properties(test_n9_poiseuille_body_force_production PROPERTIES
            LABELS "n9;qualification;integration;validation;validation-total"
            ENVIRONMENT "CFDX_PRODUCTION_SOLVER=\$<TARGET_FILE:cfdx_production_solver>;PYTHONPATH=${CMAKE_SOURCE_DIR}/src/cfdx/python"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            TIMEOUT 900)
    endif()

    if(Python3_Interpreter_FOUND)
        add_test(NAME test_missing_qualification_case_contracts COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/test_missing_qualification_case_contracts.py)
        set_tests_properties(test_missing_qualification_case_contracts PROPERTIES LABELS "qualification;validation")
    endif()

    if(Python3_Interpreter_FOUND)
        if(CFDX_HAS_HDF5)
            # End-to-end OpenFOAM result import through the production executable.
            # This is intentionally bounded: each subprocess has a 30 s timeout.
            add_test(NAME test_openfoam_results_cli_production_e2e
                     COMMAND ${Python3_EXECUTABLE} -m pytest -q
                             ${CMAKE_SOURCE_DIR}/tests/python/test_openfoam_results_cli.py)
            set_tests_properties(test_openfoam_results_cli_production_e2e PROPERTIES
                LABELS "openfoam;integration;qualification"
                ENVIRONMENT "CFDX_PRODUCTION_SOLVER=$<TARGET_FILE:cfdx_production_solver>;PYTHONPATH=${CMAKE_SOURCE_DIR}/src/cfdx/python"
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                TIMEOUT 120)
            add_test(NAME test_production_solver_e2e
                     COMMAND ${Python3_EXECUTABLE}
                             ${CMAKE_SOURCE_DIR}/tests/python/test_production_solver_e2e.py)
            set_tests_properties(test_production_solver_e2e PROPERTIES
                ENVIRONMENT "CFDX_PRODUCTION_SOLVER=$<TARGET_FILE:cfdx_production_solver>;CFDX_PRODUCTION_MESH=${CMAKE_SOURCE_DIR}/tests/data/cavity_with_patches.msh;PYTHONPATH=${CMAKE_SOURCE_DIR}/src/cfdx/python")
            add_test(NAME test_production_solver_numerics_e2e
                     COMMAND ${Python3_EXECUTABLE}
                             ${CMAKE_SOURCE_DIR}/tests/python/test_production_solver_numerics_e2e.py)
            set_tests_properties(test_production_solver_numerics_e2e PROPERTIES
                ENVIRONMENT "CFDX_PRODUCTION_SOLVER=$<TARGET_FILE:cfdx_production_solver>;PYTHONPATH=${CMAKE_SOURCE_DIR}/src/cfdx/python")
        endif()
        add_test(NAME test_meshio_bridge
                 COMMAND ${Python3_EXECUTABLE}
                         ${CMAKE_SOURCE_DIR}/tests/mesh/test_meshio_bridge.py
                         ${CMAKE_SOURCE_DIR}/scripts/meshio_import.py)
        add_test(NAME test_mesh_geometry_matrix
                 COMMAND ${Python3_EXECUTABLE}
                         ${CMAKE_SOURCE_DIR}/tests/mesh/test_mesh_geometry_matrix.py
                         ${CMAKE_SOURCE_DIR}/scripts/meshio_import.py)
    endif()

    # I/O adapter tests (M0.10-T03)
    # C++ solver adapters are now thin stubs (format parsing moved to Python).
    add_cfdx_test(test_gap_analysis tests/unit/test_gap_analysis.cpp)

    # Thermal production-path qualification label (T1 / issue #880).
    set_tests_properties(test_thermal_zones_vv PROPERTIES LABELS "qualification;validation;phase13")
    set_tests_properties(test_production_energy_stage_vv PROPERTIES LABELS "qualification;validation;phase13")
    set_tests_properties(test_production_s2s_energy_stage_vv PROPERTIES LABELS "qualification;validation;phase13")
    set_tests_properties(test_production_s2s_energy_outer_loop_vv PROPERTIES LABELS "qualification;validation;phase13")
    set_tests_properties(test_production_buoyancy_outer_loop_vv PROPERTIES LABELS "qualification;validation;phase13")

endif()


```

### `apps/verification/naca0012_wall_distance_comparison.cpp`

SHA256 `1a1aa5c83fcc21869eac050374973fa931110015960de837b981ee2e8ad1a80e`

```cpp
#include "cfdx/physics/wall_distance.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace cfdx::physics;

namespace {
struct Segment { std::size_t a; std::size_t b; };

std::string trim(const std::string& value) {
    const auto first=value.find_first_not_of(" \t\r\n");
    if(first==std::string::npos) return {};
    const auto last=value.find_last_not_of(" \t\r\n");
    return value.substr(first,last-first+1);
}

std::pair<std::size_t,std::size_t> section(const std::vector<std::string>& lines,
                                           const std::string& key) {
    const std::string prefix=key+"=";
    for(std::size_t i=0;i<lines.size();++i) {
        const std::string line=trim(lines[i]);
        if(line.rfind(prefix,0)!=0) continue;
        try {
            return {i,static_cast<std::size_t>(std::stoull(trim(line.substr(prefix.size()))))};
        } catch(const std::exception&) {
            throw std::runtime_error("invalid SU2 count for "+key);
        }
    }
    throw std::runtime_error("missing SU2 section "+key);
}

std::size_t next_data_line(const std::vector<std::string>& lines,std::size_t index) {
    while(index<lines.size()) {
        const std::string line=trim(lines[index]);
        if(!line.empty() && line.front()!='%') break;
        ++index;
    }
    if(index>=lines.size()) throw std::runtime_error("unexpected end of SU2 mesh");
    return index;
}

struct NacaMesh {
    std::vector<WallDistanceVec3> points;
    std::vector<Segment> airfoil_edges;
    double xmin{std::numeric_limits<double>::infinity()};
    double xmax{-std::numeric_limits<double>::infinity()};
    double ymin{std::numeric_limits<double>::infinity()};
    double ymax{-std::numeric_limits<double>::infinity()};
};

NacaMesh read_su2_airfoil_mesh(const std::string& path) {
    std::ifstream input(path);
    if(!input) throw std::runtime_error("cannot open SU2 mesh: "+path);
    std::vector<std::string> lines;
    for(std::string line;std::getline(input,line);) lines.push_back(line);

    const auto [point_header,npoints]=section(lines,"NPOIN");
    NacaMesh mesh;
    mesh.points.reserve(npoints);
    std::size_t cursor=point_header+1;
    for(std::size_t i=0;i<npoints;++i) {
        cursor=next_data_line(lines,cursor);
        std::istringstream row(lines[cursor++]);
        double x=0.0,y=0.0;
        if(!(row>>x>>y)) throw std::runtime_error("invalid SU2 point record");
        mesh.points.push_back({x,y,0.0});
    }

    bool found_airfoil=false;
    for(std::size_t i=0;i<lines.size();++i) {
        const std::string line=trim(lines[i]);
        const std::string tag_prefix="MARKER_TAG=";
        if(line.rfind(tag_prefix,0)!=0 || trim(line.substr(tag_prefix.size()))!="airfoil") continue;
        found_airfoil=true;
        std::size_t elems_header=i+1;
        while(elems_header<lines.size() &&
              trim(lines[elems_header]).rfind("MARKER_ELEMS=",0)!=0) ++elems_header;
        if(elems_header>=lines.size())
            throw std::runtime_error("airfoil marker has no MARKER_ELEMS section");
        std::size_t nelems=0;
        try {
            nelems=static_cast<std::size_t>(std::stoull(
                trim(trim(lines[elems_header]).substr(std::string("MARKER_ELEMS=").size()))));
        } catch(const std::exception&) {
            throw std::runtime_error("invalid airfoil MARKER_ELEMS count");
        }
        cursor=elems_header+1;
        for(std::size_t j=0;j<nelems;++j) {
            cursor=next_data_line(lines,cursor);
            std::istringstream row(lines[cursor++]);
            int element_type=0;
            long long a=-1,b=-1;
            if(!(row>>element_type>>a>>b)) throw std::runtime_error("invalid airfoil marker element");
            if(element_type!=3) throw std::runtime_error("expected SU2 line element type 3 on airfoil");
            if(a<0 || b<0 || static_cast<std::size_t>(a)>=mesh.points.size() ||
               static_cast<std::size_t>(b)>=mesh.points.size())
                throw std::runtime_error("airfoil marker node index out of range");
            mesh.airfoil_edges.push_back({static_cast<std::size_t>(a),static_cast<std::size_t>(b)});
        }
        break;
    }
    if(!found_airfoil || mesh.airfoil_edges.empty())
        throw std::runtime_error("SU2 mesh has no non-empty airfoil marker");

    for(const auto& edge:mesh.airfoil_edges) {
        const auto& a=mesh.points[edge.a];
        const auto& b=mesh.points[edge.b];
        mesh.xmin=std::min({mesh.xmin,a.x,b.x});
        mesh.xmax=std::max({mesh.xmax,a.x,b.x});
        mesh.ymin=std::min({mesh.ymin,a.y,b.y});
        mesh.ymax=std::max({mesh.ymax,a.y,b.y});
    }
    if(!(mesh.xmax>mesh.xmin) || !(mesh.ymax>mesh.ymin))
        throw std::runtime_error("airfoil marker has degenerate geometric bounds");
    return mesh;
}

// Extrude each boundary segment to a rectangular panel spanning z. For points
// at z=0, point-to-panel distance equals point-to-segment distance in x-y.
WallSurface make_extruded_wall_surface(const NacaMesh& mesh,double half_depth) {
    WallSurface surface;
    surface.points.reserve(mesh.airfoil_edges.size()*4);
    surface.triangles.reserve(mesh.airfoil_edges.size()*2);
    for(const auto& edge:mesh.airfoil_edges) {
        const auto& a=mesh.points[edge.a];
        const auto& b=mesh.points[edge.b];
        const std::size_t base=surface.points.size();
        surface.points.push_back({a.x,a.y,-half_depth});
        surface.points.push_back({b.x,b.y,-half_depth});
        surface.points.push_back({b.x,b.y, half_depth});
        surface.points.push_back({a.x,a.y, half_depth});
        surface.triangles.push_back({{base,base+1,base+2}});
        surface.triangles.push_back({{base,base+2,base+3}});
    }
    return surface;
}

bool inside_airfoil(const NacaMesh& mesh,const WallDistanceVec3& p) {
    bool inside=false;
    for(const auto& edge:mesh.airfoil_edges) {
        const auto& a=mesh.points[edge.a];
        const auto& b=mesh.points[edge.b];
        if((a.y>p.y)!=(b.y>p.y)) {
            const double x_cross=a.x+(p.y-a.y)*(b.x-a.x)/(b.y-a.y);
            if(p.x<x_cross) inside=!inside;
        }
    }
    return inside;
}

std::vector<double> independent_segment_reference(const NacaMesh& mesh,
                                                  const WallDistanceGrid& grid) {
    // Independent 2-D point-to-segment oracle. Deliberately does not call
    // point_triangle_distance2(), exact_point_distance(), or the wall BVH.
    std::vector<double> reference(grid.points.size(),std::numeric_limits<double>::infinity());
    for(std::size_t i=0;i<grid.points.size();++i) {
        if(grid.solid[i]) continue;
        const auto& p=grid.points[i];
        double best2=std::numeric_limits<double>::infinity();
        for(const auto& edge:mesh.airfoil_edges) {
            const auto& a=mesh.points[edge.a];
            const auto& b=mesh.points[edge.b];
            const double dx=b.x-a.x,dy=b.y-a.y;
            const double length2=dx*dx+dy*dy;
            if(length2<=0.0) throw std::runtime_error("degenerate airfoil marker segment");
            const double t=std::clamp(((p.x-a.x)*dx+(p.y-a.y)*dy)/length2,0.0,1.0);
            const double rx=p.x-(a.x+t*dx),ry=p.y-(a.y+t*dy);
            best2=std::min(best2,rx*rx+ry*ry);
        }
        reference[i]=std::sqrt(best2);
    }
    return reference;
}

struct ResultRow {
    std::string method;
    double l2{0.0},linf{0.0},near_l2{0.0},residual{0.0};
    std::size_t iterations{0},valid_samples{0},invalid_samples{0};
    double elapsed_ms{0.0};
    bool converged{false};
    std::string stopping_reason;
};
} // namespace

int main(int argc,char** argv) {
    try {
        if(argc<2 || argc>4) {
            std::cerr << "Usage: cfdx_naca0012_wall_distance_comparison <mesh.su2> [output.csv] [nx>=9]\n";
            return 2;
        }
        const std::string mesh_path=argv[1];
        const std::string output_path=argc>=3 ? argv[2] : "naca0012_wall_distance_comparison.csv";
        const std::size_t nx=argc>=4 ? static_cast<std::size_t>(std::stoull(argv[3])) : 61;
        if(nx<9) throw std::invalid_argument("nx must be >= 9");

        const NacaMesh mesh=read_su2_airfoil_mesh(mesh_path);
        const double chord=mesh.xmax-mesh.xmin;
        const double padding=0.20*chord;
        const double xlo=mesh.xmin-padding, xhi=mesh.xmax+padding;
        const double ylo=mesh.ymin-padding, yhi=mesh.ymax+padding;
        const double xspan=xhi-xlo, yspan=yhi-ylo;
        const auto ny=std::max<std::size_t>(
            9,static_cast<std::size_t>(std::lround(
                static_cast<double>(nx-1)*yspan/xspan))+1);
        const double hx=xspan/static_cast<double>(nx-1);
        const double hy=yspan/static_cast<double>(ny-1);
        const double depth=2.0*chord;
        const WallSurface surface=make_extruded_wall_surface(mesh,2.0*depth);
        const WallDistanceGrid grid=make_wall_distance_grid(
            nx,ny,3,{xlo,ylo,-depth},{hx,hy,depth},
            [&mesh](const WallDistanceVec3& p){ return inside_airfoil(mesh,p); });
        const auto reference=independent_segment_reference(mesh,grid);
        const std::array<WallDistanceMethod,9> methods={{
            WallDistanceMethod::EXACT_GEOMETRIC,WallDistanceMethod::SEARCH_BASED,
            WallDistanceMethod::MESH_WAVE,WallDistanceMethod::DIRECTIONAL_MESH_WAVE,
            WallDistanceMethod::POISSON,WallDistanceMethod::EIKONAL,
            WallDistanceMethod::HAMILTON_JACOBI,WallDistanceMethod::ADVECTION_DIFFUSION,
            WallDistanceMethod::HYBRID_POISSON_EIKONAL}};

        std::cout << std::setprecision(10) << "mesh=" << mesh_path
                  << " airfoil_segments=" << mesh.airfoil_edges.size()
                  << " samples=" << grid.points.size() << " nx=" << nx
                  << " ny=" << ny << " nz=3 chord=" << chord
                  << " hx=" << hx << " hy=" << hy << "\n";
        std::ofstream csv(output_path);
        if(!csv) throw std::runtime_error("cannot open CSV output: "+output_path);
        csv << "method,l2_relative,linf_relative,near_wall_l2_relative,residual_inf,iterations,"
               "valid_samples,invalid_samples,elapsed_ms,converged,stopping_reason\n";
        csv << std::setprecision(17);

        std::vector<ResultRow> rows;
        rows.reserve(methods.size());
        const std::size_t effective_n=std::max(nx,ny);
        for(const auto method:methods) {
            const std::size_t budget=wall_distance_iteration_budget(method,effective_n);
            const auto t0=std::chrono::steady_clock::now();
            const auto result=compute_wall_distance(method,surface,grid,budget);
            const auto t1=std::chrono::steady_clock::now();
            const double elapsed=std::chrono::duration<double,std::milli>(t1-t0).count();
            const auto metrics=compare_wall_distance(grid,reference,result.distance,
                                                     2.0*std::min(hx,hy));
            std::size_t valid=0,invalid=0;
            for(std::size_t i=0;i<grid.points.size();++i) {
                if(grid.solid[i]) continue;
                if(i>=result.valid.size() || !result.valid[i] ||
                   !std::isfinite(result.distance[i])) ++invalid;
                else ++valid;
            }
            ResultRow row{result.method,metrics.l2_relative,metrics.linf_relative,
                metrics.near_wall_l2_relative,result.residual_inf,result.iterations,
                valid,invalid,elapsed,result.converged,result.stopping_reason};
            if((method==WallDistanceMethod::EXACT_GEOMETRIC ||
                method==WallDistanceMethod::SEARCH_BASED) &&
               (row.l2>1e-12 || row.linf>1e-12 || row.invalid_samples!=0))
                throw std::runtime_error("geometric oracle self-consistency gate failed for "+row.method);
            if(row.invalid_samples!=0)
                std::cerr << "WARNING: " << row.method << " has " << row.invalid_samples
                          << " invalid fluid samples\n";
            csv << row.method << ',' << row.l2 << ',' << row.linf << ',' << row.near_l2 << ','
                << row.residual << ',' << row.iterations << ',' << row.valid_samples << ','
                << row.invalid_samples << ',' << row.elapsed_ms << ','
                << (row.converged ? "true" : "false") << ',' << row.stopping_reason << '\n';
            rows.push_back(std::move(row));
        }
        csv.close();
        if(!csv) throw std::runtime_error("failed while writing CSV output: "+output_path);

        std::cout << "method,l2_relative,linf_relative,near_wall_l2_relative,iterations,residual_inf,"
                     "elapsed_ms,converged,stopping_reason\n";
        for(const auto& row:rows)
            std::cout << row.method << ',' << row.l2 << ',' << row.linf << ','
                      << row.near_l2 << ',' << row.iterations << ',' << row.residual << ','
                      << row.elapsed_ms << ',' << (row.converged ? "true" : "false") << ','
                      << row.stopping_reason << '\n';
        std::cout << "csv=" << output_path << '\n';
        for(const auto& row:rows)
            if(row.invalid_samples!=0)
                throw std::runtime_error("invalid fluid distance field for "+row.method);
        return 0;
    } catch(const std::exception& error) {
        std::cerr << "NACA0012 wall-distance comparison failed: " << error.what() << '\n';
        return 1;
    }
}

```

### `apps/verification/wall_distance_complex_benchmark.cpp`

SHA256 `5d9c4d97a55ba0062c88d5fb5fe4a935937a94d55905e9c45012aa782c1c282d`

```cpp
#include "cfdx/physics/wall_distance.h"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace cfdx::physics;

namespace {

void add_triangle(WallSurface& s,std::size_t a,std::size_t b,std::size_t c) {
    s.triangles.push_back({{a,b,c}});
}

void add_box(WallSurface& s,const WallDistanceVec3& lo,const WallDistanceVec3& hi) {
    const std::size_t b=s.points.size();
    s.points.insert(s.points.end(),{
        {lo.x,lo.y,lo.z},{hi.x,lo.y,lo.z},{hi.x,hi.y,lo.z},{lo.x,hi.y,lo.z},
        {lo.x,lo.y,hi.z},{hi.x,lo.y,hi.z},{hi.x,hi.y,hi.z},{lo.x,hi.y,hi.z}});
    const std::array<std::array<std::size_t,4>,6> q={{
        {{0,3,2,1}},{{4,5,6,7}},{{0,1,5,4}},{{1,2,6,5}},{{2,3,7,6}},{{3,0,4,7}}}};
    for(const auto& f:q) { add_triangle(s,b+f[0],b+f[1],b+f[2]); add_triangle(s,b+f[0],b+f[2],b+f[3]); }
}

void add_cylinder_x(WallSurface& s,double x0,double x1,double cy,double cz,double r,std::size_t n) {
    const std::size_t b=s.points.size();
    for(std::size_t i=0;i<n;++i) {
        const double a=2.0*3.14159265358979323846*static_cast<double>(i)/static_cast<double>(n);
        s.points.push_back({x0,cy+r*std::cos(a),cz+r*std::sin(a)});
        s.points.push_back({x1,cy+r*std::cos(a),cz+r*std::sin(a)});
    }
    const std::size_t c0=s.points.size(); s.points.push_back({x0,cy,cz});
    const std::size_t c1=s.points.size(); s.points.push_back({x1,cy,cz});
    for(std::size_t i=0;i<n;++i) {
        const std::size_t j=(i+1)%n;
        add_triangle(s,b+2*i,b+2*j,b+2*j+1); add_triangle(s,b+2*i,b+2*j+1,b+2*i+1);
        add_triangle(s,c0,b+2*j,b+2*i);
        add_triangle(s,c1,b+2*i+1,b+2*j+1);
    }
}

WallSurface make_complex_wing_body_tail() {
    WallSurface s;
    // Deliberately non-trivial multi-body-like surface: fuselage, swept-ish
    // wing slabs, tail planes and a vertical fin. The reference is the
    // triangulated surface distance, not an analytic approximation.
    add_cylinder_x(s,0.0,9.0,0.0,0.0,0.72,64);
    add_box(s,{3.0,-3.8,-0.10},{6.4,3.8,0.10});
    add_box(s,{7.0,-1.9,-0.08},{8.7,1.9,0.08});
    add_box(s,{7.2,-0.12,0.0},{8.8,0.12,1.45});
    add_box(s,{-0.15,-1.35,-0.55},{1.9,1.35,0.0});
    return s;
}

bool inside_complex(const WallDistanceVec3& p) {
    // Keep the computational solid exactly consistent with add_cylinder_x():
    // the surface is a finite cylinder x in [0,9] with circular radius 0.72.
    // The previous ellipsoidal test silently changed the solid topology near
    // both cylinder end caps, invalidating the complex-geometry benchmark.
    if(p.x>=0.0 && p.x<=9.0 &&
       p.y*p.y+p.z*p.z<=0.72*0.72) return true;
    if(p.x>=3.0 && p.x<=6.4 && std::abs(p.y)<=3.8 && std::abs(p.z)<=0.10) return true;
    if(p.x>=7.0 && p.x<=8.7 && std::abs(p.y)<=1.9 && std::abs(p.z)<=0.08) return true;
    if(p.x>=7.2 && p.x<=8.8 && std::abs(p.y)<=0.12 && p.z>=0.0 && p.z<=1.45) return true;
    if(p.x>=-0.15 && p.x<=1.9 && std::abs(p.y)<=1.35 && p.z>=-0.55 && p.z<=0.0) return true;
    return false;
}

struct Row {
    std::string method; double l2,linf,near_l2,violations,residual;
    std::size_t invalid,iterations; double ms,init_ms,poisson_ms;
    std::size_t poisson_iterations; double poisson_residual; bool converged;
    double min_distance,max_distance;
    std::size_t wall_ray_hits,wall_ray_misses,wall_fallbacks,wall_bad_alignment;
    double wall_min_alignment,poisson_phi_min,poisson_phi_max,poisson_grad_min,poisson_grad_max;
    double poisson_distance_l2_error,poisson_distance_linf_error;
};

} // namespace

int main(int argc,char** argv) {
    std::string output = argc>1 ? argv[1] : "wall_distance_complex_benchmark.csv";
    const WallSurface surface=make_complex_wing_body_tail();
    const WallDistanceGrid grid=make_wall_distance_grid(
        40,30,22,{-1.5,-5.0,-2.2},{0.30,0.34,0.21},inside_complex);
    const auto reference=exact_reference(surface,grid);
    const double h=std::min({grid.spacing.x,grid.spacing.y,grid.spacing.z});
    const WallDistanceBvh poisson_bvh(surface);
    const auto poisson_operator_audit=audit_poisson_operator(poisson_bvh,grid);
    const auto poisson_offset_audit=audit_poisson_wall_offsets(poisson_bvh,grid);
    const std::size_t max_id=poisson_offset_audit.max_wall_coefficient_cell;
    const std::size_t max_k=max_id/(grid.nx*grid.ny);
    const std::size_t max_rem=max_id%(grid.nx*grid.ny);
    const std::size_t max_j=max_rem/grid.nx;
    const std::size_t max_i=max_rem%grid.nx;
    const auto max_p=grid.points[max_id];

    const std::array<WallDistanceMethod,9> methods={{
        WallDistanceMethod::EXACT_GEOMETRIC,
        WallDistanceMethod::SEARCH_BASED,
        WallDistanceMethod::MESH_WAVE,
        WallDistanceMethod::DIRECTIONAL_MESH_WAVE,
        WallDistanceMethod::POISSON,
        WallDistanceMethod::EIKONAL,
        WallDistanceMethod::HAMILTON_JACOBI,
        WallDistanceMethod::ADVECTION_DIFFUSION,
        WallDistanceMethod::HYBRID_POISSON_EIKONAL}};

    // Iteration budget must scale with resolution: the pseudo-time methods
    // need O(n^2) steps and Poisson needs O(n). A single fixed budget across
    // resolutions silently turns a convergence gate into a resolution gate.
    // Use the longest grid axis as the effective n because it bounds the
    // wave-propagation distance across the domain.
    const std::size_t n_effective=std::max({grid.nx,grid.ny,grid.nz});

    std::vector<Row> rows;
    std::cout << "CFDX wall-distance complex-geometry benchmark\n";
    std::cout << "surface_vertices=" << surface.points.size()
              << " surface_triangles=" << surface.triangles.size()
              << " samples=" << grid.points.size() << " h=" << h
              << " n_effective=" << n_effective << "\n";
    std::cout << "poisson_operator,fluid_nodes=" << poisson_operator_audit.fluid_nodes
              << ",fluid_fluid_faces=" << poisson_operator_audit.fluid_fluid_faces
              << ",solid_faces=" << poisson_operator_audit.solid_faces
              << ",outer_faces=" << poisson_operator_audit.outer_faces
              << ",min_diagonal=" << poisson_operator_audit.min_diagonal
              << ",max_diagonal=" << poisson_operator_audit.max_diagonal
              << ",min_diagonal_dominance=" << poisson_operator_audit.min_diagonal_dominance
              << ",symmetry_error=" << poisson_operator_audit.symmetry_error
              << ",min_delta_over_h=" << poisson_offset_audit.min_delta_over_h
              << ",max_delta_over_h=" << poisson_offset_audit.max_delta_over_h
              << ",max_wall_coefficient=" << poisson_offset_audit.max_wall_coefficient
              << ",max_wall_coefficient_cell=" << poisson_offset_audit.max_wall_coefficient_cell
              << ",max_wall_coefficient_ijk=" << max_i << ":" << max_j << ":" << max_k
              << ",max_wall_coefficient_point=" << max_p.x << ":" << max_p.y << ":" << max_p.z
              << ",max_wall_coefficient_delta=" << poisson_offset_audit.max_wall_coefficient_delta
              << ",max_wall_coefficient_h=" << poisson_offset_audit.max_wall_coefficient_h
              << ",degenerate_wall_offsets=" << poisson_offset_audit.degenerate_count << "\n";
    std::cout << "method,l2_relative,linf_relative,near_wall_l2_relative,monotonicity_violations,invalid,iterations,residual_inf,time_ms,eikonal_init_ms,poisson_stage_ms,poisson_iterations,poisson_residual_inf,converged,min_distance,max_distance,wall_ray_hits,wall_ray_misses,wall_fallbacks,wall_bad_alignment,wall_min_alignment,poisson_phi_min,poisson_phi_max,poisson_grad_min,poisson_grad_max,poisson_distance_l2_error,poisson_distance_linf_error\n";

    for(const auto method:methods) {
        const std::size_t budget=wall_distance_iteration_budget(method,n_effective);
        double init_ms=0.0, poisson_ms=0.0;
        if(method==WallDistanceMethod::HAMILTON_JACOBI) {
            const auto p0=std::chrono::steady_clock::now();
            std::size_t init_iter=0;
            (void)eikonal_fast_sweep(surface,grid,budget,&init_iter);
            const auto p1=std::chrono::steady_clock::now();
            init_ms=std::chrono::duration<double,std::milli>(p1-p0).count();
        } else if(method==WallDistanceMethod::HYBRID_POISSON_EIKONAL) {
            const auto p0=std::chrono::steady_clock::now();
            std::size_t poisson_iter=0; double poisson_residual=0.0;
            (void)poisson_distance(surface,grid,budget,1.5,&poisson_iter,&poisson_residual);
            const auto p1=std::chrono::steady_clock::now();
            poisson_ms=std::chrono::duration<double,std::milli>(p1-p0).count();
        }
        const auto t0=std::chrono::steady_clock::now();
        const auto result=compute_wall_distance(method,surface,grid,budget);
        const auto t1=std::chrono::steady_clock::now();
        const double ms=std::chrono::duration<double,std::milli>(t1-t0).count();
        const auto m=compare_wall_distance(grid,reference,result.distance,2.0*h);
        std::size_t invalid=0;
        double min_distance=std::numeric_limits<double>::infinity();
        double max_distance=0.0;
        for(std::size_t i=0;i<result.distance.size();++i) {
            if(grid.solid[i]) continue;
            if(!result.valid[i] || !std::isfinite(result.distance[i])) { ++invalid; continue; }
            min_distance=std::min(min_distance,result.distance[i]);
            max_distance=std::max(max_distance,result.distance[i]);
        }
        if(!std::isfinite(min_distance)) min_distance=0.0;
        rows.push_back({result.method,m.l2_relative,m.linf_relative,m.near_wall_l2_relative,m.monotonicity_violations,
            result.residual_inf,invalid,result.iterations,ms,init_ms,poisson_ms,result.auxiliary_iterations,
            result.auxiliary_residual_inf,result.converged,min_distance,max_distance,
            result.poisson_wall_ray_hits,result.poisson_wall_ray_misses,result.poisson_wall_fallbacks,
            result.poisson_wall_bad_alignment,result.poisson_wall_min_alignment,result.poisson_phi_min,
            result.poisson_phi_max,result.poisson_grad_min,result.poisson_grad_max,
            result.poisson_distance_l2_error,result.poisson_distance_linf_error});
        std::cout << result.method << "," << std::setprecision(8)
                  << m.l2_relative << "," << m.linf_relative << ","
                  << m.near_wall_l2_relative << "," << m.monotonicity_violations << ","
                  << invalid << "," << result.iterations << "," << result.residual_inf << ","
                  << ms << "," << result.auxiliary_iterations << "," << result.auxiliary_residual_inf
                  << "," << (result.converged ? "true" : "false") << "\n";
    }

    std::ofstream csv(output);
    if(!csv) throw std::runtime_error("cannot open benchmark output: "+output);
    csv << "method,l2_relative,linf_relative,near_wall_l2_relative,monotonicity_violations,invalid,iterations,residual_inf,time_ms,eikonal_init_ms,poisson_stage_ms,poisson_iterations,poisson_residual_inf,converged,min_distance,max_distance,wall_ray_hits,wall_ray_misses,wall_fallbacks,wall_bad_alignment,wall_min_alignment,poisson_phi_min,poisson_phi_max,poisson_grad_min,poisson_grad_max,poisson_distance_l2_error,poisson_distance_linf_error\n";
    for(const auto& r:rows)
        csv << r.method << "," << r.l2 << "," << r.linf << "," << r.near_l2 << ","
            << r.violations << "," << r.invalid << "," << r.iterations << "," << r.residual << "," << r.ms
            << "," << r.init_ms << "," << r.poisson_ms << "," << r.poisson_iterations
            << "," << r.poisson_residual << "," << (r.converged ? "true" : "false") << "," << r.min_distance << "," << r.max_distance
            << "," << r.wall_ray_hits << "," << r.wall_ray_misses << "," << r.wall_fallbacks << ","
            << r.wall_bad_alignment << "," << r.wall_min_alignment << "," << r.poisson_phi_min << ","
            << r.poisson_phi_max << "," << r.poisson_grad_min << "," << r.poisson_grad_max << ","
            << r.poisson_distance_l2_error << "," << r.poisson_distance_linf_error << "\n";
    csv.close();

    // The exact method must be an exact self-reference. This is a regression
    // gate for the benchmark infrastructure, not a physical validation claim.
    if(rows.front().l2>1e-14 || rows.front().near_l2>1e-14)
        throw std::runtime_error("exact wall-distance reference is not self-consistent");
    for(const auto& row:rows)
        if(row.invalid!=0)
            throw std::runtime_error("invalid fluid distance field for "+row.method);
    return 0;
}

```

### `docs/validation/NACA0012_IMPORT_REPAIR_2026-10-25.md`

SHA256 `8b9f7725fdafc564f2c92469995d28fbec2ac4f87ae5a5751cbeb23d8f54e117`

````markdown
# NACA0012: importer repair and verification boundary, 2026-10-25

Baseline: `71858a20fbfc67006824936d0646b3dd1cb79394`. This record distinguishes
verified input/topology repairs from physical flow qualification. No NACA force
or turbulence status is promoted by this change.

## The different NACA cases

| Path | Declared problem | Evidence supplied by that path |
|---|---|---|
| `tests/data/su2/inv_NACA0012_basic.cfg` | Euler, Mach 0.8, incidence 1.25 degrees | SU2 format-reader fixture; no CFDX compressible production solver |
| `tests/data/su2/incomp_NACA0012.cfg` | Incompressible RANS, SA, incidence 10 degrees; source values imply Re approximately 6 million | Imported case metadata; production executable rejects RANS |
| `test_naca0012_qualification` | Incompressible laminar Re=1000, incidence zero, generated O-grid | SIMPLE + upwind solve in the long test; the quick test does not solve the flow |
| N9-S8 | Laminar Re=1000, incidence zero on imported NASA/TMR meshes | Contract and planned six-algorithm/three-mesh campaign; not a turbulent experimental comparison |
| TUR-NACA0012 | Turbulent NACA validation against independent data | Still requires a production RANS path and retained physical evidence |

The NASA turbulent validation problem uses Re=6 million and a modified sharp
trailing edge. Its grid may be reused in a declared laminar study, but the
experimental turbulent drag is not the oracle for a Re=1000 laminar solve.
The generated O-grid uses the finite-thickness original trailing edge; those
geometry definitions must also be recorded before interpreting a drag difference.

Primary references:

- https://tmbwg.github.io/turbmodels/naca0012_grids.html
- https://tmbwg.github.io/turbmodels/naca0012_val.html
- https://tmbwg.github.io/turbmodels/NACA0012_grids/n0012_225-65.nmf

## Reproduced defects in the original NASA/TMR importer

The official `n0012_225-65.p2dfmt` starts with three integers: `1`, `225`, `65`.
The first integer is the number of blocks. The old parser instead interpreted
`1` and `225` as the grid dimensions and rejected the actual file with:

```text
ValueError: expected 452 numeric tokens for 1x225, got 29253
```

The advertised `.gz` files were also read as plain ASCII.

After bypassing only that header defect for diagnosis, the original `convert()`
assigned every j=0 edge to `airfoil`. On the official 225x65 grid this produced
224 wall faces instead of 128. Ninety-six of these faces are the two wake cuts;
their centres extend to x=453.135486 chords. The neutral map explicitly assigns
wall points 49..177 and a reversed one-to-one connection 1..49 to 225..177.
These are **48 pairs of internal wake faces**, not solid walls. The original
converter kept their 49 coincident vertex pairs topologically disconnected.

The i=0 and i=ni-1 sides were named `wake_lower` and `wake_upper`, although the
neutral map identifies them as downstream farfield. Extrusion patches were
named `_front_cells` and `_back_cells`, receiving generic patch type 7 from the
writer instead of EMPTY type 6.

The previous `test_n9_s8_naca0012_contract.py` checked filename dimensions and
words in a chapter. It did not open a PLOT3D file or inspect its converted
connectivity, so none of these defects were exercised.

## Repair

`import_tmr_naca0012.py` now:

1. reads the single-block formatted PLOT3D header, validates token count and
   finite coordinates, and supports gzip;
2. identifies the reversed coincident j=0 wake endpoint interval rather than
   inventing a wall-index ratio from the filename;
3. validates that the detected interval terminates at the unit-chord trailing
   edge, so a partly corrupted wake cannot silently become extra solid wall;
4. joins the paired vertices in both extrusion planes, creating ordinary
   internal owner/neighbour faces along the wake;
5. exports only the airfoil as WALL, combines all outer sides as farfield, and
   exports `front` and `back` as EMPTY;
6. retains the source SHA256 and topology counts in the HDF5 provenance. The
   status remains `imported-provenance-only`.

The implementation is deliberately scoped to the advertised single-block,
unit-chord NASA NACA C-grid family. Other topology is rejected explicitly.
It does not claim to be a generic PLOT3D multiblock connector.

## Verification actually executed

Nine added importer regression cases fail **8/9 on the baseline** and pass
**9/9 after repair**. They cover the real header form, compressed input,
malformed input, nonfinite coordinates, an independent miniature C-grid with
exact face-count oracles, missing wake matching, partial wake corruption and
provenance. They are collected by the existing full `tests/python` CI job.

The official 225x65 grid was then imported and read through the actual native
CFDX HDF5 reader, `validate_mesh` and `build_fv_geometry`, linked against the
baseline production core. The following are measured values, not expected
values copied into a numerical solution:

| Native check | Measured result |
|---|---:|
| Cells | 14,336 |
| Points after wake joining | 29,152 |
| Faces | 57,584 |
| Internal faces | 28,432 |
| Internal wake faces | 48 |
| Airfoil wall faces | 128 |
| Farfield faces | 352 |
| Front / back EMPTY faces | 14,336 each |
| Unassigned physical boundary faces | 0 |
| Minimum cell volume | 3.3070392178011831e-9 |
| Minimum oriented face dot product Sf·d | 1.6535196089005912e-9 |
| Geometry errors | 0 |
| Geometry warnings | 1, large aspect ratio |
| Maximum nonorthogonality | 12.06419507085366 degrees |
| Maximum skewness | 0.15920642368411458 |
| Maximum absolute surface closure | 3.5527136788005009e-15 |

Every cell volume and every face-normal orientation was checked for finiteness
and strict positivity. All boundary faces belonged to exactly one patch.

Official source hashes:

- `n0012_225-65.p2dfmt`: `66817d6bbc9353cb0ca9867f958da1bcce6485f7d68ec7038362a341bdb0fa8c`
- `n0012_225-65.nmf`: `c6c32668f19bbd1a89510dfeb3a43c924137ae90afa7772eeb83875f9d9b34b0`

Reproduction, with assets downloaded deliberately before the command:

```bash
python scripts/import_tmr_naca0012.py \
  --input /path/to/n0012_225-65.p2dfmt \
  --output /path/to/n0012_225-65.h5
PYTHONPATH=src/cfdx/python python -m pytest -q \
  tests/python/test_tmr_naca0012_import.py \
  tests/python/test_su2_naca0012_tutorial.py \
  tests/python/test_su2_naca0012_validation.py
```

The three selected Python files contain **12 passing tests** after the complete
change. A passed mesh import is not a passed flow simulation.

## Honest SU2 tutorial scope

The original Python command stopped at canonical conversion: JST has no faithful
registered CFDX convection counterpart. Its error printer also classified
`unsupported_nonblocking` entries as blocking because it searched for the
substring `block`. Its unused "theory" comparison treated a freestream static
pressure as stagnation pressure and did not provide an aerodynamic oracle.
The shell wrapper used unsupported options `--case`, `--output`, `--max-iter`.
The production executable accepts `--mesh`, `--output-dir`, `--iterations` and
explicitly rejects the compressible energy model.

The Python command now retains distinct reader, boundary/initial-condition,
canonical-conversion and physical-run verdicts in `reader_report.json`; the
shell command delegates to it. Its default exit code remains nonzero (2) when
full aerodynamic validation is unavailable. An explicit `--reader-only` permits
success for the restricted input-reader scope:

```bash
python validation/run_su2_naca0012.py --reader-only --output-dir /path/to/evidence
```

No JST substitution is introduced and no solver result is manufactured.

The current SU2 adapter still has these independently observed metadata defects,
which this importer-focused patch does not repair:

- `MARKER_EULER=(airfoil)` becomes `wall_no_slip` through a name-based fallback;
- `farfield` becomes `outlet`/`fixed` with zero velocity;
- initial scalar speed is 272 m/s, but the vector contains the unit incidence
  direction `[0.9997620270799091, 0.02181488503456112, 0]`.

The updated report exposes all three mismatches. Existing format-reader tests
check patch names and cardinalities, not faithful boundary/initial-condition
semantics. These defects must be fixed before that importer can be described as
physically equivalent to the original SU2 case.

## Remaining numerical limitations relevant to the laminar NACA solve

The existing generated O-grid n64 is topologically valid, but native geometry
inspection gives nonorthogonality 80.550987 degrees, maximum skewness 0.986432
and 82 warnings. Its wall-adjacent centre distance has median 0.0268116 chord,
minimum 0.00511106 and maximum 0.0326396. It is a coarse, highly skewed starting
mesh. The indexing of its circular outer contour also breaks the reflection
pairing of the airfoil surface: the maximum corresponding reflection difference
in the first layer is 0.00319299 chord. These are measured quality limitations;
their contribution to Cd/Cl must be quantified by controlled flow comparisons.

On both actual meshes, the exact manufactured affine field phi=x+2y was supplied
to the production `gauss_gradient_with_boundary`. Only cells with no physical
boundary face were included, so the result is independent of boundary-value
approximations. WLS received explicit exact Dirichlet values at physical faces
and the correct zero spanwise derivative.

| Mesh | Interior cells | Production Gauss gradient max absolute error | WLS max absolute error |
|---|---:|---:|---:|
| Generated O-grid n64 | 2,838 | 0.10987995423605813 | 1.1082784313214814e-14 |
| Repaired TMR 225x65 | 13,858 | 0.17283073646457187 | 5.4393680889484996e-9 |
| Repaired TMR 449x129 | 56,386 | 0.10357401848014633 | 1.3019123605904252e-8 |
| Repaired TMR 897x257 | 227,458 | 0.07351694706579924 | 2.9871163766996353e-8 |

This demonstrates a reconstruction error before any nonlinear flow solve. The
production helper still uses Euclidean interpolation distances, separately from
`core/numerics/gradient.h`, and has no face-skew correction. Merely correcting a
standalone gradient kernel cannot certify this production path. These numbers
are **not** drag errors and do not by themselves establish the whole cause of
an aerodynamic discrepancy.

Further requirements remain unchanged: fix and independently verify the
pressure/velocity face-flux contract; validate the actual production gradient
and nonorthogonal diffusion; freeze consistent farfield conditions; compare
pressure and viscous force parts and surface Cp/Cf; execute the declared three
meshes and six algorithms; retain true residuals, independent balances and an
appropriate physical reference. A unit test, quick contract, converged linear
solve or this mesh repair cannot replace those results.

## Refined official grids and boundary-membership cost

The same native reader/geometry checks were also completed on the two finer
official files. Their gzip payloads were retrieved from the NASA archive,
`https://www.nasa.gov/wp-content/uploads/2026/02/naca0012-grids.zip`, and checked
against its ZIP member CRC and length before decompression.

| Check | 449x129 | 897x257 |
|---|---:|---:|
| Cells | 57,344 | 229,376 |
| Points after joining | 115,648 | 460,672 |
| Faces | 229,856 | 918,464 |
| Internal faces | 114,208 | 457,792 |
| Internal wake faces | 96 | 192 |
| Airfoil wall faces | 256 | 512 |
| Farfield faces | 704 | 1,408 |
| Front/back EMPTY faces each | 57,344 | 229,376 |
| Unassigned boundary faces | 0 | 0 |
| Minimum cell volume | 5.1820719983179846e-10 | 9.579601570460766e-11 |
| Minimum oriented face dot product | 2.5910359991589918e-10 | 4.7897995986945982e-11 |
| Geometry errors | 0 | 0 |
| Geometry warnings | 1, aspect ratio | 1, aspect ratio |
| Maximum nonorthogonality, degrees | 9.2209192083970688 | 6.1731507752808854 |
| Maximum skewness | 0.10395424486342719 | 0.055277414197076311 |
| Maximum absolute surface closure | 1.7763568394002505e-15 | 4.4408920985006262e-16 |

Source SHA256:

- `n0012_449-129.p2dfmt.gz`: `6d60dfa6ef9c7544ee156276fcc5487d7eacddadbdf6e81cc6f10e78fb03eff8`
- `n0012_897-257.p2dfmt.gz`: `c603ea77017ae8b97af581237b81a8251fd2d3af1c49c8c6f2cd2687aafcb28e`

During these checks, `BoundaryPatches::is_consistent_with_mesh` was found to
rescan every patch face for each physical boundary face. With the two large
EMPTY patches, this is quadratic in the number of cells. The native HDF5
reader invokes this validation, and subsequent `validate_mesh` invokes it again.
The repair retains the duplicate/range and internal-face checks, records listed
faces in a bitmap, and then scans the boundary once. It is linear in the number
of faces plus patch entries.

A benchmark on the three actual imported grids called the exact original
function body and its replacement on the same loaded mesh. Both returned true.
These are single-process wall-clock observations for that function alone,
not an end-to-end CFD speedup guarantee:

| Grid | Boundary faces | Original, seconds | Replacement, seconds |
|---|---:|---:|---:|
| 225x65 | 29,152 | 0.121096512 | 0.000140853 |
| 449x129 | 115,648 | 2.34260247 | 0.000537623 |
| 897x257 | 460,672 | 31.2930124 | 0.002588804 |

Five membership regression cases explicitly check a complete patch cover,
an absent face, a repeated face, an internal face in a patch, and an out-of-range
face. The expanded `test_mesh` passes 16/16 both with the original and repaired
function, demonstrating preserved acceptance/rejection decisions. The existing
`test_boundary` passes 11/11 with the repair. The complete native geometry and
manufactured-gradient diagnostics on all three actual grids are byte-for-byte
identical before and after this membership-only change. No timing threshold is
used as a unit-test oracle.

````

### `docs/validation/NACA0012_RUNTIME_AUDIT_2026-10-25.md`

SHA256 `bf3b92d56cde540272768bf9d013db16d4aab849165c96c706c0bc80ea548167`

````markdown
# NACA0012 runtime audit, 2026-10-25

Baseline: `71858a20fbfc67006824936d0646b3dd1cb79394`.
This is an execution record for the existing laminar qualification, separate
from the NASA mesh importer repair. It does not promote the NACA case to
physically validated status.

## Actual baseline outcome

The complete baseline CTest campaign ran `test_naca0012_qualification` on the
three meshes generated by the repository. It aborted after 20.01 seconds.
`test_naca0012_wall_shear_audit` subsequently failed by the same mechanism after
19.78 seconds. Their quick counterparts passed because they do not execute the
full three-level flow campaign.

The controls in `tests/validation/test_naca0012_geometry.cpp` were unmodified:
laminar incompressible Re=1000, zero incidence, SIMPLE, upwind convection,
linear tolerance 1e-10 and limit 4000, nonlinear tolerance 1e-8 and limit 8000,
velocity relaxation 0.6 and pressure relaxation 0.25. The source prescribes
U=(1,0,0) and p=0 at the farfield, no slip on the airfoil, and spanwise EMPTY
faces. This is not the turbulent NASA Re=6-million validation case.

| Level | Actual outcome |
|---|---|
| n64, 3,096 cells | Solver converged in 435 iterations |
| n128, 12,336 cells | Pressure solve returned NOT_APPLICABLE at nonlinear iteration 20; three lower-relaxation retries also failed |
| n256, 49,248 cells | Not reached by either full test |

The converged n64 output was Cd=0.0903395, Cl=1.55115e-5,
Cm=-0.00094537, Cd_pressure=0.0345473, Cd_viscous=0.0557922 and
continuity_linf=1.50029e-10. A converged coarse solution does not establish
mesh convergence or agreement with a physical reference. The final test itself
labels its prospective successful outcome DIAGNOSTIC, with the physical
reference gate still pending.

The wall-shear audit applied three postprocessors to the same n64 solution:

| Reconstruction | Viscous Cd |
|---|---:|
| Production Gauss | 0.0557922 |
| Weighted least squares | 0.0581681 |
| One-dimensional wall-normal estimate | 0.106364 |

Its reported relative spread was 0.688604. This is a measured sensitivity of
postprocessing on one coarse field, not three independent flow solutions and
not proof that the largest value is the correct one.

## The default failure message hides the actual residual

The baseline exception prints:

```text
pressure-correction solve did not converge
(status=3, iterations=3, residual=0.000000, relative=0.000000)
```

Status 3 is `SolverStatus::NOT_APPLICABLE`, not exhaustion of the 4000-iteration
linear budget. The automatic pressure plan resolves to CG with
`NativeBoomerAMGPreconditioner`, using `ReusableCgContext`.

A diagnostic-only local include overlay captured the real matrix, RHS and CG
iterate on that failure. It changed neither the flow formulas nor the controls.
A second diagnostic run printed every nonlinear iteration and captured the
first pressure matrix as well. These runs reproduced the original failure;
they are follow-up observations of the same defect, not additional independent
qualification cases.

At the first failed pressure solve:

- RHS norm: 5.5156653682823174e12;
- true residual norm, recomputed as ||b-Ax||: 1.7332238406933466e11;
- true relative residual: 0.031423658343383246;
- requested relative tolerance: 1e-10;
- stored `SolverResult.residual` and `residual_relative`: both zero;
- initial preconditioned product r^T M^-1 r: -4.3863563470325308e37.

The early NOT_APPLICABLE exits in `cg_solver.h` do not consistently populate
the residual fields. An independently instrumented replay identified the
specific exit here: the nonpositive `replaced_rsnew` check after residual
replacement, at baseline `cg_solver.h:361`. The initial r^T M^-1 r had already
been negative; CG only checked its finiteness and took a negative first step
coefficient. Both the missing true-residual reporting and the missing initial
positivity check deserve explicit regression coverage.

## Nonlinear divergence precedes the linear refusal

Printing every iteration instead of only 1 and 200 exposed the earlier growth:

| Nonlinear iteration | Momentum RHS-relative residual | Momentum raw residual | Pressure linear relative residual |
|---|---:|---:|---:|
| 1 | 0.02170991661 | 0.01063918177 | 5.796087420e-11 |
| 5 | 0.2776201489 | 0.1360274329 | 2.989268378e-11 |
| 10 | 1.377742198 | 11.22542016 | 1.022865973e-11 |
| 13 | 51.84263613 | 6062.070165 | 2.017762218e-11 |
| 16 | 8479.466128 | 2.491988208e9 | 2.603139310e-11 |
| 19 | 2.760647353e9 | 5.450210718e24 | 5.888560796e-11 |

At iteration 19 the flux/velocity mismatch was 1.8018721310623572e14.
The linear pressure solves still met their requested tolerance up to that
iteration. Therefore the final AMG/CG refusal is a later failure in a nonlinear
calculation that was already diverging. Replacing the last pressure solver
alone cannot certify a physically meaningful flow.

The earliest largest momentum errors are adjacent to the airfoil. This narrows
where to investigate, but does not isolate which production discretization
term initiates the instability. The existing generated O-grids have the
following independently measured geometry/reconstruction limitations:

| O-grid level | Max nonorthogonality | Max skewness | Geometry warnings | Affine production-gradient maximum error |
|---|---:|---:|---:|---:|
| n64 | 80.550987 degrees | 0.986432 | 82 | 0.1098799542 |
| n128 | 81.280277 degrees | 0.988442 | 163 | 0.1671984660 |
| n256 | 81.647616 degrees | 0.989393 | 351 | 0.1594223136 |

All three meshes have positive finite volumes, correctly oriented faces and
zero topological/geometry errors. The affine test uses phi=x+2y, includes only
cells without a physical boundary face, and precedes any flow solve. Its
corresponding WLS errors are 1.11e-14, 6.14e-14 and 4.79e-13. These are
reconstruction errors, not measured drag errors. A controlled production
discretization comparison is still needed to assign the primary nonlinear
instability to a particular term.

### One controlled relaxation diagnostic

A separate n128 A/B experiment changed only the initial velocity relaxation
from 0.6 to 0.3. Mesh, boundary conditions, schemes, alpha_p=0.25, residual
tolerances and the configured 8000/4000 nonlinear/linear iteration budgets
were preserved. The wrapper declared an external 300-second time budget and
an early-stop threshold of momentum RHS-relative residual greater than 1e6.
Printing each iteration and collecting the initial pressure snapshot were
observational changes only.

This run was stopped at iteration 48, after 8.019055 seconds, for manifest
nonlinear divergence. The last measured momentum RHS-relative residual was
3.2768351750380591e7, raw momentum residual 1.1839691475182645e19,
normalized continuity 1.9507736565740657, and flux/velocity mismatch
2.2336081815654251e10. The pressure linear solve still converged to
2.3469556868328088e-11 at that iteration, below the unchanged 1e-10 tolerance.

Thus halving this relaxation delays but does not remove the nonlinear
instability. This single comparison does not rule out other stable controls;
no parameter sweep was used and no force qualification is inferred. The
independent pressure-solver failure remains a later issue, rather than the
only reason the original flow calculation fails.

## The captured pressure matrix remains grounded and SPD

The matrix has 12,336 rows and 61,166 entries. It is exactly symmetric, has a
positive diagonal, strictly negative off-diagonal entries, one connected
component and 257 strictly dominant rows attached to prescribed farfield
pressure. There is no missing pressure gauge in this snapshot.

The diagonal ranges from 1.6798584831606658e-16 to
0.011604898446567188, a ratio of about 6.91e13. The independent normalized
analysis of B=D^-1/2 A D^-1/2 found:

- lambda_min = 1.21312018e-9;
- lambda_max approximately 1.999999999;
- normalized condition number approximately 1.649e9;
- positive LDL^T pivots, smallest 7.707e-8;
- factor symmetry defect 1.67e-16.

The fine matrix's SPD character is strongly supported, while its scaling and
conditioning have become extreme. A separate SciPy sparse direct solve gave
relative residual 2.41565e-10 and a maximum pressure correction 2.80923e24;
this is another observation of the damaged nonlinear state, not an accepted
flow solution.

## Reusing the AMG hierarchy breaks its CG-compatible action

At the failure, the reusable context had one hierarchy setup and 19 numeric
updates. An isolated native replay built AMG from the first pressure matrix
and then updated it directly to the captured iteration-20 matrix. It reproduced
the original negative r^T M^-1 r and the complete failed CG result bit for bit.
Thus intermediate flow states are not needed to reproduce this particular
AMG reuse defect.

All following results use the same captured A and b, zero initial correction,
linear tolerance 1e-10 and budget 4000:

| Linear path | Initial r^T M^-1 r | Outcome | Iterations | Recomputed relative residual |
|---|---:|---|---:|---:|
| Original first-step AMG hierarchy + numeric update, CG | -4.386356347e37 | NOT_APPLICABLE | 3 | 0.03142365834 |
| Newly built AMG, CG | 2.113141873e33 | MAX_ITER_REACHED / exhausted useful recurrence | 126 | 5.482475639e-10 |
| Newly built AMG, FGMRES | 2.113141873e33 | CONVERGED | 517 | 9.873355974e-11 |
| Jacobi, CG | 1.530768460e33 | MAX_ITER_REACHED | 4000 | 6.824463987e-9 |

FGMRES convergence is a real linear snapshot result. No force, conservation or
nonlinear-convergence claim is inferred from it.

### Measured mechanism inside the reused hierarchy

The native AMG implementation uses six pre- and six post-Jacobi sweeps with
fixed omega=0.7. `update_values` retains the old interpolation P and rebuilds
the Galerkin operators. It does not requalify the stability of that fixed
smoother on the refreshed levels.

A read-only diagnostic getter exported the actual stored level matrices.
Its direct `MatrixFreeVcyclePreconditioner` replay first reproduced the native
wrapper's r^T M^-1 r exactly, so the level analysis concerns the same hierarchy.
For a normalized symmetric level operator, Jacobi's eigenmode factor is
1-omega*lambda. Stability requires lambda_max < 2/omega = 2.857142857.

| Reused hierarchy level | Unknowns | Measured lambda_max | High-mode factor 1-0.7*lambda_max |
|---|---:|---:|---:|
| 3 | 1297 | 2.976252949 | -1.083377064 |
| 4 | 566 | 3.518576285 | -1.463003399 |
| 5 | 288 | 3.987227820 | -1.791059474 |
| 6 | 163 | 4.455448843 | -2.118814190 |
| 7 | 116 | 4.173164525 | -1.921215168 |
| 8 | 102 | 3.337720948 | -1.336404663 |

The level-6 eigenpair residual is 2.20e-15. These modes are amplified on each
smoothing sweep. In contrast, every coarse level of newly built AMG is below
the stability bound; its largest observed coarse lambda_max is 1.988241902.
The issue is not an indefinite fine pressure matrix. It is an old hierarchy
whose refreshed levels no longer support the configured fixed Jacobi damping.

A repair must preserve the explicit numeric-only contract of `update_values`:
either requalify/adapt smoothing at refresh, or report lost applicability so
the caller can deliberately rebuild a hierarchy. It must also retain true
residual checks. Simply increasing linear budgets, weakening tolerances, or
declaring the n64 force validated would not address these observations.

## Retained evidence and reproduction

Workspace evidence is in `audit_naca/flow/` next to the baseline checkout:

- `naca128_diagnostic.cpp` and `include/`: diagnostic-only wrapper and overlay;
- `naca128_diagnostic.log`: original failure and four true residual snapshots;
- `pressure_failure_1.mtx`: failing iteration-20 matrix in Matrix Market format;
- `pressure_failure_1.vectors`: columns RHS and failed CG correction, 17 digits;
- `trace/naca128_trace.log`: every accepted nonlinear iteration;
- `trace/pressure_failure_1.mtx`: initial iteration-1 matrix;
- `relaxed/`: source fixture with the single alpha_u change, bounded runner,
  and complete `naca128_alpha03.log`;
- `replay_pressure_snapshot.cpp` and `.log`: identical-system linear comparison;
- `export_amg_levels.cpp`, `amg_include/`, `amg_levels/`: actual refreshed/fresh
  level matrices and spectral analysis script, JSON and log.

The separate solver audit retains `naca_pressure_matrix_analysis.*` and the
CG exit trace. The original baseline CTest log is
`audit_output/ctest_release_all.log`; the full NACA baseline was not rerun as a
second qualification campaign.

To replay the isolated linear failure, from the retained evidence directory:

```bash
./replay_pressure_snapshot \
  pressure_failure_1.mtx pressure_failure_1.vectors \
  trace/pressure_failure_1.mtx

cd amg_levels
../export_amg_levels ../trace/pressure_failure_1.mtx \
  ../pressure_failure_1.mtx ../pressure_failure_1.vectors
python analyze_levels.py
```

The executable sources compile with C++17 against the unchanged baseline
headers and `libcfdx_core.a`. The flow diagnostic additionally needs the native
HDF5 build. Neither diagnostic overlay is a production change. The importer
and boundary-membership fixes are recorded separately in
`NACA0012_IMPORT_REPAIR_2026-10-25.md`.

````

### `docs/validation/SU2_NACA0012_READER.md`

SHA256 `fb44dc17627997b9666cdddb19e81b34ef22db780550bce20177901bf3311dba`

````markdown
# SU2 NACA0012 validation case — reader qualification

## Purpose

This case is the external-solver I/O qualification case for CFDX.

It replaces the axisymmetric sphere case as the **format-reader validation
target** for the unified solver-I/O path. It does **not** remove or invalidate
the native CFDX axisymmetric VMFL036 capability; that remains a separate
physics/solver qualification.

## Source case

The case is the SU2 NACA0012 Euler validation/tutorial family:

- 2-D NACA0012 external flow;
- Euler solver;
- Mach number 0.8;
- angle of attack 1.25 degrees;
- freestream pressure 101325 Pa;
- freestream temperature 288.15 K;
- SU2 mesh format;
- two boundary markers: `airfoil` and `farfield`.

The repository fixture is the existing SU2 NACA0012 mesh/config pair under
`tests/data/su2/`. The mesh contains 5233 vertices and 10216 triangular
cells.

SU2's current regression suite also contains a NACA0012 Euler regression
case, with documented iteration/QoI reference values. The CFDX PR intentionally
uses the SU2 input files as an **input-format/readability qualification**, not
as evidence that CFDX reproduces SU2's solver result.

## What the PR verifies

The test reconstructs a conventional SU2 case directory:

```text
naca0012/
├── mesh.su2
└── config.cfg
```

and calls the same directory-level `Su2Adapter.convert()` entry point used
by the unified importer.

It verifies:

1. SU2 mesh dimensionality;
2. vertex and cell cardinalities;
3. triangular cell connectivity;
4. boundary-marker discovery;
5. SU2 configuration parsing;
6. Euler physics identification;
7. freestream pressure/temperature;
8. iteration count;
9. conversion without blocking Gap Analysis findings.

## Results are deliberately not used

The repository's small `mismatched_solution.csv` fixture contains only 10 rows
while the NACA0012 mesh contains 10216 cells. It is therefore not a valid
solution for this mesh, and it is deliberately kept out of the mesh directory
so the adapter does not auto-discover it and block every conversion of this
case.

The test does **not** weaken the cardinality check and does not use that
fixture as validation evidence. Missing solution output is recorded as a
non-blocking diagnostic because this PR validates case reading, not result
import.

A future result-import qualification should use an actual SU2 solution whose
field localisation/cardinality matches the mesh.

## Acceptance

The reader qualification passes only when:

- the real SU2 mesh is parsed;
- the real SU2 configuration is parsed;
- mesh/config metadata are internally coherent;
- no blocking Gap Analysis finding is produced.

No CFDX numerical result is compared here. Numerical validation remains a
separate solver-level campaign.

## Reproducible command and remaining fidelity limits (2026-10-25)

Run the explicitly limited reader scope with:

```bash
python validation/run_su2_naca0012.py --reader-only --output-dir reader-evidence
```

The JSON report separates input parsing from boundary/initial-condition
fidelity, canonical numerical conversion and physical simulation. Without
`--reader-only`, the legacy command exits nonzero when the full aerodynamic
validation cannot run. It never substitutes another scheme for unsupported JST
and never presents an isentropic static-pressure expression as a NACA Cp oracle.

The current adapter's Euler wall, farfield velocity and velocity-vector scaling
still require correction; a reader PASS only covers the explicitly listed
format checks. See [the importer repair and audit record](NACA0012_IMPORT_REPAIR_2026-10-25.md)
for measured counterexamples and the separate NASA/TMR topology repair.

````

### `scripts/check_cpp_test_execution.py`

SHA256 `b5fca3597d9ead441428cc291aeba12ef55e54f707d4a68e1138f6333b7d41d8`

```python
#!/usr/bin/env python3
"""Reject deferred CFDX harness registrations without a run_all() call.

This narrow static guard recognizes literal-name/lambda run_case registrations
from common/test_harness.h. Local numerical helpers named run_case, comments,
and string literals are not test registrations. It complements execution; it
does not prove that an arbitrary call is reachable or that an oracle is sound.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re


def cpp_tokens(source: str) -> list[tuple[str, str]]:
    tokens: list[tuple[str, str]] = []
    i = 0
    while i < len(source):
        if source[i].isspace():
            i += 1
        elif source.startswith("//", i):
            end = source.find("\n", i + 2)
            i = len(source) if end < 0 else end + 1
        elif source.startswith("/*", i):
            end = source.find("*/", i + 2)
            i = len(source) if end < 0 else end + 2
        elif source.startswith('R"', i):
            opening = source.find("(", i + 2, i + 19)
            if opening < 0:
                i += 1
                continue
            closing = ")" + source[i + 2 : opening] + '"'
            end = source.find(closing, opening + 1)
            end = len(source) if end < 0 else end + len(closing)
            tokens.append(("string", source[i:end]))
            i = end
        elif source[i] in "\"'":
            quote, start = source[i], i
            i += 1
            while i < len(source):
                if source[i] == "\\":
                    i += 2
                elif source[i] == quote:
                    i += 1
                    break
                else:
                    i += 1
            tokens.append(("string", source[start:i]))
        elif source[i].isalpha() or source[i] == "_":
            match = re.match(r"[A-Za-z_][A-Za-z_0-9]*", source[i:])
            assert match is not None
            value = match.group(0)
            tokens.append(("identifier", value))
            i += len(value)
        else:
            tokens.append(("punctuation", source[i]))
            i += 1
    return tokens


def inspect_source(source: str) -> tuple[int, bool]:
    tokens = cpp_tokens(source)
    values = [token[1] for token in tokens]
    uses_harness = any(
        values[i : i + 3] == ["#", "include", '"common/test_harness.h"']
        for i in range(len(tokens) - 2)
    )
    if not uses_harness:
        return 0, False
    registrations = sum(
        values[i] == "run_case"
        and values[i + 1] == "("
        and tokens[i + 2][0] == "string"
        and values[i + 3 : i + 5] == [",", "["]
        for i in range(len(tokens) - 4)
    )
    executes = any(
        values[i : i + 3] == ["run_all", "(", ")"]
        for i in range(len(tokens) - 2)
    )
    return registrations, registrations > 0 and not executes


def self_test() -> None:
    prefix = '#include "common/test_harness.h"\n'
    registration = 'int main() { run_case("example", [] { fail(); }); return 0; }'
    assert inspect_source(prefix + registration) == (1, True)
    assert inspect_source(prefix + registration.replace("return 0", "return run_all()")) == (1, False)
    assert inspect_source(prefix + '// run_case("ignored", [] {});\n') == (0, False)
    assert inspect_source(prefix + 'const char *x = R"tag(run_case("ignored", [] {});)tag";') == (0, False)
    assert inspect_source(prefix + 'double run_case(int n) { return n; } int main() { run_case(16); }') == (0, False)
    assert inspect_source(prefix + registration + '// run_all()\nconst char *x="run_all()";') == (1, True)
    print("Harness guard self-tests: 6 passed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    files, cases, missing = 0, 0, []
    for path in sorted((args.root / "tests").rglob("*.cpp")):
        count, absent = inspect_source(path.read_text(encoding="utf-8"))
        if count:
            files += 1
            cases += count
        if absent:
            missing.append((path.relative_to(args.root), count))
    for path, count in missing:
        print(f"FAIL: {path}: {count} deferred case(s) but no run_all() call")
    print(f"Harness execution guard: {files} files, {cases} registered cases, {len(missing)} missing executions")
    return 1 if missing else 0


if __name__ == "__main__":
    raise SystemExit(main())

```

### `scripts/import_tmr_naca0012.py`

SHA256 `e7527a1410f35d8689b75d9ee069b067f360b4c6c998bfcf28e0ef8a95deee81`

```python
#!/usr/bin/env python3
"""Import the single-block NASA/TMR NACA0012 PLOT3D C-grid.

The j=0 boundary contains the airfoil AND two coincident wake cuts.  The
wake cuts are joined as internal faces; only the closed airfoil is a wall.
This is an input/topology conversion, not a physical qualification.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
from pathlib import Path
import re

import h5py
import numpy as np

from generate_2d_validation_meshes import write_mesh

TMR_URL = "https://tmbwg.github.io/turbmodels/naca0012_grids.html"
TMR_ARCHIVE_URL = "https://www.nasa.gov/wp-content/uploads/2026/02/naca0012-grids.zip"


def read_plot3d_2d(path: Path):
    """Read NASA's formatted multigrid PLOT3D, with one 2-D block.

    The first record is nbl, not ni.  Compressed files in the advertised
    grid family use the same format and are read without manual unpacking.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="ascii") as stream:
        tokens = stream.read().replace(",", " ").split()
    if len(tokens) < 3:
        raise ValueError(f"PLOT3D block header is missing: {path}")
    nblocks = int(tokens[0])
    if nblocks != 1:
        raise ValueError(f"only single-block NASA/TMR PLOT3D is supported, got {nblocks}")
    ni, nj = int(tokens[1]), int(tokens[2])
    if ni < 5 or nj < 2:
        raise ValueError(f"invalid NASA/TMR C-grid dimensions: {ni}x{nj}")
    expected = 3 + 2 * ni * nj
    if len(tokens) != expected:
        raise ValueError(f"expected {expected} numeric tokens for {ni}x{nj}, got {len(tokens)}")
    values = np.asarray(
        [float(value.replace("D", "E").replace("d", "e")) for value in tokens[3:]],
        dtype=np.float64,
    )
    if not np.isfinite(values).all():
        raise ValueError("NASA/TMR grid coordinates must be finite")
    return values[:ni * nj].reshape((nj, ni)), values[ni * nj:].reshape((nj, ni))


def convert(x, y):
    """Join the coincident j=0 wake and extrude with EMPTY front/back.

    NASA's neutral map declares a reversed one-to-one connection between
    the leading and trailing index intervals on j=0.  Determine that
    interval from coincident endpoint pairs, including the trailing-edge
    pair.  No wall-index ratio is guessed from the filename.  Reject a
    different topology instead of manufacturing a wall in its wake.
    """
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    if x.ndim != 2 or x.shape != y.shape or min(x.shape) < 2:
        raise ValueError("C-grid x/y arrays must have the same two-dimensional shape")
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("C-grid coordinates must be finite")
    nj, ni = x.shape
    # Absolute tolerance is in the unit-chord NASA grid coordinates.  It
    # is deliberately independent of the distant farfield radius.
    paired_points = 0
    for i in range(ni // 2):
        opposite = ni - 1 - i
        if max(abs(x[0, i] - x[0, opposite]),
               abs(y[0, i] - y[0, opposite])) > 1e-12:
            break
        paired_points += 1
    first_wall = paired_points - 1
    last_wall = ni - paired_points
    if paired_points < 2 or last_wall - first_wall < 3:
        raise ValueError("expected a reversed coincident j=0 wake and a nondegenerate airfoil")

    # The supplied NASA family is normalized to unit chord.  Check the
    # independently known geometry before accepting the detected seam: a
    # partly corrupted wake must not be reclassified as extra solid wall.
    wall_x = x[0, first_wall:last_wall + 1]
    wall_y = y[0, first_wall:last_wall + 1]
    geometry_tolerance = 1e-6  # admits the documented original TMR TE typo.
    if (abs(wall_x[0] - 1.0) > geometry_tolerance or
            abs(wall_y[0]) > geometry_tolerance or
            abs(wall_x.min()) > geometry_tolerance or
            wall_x.max() > 1.0 + geometry_tolerance or
            (np.diff(x[0, :paired_points]) >= 0.0).any()):
        raise ValueError("matched wake must terminate at the unit-chord NACA trailing edge")

    # Merge only the explicitly matched wake vertices.  This gives each
    # wake face two owners through write_mesh's ordinary connectivity
    # matching; it is not a periodic/no-slip physical boundary.
    representatives = np.arange(ni * nj)
    for i in range(paired_points):
        representatives[ni - 1 - i] = i
    unique, inverse = np.unique(representatives, return_inverse=True)
    nplane = len(unique)
    flat_x, flat_y = x.ravel(), y.ravel()
    points = [(float(flat_x[i]), float(flat_y[i]), z)
              for z in (0.0, 1.0) for i in unique]

    def p(j, i, layer=0):
        return int(inverse[j * ni + i]) + layer * nplane

    cells = []
    for j in range(nj - 1):
        for i in range(ni - 1):
            lower = (p(j, i), p(j, i + 1), p(j + 1, i + 1), p(j + 1, i))
            cell = lower + tuple(v + nplane for v in lower)
            if len(set(cell)) != 8:
                raise ValueError("wake stitching created a degenerate hexahedron")
            cells.append(cell)
    airfoil = [(p(0, i + 1), p(0, i), p(0, i, 1), p(0, i + 1, 1))
               for i in range(first_wall, last_wall)]
    farfield = [(p(nj - 1, i), p(nj - 1, i + 1),
                 p(nj - 1, i + 1, 1), p(nj - 1, i, 1))
                for i in range(ni - 1)]
    # i=0 and i=ni-1 are downstream farfield boundaries in the neutral
    # map.  The actual wake connection is along j=0, handled above.
    farfield += [(p(j, 0), p(j + 1, 0), p(j + 1, 0, 1), p(j, 0, 1))
                 for j in range(nj - 1)]
    farfield += [(p(j + 1, ni - 1), p(j, ni - 1),
                  p(j, ni - 1, 1), p(j + 1, ni - 1, 1))
                 for j in range(nj - 1)]
    front = [tuple(cell[:4]) for cell in cells]
    back = [tuple(reversed(cell[4:])) for cell in cells]
    return points, cells, {"airfoil": airfoil, "farfield": farfield,
                           "front": front, "back": back}


def infer_dimensions(path: Path):
    match = re.search(r"_(\d+)-(\d+)\.p2dfmt(?:\.gz)?$", path.name)
    if not match:
        raise ValueError(f"cannot infer TMR grid dimensions from {path.name}")
    return int(match.group(1)), int(match.group(2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    ni_name, nj_name = infer_dimensions(args.input)
    x, y = read_plot3d_2d(args.input)
    nj, ni = x.shape
    if (ni, nj) != (ni_name, nj_name):
        raise ValueError(f"filename declares {ni_name}x{nj_name}, file contains {ni}x{nj}")
    points, cells, boundaries = convert(x, y)
    write_mesh(args.output, points, cells, boundaries)
    with h5py.File(args.output, "r+") as h5:
        h5.attrs["validation_source"] = "NASA TMR NACA0012"
        h5.attrs["validation_source_url"] = TMR_URL
        h5.attrs["validation_archive_url"] = TMR_ARCHIVE_URL
        h5.attrs["validation_source_sha256"] = hashlib.sha256(args.input.read_bytes()).hexdigest()
        h5.attrs["validation_grid"] = f"{ni}x{nj}"
        h5.attrs["validation_topology"] = "C-grid; j=0 reversed wake stitched"
        h5.attrs["validation_wall_faces"] = len(boundaries["airfoil"])
        h5.attrs["validation_wake_pairs"] = (ni - 1 - len(boundaries["airfoil"])) // 2
        h5.attrs["validation_physics"] = "N9-S8 laminar incompressible campaign"
        h5.attrs["validation_status"] = "imported-provenance-only"
    print(f"Imported {args.input} -> {args.output} ({ni}x{nj}); "
          f"airfoil_faces={len(boundaries['airfoil'])}; physical qualification NOT RUN")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

```

### `src/cfdx/core/linalg/cg_solver.h`

SHA256 `e465e903e239fe9a29ec318cb7d665785b59f4dbe7dfae8ef7ae33c0322b0d77`

```cpp
// M0.8-T04 — Conjugate Gradient solver
//
// Spécification CFDX v0.7 §35 :
//   Premiers solveurs : CG, BiCGStab, GMRES.
//
// CG (Hestenes-Stiefel) :
//   Résout A x = b pour A symétrique définie positive.
//  avec préconditionneur diagonal (Jacobi) par défaut.

#pragma once

#include "sparse_matrix.h"
#include "vector.h"
#include "mixed_precision.h"
#include "solver_status.h"
#include "krylov_reductions.h"
#include "null_space.h"
#include "preconditioner.h"
#include <cstddef>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <vector>

namespace cfdx {
namespace core {

enum class CGResidualReplacementPolicy : std::uint8_t {
    Adaptive,
    Disabled,
    Periodic
};

// Résout A x = b par la méthode du gradient conjugué (CG).
//
// Args:
//   A          : matrice symétrique définie positive.
//   b          : vecteur second membre.
//   x          : solution (initialisé à l'entrée, modifié en sortie).
//   max_iter   : nombre maximum d'itérations.
//   tolerance  : tolérance sur le résidu relatif.
//
// Retourne un SolverResult.
namespace detail {
inline SolverResult solve_cg_impl(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
    std::size_t max_iter,
    double tolerance,
    PrecisionPolicy precision,
    KrylovReductionPolicy reduction,
    Preconditioner* preconditioner,
    const NullSpaceProjector* null_space,
    bool setup_preconditioner,
    CGResidualReplacementPolicy residual_policy = CGResidualReplacementPolicy::Adaptive,
    std::size_t periodic_replacement_interval = 0)
{
    SolverResult result;

    if (A.n_rows() != A.n_cols()) {
        result.status = SolverStatus::NOT_APPLICABLE;
        return result;
    }
    if (b.size() != A.n_rows()) {
        result.status = SolverStatus::NOT_APPLICABLE;
        return result;
    }
    if (x.size() != A.n_cols()) {
        result.status = SolverStatus::NOT_APPLICABLE;
        return result;
    }

    if (null_space &&
        (null_space->dimension() != A.n_rows() ||
         !null_space->is_compatible(b) ||
         !null_space->is_null_space(A))) {
        result.status = SolverStatus::NOT_APPLICABLE;
        return result;
    }

    const std::size_t n = A.n_rows();
    for (std::size_t k = 0; k < A.nnz(); ++k)
        if (!std::isfinite(A.values_data()[k])) {
            result.status = SolverStatus::DIVERGED;
            return result;
        }
    for (std::size_t i = 0; i < b.size(); ++i)
        if (!std::isfinite(b(i)) || !std::isfinite(x(i))) {
            result.status = SolverStatus::DIVERGED;
            return result;
        }
    if (null_space) null_space->remove(x);
    // A valid system can fail before the usual convergence/final-residual
    // reporting path (preconditioner failure, non-positive curvature, or
    // rejected residual replacement). Preserve its current b-A*x diagnostic
    // instead of returning SolverResult's default zero residuals. This helper
    // is used only on failure; it does not change an iterate or a stopping gate.
    const auto not_applicable = [&](std::size_t iterations) {
        result.status = SolverStatus::NOT_APPLICABLE;
        result.iterations = iterations;
        Vector failure_residual(n);
        mixed_precision_true_residual(A, b, x, failure_residual);
        result.residual =
            krylov_norm2(failure_residual, SolverPrecision::FP64, reduction);
        const double rhs_norm =
            krylov_norm2(b, SolverPrecision::FP64, reduction);
        result.residual_relative = rhs_norm > 0.0
            ? result.residual / rhs_norm
            : (result.residual == 0.0
                ? 0.0 : std::numeric_limits<double>::infinity());
        return result;
    };
    const double* Av = A.values_data();
    const bool mp32 = precision.enabled && precision.operator_precision == SolverPrecision::FP32;
    const SolverPrecision redp = precision.enabled ? precision.reduction_precision : SolverPrecision::FP64;
    const auto* Ac = A.columns_data();
    const auto* Ar = A.row_offsets_data();

    // Preserve the historical Jacobi path when no external preconditioner is
    // supplied. This keeps the native solver behavior and API unchanged.
    std::vector<double> M(n, 0.0);
    if (preconditioner) {
        if (setup_preconditioner && !preconditioner->setup(A)) {
            return not_applicable(0);
        }
    } else {
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t k = Ar[i]; k < Ar[i + 1]; ++k) {
                if (Ac[k] == static_cast<std::uint32_t>(i)) {
                    M[i] = Av[k];
                    break;
                }
            }
            if (!(M[i] > 0.0) || !std::isfinite(M[i])) {
                return not_applicable(0);
            }
        }
    }

    // r = b - A x
    std::vector<double> r(n);
    Vector rv(n); mixed_precision_true_residual(A,b,x,rv);
    if (null_space) null_space->remove(rv);
    for (std::size_t i=0;i<n;++i) r[i]=rv(i);

    // z = M^{-1} r
    std::vector<double> z(n);
    Vector z_vector(n);
    Vector r_vector(n);
    const auto apply_preconditioner = [&]() {
        if (preconditioner) {
            for (std::size_t i = 0; i < n; ++i) r_vector(i) = r[i];
            if (!preconditioner->apply(r_vector, z_vector)) return false;
            for (std::size_t i = 0; i < n; ++i) z[i] = z_vector(i);
        } else {
            for (std::size_t i = 0; i < n; ++i) z[i] = r[i] / M[i];
        }
        if (null_space) null_space->remove(z);
        for (const double value : z)
            if (!std::isfinite(value)) return false;
        return true;
    };
    if (!apply_preconditioner()) {
        return not_applicable(0);
    }

    // p = z
    std::vector<double> p = z;

    for (std::size_t i = 0; i < n; ++i) z_vector(i) = z[i];
    for (std::size_t i = 0; i < n; ++i) r_vector(i) = r[i];
    double rsold = krylov_dot(r_vector, z_vector, redp, reduction);

    const double b_norm = krylov_norm2(b, SolverPrecision::FP64, reduction);
    const double tol_abs = tolerance * std::max(b_norm, 1e-15);

    if (!std::isfinite(rsold)) { result.status=SolverStatus::DIVERGED; return result; }

    const double initial_true_residual =
        krylov_norm2(r_vector, SolverPrecision::FP64, reduction);
    if (preconditioner) {
        result.diagnostics.push_back(
            {0, initial_true_residual, initial_true_residual,
             rsold, 0.0, 0.0, 0.0});
    }
    if ((preconditioner && initial_true_residual <= tol_abs) ||
        (!preconditioner && rsold < tol_abs * tol_abs)) {
        result.status = SolverStatus::CONVERGED;
        result.iterations = 0;
        result.residual = preconditioner
            ? initial_true_residual
            : krylov_norm2(rv, redp, reduction);
        result.residual_relative = (b_norm > 0.0) ? result.residual / b_norm : 0.0;
        return result;
    }

    // Residual replacement only pays for itself while b-Ax is still
    // decreasing. Once the iteration reaches the attainable accuracy of the
    // system -- about eps*cond(A) in relative terms, which no solver, direct or
    // iterative, can improve upon in FP64 -- the recursive residual keeps
    // shrinking while b-Ax does not. The gap criterion below is then satisfied
    // at every remaining iteration, and the restart that follows each
    // replacement discards the Krylov direction, degrading CG into
    // preconditioned steepest descent for the rest of the budget. Requiring
    // measurable progress since the previous replacement bounds the number of
    // replacements to O(log(r0/attainable_floor)).
    double residual_at_last_replacement = std::numeric_limits<double>::infinity();

    for (std::size_t iter = 1; iter <= max_iter; ++iter) {
        // Ap = A p
        std::vector<double> Ap(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (std::size_t k = Ar[i]; k < Ar[i + 1]; ++k) {
                sum += mp32 ? static_cast<double>(static_cast<float>(Av[k])*static_cast<float>(p[Ac[k]])) : Av[k]*p[Ac[k]];
            }
            Ap[i] = sum;
        }
        if (null_space) null_space->remove(Ap);

        // alpha = rsold / (p · Ap)
        Vector p_vector(n), Ap_vector(n);
        for (std::size_t i = 0; i < n; ++i) {
            p_vector(i) = p[i];
            Ap_vector(i) = Ap[i];
        }
        const double pAp = krylov_dot(p_vector, Ap_vector, redp, reduction);
        if (!(pAp > 0.0) || !std::isfinite(pAp)) {
            return not_applicable(iter);
        }
        const double alpha = rsold / pAp;
        if (!std::isfinite(alpha)) { result.status=SolverStatus::DIVERGED; result.iterations=iter; return result; }

        // x = x + alpha p
        // r = r - alpha Ap
        for (std::size_t i = 0; i < n; ++i) {
            x(i) += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }
        if (null_space) {
            null_space->remove(x);
            null_space->remove(r);
        }

        for (std::size_t i = 0; i < n; ++i)
            if (!std::isfinite(x(i)) || !std::isfinite(r[i])) { result.status=SolverStatus::DIVERGED; result.iterations=iter; return result; }

        // z = M^{-1} r
        if (!apply_preconditioner()) {
            return not_applicable(iter);
        }

        // rsnew = r · z
        for (std::size_t i = 0; i < n; ++i) {
            r_vector(i) = r[i];
            z_vector(i) = z[i];
        }
        const double rsnew = krylov_dot(r_vector, z_vector, redp, reduction);

        if (!std::isfinite(rsnew)) { result.status=SolverStatus::DIVERGED; result.iterations=iter; return result; }
        // A preconditioned norm is not a true residual norm. Re-evaluate
        // b-A*x in FP64 for the external-preconditioner path before deciding
        // convergence, just as the reported final residual is evaluated.
        const double res = preconditioner
            ? mixed_precision_true_residual(A, b, x, rv)
            : std::sqrt(std::abs(rsnew));
        if (preconditioner &&
            (iter == 1 || iter % 100 == 0 || iter == max_iter)) {
            const double recursive_residual =
                std::sqrt(std::max(0.0, krylov_dot(r_vector, r_vector,
                                                   SolverPrecision::FP64,
                                                   reduction)));
            result.diagnostics.push_back(
                {iter, res, recursive_residual, rsnew, pAp, alpha,
                 rsold != 0.0 ? rsnew / rsold : 0.0});
        }
        if (res < tol_abs) {
            result.status = SolverStatus::CONVERGED;
            result.iterations = iter;
            result.residual = res;
            result.residual_relative = (b_norm > 0.0) ? res / b_norm : 0.0;
            return result;
        }

        // The recursive CG residual can drift away from b-Ax in finite
        // precision, especially once an effective AMG preconditioner has
        // reduced the true residual by many orders of magnitude. Continuing
        // with the stale recurrence can eventually destroy conjugacy and make
        // p^T A p appear non-positive even though A is SPD. Replacing the
        // recursive residual when the gap becomes material is a standard
        // residual-replacement safeguard; the true residual is already
        // available here, so this adds no extra fine-grid matvec.
        if (preconditioner) {
            double recursive_r2 = 0.0;
            for (std::size_t i = 0; i < n; ++i)
                recursive_r2 += r[i] * r[i];
            const double recursive_res = std::sqrt(recursive_r2);
            const double gap = std::abs(res - recursive_res);
            const double scale = std::max(res, recursive_res);
            result.max_true_recursive_gap =
                std::max(result.max_true_recursive_gap, gap);
            if (recursive_res > 0.0 && std::isfinite(recursive_res)) {
                const double ratio = res / recursive_res;
                if (std::isfinite(ratio)) {
                    result.min_true_recursive_ratio =
                        std::min(result.min_true_recursive_ratio, ratio);
                    result.max_true_recursive_ratio =
                        std::max(result.max_true_recursive_ratio, ratio);
                }
            }
            result.max_true_residual = std::max(result.max_true_residual, res);
            result.min_true_residual = std::min(result.min_true_residual, res);
            constexpr double residual_gap_fraction = 0.25;
            const bool made_progress =
                res < 0.5 * residual_at_last_replacement;
            // A non-positive r^T z means the recursive residual has been
            // annihilated by cancellation. z = M^{-1}r is then numerically
            // zero, so the next search direction is zero and p^T A p vanishes.
            // That is a dead recurrence, not an unusable preconditioner, and it
            // must be recognized here: otherwise it surfaces one iteration
            // later as a spurious NOT_APPLICABLE from the p^T A p guard, which
            // wrongly tells the caller to abandon a perfectly valid
            // preconditioner.
            const bool recurrence_collapsed = !(rsnew > 0.0);
            const bool adaptive_replacement =
                residual_policy == CGResidualReplacementPolicy::Adaptive &&
                std::isfinite(recursive_res) && scale > 0.0 &&
                gap > residual_gap_fraction * scale && made_progress;
            const bool periodic_replacement =
                residual_policy == CGResidualReplacementPolicy::Periodic &&
                periodic_replacement_interval > 0 &&
                (iter % periodic_replacement_interval == 0);
            // Repairing a collapsed recurrence is itself a replacement, so the
            // Disabled policy must not do it, and it is pointless once b-Ax has
            // stopped improving.
            const bool repair_collapse =
                recurrence_collapsed && made_progress &&
                residual_policy != CGResidualReplacementPolicy::Disabled;
            if (recurrence_collapsed && !repair_collapse &&
                !periodic_replacement) {
                // The recurrence is dead and replacing it cannot help, so the
                // iteration has extracted all the accuracy this system allows.
                // Report the accuracy actually achieved, at the iteration where
                // it was reached, instead of spending the remaining budget on a
                // zero search direction.
                result.status = SolverStatus::MAX_ITER_REACHED;
                result.iterations = iter;
                result.residual = res;
                result.residual_relative = (b_norm > 0.0) ? res / b_norm : 0.0;
                return result;
            }
            if (adaptive_replacement || periodic_replacement ||
                repair_collapse) {
                ++result.residual_replacements;
                residual_at_last_replacement = res;
                if (result.first_residual_replacement == 0)
                    result.first_residual_replacement = iter;
                result.last_residual_replacement = iter;
                for (std::size_t i = 0; i < n; ++i) r[i] = rv(i);
                if (null_space) null_space->remove(r);
                if (!apply_preconditioner()) {
                    return not_applicable(iter);
                }
                for (std::size_t i = 0; i < n; ++i) {
                    r_vector(i) = r[i];
                    z_vector(i) = z[i];
                }
                const double replaced_rsnew =
                    krylov_dot(r_vector, z_vector, redp, reduction);
                if (!std::isfinite(replaced_rsnew) || !(replaced_rsnew > 0.0)) {
                    return not_applicable(iter);
                }
                // The replacement recomputes b-Ax after the normal convergence
                // check. Re-test the true residual here so a converged iterate
                // produced exactly at a replacement is not reported as
                // MAX_ITER_REACHED merely because the replacement was the last
                // operation of the iteration.
                if (res <= tol_abs) {
                    result.status = SolverStatus::CONVERGED;
                    result.iterations = iter;
                    result.residual = res;
                    result.residual_relative = (b_norm > 0.0) ? res / b_norm : 0.0;
                    return result;
                }
                // Restart the search direction after replacement. This avoids
                // carrying a direction built from a residual that no longer
                // represents b-Ax while retaining the current iterate.
                p = z;
                rsold = replaced_rsnew;
                continue;
            }
        }

        // beta = rsnew / rsold
        const double beta = rsnew / rsold;

        // p = z + beta p
        for (std::size_t i = 0; i < n; ++i) {
            p[i] = z[i] + beta * p[i];
        }
        if (null_space) null_space->remove(p);

        rsold = rsnew;
    }

    result.status = SolverStatus::MAX_ITER_REACHED;
    result.iterations = max_iter;
    result.residual = mixed_precision_true_residual(A,b,x,rv);
    result.residual_relative = (b_norm > 0.0) ? result.residual / b_norm : 0.0;
    return result;
}
} // namespace detail

// Diagnostic/qualification entry point: keeps the production solve_cg API
// unchanged while allowing V&V campaigns to isolate residual-replacement
// behavior without changing the stopping criterion.
inline SolverResult solve_cg_controlled(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
    Preconditioner& preconditioner,
    std::size_t max_iter,
    double tolerance,
    CGResidualReplacementPolicy residual_policy,
    std::size_t periodic_replacement_interval = 0)
{
    return detail::solve_cg_impl(
        A, b, x, max_iter, tolerance, {}, {},
        &preconditioner, nullptr, true,
        residual_policy, periodic_replacement_interval);
}

inline SolverResult solve_cg(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
    std::size_t max_iter = 1000,
    double tolerance = 1e-12,
    PrecisionPolicy precision = {},
    KrylovReductionPolicy reduction = {})
{
    return detail::solve_cg_impl(
        A, b, x, max_iter, tolerance, precision, reduction,
        nullptr, nullptr, true);
}

inline SolverResult solve_cg(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
    Preconditioner& preconditioner,
    std::size_t max_iter = 1000,
    double tolerance = 1e-12,
    PrecisionPolicy precision = {},
    KrylovReductionPolicy reduction = {})
{
    return detail::solve_cg_impl(
        A, b, x, max_iter, tolerance, precision, reduction,
        &preconditioner, nullptr, true);
}

inline SolverResult solve_cg(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
    const NullSpaceProjector& null_space,
    std::size_t max_iter = 1000,
    double tolerance = 1e-12,
    PrecisionPolicy precision = {},
    KrylovReductionPolicy reduction = {})
{
    return detail::solve_cg_impl(
        A, b, x, max_iter, tolerance, precision, reduction,
        nullptr, &null_space, true);
}

inline SolverResult solve_cg(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
    Preconditioner& preconditioner,
    const NullSpaceProjector& null_space,
    std::size_t max_iter = 1000,
    double tolerance = 1e-12,
    PrecisionPolicy precision = {},
    KrylovReductionPolicy reduction = {})
{
    return detail::solve_cg_impl(
        A, b, x, max_iter, tolerance, precision, reduction,
        &preconditioner, &null_space, true);
}

// Solve with a preconditioner that has already been prepared for A. This is
// the low-level path used by reusable nonlinear/pressure-solver contexts.
inline SolverResult solve_cg_prepared(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
    Preconditioner& preconditioner,
    std::size_t max_iter = 1000,
    double tolerance = 1e-12,
    PrecisionPolicy precision = {},
    KrylovReductionPolicy reduction = {},
    const NullSpaceProjector* null_space = nullptr)
{
    return detail::solve_cg_impl(
        A, b, x, max_iter, tolerance, precision, reduction,
        &preconditioner, null_space, false);
}

}  // namespace core
}  // namespace cfdx

```

### `src/cfdx/core/mesh/boundary.cpp`

SHA256 `58dbf7df4109dcb44f7c2c55a5b6f1bc14e78b076243cb486dca4a509f7584bd`

```cpp
// M0.1-T05 — Boundary patches implementation

#include "cfdx/core/mesh/boundary.h"
#include "cfdx/core/mesh/mesh.h"

namespace cfdx {
namespace core {

bool BoundaryPatches::is_consistent_with_mesh(const Mesh& m) const {
    if (!is_consistent(m.n_faces())) return false;

    const FaceOwnership& own = m.ownership();
    std::vector<bool> covered(m.n_faces(), false);

    // is_consistent() has already rejected duplicate/out-of-range ids.
    // Mark each listed boundary face once, retaining the internal-face gate.
    for (const auto& p : patches_) {
        for (const auto fid : p.face_ids) {
            if (own.neighbour(fid) != FaceOwnership::BOUNDARY) {
                return false;
            }
            covered[fid] = true;
        }
    }

    // Membership must be tested in O(1), not by rescanning every patch for
    // every boundary face. In a one-cell-thick 2-D extrusion, the two EMPTY
    // patches contain O(n_cells) faces and the former search was quadratic.
    for (std::size_t f = 0; f < m.n_faces(); ++f) {
        if (own.neighbour(f) == FaceOwnership::BOUNDARY && !covered[f]) {
            return false;
        }
    }

    return true;
}

}  // namespace core
}  // namespace cfdx

```

### `src/cfdx/core/numerics/temporal.h`

SHA256 `c59e0a7169622fe0289b74da08806168865869dfe22db7b0f64c7e646f1f59ea`

```cpp
// M0.13 — Temporal Discretization
//
// Spécification CFDX v0.7 §33-35 :
//   Schémas d'intégration temporelle pour ∂φ/∂t + ∇·F = S
//
//   Schémas supportés :
//     - Euler explicite : φ^{n+1} = φ^n + Δt * RHS(φ^n)
//     - Euler implicite : φ^{n+1} = φ^n + Δt * RHS(φ^{n+1})
//     - Crank-Nicolson : φ^{n+1} = φ^n + 0.5*Δt * (RHS(φ^n) + RHS(φ^{n+1}))
//     - BDF2 : (3φ^{n+1} - 4φ^n + φ^{n-1}) / (2Δt) = RHS(φ^{n+1})
//
//   Local time stepping : Δt_local = CFL * h / |u| pour convergence régime permanent
//   Contrôle adaptatif : Δt^{n+1} = Δt^n * min(max_factor, max(min_factor, CFL_target / CFL_current))

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/index_types.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include <cstddef>
#include <functional>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>

namespace cfdx {
namespace core {

enum class TimeScheme : std::uint8_t {
    EULER_EXPLICIT = 0,
    EULER_IMPLICIT,
    CRANK_NICOLSON,
    BDF2
};

inline const char* to_string(TimeScheme s) {
    switch (s) {
        case TimeScheme::EULER_EXPLICIT:  return "euler_explicit";
        case TimeScheme::EULER_IMPLICIT:  return "euler_implicit";
        case TimeScheme::CRANK_NICOLSON:  return "crank_nicolson";
        case TimeScheme::BDF2:            return "bdf2";
        default:                           return "unknown";
    }
}

inline TimeScheme time_scheme_from_string(const std::string& s) {
    if (s == "euler_explicit")  return TimeScheme::EULER_EXPLICIT;
    if (s == "euler_implicit")  return TimeScheme::EULER_IMPLICIT;
    if (s == "crank_nicolson")  return TimeScheme::CRANK_NICOLSON;
    if (s == "bdf2")            return TimeScheme::BDF2;
    throw std::runtime_error("time_scheme_from_string: unknown scheme '" + s + "'");
}

// RHS function type: computes spatial residual RHS = -∇·F + S
// Input: current field φ, Output: residual (size = n_cells, dim = field dimension)
using RhsFunction = std::function<void(const Field<double, Location::CELL>&, Field<double, Location::CELL>&)>;

// Time integration context - holds history for multi-step schemes
struct TimeIntegrationContext {
    std::size_t n_cells = 0;
    std::size_t field_dim = 0;

    // For BDF2: stores φ^{n-1} and φ^n
    Field<double, Location::CELL> phi_prev;  // φ^{n-1}
    Field<double, Location::CELL> phi_curr;  // φ^n

    bool has_prev = false;  // true after first step

    // Previous time-step size used by BDF2 for variable-step integration.
    // 0 means "not yet set": the constant-step formula (equal to ω = 1) applies.
    double dt_prev = 0.0;

    TimeIntegrationContext() = default;

    TimeIntegrationContext(std::size_t n, std::size_t dim, const std::string& name)
        : n_cells(n), field_dim(dim),
          phi_prev(n, name + "_prev", "1", dim),
          phi_curr(n, name + "_curr", "1", dim)
    {}

    void initialize(const Field<double, Location::CELL>& phi0) {
        phi_curr = phi0;  // copy
        has_prev = false;
    }

    // Shift history: φ^{n-1} = φ^n, φ^n = φ^{n+1}
    void shift(const Field<double, Location::CELL>& phi_new) {
        phi_prev = phi_curr;
        phi_curr = phi_new;
        has_prev = true;
    }
};

// Variable-step BDF2 coefficients. With ω = dt / dt_prev:
//   φ^{n+1} = ( a0 φ^n - a1 φ^{n-1} + b dt f^{n+1} ) / (1 + 2ω)
//   a0 = (1+ω)^2,  a1 = ω^2,  b = (1+ω)
// For ω = 1 this reduces to the constant-step form (4/3, 1/3, 2/3).
inline void bdf2_coefficients(double dt_prev, double dt,
                              double& a0, double& a1, double& denom, double& b)
{
    const double dp = (dt_prev > 0.0) ? dt_prev : dt;   // unset history -> uniform
    const double w = dt / dp;
    a0 = (1.0 + w) * (1.0 + w);
    a1 = w * w;
    denom = 1.0 + 2.0 * w;
    b = (1.0 + w) * dt;
}

// Advance field in time by one step
// Returns the new field φ^{n+1}
inline Field<double, Location::CELL> advance_time(
    const Field<double, Location::CELL>& phi,
    double dt,
    const RhsFunction& rhs_func,
    TimeScheme scheme,
    TimeIntegrationContext* ctx = nullptr)
{
    if (!(dt > 0.0) || !std::isfinite(dt)) {
        throw std::invalid_argument("advance_time: dt must be finite and strictly positive");
    }
    if (!rhs_func) {
        throw std::invalid_argument("advance_time: rhs function must be valid");
    }

    const std::size_t n_cells = phi.size();
    const std::size_t dim = phi.dimension();
    if (n_cells == 0 || dim == 0) {
        throw std::invalid_argument("advance_time: field must have non-zero size and dimension");
    }
    if (ctx && (ctx->n_cells != n_cells || ctx->field_dim != dim)) {
        throw std::invalid_argument("advance_time: time-integration context does not match field shape");
    }
    const auto require_finite = [&](const Field<double, Location::CELL>& field,
                                    const char* message) {
        for (std::size_t d = 0; d < dim; ++d) {
            for (std::size_t cell = 0; cell < n_cells; ++cell) {
                if (!std::isfinite(field.component_data(d)[cell])) {
                    throw std::runtime_error(message);
                }
            }
        }
    };
    require_finite(phi, "advance_time: input field contains a non-finite value");
    Field<double, Location::CELL> phi_new(
        n_cells, "phi_new", phi.metadata().unit, dim);

    Field<double, Location::CELL> rhs(
        n_cells, "rhs", phi.metadata().unit + "/s", dim);
    Field<double, Location::CELL> rhs_iter(
        n_cells, "rhs_iter", phi.metadata().unit + "/s", dim);

    rhs_func(phi, rhs);
    for (std::size_t d = 0; d < dim; ++d) {
        const double* rhs_data = rhs.component_data(d);
        for (std::size_t cell = 0; cell < n_cells; ++cell) {
            if (!std::isfinite(rhs_data[cell])) {
                throw std::runtime_error("advance_time: RHS contains a non-finite value");
            }
        }
    }

    constexpr int max_iterations = 100;
    constexpr double tolerance = 1e-12;

    auto solve_fixed_point = [&](auto&& predictor) {
        predictor();
        for (int iteration = 0; iteration < max_iterations; ++iteration) {
            rhs_func(phi_new, rhs_iter);
            require_finite(rhs_iter, "advance_time: implicit RHS contains a non-finite value");
            double max_delta = 0.0;
            for (std::size_t d = 0; d < dim; ++d) {
                const double* old_data = phi_new.component_data(d);
                const double* phi_data = phi.component_data(d);
                const double* rhs_data = rhs_iter.component_data(d);
                double* new_data = phi_new.component_data(d);
                (void)old_data;
                for (std::size_t cell = 0; cell < n_cells; ++cell) {
                    double candidate = 0.0;
                    if (scheme == TimeScheme::EULER_IMPLICIT) {
                        candidate = phi_data[cell] + dt * rhs_data[cell];
                    } else if (scheme == TimeScheme::CRANK_NICOLSON) {
                        const double* rhs_old = rhs.component_data(d);
                        candidate = phi_data[cell] +
                            0.5 * dt * (rhs_old[cell] + rhs_data[cell]);
                    } else {
                        const double* prev = ctx->phi_prev.component_data(d);
                        double a0 = 0.0, a1 = 0.0, denom = 0.0, b = 0.0;
                        bdf2_coefficients(ctx->dt_prev, dt, a0, a1, denom, b);
                        candidate = (a0 * phi_data[cell] - a1 * prev[cell] +
                             b * rhs_data[cell]) / denom;
                    }
                    if (!std::isfinite(candidate)) {
                        throw std::runtime_error("advance_time: implicit candidate contains a non-finite value");
                    }
                    max_delta = std::max(max_delta, std::abs(candidate - new_data[cell]));
                    new_data[cell] = candidate;
                }
            }
            if (max_delta <= tolerance) {
                return;
            }
        }
        throw std::runtime_error("advance_time: implicit iteration did not converge");
    };

    switch (scheme) {
        case TimeScheme::EULER_EXPLICIT:
            for (std::size_t d = 0; d < dim; ++d) {
                const double* phi_data = phi.component_data(d);
                const double* rhs_data = rhs.component_data(d);
                double* new_data = phi_new.component_data(d);
                for (std::size_t cell = 0; cell < n_cells; ++cell) {
                    new_data[cell] = phi_data[cell] + dt * rhs_data[cell];
                }
            }
            break;

        case TimeScheme::EULER_IMPLICIT:
            solve_fixed_point([&] {
                for (std::size_t d = 0; d < dim; ++d) {
                    const double* phi_data = phi.component_data(d);
                    double* new_data = phi_new.component_data(d);
                    for (std::size_t cell = 0; cell < n_cells; ++cell) {
                        new_data[cell] = phi_data[cell];
                    }
                }
            });
            break;

        case TimeScheme::CRANK_NICOLSON:
            solve_fixed_point([&] {
                for (std::size_t d = 0; d < dim; ++d) {
                    const double* phi_data = phi.component_data(d);
                    const double* rhs_data = rhs.component_data(d);
                    double* new_data = phi_new.component_data(d);
                    for (std::size_t cell = 0; cell < n_cells; ++cell) {
                        new_data[cell] = phi_data[cell] + dt * rhs_data[cell];
                    }
                }
            });
            break;

        case TimeScheme::BDF2:
            if (!ctx || !ctx->has_prev) {
                // Bootstrap BDF2 with a genuinely implicit Euler step.
                for (std::size_t d = 0; d < dim; ++d) {
                    const double* phi_data = phi.component_data(d);
                    double* new_data = phi_new.component_data(d);
                    for (std::size_t cell = 0; cell < n_cells; ++cell) {
                        new_data[cell] = phi_data[cell];
                    }
                }
                for (int iteration = 0; iteration < max_iterations; ++iteration) {
                    rhs_func(phi_new, rhs_iter);
                    require_finite(rhs_iter, "advance_time: BDF2 bootstrap RHS contains a non-finite value");
                    double max_delta = 0.0;
                    for (std::size_t d = 0; d < dim; ++d) {
                        const double* phi_data = phi.component_data(d);
                        const double* rhs_data = rhs_iter.component_data(d);
                        double* new_data = phi_new.component_data(d);
                        for (std::size_t cell = 0; cell < n_cells; ++cell) {
                            const double candidate = phi_data[cell] + dt * rhs_data[cell];
                            if (!std::isfinite(candidate)) {
                                throw std::runtime_error("advance_time: BDF2 bootstrap candidate contains a non-finite value");
                            }
                            max_delta = std::max(max_delta, std::abs(candidate - new_data[cell]));
                            new_data[cell] = candidate;
                        }
                    }
                    if (max_delta <= tolerance) break;
                    if (iteration == max_iterations - 1) {
                        throw std::runtime_error("advance_time: implicit Euler bootstrap did not converge");
                    }
                }
            } else {
                solve_fixed_point([&] {});
            }
            break;
    }

    // Reject arithmetic overflow and invalid iterates before committing history.
    require_finite(phi_new, "advance_time: output field contains a non-finite value");
    if (ctx) {
        ctx->dt_prev = dt;
        ctx->shift(phi_new);
    }
    return phi_new;
}


// Compute stable time step based on CFL condition
// CFL = |u| * Δt / h ≤ CFL_max
// For convection: Δt = CFL_max * h / (|u| + c) where c = speed of sound (if compressible)
inline double compute_cfl_time_step(
    const Mesh& mesh,
    const Field<double, Location::CELL>& velocity,  // dim=3
    double cfl_max = 0.5,
    double speed_of_sound = 0.0)
{
    const std::size_t n_cells = mesh.n_cells();
    const double* ux = velocity.component_data(0);
    const double* uy = velocity.component_data(1);
    const double* uz = velocity.component_data(2);

    // Compute cell characteristic length h = V^(1/3) (or similar)
    // Need cell volumes
    std::vector<Vec3> face_centres(mesh.n_faces());
    std::vector<Vec3> face_Sf(mesh.n_faces());

    const PointCloud& pts = mesh.points();
    const double* px = pts.x_data();
    const double* py = pts.y_data();
    const double* pz = pts.z_data();
    const auto* face_verts = mesh.faces().vertices_data();
    const auto* face_offsets = mesh.faces().offsets_data();

    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const Offset off = face_offsets[f];
        const Offset n = face_offsets[f + 1] - off;
        const FaceGeometry fg = compute_face_geometry(px, py, pz, face_verts, off, n);
        face_centres[f] = fg.centre;
        face_Sf[f] = fg.Sf;
    }

    std::vector<double> cell_volume(n_cells, 0.0);
    const CellConnectivity& cells = mesh.cells();
    const auto* cell_faces = cells.faces_data();
    const auto* cell_offsets = cells.offsets_data();

    for (std::size_t c = 0; c < n_cells; ++c) {
        const Offset off = cell_offsets[c];
        const Offset n = cell_offsets[c + 1] - off;
        const CellGeometry cg = compute_cell_geometry_oriented(
            face_centres.data(), face_Sf.data(), cell_faces + off, n,
            static_cast<CellIndex>(c), mesh.ownership());
        cell_volume[c] = cg.volume;
    }

    double dt_min = 1e30;
    for (std::size_t c = 0; c < n_cells; ++c) {
        if (cell_volume[c] <= 0) continue;

        const double h = std::cbrt(cell_volume[c]);  // characteristic length
        const double vel_mag = std::sqrt(ux[c]*ux[c] + uy[c]*uy[c] + uz[c]*uz[c]);
        const double wave_speed = vel_mag + speed_of_sound;

        if (wave_speed > 1e-12) {
            double dt_local = cfl_max * h / wave_speed;
            dt_min = std::min(dt_min, dt_local);
        }
    }

    return dt_min;
}

// Adaptive time step controller
struct AdaptiveTimeStepper {
    double dt = 1e-3;          // current time step
    double dt_min = 1e-12;     // minimum allowed
    double dt_max = 1.0;       // maximum allowed
    double cfl_target = 0.5;   // target CFL
    double cfl_current = 0.0;  // current CFL (computed externally)
    double safety_factor = 0.9;
    double growth_factor = 1.2;
    double shrink_factor = 0.5;

    // Update time step based on current CFL
    void update(double cfl_achieved) {
        cfl_current = cfl_achieved;
        if (cfl_achieved > 0) {
            double ratio = cfl_target / cfl_achieved;
            double new_dt = dt * safety_factor * ratio;
            // Limit growth/shrink
            if (new_dt > dt * growth_factor) new_dt = dt * growth_factor;
            if (new_dt < dt * shrink_factor) new_dt = dt * shrink_factor;
            dt = std::clamp(new_dt, dt_min, dt_max);
        }
    }

    // Called when step is rejected (e.g., non-convergence)
    void reject() {
        dt = std::max(dt * shrink_factor, dt_min);
    }
};

// Local time stepping for steady-state acceleration
// Δt_local = CFL * h / |u| for each cell
inline std::vector<double> compute_local_time_steps(
    const Mesh& mesh,
    const Field<double, Location::CELL>& velocity,
    double cfl = 0.5)
{
    const std::size_t n_cells = mesh.n_cells();
    const double* ux = velocity.component_data(0);
    const double* uy = velocity.component_data(1);
    const double* uz = velocity.component_data(2);

    // Compute cell volumes
    std::vector<Vec3> face_centres(mesh.n_faces());
    std::vector<Vec3> face_Sf(mesh.n_faces());

    const auto* face_verts = mesh.faces().vertices_data();
    const auto* face_offsets = mesh.faces().offsets_data();
    (void)face_verts; // used in compute_face_geometry

    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const Offset off = face_offsets[f];
        const Offset n = face_offsets[f + 1] - off;
        const FaceGeometry fg = compute_face_geometry(
            mesh.points().x_data(), mesh.points().y_data(), mesh.points().z_data(), face_verts, off, n);
        face_centres[f] = fg.centre;
        face_Sf[f] = fg.Sf;
    }

    std::vector<double> cell_volume(n_cells, 0.0);
    const CellConnectivity& cells = mesh.cells();
    const auto* cell_faces = cells.faces_data();
    const auto* cell_offsets = cells.offsets_data();
    (void)cell_faces; // used in compute_cell_geometry

    for (std::size_t c = 0; c < n_cells; ++c) {
        const Offset off = cell_offsets[c];
        const Offset n = cell_offsets[c + 1] - off;
        (void)n; // used in compute_cell_geometry
        const CellGeometry cg = compute_cell_geometry_oriented(
            face_centres.data(), face_Sf.data(), cell_faces + off, n,
            static_cast<CellIndex>(c), mesh.ownership());
        cell_volume[c] = cg.volume;
    }

    // Large initial dt for cells with zero velocity
    const double dt_max = 1e30;
    std::vector<double> dt_local(n_cells, dt_max);
    for (std::size_t c = 0; c < n_cells; ++c) {
        if (cell_volume[c] <= 0) continue;
        const double h = std::cbrt(cell_volume[c]);
        const double vel_mag = std::sqrt(ux[c]*ux[c] + uy[c]*uy[c] + uz[c]*uz[c]);
        if (vel_mag > 1e-12) {
            dt_local[c] = cfl * h / vel_mag;
        }
    }
    return dt_local;
}

// Apply local time stepping: φ^{n+1}_c = φ^n_c + Δt_c * RHS_c
inline void apply_local_time_stepping(
    const Field<double, Location::CELL>& phi,
    const std::vector<double>& dt_local,
    const RhsFunction& rhs_func,
    Field<double, Location::CELL>& phi_new)
{
    const std::size_t n_cells = phi.size();
    const std::size_t dim = phi.dimension();

    Field<double, Location::CELL> rhs(n_cells, "rhs", phi.metadata().unit + "/s", dim);
    rhs_func(phi, rhs);

    for (std::size_t c = 0; c < n_cells; ++c) {
        for (std::size_t d = 0; d < dim; ++d) {
            const double* phi_data = phi.component_data(d);
            const double* rhs_data = rhs.component_data(d);
            double* new_data = phi_new.component_data(d);
            new_data[c] = phi_data[c] + dt_local[c] * rhs_data[c];
        }
    }
}

}  // namespace core
}  // namespace cfdx

```

### `src/cfdx/physics/steady_incompressible_solver.h`

SHA256 `a9f2df548482d337c653a9ee0b0ee0d49155ab38c0786f203b113a29ad6f5d45`

```cpp
#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/linalg/bicgstab_solver.h"
#include "cfdx/core/linalg/coupled_amg_schur.h"
#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/linear_solver_context.h"
#include "cfdx/core/linalg/mgr_preconditioner.h"
#include "cfdx/core/linalg/coupled_amg_schur.h"
#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "cfdx/core/linalg/preconditioner.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/numerics/gradient.h"
#include "cfdx/io/restart/dat_restart.h"
#include "cfdx/physics/advanced_convergence.h"
#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/pressure_velocity_algorithms.h"
#include "cfdx/physics/pressure_velocity_system.h"
#include "cfdx/physics/solver_control.h"
#include "cfdx/physics/sst_production_coupling.h"
#include "cfdx/physics/timestep_control.h"
#include "cfdx/core/numerics/conservation.h"
#include "cfdx/core/numerics/temporal.h"
#include "cfdx/cfdx/core/utils/convergence_monitor.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <cstddef>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <numeric>
#include <stdexcept>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <functional>

namespace cfdx::physics {

struct VelocityBoundaryCondition {
    enum class Type { FIXED_VALUE, ZERO_GRADIENT };
    Type type = Type::ZERO_GRADIENT;
    cfdx::core::Vec3 value{0.0, 0.0, 0.0};
};

using VelocityBoundaryConditions = std::map<std::string, VelocityBoundaryCondition>;

enum class IncompressibleProbeField {
    PRESSURE,
    U_X,
    U_Y,
    U_Z,
    U_MAGNITUDE,
};

struct IncompressiblePointProbe {
    std::string name;
    cfdx::core::Vec3 location{0.0, 0.0, 0.0};
    IncompressibleProbeField field = IncompressibleProbeField::PRESSURE;
};

struct IncompressibleProbeSample {
    std::string name;
    std::size_t iteration = 0;
    double time = 0.0;
    double value = 0.0;
};

struct DiagnosticsControls {
    bool coupled_matrix_summary = false;
    bool freeze_state_probe = false;
    long debug_cell = -1;
    bool debug_cell_auto = false;
    // Bounded progress trace for long validation solves. Disabled by default.
    bool iteration_trace = false;
    std::size_t iteration_trace_frequency = 50;
};

struct IncompressibleTransientControls {
    bool enabled = false;
    cfdx::core::TimeScheme scheme = cfdx::core::TimeScheme::BDF2;
    double dt = 0.0;
    double dt_previous = 0.0;
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* previous = nullptr;
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* previous_previous = nullptr;
    bool history_valid = false;
};

// ---------------------------------------------------------------------------
// Initialization contract (N7)
// ---------------------------------------------------------------------------
// The initial state is part of the numerical method: it fixes the denominator
// of every relative residual reported by the solve, because
// result.reference_momentum_residual is taken from the first iteration. The
// strategy is therefore declared explicitly and validated rather than inferred
// from whether a restart path happens to be non-empty.
//
// Automatic selection is deliberately NOT implemented. Choosing a strategy by
// probing the filesystem or mesh quality would change the initial state, and
// therefore the relative-residual reference, without the resulting convergence
// being attributable to a declared numerical method.
enum class InitializationMode {
    // Caller-provided fields are used unchanged. The caller is responsible for
    // having produced a meaningful state.
    Provided,
    // Fields are reset to a uniform state before the solve.
    Uniform,
    // Fields are overwritten from a CFDX-DAT checkpoint. Requires
    // controls.initialization.restart_path to be non-empty.
    Restart
};

inline const char* to_string(InitializationMode mode) {
    switch (mode) {
        case InitializationMode::Provided: return "provided";
        case InitializationMode::Uniform: return "uniform";
        case InitializationMode::Restart: return "restart";
    }
    return "unknown";
}

// Parses an initialization mode name. Returns false for an unrecognised name so
// that callers reject it instead of silently defaulting.
inline bool parse_initialization_mode(const std::string& name, InitializationMode& out) {
    if (name == "provided") { out = InitializationMode::Provided; return true; }
    if (name == "uniform") { out = InitializationMode::Uniform; return true; }
    if (name == "restart") { out = InitializationMode::Restart; return true; }
    return false;
}

struct InitializationControls {
    InitializationMode mode = InitializationMode::Provided;
    // Uniform state applied when mode == Uniform.
    cfdx::core::Vec3 uniform_velocity{0.0, 0.0, 0.0};
    double uniform_pressure = 0.0;
    // Checkpoint path applied when mode == Restart.
    std::string restart_path;
};

// Read-only view of the exact coupled saddle-point system handed to the
// Krylov solver, exposed to the optional evidence hook in
// IncompressibleSolverControls. The hook fires at the Krylov call point -
// after assembly, gauge treatment, finite checks and preconditioner setup -
// so the aliased matrix, RHS and initial guess ARE the exact solve arguments
// (same objects, no copy, no further mutation). The references are only
// valid for the duration of the callback: a consumer that needs the data
// beyond the call must copy it.
struct CoupledAssemblySnapshot {
    const cfdx::core::SparseMatrix* matrix = nullptr;
    const cfdx::core::Vector* rhs = nullptr;
    // Production initial guess ([U_old; p_old] with the gauge row set to the
    // reference value), i.e. the x passed to the Krylov solver.
    const cfdx::core::Vector* initial_guess = nullptr;
    // Momentum-weighted Rhie-Chow pressure Schur diagonal, in both forms:
    // pinned (reference cell = 1.0, exactly as production leaves it at solve
    // time) and raw (as assembled, before the pin). Diagnostic evidence
    // channels; production preconditioning does not consume them.
    const std::vector<double>* pressure_schur_diagonal = nullptr;
    const std::vector<double>* pressure_schur_diagonal_raw = nullptr;
    // Production PCD pressure operators of this assembly, exactly as passed
    // to PcdSchurApproximation (reference row already imposed by the
    // assembler). Null when the resolved coupled plan is not PCD.
    const cfdx::core::SparseMatrix* pcd_pressure_mass = nullptr;
    const cfdx::core::SparseMatrix* pcd_pressure_laplacian = nullptr;
    const cfdx::core::SparseMatrix* pcd_pressure_convection_diffusion = nullptr;
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* velocity_old = nullptr;
    std::size_t reference_cell = 0;
    bool has_fixed_pressure = false;
};

struct SSTProductionRuntime {
    // State must be explicitly initialized from case data or restored from a
    // DAT checkpoint with the helpers in sst_production_coupling.h. The
    // production solver advances this state transactionally when non-null.
    SSTProductionState* state = nullptr;
    TurbulenceTransportControls transport;
    ScalarBoundaryConditions k_boundary_conditions;
    ScalarBoundaryConditions omega_boundary_conditions;
    // Optional per-face omega data, including the existing geometric wall-omega
    // contract. Kept owned by the runtime so the pointer passed to transport
    // remains valid for the duration of the update.
    ScalarBoundaryFaceValues omega_face_values;
    std::vector<std::string> wall_patch_names;
    std::size_t max_transport_iterations = 100;
    double transport_tolerance = 1e-8;
};

struct IncompressibleSolverControls {
    InitializationControls initialization;
    // Opt-in only. The production RANS selector remains fail-closed; supplying
    // this runtime explicitly opts the caller into the SST-coupled path.
    SSTProductionRuntime* sst_runtime = nullptr;
    PressureVelocityAlgorithm algorithm = PressureVelocityAlgorithm::SIMPLE;
    IncompressibleTransientControls transient;
    CouplingControls coupling;
    ConvergenceCriteria convergence;
    std::size_t linear_max_iterations = 1000;
    double linear_tolerance = 1e-10;
    ConvergenceAccelerationControls acceleration;
    double density = 1.0;
    double kinematic_viscosity = 1.0e-3;
    double turbulent_viscosity = 0.0;
    // Optional cellwise kinematic eddy viscosity (m^2/s). Replaces the scalar
    // value when supplied; the production RANS guard remains active.
    std::vector<double> cell_turbulent_viscosity;
    cfdx::core::Vec3 body_force{0.0, 0.0, 0.0};
    // Optional cellwise body force per unit mass (m/s^2). Replaces the
    // uniform controls.body_force vector when supplied (Boussinesq
    // buoyancy coupling, issue #885 step 4, fiche methode 10): the
    // assembly multiplies it by the constant density exactly like the
    // uniform channel, in both the momentum predictor and the
    // independent momentum conservation audit. Supplying a nonzero
    // uniform body_force together with the cellwise channel is rejected,
    // mirroring the cellwise turbulent-viscosity contract.
    std::vector<cfdx::core::Vec3> cell_body_force;
    PressureGaugePolicy pressure_gauge_policy = PressureGaugePolicy::REFERENCE_CELL;
    std::size_t pressure_reference_cell = 0;
    double pressure_reference_value = 0.0;
    cfdx::core::LinearSolverRequest momentum_linear_solver;
    cfdx::core::LinearSolverRequest pressure_linear_solver;
    cfdx::core::LinearSolverRequest coupled_linear_solver;
    // PCD pressure-side solves are nested inside the coupled preconditioner.
    // They deliberately have independent solver controls so the outer coupled
    // iteration budget can never become an accidental inner Krylov budget.
    cfdx::core::LinearSolverRequest pcd_laplacian_linear_solver{
        cfdx::core::KrylovModel::GMRES,
        cfdx::core::PreconditionerModel::NativeAMG,
        40,
        false,
        cfdx::core::NullSpaceModel::None};
    // Fp is the nested PCD pressure convection-diffusion operator. The
    // Ghia Re=100 qualification showed that FGMRES(40) repeatedly hit the
    // 100-iteration cap with true relative residuals above the unchanged
    // 1e-8 gate; an isolated FP64 experiment converged in roughly 53-56
    // iterations with restart 512. Give this inner solve enough Krylov space
    // to avoid the measured restart plateau. The iteration budget and
    // tolerance remain unchanged; the effective restart is clamped to n_rows
    // at the call site.
    cfdx::core::LinearSolverRequest pcd_convection_diffusion_linear_solver{
        cfdx::core::KrylovModel::FGMRES,
        cfdx::core::PreconditionerModel::ILU0,
        512,
        false,
        cfdx::core::NullSpaceModel::None};
    std::size_t pcd_laplacian_max_iterations = 100;
    std::size_t pcd_convection_diffusion_max_iterations = 100;
    double pcd_laplacian_tolerance = 1e-8;
    double pcd_convection_diffusion_tolerance = 1e-8;
    bool use_bounded_convection = true;
    ConvectionScheme convection_scheme = ConvectionScheme::UPWIND;
    DiagnosticsControls diagnostics;
    NonlinearRetryControls nonlinear_retry;
    std::vector<IncompressiblePointProbe> probes;
    // Optional quantities of interest gated by the convergence monitor. Each
    // gate must name a configured probe; a gate that names nothing is a
    // misconfiguration and is rejected rather than silently never satisfied.
    std::vector<cfdx::core::QoIGate> qoi_gates;
    AdaptiveRelaxationControls adaptive_relaxation;
    std::function<void(const IncompressibleProbeSample&)> probe_callback;
    // Called after each completed nonlinear iteration with the authoritative solver state.
    // Return false to stop the solve after the current iteration.
    using IterationOutputCallback = std::function<bool(
        std::size_t,
        double,
        const cfdx::core::Mesh&,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&)>;
    IterationOutputCallback iteration_output_callback;
    // Evidence hook called after the declared initialization strategy has been
    // applied and before geometry assembly or the first residual is formed.
    // References are valid only for the duration of the callback.
    using InitializedStateCallback = std::function<void(
        const cfdx::core::Mesh&,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&)>;
    InitializedStateCallback initialized_state_callback;
    // Evidence-only hook: invoked at the exact Krylov call point of the
    // coupled solve - after assembly, gauge treatment, finite checks and
    // preconditioner setup - so the snapshot matrix, RHS and initial guess
    // are the very objects handed to the solver. Null by default: no
    // capture, no overhead, no behavior change. It exists so an isolated
    // benchmark can replay the REAL solved system instead of a divergent
    // replica - see tests/validation/test_pcd_schur_action_benchmark.cpp.
    std::function<void(const CoupledAssemblySnapshot&)> coupled_assembly_callback;
};

struct IncompressibleIteration {
    std::size_t iteration = 0;
    double momentum_residual = std::numeric_limits<double>::infinity();
    double pressure_residual = std::numeric_limits<double>::infinity();
    double continuity_l1 = std::numeric_limits<double>::infinity();
    double continuity_linf = std::numeric_limits<double>::infinity();
    double continuity_normalized = std::numeric_limits<double>::infinity();
    // Explicit normalization evidence: raw residual / max(abs(momentum RHS))
    // and raw residual / first-iteration raw residual are different metrics.
    double momentum_rhs_scale = std::numeric_limits<double>::infinity();
    double momentum_equation_residual_reference = std::numeric_limits<double>::infinity();
    double momentum_equation_residual_rhs_relative = std::numeric_limits<double>::infinity();
    std::array<double, 3> momentum_rhs_scale_components{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    std::array<double, 3> momentum_equation_residual_rhs_relative_components{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    double momentum_equation_residual = std::numeric_limits<double>::infinity();
    double momentum_equation_residual_relative = std::numeric_limits<double>::infinity();
    double velocity_change_inf = std::numeric_limits<double>::infinity();
    double pressure_change_inf = std::numeric_limits<double>::infinity();
    double effective_alpha_u = 0.0;
    double effective_alpha_p = 0.0;
    double nonlinear_convergence_metric = std::numeric_limits<double>::infinity();
    std::size_t momentum_linear_iterations = 0;
    std::size_t pressure_linear_iterations = 0;
    std::size_t pressure_correctors_used = 0;
    double linear_tolerance_used = std::numeric_limits<double>::quiet_NaN();
    // Conservative pressure-corrected flux continuity at the end of the iteration.
    double corrected_flux_continuity_linf = std::numeric_limits<double>::infinity();
    // Continuity obtained by reconstructing the face flux from the corrected cell velocity.
    double reconstructed_velocity_continuity_linf = std::numeric_limits<double>::infinity();
    // Maximum face-flux discrepancy between the conservative flux and reconstructed U flux.
    double flux_velocity_mismatch_linf = std::numeric_limits<double>::infinity();
    // Boundary-separated evidence channel: same discrepancy restricted to
    // fixed-pressure boundary faces, where the two reconstruction paths
    // differ analytically (Rhie-Chow boundary correction plus the
    // predictor/corrector velocity response). Diagnostic only.
    double flux_velocity_mismatch_boundary_linf = std::numeric_limits<double>::infinity();
    // Normalized internal-face mismatch
    // E_phi = ||phi_RC - phi_U||_inf / max(||phi_RC||_inf, rho*max|U|*max|Sf|),
    // with an explicit physical floor so the value is comparable across
    // meshes and never divides by a vanishing flux scale.
    double flux_velocity_mismatch_normalized = std::numeric_limits<double>::infinity();
    // Operator identity: max |phi_stored - phi_rederived| where phi_rederived
    // is the same Rhie-Chow operator re-derived on the current state.
    // Expected zero up to round-off; nonzero flags a stale or mutated
    // authoritative flux, or an omitted/duplicated SIMPLEC response term.
    double rhie_chow_identity_linf = std::numeric_limits<double>::infinity();
    // SIMPLEC coefficient audit: compare the signed row sum used by the
    // production denominator with the absolute-value alternative. These are
    // diagnostics only; validation gates remain independent.
    std::array<double, 3> simplec_signed_offdiag_sum_linf{0.0, 0.0, 0.0};
    std::array<double, 3> simplec_abs_offdiag_sum_linf{0.0, 0.0, 0.0};
    std::array<double, 3> simplec_denominator_signed_linf{0.0, 0.0, 0.0};
    std::array<double, 3> simplec_denominator_abs_linf{0.0, 0.0, 0.0};
    std::array<double, 3> simplec_denominator_gap_linf{0.0, 0.0, 0.0};
    std::array<std::size_t, 3> simplec_denominator_worst_cell{0, 0, 0};

    // Component-wise and location-aware diagnostics for independently
    // reassembled momentum equations. These are diagnostic signals only;
    // acceptance thresholds remain owned by the validation tests.
    std::array<double, 3> momentum_equation_residual_components{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    double momentum_equation_residual_internal = std::numeric_limits<double>::infinity();
    double momentum_equation_residual_boundary = std::numeric_limits<double>::infinity();
    double pressure_gradient_linf = std::numeric_limits<double>::infinity();
    double pressure_gradient_l2 = std::numeric_limits<double>::infinity();
    std::size_t momentum_residual_cell = 0;
    double momentum_residual_no_pressure = std::numeric_limits<double>::infinity();
    double momentum_pressure_contribution = std::numeric_limits<double>::infinity();
    std::string momentum_residual_patch;
    // Independent finite-volume conservation diagnostics for the authoritative
    // pressure-corrected mass flux. These are distinct from linear residuals.
    double mass_boundary_flux = std::numeric_limits<double>::infinity();
    double mass_global_cell_balance = std::numeric_limits<double>::infinity();
    double mass_local_l1 = std::numeric_limits<double>::infinity();
    double mass_local_linf = std::numeric_limits<double>::infinity();
    double mass_local_l2 = std::numeric_limits<double>::infinity();
    double mass_normalized_imbalance = std::numeric_limits<double>::infinity();
    std::size_t mass_worst_cell = 0;
    std::size_t mass_nonfinite_faces = 0;
    std::size_t boundedness_nonfinite_velocity = 0;
    std::array<double, 3> momentum_conservation_residual{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    std::array<double, 3> momentum_conservation_normalized{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    std::array<std::size_t, 3> momentum_conservation_worst_cell{0,0,0};
};

struct IncompressibleSolveResult {
    bool converged = false;
    std::size_t iterations = 0;
    std::vector<IncompressibleIteration> history;
    // Retained point-probe samples, one per configured probe per iteration.
    // This is the artifact a probe export is written from.
    std::vector<IncompressibleProbeSample> probe_samples;
    // Final conservative face mass flux used by the pressure-velocity
    // continuity operator. This is exposed as evidence so validation can
    // compare algorithms on the same authoritative face operator rather than
    // reconstructing phi from the converged cell velocity.
    cfdx::core::Field<double, cfdx::core::Location::FACE> authoritative_mass_flux;
    double reference_momentum_residual = 0.0;
    cfdx::core::LinearSolverContextStats pressure_linear_context;
    cfdx::core::ConvergenceStatus convergence_status = cfdx::core::ConvergenceStatus::CONTINUE;
    std::string convergence_reason;
    // The plans the dispatcher actually resolved, not just the request. A
    // requested method can legitimately resolve to another one, and a report
    // cannot show that unless the resolved plan is recorded here.
    cfdx::core::LinearSolverPlan pressure_linear_plan;
    bool pressure_linear_plan_resolved = false;
    cfdx::core::LinearSolverPlan coupled_linear_plan;
    bool coupled_linear_plan_resolved = false;
    CoupledSchurModel requested_coupled_schur_model = CoupledSchurModel::BlockLocal;
    CoupledSchurModel resolved_coupled_schur_model = CoupledSchurModel::BlockLocal;
    bool coupled_schur_model_resolved = false;
};

// Applies the declared initialization strategy to the caller's fields.
//
// legacy_restart_path is the positional restart argument of
// solve_steady_incompressible. It supplies the checkpoint for a declared
// Restart mode, so existing callers that only pass the positional argument keep
// working. A caller that declares Provided or Uniform while also passing a
// legacy path gets the path ignored rather than silently applied, because that
// would make the initial state depend on an undeclared input.
inline void apply_initialization(
    const InitializationControls& init,
    const cfdx::core::Mesh& mesh,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    cfdx::io::DatRestartFields& restart_fields,
    const std::string& legacy_restart_path)
{
    using namespace cfdx::core;
    const std::string& restart_source =
        init.restart_path.empty() ? legacy_restart_path : init.restart_path;
    switch (init.mode) {
        case InitializationMode::Provided:
            // Caller-owned state is used as-is.
            break;
        case InitializationMode::Uniform:
            for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
                U.component_data(0)[c] = init.uniform_velocity.x;
                U.component_data(1)[c] = init.uniform_velocity.y;
                U.component_data(2)[c] = init.uniform_velocity.z;
                p(c) = init.uniform_pressure;
            }
            break;
        case InitializationMode::Restart: {
            if (restart_source.empty())
                throw std::invalid_argument(
                    "initialization mode 'restart' requires a checkpoint path");
            (void)cfdx::io::read_dat_restart_fields(
                restart_source, mesh, U, p, restart_fields);
            break;
        }
    }
}

inline void validate_initialization_controls(const InitializationControls& init) {
    switch (init.mode) {
        case InitializationMode::Provided:
            break;
        case InitializationMode::Uniform:
            if (!std::isfinite(init.uniform_pressure) ||
                !std::isfinite(init.uniform_velocity.x) ||
                !std::isfinite(init.uniform_velocity.y) ||
                !std::isfinite(init.uniform_velocity.z))
                throw std::invalid_argument(
                    "uniform initialization state must be finite");
            break;
        case InitializationMode::Restart:
            // The checkpoint may come from init.restart_path or from the
            // positional argument; apply_initialization resolves which.
            break;
    }
}

inline void validate_incompressible_controls(
    const IncompressibleSolverControls& c, std::size_t n_cells)
{
    validate_initialization_controls(c.initialization);
    validate_coupling_controls(c.coupling);
    validate_convergence_criteria(c.convergence);
    validate_convergence_acceleration_controls(c.acceleration);
    validate_adaptive_relaxation_controls(c.adaptive_relaxation);
    validate_nonlinear_retry_controls(c.nonlinear_retry);
    if (c.sst_runtime) {
        if (!c.sst_runtime->state)
            throw std::invalid_argument("SST production runtime requires an initialized state");
        if (c.sst_runtime->transport.model != TurbulenceModel::SST)
            throw std::invalid_argument("SST production runtime requires SST transport controls");
        if (c.sst_runtime->transport.density != c.density ||
            c.sst_runtime->transport.molecular_viscosity != c.kinematic_viscosity)
            throw std::invalid_argument(
                "SST transport density and molecular viscosity must match the carrier material properties");
        if (c.sst_runtime->max_transport_iterations == 0 ||
            !(c.sst_runtime->transport_tolerance > 0.0) ||
            !std::isfinite(c.sst_runtime->transport_tolerance))
            throw std::invalid_argument("SST production runtime has invalid transport controls");
        validate_sst_production_state(*c.sst_runtime->state, n_cells);
        if (!c.cell_turbulent_viscosity.empty() || c.turbulent_viscosity != 0.0)
            throw std::invalid_argument(
                "SST production runtime cannot be combined with prescribed turbulent viscosity");
    }
    if (!c.qoi_gates.empty()) {
        validate_convergence_monitor_controls(cfdx::core::ConvergenceMonitorControls{
            cfdx::core::ConvergenceGate{}, 1.0, 1, 1000, 1, 1.0, 1, 10.0, c.qoi_gates});
        for (const auto& gate : c.qoi_gates) {
            const bool declared =
                std::any_of(c.probes.begin(), c.probes.end(),
                            [&](const IncompressiblePointProbe& probe) {
                                return probe.name == gate.name;
                            });
            if (!declared)
                throw std::invalid_argument(
                    "convergence QoI gate '" + gate.name +
                    "' does not name a configured probe");
        }
    }
    if (c.transient.enabled) {
        if (!(c.transient.dt > 0.0) || !std::isfinite(c.transient.dt))
            throw std::invalid_argument("transient dual-time dt must be finite and positive");
        if (c.transient.scheme == cfdx::core::TimeScheme::BDF2 &&
            c.transient.history_valid &&
            (!(c.transient.dt_previous > 0.0) || !std::isfinite(c.transient.dt_previous)))
            throw std::invalid_argument("transient BDF2 previous dt must be finite and positive");
        if (!c.transient.previous || c.transient.previous->size() != n_cells ||
            c.transient.previous->dimension() != 3)
            throw std::invalid_argument("transient dual-time previous velocity history is missing or invalid");
        if (c.transient.scheme == cfdx::core::TimeScheme::BDF2 && c.transient.history_valid &&
            (!c.transient.previous_previous || c.transient.previous_previous->size() != n_cells ||
             c.transient.previous_previous->dimension() != 3))
            throw std::invalid_argument("transient BDF2 second velocity history is missing or invalid");
    }
    if (!(c.density > 0.0) || !std::isfinite(c.density) ||
        !(c.kinematic_viscosity >= 0.0) || !std::isfinite(c.kinematic_viscosity) ||
        !(c.turbulent_viscosity >= 0.0) || !std::isfinite(c.turbulent_viscosity))
        throw std::invalid_argument("invalid incompressible material properties");
    // Reject overflow in the effective dynamic viscosity before applying
    // initialization. Otherwise an invalid material combination can mutate
    // caller-owned U/p fields before assembly discovers non-finite coefficients.
    const auto validate_effective_viscosity = [&](double nu_t) {
        const double total_kinematic_viscosity = c.kinematic_viscosity + nu_t;
        const double effective_dynamic_viscosity =
            c.density * total_kinematic_viscosity;
        if (!std::isfinite(total_kinematic_viscosity) ||
            !std::isfinite(effective_dynamic_viscosity))
            throw std::invalid_argument(
                "effective dynamic viscosity must be finite");
    };
    validate_effective_viscosity(c.turbulent_viscosity);
    if (!c.cell_turbulent_viscosity.empty()) {
        if (c.cell_turbulent_viscosity.size() != n_cells)
            throw std::invalid_argument(
                "cellwise turbulent viscosity size does not match mesh");
        if (c.turbulent_viscosity != 0.0)
            throw std::invalid_argument(
                "specify either scalar or cellwise turbulent viscosity, not both");
        for (const double nu_t : c.cell_turbulent_viscosity) {
            if (!(nu_t >= 0.0) || !std::isfinite(nu_t))
                throw std::invalid_argument(
                    "cellwise turbulent viscosity must be finite and non-negative");
            validate_effective_viscosity(nu_t);
        }
    }
    if (!c.cell_body_force.empty()) {
        if (c.cell_body_force.size() != n_cells)
            throw std::invalid_argument(
                "cellwise body force size does not match mesh");
        if (c.body_force.x != 0.0 || c.body_force.y != 0.0 ||
            c.body_force.z != 0.0)
            throw std::invalid_argument(
                "specify either a uniform or a cellwise body force, not both");
        for (const auto& force : c.cell_body_force) {
            if (!std::isfinite(force.x) || !std::isfinite(force.y) ||
                !std::isfinite(force.z))
                throw std::invalid_argument(
                    "cellwise body force components must be finite");
        }
    }
    if (c.pressure_reference_cell >= n_cells)
        throw std::invalid_argument("pressure reference cell out of range");
    if (c.linear_max_iterations == 0 || !(c.linear_tolerance > 0.0))
        throw std::invalid_argument("invalid linear solver controls");
}

inline void apply_pressure_gauge_policy(
    cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    const FvGeometry& geometry,
    PressureGaugePolicy policy,
    std::size_t reference_cell,
    double reference_value,
    bool has_fixed_pressure_boundary)
{
    if (has_fixed_pressure_boundary || policy == PressureGaugePolicy::NONE) return;
    if (policy == PressureGaugePolicy::REFERENCE_CELL) {
        const double shift = reference_value - p(reference_cell);
        for (std::size_t c = 0; c < p.size(); ++c) p(c) += shift;
        return;
    }
    if (policy == PressureGaugePolicy::ZERO_MEAN) {
        double volume = 0.0, weighted_sum = 0.0;
        for (std::size_t c = 0; c < p.size(); ++c) {
            volume += geometry.cell_volumes[c];
            weighted_sum += geometry.cell_volumes[c] * p(c);
        }
        if (!(volume > 0.0) || !std::isfinite(volume) || !std::isfinite(weighted_sum))
            throw std::runtime_error("pressure gauge: invalid volume-weighted mean");
        const double mean = weighted_sum / volume;
        for (std::size_t c = 0; c < p.size(); ++c) p(c) -= mean;
    }
}

inline cfdx::core::Field<double, cfdx::core::Location::CELL>
gauss_gradient_with_boundary(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& field,
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const ScalarBoundaryConditions& bcs)
{
    using namespace cfdx::core;
    if (field.dimension() != 1 || field.size() != mesh.n_cells())
        throw std::invalid_argument("gauss_gradient_with_boundary: invalid field");

    Field<double, Location::CELL> grad(
        mesh.n_cells(), field.name() + "_grad_bc", field.metadata().unit + "/m", 3);
    const auto& own = mesh.ownership();

    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        Vec3 sum{0.0, 0.0, 0.0};
        const Offset off = mesh.cells().offsets_data()[c];
        const Offset count = mesh.cells().offsets_data()[c + 1] - off;
        for (Offset k = 0; k < count; ++k) {
            const std::size_t f = mesh.cells().faces_data()[off + k];
            const bool owner = own.owner(f) == c;
            Vec3 Sf = owner ? geometry.face_area_vectors[f] : geometry.face_area_vectors[f] * -1.0;
            double vf = field(c);

            if (own.neighbour(f) >= 0) {
                // The mesh stores owner/neighbour independent of the cell being
                // visited.  When c is the neighbour, own.neighbour(f) == c;
                // using it again would therefore collapse the face interpolation
                // to field(c) and halve the Green-Gauss pressure gradient.
                const std::size_t other = owner
                    ? static_cast<std::size_t>(own.neighbour(f))
                    : static_cast<std::size_t>(own.owner(f));
                const double dc = (geometry.face_centres[f] -
                                   geometry.cell_centres[c]).mag();
                const double dn = (geometry.cell_centres[other] -
                                   geometry.face_centres[f]).mag();
                if (!(dc + dn > 0.0))
                    throw std::runtime_error(
                        "gauss_gradient_with_boundary: degenerate face distance");
                const double w = dn / (dc + dn);
                vf = w * field(c) + (1.0 - w) * field(other);
            } else {
                const std::size_t p = geometry.face_patch[f];
                if (p < mesh.boundary().n_patches() &&
                    mesh.boundary().patch(p).type == PatchType::EMPTY)
                    continue;
                if (p < mesh.boundary().n_patches()) {
                    const auto& name = mesh.boundary().patch(p).name;
                    const auto it = bcs.find(name);
                    if (it != bcs.end() && it->second.type == ScalarBoundaryType::FIXED_VALUE)
                        vf = it->second.value;
                }
            }
            sum += Sf * vf;
        }
        const double invV = 1.0 / geometry.cell_volumes[c];
        grad.component_data(0)[c] = sum.x * invV;
        grad.component_data(1)[c] = sum.y * invV;
        grad.component_data(2)[c] = sum.z * invV;
    }
    return grad;
}

inline cfdx::core::Field<double, cfdx::core::Location::FACE>
make_mass_flux(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    double rho,
    const VelocityBoundaryConditions& bcs,
    const ScalarBoundaryConditions& pressure_bcs = {})
{
    using namespace cfdx::core;
    (void)pressure_bcs;
    Field<double, Location::FACE> flux(mesh.n_faces(), "phi", "kg/s", 1);
    const auto& own = mesh.ownership();

    for (std::size_t c = 0; c < U.size(); ++c) {
        const double ux = U.component_data(0)[c];
        const double uy = U.component_data(1)[c];
        const double uz = U.component_data(2)[c];
        if (!std::isfinite(ux) || !std::isfinite(uy) || !std::isfinite(uz))
            throw std::runtime_error("make_mass_flux: non-finite cell velocity");
    }
    for (std::size_t f=0; f<mesh.n_faces(); ++f) {
        const std::size_t o=own.owner(f);
        Vec3 Uf;
        U.get(o,Uf.x,Uf.y,Uf.z);
        if (own.neighbour(f)>=0) {
            const std::size_t n=static_cast<std::size_t>(own.neighbour(f));
            Vec3 Un;
            U.get(n,Un.x,Un.y,Un.z);
            Uf=(Uf+Un)*0.5;
        } else {
            const std::size_t p=geometry.face_patch[f];
            if (p < mesh.boundary().n_patches() &&
                mesh.boundary().patch(p).type == PatchType::EMPTY) {
                flux(f) = 0.0;
                continue;
            }
            if(p<mesh.boundary().n_patches()) {
                const auto& name=mesh.boundary().patch(p).name;
                const auto it=bcs.find(name);
                if(it!=bcs.end() && it->second.type==VelocityBoundaryCondition::Type::FIXED_VALUE)
                    Uf=it->second.value;
            }
        }
        flux(f)=rho*Uf.dot(geometry.face_area_vectors[f]);
    }
    return flux;
}

// Orthogonal implicit coefficient shared by the Rhie-Chow face flux and the
// pressure-correction matrix. dAU_diagonal contains the volume-weighted
// directional inverse momentum diagonal (V*rAU), interpolated to the face.
// The implicit pressure difference is along the owner-to-neighbour line, so
// the matching surface measure is Sf.e, not |Sf|. They coincide only on an
// orthogonal face. The deferred non-orthogonal part is handled separately by
// the Rhie-Chow flux reconstruction.
inline double pressure_face_orthogonal_coefficient(
    const cfdx::core::Vec3& Sf,
    const cfdx::core::Vec3& dvec,
    const cfdx::core::Vec3& dAU_diagonal)
{
    const double d = dvec.mag();
    const double area = Sf.mag();
    if (!(d > 0.0) || !std::isfinite(d) || !(area > 0.0) || !std::isfinite(area))
        throw std::runtime_error("pressure_face_orthogonal_coefficient: invalid face geometry");
    const cfdx::core::Vec3 e{dvec.x / d, dvec.y / d, dvec.z / d};
    const double orthogonal_area = Sf.dot(e);
    if (!(orthogonal_area > 0.0) || !std::isfinite(orthogonal_area))
        throw std::runtime_error(
            "pressure_face_orthogonal_coefficient: face orientation is inconsistent with cell-centre direction");
    const double rfn = dAU_diagonal.x * e.x * e.x +
                       dAU_diagonal.y * e.y * e.y +
                       dAU_diagonal.z * e.z * e.z;
    const double coefficient = rfn * orthogonal_area / d;
    if (!std::isfinite(coefficient))
        throw std::runtime_error("pressure_face_orthogonal_coefficient: non-finite coefficient");
    return coefficient;
}

// Deferred part of the pressure-response face flux. The implicit term uses
// only the centre-line pressure jump; this remainder makes the sum equal to
// (diag(dAU) * grad(p)) dot Sf for a linear pressure field, including skewed
// faces and anisotropic momentum response.
inline double pressure_face_nonorthogonal_flux(
    const cfdx::core::Vec3& Sf,
    const cfdx::core::Vec3& dvec,
    const cfdx::core::Vec3& dAU_diagonal,
    const cfdx::core::Vec3& grad_p)
{
    const double d = dvec.mag();
    if (!(d > 0.0) || !std::isfinite(d))
        throw std::runtime_error("pressure_face_nonorthogonal_flux: invalid centre distance");
    const cfdx::core::Vec3 e{dvec.x / d, dvec.y / d, dvec.z / d};
    const double orthogonal_area = Sf.dot(e);
    const double rfn = dAU_diagonal.x * e.x * e.x +
                       dAU_diagonal.y * e.y * e.y +
                       dAU_diagonal.z * e.z * e.z;
    const double full_response_flux =
        dAU_diagonal.x * grad_p.x * Sf.x +
        dAU_diagonal.y * grad_p.y * Sf.y +
        dAU_diagonal.z * grad_p.z * Sf.z;
    const double implicit_gradient_flux =
        rfn * grad_p.dot(e) * orthogonal_area;
    const double deferred = full_response_flux - implicit_gradient_flux;
    if (!std::isfinite(deferred))
        throw std::runtime_error("pressure_face_nonorthogonal_flux: non-finite deferred flux");
    return deferred;
}

inline double rhie_chow_pressure_flux_internal(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    std::size_t face,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& grad_p,
    const std::array<std::vector<double>, 3>& rAU,
    double rho)
{
    using namespace cfdx::core;
    const auto nr = mesh.ownership().neighbour(face);
    if (nr < 0) return 0.0;
    const std::size_t o = mesh.ownership().owner(face);
    const std::size_t n = static_cast<std::size_t>(nr);
    const Vec3 Sf = geometry.face_area_vectors[face];
    const double area = Sf.mag();
    if (!(area > 0.0) || !std::isfinite(area))
        throw std::runtime_error("rhie_chow_pressure_flux_internal: degenerate face");
    const Vec3 dvec = geometry.cell_centres[n] - geometry.cell_centres[o];
    const double d = dvec.mag();
    if (!(d > 0.0) || !std::isfinite(d))
        throw std::runtime_error("rhie_chow_pressure_flux_internal: degenerate centre distance");
    // rAU is the inverse integrated momentum diagonal (1/A_P). The volume
    // factor is introduced here, at the velocity/pressure projection stage,
    // because the pressure gradient is a volumetric momentum source. Thus
    // cx/cy/cz are dAU = V*rAU, not a redefinition of rAU itself.
    const double cx = 0.5*(geometry.cell_volumes[o] * rAU[0][o] +
                           geometry.cell_volumes[n] * rAU[0][n]);
    const double cy = 0.5*(geometry.cell_volumes[o] * rAU[1][o] +
                           geometry.cell_volumes[n] * rAU[1][n]);
    const double cz = 0.5*(geometry.cell_volumes[o] * rAU[2][o] +
                           geometry.cell_volumes[n] * rAU[2][n]);
    const Vec3 gp{
        0.5*(grad_p.component_data(0)[o] + grad_p.component_data(0)[n]),
        0.5*(grad_p.component_data(1)[o] + grad_p.component_data(1)[n]),
        0.5*(grad_p.component_data(2)[o] + grad_p.component_data(2)[n])};
    // The face response uses the same dAU = V/A_P operator as the cell
    // velocity reconstruction. Keeping the orthogonal pressure difference
    // implicit and the non-orthogonal remainder deferred is important: the
    // continuity matrix and the conservative Rhie-Chow flux must represent
    // the same discrete operator, otherwise a residual floor can survive
    // even when the momentum equations appear converged.
    // Orthogonal pressure difference is implicit in the pressure matrix.
    // The non-orthogonal remainder is deferred explicitly using the same
    // Gauss gradient. This keeps the face flux and pressure equation on the
    // same discrete operator on arbitrary meshes.
    const double implicit_coefficient = pressure_face_orthogonal_coefficient(
        Sf, dvec, {cx, cy, cz});
    const double non_orthogonal_flux = pressure_face_nonorthogonal_flux(
        Sf, dvec, {cx, cy, cz}, gp);
    return rho * (
        implicit_coefficient * (p(n) - p(o)) +
        non_orthogonal_flux);
}

inline cfdx::core::Field<double, cfdx::core::Location::FACE>
make_rhie_chow_mass_flux(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& U,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& p,
    const std::array<std::vector<double>, 3>& rAU,
    double rho,
    const VelocityBoundaryConditions& bcs,
    const ScalarBoundaryConditions& pressure_bcs = {})
;

inline cfdx::core::Field<double, cfdx::core::Location::FACE>
make_rhie_chow_mass_flux(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& U,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& p,
    const std::vector<double>& rAU,
    double rho,
    const VelocityBoundaryConditions& bcs)
{
    std::array<std::vector<double>, 3> directional_rAU{rAU, rAU, rAU};
    return make_rhie_chow_mass_flux(
        mesh, geometry, U, p, directional_rAU, rho, bcs, {});
}

inline cfdx::core::Field<double, cfdx::core::Location::FACE>
make_rhie_chow_mass_flux(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& U,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& p,
    const std::array<std::vector<double>, 3>& rAU,
    double rho,
    const VelocityBoundaryConditions& bcs,
    const ScalarBoundaryConditions& pressure_bcs)
{
    using namespace cfdx::core;
    const auto& own = mesh.ownership();
    if (p.size() != mesh.n_cells())
        throw std::invalid_argument("make_rhie_chow_mass_flux: field size mismatch");
    // Boundary pressure is passed explicitly to the gradient operator. A
    // boundary face has no neighbour coefficient: any dAU/rAU interpolation
    // must therefore use the owner-cell value rather than a ghost/zero entry.
    const auto grad_p = gauss_gradient_with_boundary(
        p, mesh, geometry, pressure_bcs);
    for (const auto& component : rAU) {
        if (component.size() != mesh.n_cells())
            throw std::invalid_argument("make_rhie_chow_mass_flux: inverse momentum diagonal size mismatch");
        for (const double value : component) {
            if (!std::isfinite(value) || value <= 0.0)
                throw std::invalid_argument(
                    "make_rhie_chow_mass_flux: inverse momentum diagonal must be finite and positive");
        }
    }

    // Momentum-weighted (Rhie-Chow) face flux:
    //
    //     phi_f = rho [ H/A |_f . Sf - rAU_f * grad(p)_f . Sf ]
    //
    // The first term is linearly interpolated from the cell-centred predictor.
    // The implicit pressure jump is taken along the owner-neighbour centre
    // line and multiplied by Sf.e, not |Sf|. The tangential remainder is a
    // deferred non-orthogonal correction. The same implicit coefficient is
    // assembled in the pressure-correction matrix below; it must never be
    // replaced by an unrelated least-squares fit of the cell velocity.
    //
    // For a fixed-pressure boundary the boundary pressure is known. For a
    // zero-gradient pressure boundary there is no normal pressure correction.
    Field<double, Location::FACE> flux(mesh.n_faces(), "phi", "kg/s", 1);
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const std::size_t o = own.owner(f);
        Vec3 Hf;
        U.get(o, Hf.x, Hf.y, Hf.z);
        const auto nr = own.neighbour(f);
        double correction = 0.0;
        if (nr >= 0) {
            const std::size_t n = static_cast<std::size_t>(nr);
            Vec3 Hn;
            U.get(n, Hn.x, Hn.y, Hn.z);
            Hf = (Hf + Hn) * 0.5;
            correction = rhie_chow_pressure_flux_internal(
                mesh, geometry, f, p, grad_p, rAU, rho);
        } else {
            const std::size_t patch = geometry.face_patch[f];
            if (patch < mesh.boundary().n_patches() &&
                mesh.boundary().patch(patch).type == PatchType::EMPTY) {
                flux(f) = 0.0;
                continue;
            }
            if (patch < mesh.boundary().n_patches()) {
                const auto& name = mesh.boundary().patch(patch).name;
                const auto ubc = bcs.find(name);
                if (ubc != bcs.end() &&
                    ubc->second.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                    Hf = ubc->second.value;

                const auto pbc = pressure_bcs.find(name);
                if (pbc != pressure_bcs.end() &&
                    pbc->second.type == ScalarBoundaryType::FIXED_VALUE) {
                    const std::size_t o = own.owner(f);
                    const Vec3 Sf = geometry.face_area_vectors[f];
                    const Vec3 dvec =
                        geometry.face_centres[f] - geometry.cell_centres[o];
                    const double d = dvec.mag();
                    if (!(d > 0.0) || !std::isfinite(d))
                        throw std::runtime_error("make_rhie_chow_mass_flux: degenerate boundary distance");
                    const Vec3 gp{
                        grad_p.component_data(0)[o],
                        grad_p.component_data(1)[o],
                        grad_p.component_data(2)[o]};
                    const Vec3 dAU_boundary{
                        geometry.cell_volumes[o] * rAU[0][o],
                        geometry.cell_volumes[o] * rAU[1][o],
                        geometry.cell_volumes[o] * rAU[2][o]};
                    const double implicit_coefficient =
                        pressure_face_orthogonal_coefficient(
                            Sf, dvec, dAU_boundary);
                    const double deferred = pressure_face_nonorthogonal_flux(
                        Sf, dvec, dAU_boundary, gp);
                    correction = rho * (
                        implicit_coefficient * (pbc->second.value - p(o)) +
                        deferred);
                }
            }
        }
        flux(f) = rho * Hf.dot(geometry.face_area_vectors[f]) - correction;
    }
    return flux;
}

inline ScalarEquation assemble_momentum_component(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& pressure_gradient_component,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& body_component,
    double effective_dynamic_viscosity,
    const ScalarBoundaryConditions& velocity_bcs,
    std::size_t component,
    bool bounded,
    ConvectionScheme convection_scheme = ConvectionScheme::UPWIND,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* convected_field = nullptr)
{
    cfdx::core::Field<double, cfdx::core::Location::CELL> source(
        mesh.n_cells(), "momentum_source", "N/m3", 1);
    cfdx::core::Field<double, cfdx::core::Location::CELL> sp(
        mesh.n_cells(), "momentum_sp", "kg/m3/s", 1);
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        if (pressure_gradient_component.dimension() != 3 ||
            pressure_gradient_component.size() != mesh.n_cells())
            throw std::invalid_argument("assemble_momentum_component: pressure gradient must be a 3-component cell field");
        source(c) = body_component(c) - pressure_gradient_component.component_data(component)[c];
        sp(c) = 0.0;
    }
    if (component >= 3)
        throw std::invalid_argument("assemble_momentum_component: component index out of range");
    return assemble_scalar_equation(
        mesh, geometry, mass_flux, effective_dynamic_viscosity,
        source, sp, velocity_bcs, bounded, nullptr, nullptr, nullptr, nullptr,
        convection_scheme, convected_field);
}

inline ScalarEquation assemble_momentum_component(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& pressure_gradient_component,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& body_component,
    const std::vector<double>& cell_effective_dynamic_viscosity,
    const ScalarBoundaryConditions& velocity_bcs,
    std::size_t component,
    bool bounded,
    ConvectionScheme convection_scheme = ConvectionScheme::UPWIND,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* convected_field = nullptr)
{
    const std::size_t nc = mesh.n_cells();
    if (cell_effective_dynamic_viscosity.size() != nc)
        throw std::invalid_argument(
            "assemble_momentum_component: cell viscosity size does not match mesh");
    for (const double viscosity : cell_effective_dynamic_viscosity) {
        if (!std::isfinite(viscosity) || viscosity < 0.0)
            throw std::invalid_argument(
                "assemble_momentum_component: cell viscosity must be finite and non-negative");
    }
    if (component >= 3)
        throw std::invalid_argument("assemble_momentum_component: component index out of range");
    if (pressure_gradient_component.dimension() != 3 ||
        pressure_gradient_component.size() != nc ||
        body_component.size() != nc || body_component.dimension() != 1)
        throw std::invalid_argument(
            "assemble_momentum_component: source/pressure-gradient field shape mismatch");

    cfdx::core::Field<double, cfdx::core::Location::CELL> source(
        nc, "momentum_source", "N/m3", 1);
    cfdx::core::Field<double, cfdx::core::Location::CELL> sp(
        nc, "momentum_sp", "kg/m3/s", 1);
    for (std::size_t c = 0; c < nc; ++c) {
        source(c) = body_component(c) -
                    pressure_gradient_component.component_data(component)[c];
        sp(c) = 0.0;
    }
    return assemble_scalar_equation(
        mesh, geometry, mass_flux, 0.0, source, sp, velocity_bcs, bounded,
        nullptr, nullptr, nullptr, &cell_effective_dynamic_viscosity,
        convection_scheme, convected_field);
}

inline ScalarEquation assemble_momentum_component(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& pressure_gradient_component,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& body_component,
    double effective_dynamic_viscosity,
    const BoundaryConstraintMap& boundary_constraints,
    std::size_t component,
    bool bounded,
    double time = 0.0,
    ConvectionScheme convection_scheme = ConvectionScheme::UPWIND,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* convected_field = nullptr)
{
    if (component >= 3)
        throw std::invalid_argument("assemble_momentum_component: component index out of range");
    if (pressure_gradient_component.dimension() != 3 ||
        pressure_gradient_component.size() != mesh.n_cells())
        throw std::invalid_argument("assemble_momentum_component: pressure gradient must be a 3-component cell field");
    cfdx::core::Field<double, cfdx::core::Location::CELL> source(
        mesh.n_cells(), "momentum_source", "N/m3", 1);
    cfdx::core::Field<double, cfdx::core::Location::CELL> sp(
        mesh.n_cells(), "momentum_sp", "kg/m3/s", 1);
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        source(c) = body_component(c) - pressure_gradient_component.component_data(component)[c];
        sp(c) = 0.0;
    }
    const std::string field = "U." + std::string(1, static_cast<char>('x' + component));
    const auto face_conditions = resolve_scalar_boundary_constraints(
        mesh, geometry.face_centres, boundary_constraints, field, time, &mass_flux);
    return assemble_scalar_equation(
        mesh, geometry, mass_flux, effective_dynamic_viscosity,
        source, sp, {}, bounded, nullptr, nullptr, nullptr, nullptr,
        convection_scheme, convected_field, &face_conditions);
}

inline void add_dual_time_mass_term(
    ScalarEquation& equation,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& u_n,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& u_nm1,
    const FvGeometry& geometry,
    double rho,
    double dt,
    double dt_previous,
    bool bdf2,
    std::size_t component)
{
    if (!(dt > 0.0) || !std::isfinite(dt) || !(rho > 0.0) || !std::isfinite(rho))
        throw std::invalid_argument("add_dual_time_mass_term: invalid rho/dt");
    if (u_n.dimension() != 3 || u_nm1.dimension() != 3 ||
        u_n.size() != equation.rhs.size() || u_nm1.size() != equation.rhs.size())
        throw std::invalid_argument("add_dual_time_mass_term: history shape mismatch");
    if (component >= 3) throw std::invalid_argument("add_dual_time_mass_term: invalid component");

    double a0 = 1.0, a1 = 0.0, denom = 1.0, b = dt;
    if (bdf2) {
        if (!(dt_previous > 0.0) || !std::isfinite(dt_previous))
            throw std::invalid_argument("add_dual_time_mass_term: invalid previous dt for BDF2");
        cfdx::core::bdf2_coefficients(dt_previous, dt, a0, a1, denom, b);
    }
    for (std::size_t c = 0; c < equation.rhs.size(); ++c) {
        const double mass = rho * geometry.cell_volumes[c] * denom / b;
        const double source = rho * geometry.cell_volumes[c] *
            (a0 * u_n.component_data(component)[c] -
             a1 * u_nm1.component_data(component)[c]) / b;
        if (!(mass > 0.0) || !std::isfinite(mass) || !std::isfinite(source))
            throw std::runtime_error("add_dual_time_mass_term: non-finite temporal coefficient");
        auto* rows = equation.matrix.row_offsets_data();
        auto* cols = equation.matrix.columns_data();
        auto* values = equation.matrix.values_data();
        bool found = false;
        for (std::uint32_t k = rows[c]; k < rows[c + 1]; ++k) {
            if (cols[k] == c) {
                values[k] += mass;
                found = true;
                break;
            }
        }
        if (!found) throw std::runtime_error("add_dual_time_mass_term: missing momentum diagonal");
        equation.diagonal[c] += mass;
        equation.rhs(c) += source;
    }
}

inline void apply_velocity_boundary_conditions(
    const cfdx::core::Mesh& mesh,
    const VelocityBoundaryConditions& bcs,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& U)
{
    // Cell-centred finite-volume fields have no independent patch storage.
    // Fixed-value conditions are therefore enforced through the matrix and
    // the face flux; no cell overwrite is performed here.
    (void)mesh;
    (void)bcs;
    (void)U;
}

inline double sample_incompressible_probe(
    const IncompressiblePointProbe& probe,
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& p)
{
    if (mesh.n_cells() == 0)
        throw std::invalid_argument("sample_incompressible_probe: mesh has no cells");
    std::size_t nearest = 0;
    double best_distance = std::numeric_limits<double>::infinity();
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        const auto delta = geometry.cell_centres[c] - probe.location;
        const double distance = delta.mag2();
        if (distance < best_distance) {
            best_distance = distance;
            nearest = c;
        }
    }
    switch (probe.field) {
    case IncompressibleProbeField::PRESSURE: return p(nearest);
    case IncompressibleProbeField::U_X: { cfdx::core::Vec3 v; U.get(nearest, v.x, v.y, v.z); return v.x; }
    case IncompressibleProbeField::U_Y: { cfdx::core::Vec3 v; U.get(nearest, v.x, v.y, v.z); return v.y; }
    case IncompressibleProbeField::U_Z: { cfdx::core::Vec3 v; U.get(nearest, v.x, v.y, v.z); return v.z; }
    case IncompressibleProbeField::U_MAGNITUDE: {
        cfdx::core::Vec3 velocity;
        U.get(nearest, velocity.x, velocity.y, velocity.z);
        return velocity.mag();
    }
    }
    throw std::invalid_argument("sample_incompressible_probe: unsupported field");
}

inline void relax_momentum_equation(
    ScalarEquation& equation,
    const cfdx::core::Vector& old_field,
    double alpha)
{
    if (!(alpha > 0.0 && alpha <= 1.0) ||
        old_field.size() != equation.rhs.size())
        throw std::invalid_argument("relax_momentum_equation: invalid relaxation");

    if (alpha == 1.0) return;

    const auto* rows = equation.matrix.row_offsets_data();
    const auto* cols = equation.matrix.columns_data();
    auto* values = equation.matrix.values_data();

    // Standard equation relaxation:
    //
    //   A u = b
    //
    // becomes
    //
    //   [A_off + A_P/alpha] u
    //       = b + (1-alpha)/alpha A_P u_old.
    //
    // This preserves the exact algebraic relation used later to construct
    // HbyA and rAU. A field blend after the solve would break that relation.
    for (std::size_t c = 0; c < old_field.size(); ++c) {
        const double a = equation.diagonal[c];
        if (!(a > 0.0) || !std::isfinite(a))
            throw std::runtime_error("relax_momentum_equation: invalid diagonal");

        double* diagonal_value = nullptr;
        for (std::uint32_t k = rows[c]; k < rows[c + 1]; ++k) {
            if (cols[k] == c) {
                diagonal_value = &values[k];
                break;
            }
        }
        if (!diagonal_value)
            throw std::runtime_error("relax_momentum_equation: missing diagonal entry");

        const double added = a * (1.0 / alpha - 1.0);
        *diagonal_value += added;
        equation.diagonal[c] += added;
        equation.rhs(c) += added * old_field(c);
    }
}

struct PcdPressureOperators {
    cfdx::core::SparseMatrix mass;
    cfdx::core::SparseMatrix laplacian;
    cfdx::core::SparseMatrix convection_diffusion;
    std::unique_ptr<cfdx::core::Preconditioner> laplacian_preconditioner;
    std::unique_ptr<cfdx::core::Preconditioner> convection_diffusion_preconditioner;
};

inline cfdx::core::SparseMatrix impose_pcd_reference_row(
    const cfdx::core::SparseMatrix& input, std::size_t reference_cell)
{
    const std::size_t n = input.n_rows();
    if (n == 0 || reference_cell >= n)
        throw std::invalid_argument("PCD reference cell is outside pressure operator");
    cfdx::core::SparseMatrix result(n, n);
    const auto* ro = input.row_offsets_data();
    const auto* co = input.columns_data();
    const auto* va = input.values_data();
    for (std::size_t r = 0; r < n; ++r) {
        if (r == reference_cell) {
            result.push_back(r, r, 1.0);
            continue;
        }
        for (std::uint32_t k = ro[r]; k < ro[r + 1]; ++k)
            result.push_back(r, co[k], va[k]);
    }
    result.finalize();
    return result;
}

// PCD pressure-side operators. The convection_scheme parameter is accepted
// but deliberately NOT forwarded to assemble_scalar_equation, and this is
// not a TODO. Every ConvectionScheme except UPWIND is implemented in
// assemble_scalar_equation as a DEFERRED CORRECTION: an explicit
// right-hand-side term evaluated on the convected field, riding the SAME
// conservative first-order upwind matrix. These operators are matrix-only
// preconditioner ingredients, so the assembled matrix is identical for
// every scheme and ConvectionScheme::UPWIND below is not an approximation
// of the true pressure operator but its exact matrix for every scheme
// setting. Forwarding a non-UPWIND scheme would also throw, because a
// preconditioner assembly has no convected field to evaluate the deferred
// correction against. Audit context: issue #663 (NACA0012 convection-scheme
// sensitivity campaign).
inline PcdPressureOperators assemble_pcd_pressure_operators(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& velocity,
    double rho,
    double kinematic_viscosity,
    const VelocityBoundaryConditions& velocity_bcs,
    const ScalarBoundaryConditions& pressure_bcs,
    std::size_t reference_cell,
    ConvectionScheme /*convection_scheme — unused by design, see comment above*/)
{
    using namespace cfdx::core;
    const std::size_t nc = mesh.n_cells();
    if (velocity.size() != nc || velocity.dimension() != 3)
        throw std::invalid_argument("PCD pressure operator assembly: invalid velocity field");
    if (!(rho > 0.0) || !std::isfinite(rho) ||
        !(kinematic_viscosity >= 0.0) || !std::isfinite(kinematic_viscosity))
        throw std::invalid_argument("PCD pressure operator assembly: invalid fluid properties");

    PcdPressureOperators result;
    result.mass = SparseMatrix(nc, nc);
    for (std::size_t c = 0; c < nc; ++c)
        result.mass.push_back(c, c, geometry.cell_volumes[c]);
    result.mass.finalize();

    Field<double, Location::FACE> zero_flux(mesh.n_faces(), "pcd_zero_flux", "", 1);
    auto mass_flux = make_mass_flux(mesh, geometry, velocity, rho, velocity_bcs);

    Field<double, Location::CELL> zero_source(nc, "pcd_zero_source", "", 1);
    Field<double, Location::CELL> zero_implicit(nc, "pcd_zero_implicit", "", 1);

    const auto diffusion = assemble_scalar_equation(
        mesh, geometry, zero_flux, 1.0,
        zero_source, zero_implicit, pressure_bcs, true, nullptr, nullptr,
        nullptr, nullptr, ConvectionScheme::UPWIND, nullptr);

    const auto convection_diffusion = assemble_scalar_equation(
        mesh, geometry, mass_flux, rho * kinematic_viscosity,
        zero_source, zero_implicit, pressure_bcs, true, nullptr, nullptr,
        nullptr, nullptr, ConvectionScheme::UPWIND, nullptr);

    // A fixed pressure boundary already removes the constant null mode of
    // the pressure operators, so the gauge reference row is imposed only
    // when no such boundary exists. The fixed-pressure decision is made
    // here, from pressure_bcs, so every caller passes its true reference
    // cell and the caller-side gauge policy cannot drift from the
    // assembly. A gauge sentinel (>= nc) still leaves the operators
    // untouched for direct callers that pass one.
    bool fixed_pressure_gauge = false;
    for (const auto& [bc_name, bc] : pressure_bcs) {
        (void)bc_name;
        fixed_pressure_gauge = fixed_pressure_gauge ||
            bc.type == ScalarBoundaryType::FIXED_VALUE;
    }
    const std::size_t gauge_cell = fixed_pressure_gauge
        ? std::numeric_limits<std::size_t>::max() : reference_cell;
    if (gauge_cell < nc) {
        result.laplacian = impose_pcd_reference_row(diffusion.matrix, gauge_cell);
        result.convection_diffusion =
            impose_pcd_reference_row(convection_diffusion.matrix, gauge_cell);
    } else {
        result.laplacian = diffusion.matrix;
        result.convection_diffusion = convection_diffusion.matrix;
    }
    return result;
}

// Coherent gauge for the nested algebraic-Schur pressure operators
// (the LSC/BFBt stabilized operator). The coupled production system
// imposes the pressure reference by eliminating the reference row AND
// column and inserting a unit diagonal (see the reduced coupled
// assembly), so the nested operator mirrors exactly that structure.
//
// The reference right-hand side entry is deliberately NOT zeroed and
// the nested correction at the reference is deliberately NOT pinned to
// zero: the outer Krylov iteration needs the gauge component of the
// preconditioned residual to correct the coupled gauge row (see the
// PR #832 note on the initial gauge row, and the measured
// frozen-residual failure documented in pcd_schur.h). A fixed pressure
// boundary already removes the constant null mode, so the caller
// passes a sentinel (>= n) and this function is a no-op in that case.
inline void apply_nested_schur_pressure_gauge(
    cfdx::core::SparseMatrix& pressure_operator, std::size_t gauge_cell)
{
    const std::size_t n = pressure_operator.n_rows();
    if (n == 0 || gauge_cell >= n)
        return;
    for (std::size_t k = pressure_operator.row_offsets_data()[gauge_cell];
         k < pressure_operator.row_offsets_data()[gauge_cell + 1]; ++k) {
        pressure_operator.values_data()[k] = 0.0;
    }
    for (std::size_t row = 0; row < n; ++row) {
        if (row == gauge_cell)
            continue;
        for (std::size_t k = pressure_operator.row_offsets_data()[row];
             k < pressure_operator.row_offsets_data()[row + 1]; ++k) {
            if (pressure_operator.columns_data()[k] == gauge_cell)
                pressure_operator.values_data()[k] = 0.0;
        }
    }
    bool ref_diag_set = false;
    for (std::size_t k = pressure_operator.row_offsets_data()[gauge_cell];
         k < pressure_operator.row_offsets_data()[gauge_cell + 1]; ++k) {
        if (pressure_operator.columns_data()[k] == gauge_cell) {
            pressure_operator.values_data()[k] = 1.0;
            ref_diag_set = true;
            break;
        }
    }
    if (!ref_diag_set) {
        cfdx::core::SparseMatrix rebuilt(n, n);
        for (std::size_t row = 0; row < n; ++row) {
            if (row == gauge_cell) {
                rebuilt.push_back(row, row, 1.0);
                continue;
            }
            for (std::size_t k = pressure_operator.row_offsets_data()[row];
                 k < pressure_operator.row_offsets_data()[row + 1]; ++k) {
                rebuilt.push_back(row, pressure_operator.columns_data()[k],
                                  pressure_operator.values_data()[k]);
            }
        }
        rebuilt.finalize();
        pressure_operator = std::move(rebuilt);
    }
}

// Best-effort contract for the nested pressure solves invoked inside a
// coupled Schur approximation action (PCD Fp/Kp, LSC/BFBt H). These
// callbacks are preconditioner applications, not standalone solves: the
// Schur action itself is only an approximation of the true Schur inverse
// (pcd_schur.h, lsc_bfbt_schur.h), and the outer coupled FGMRES owns the
// final accuracy through its independently recomputed true-residual gate.
// Requiring strict CONVERGED turns a single benign inner shortfall into a
// total outer NOT_APPLICABLE abort: on the N9 Couette campaign one Fp apply
// that could not reach the last digit of its relative tolerance killed the
// whole coupled solve after 168 successful applies, because by then the
// outer residual had shrunk enough for the nested rhs direction to hit its
// own rounding floor. A nested run that exhausted its iteration budget but
// produced a finite best iterate (the Krylov drivers restore the iterate
// with the best verified true residual) remains a valid preconditioner
// action. Only a genuinely unusable result fails the action: divergence,
// breakdown, an inapplicable inner setup, or a non-finite iterate.
inline bool nested_schur_pressure_solve_usable(
    const cfdx::core::SolverResult& result,
    const cfdx::core::Vector& x,
    std::size_t expected_size)
{
    if (result.status == cfdx::core::SolverStatus::CONVERGED)
        return true;
    if (result.status != cfdx::core::SolverStatus::MAX_ITER_REACHED)
        return false;
    if (x.size() != expected_size || !x.is_valid())
        return false;
    for (std::size_t i = 0; i < x.size(); ++i)
        if (!std::isfinite(x(i)))
            return false;
    return true;
}

inline cfdx::core::SolverResult solve_coupled_momentum_continuity(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    ScalarEquation ex,
    ScalarEquation ey,
    ScalarEquation ez,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& U_old,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& p_old,
    const VelocityBoundaryConditions& velocity_bcs,
    const ScalarBoundaryConditions& pressure_bcs,
    double rho,
    double kinematic_viscosity,
    std::size_t reference_cell,
    double reference_value,
    std::size_t max_iterations,
    double tolerance,
    const DiagnosticsControls& diagnostics,
    CoupledSchurModel schur_model,
    const cfdx::core::LinearSolverRequest& solver_request,
    const IncompressibleSolverControls& controls,
    // Live momentum under-relaxation consumed only by the coupled SIMPLEC
    // Schur approximation. The coupled matrix itself stays unrelaxed;
    // Picard relaxation is applied to the field update after the solve.
    double momentum_under_relaxation,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    cfdx::core::LinearSolverPlan* resolved_linear_plan = nullptr)
{
    using namespace cfdx::core;
    const std::size_t nc = mesh.n_cells();
    const std::size_t nv = 3 * nc;
    const std::size_t np = nc;
    const std::size_t n = nv + np;

    SparseMatrix A(n, n);
    Vector b(n, 0.0);

    // The momentum equations received the current pressure explicitly during
    // segregated assembly. Undo that explicit contribution here: the coupled
    // matrix contains the pressure gradient implicitly instead.
    const auto grad_p_old =
        gauss_gradient_with_boundary(p_old, mesh, geometry, pressure_bcs);
    auto add_momentum_block = [&](const ScalarEquation& eq, std::size_t component) {
        const std::size_t row_base = component * nc;
        const auto* ro = eq.matrix.row_offsets_data();
        const auto* co = eq.matrix.columns_data();
        const auto* va = eq.matrix.values_data();
        const double* old_component = U_old.component_data(component);
        for (std::size_t r = 0; r < nc; ++r) {
            const double diagonal = eq.diagonal[r];
            if (!(diagonal > 0.0) || !std::isfinite(diagonal))
                throw std::runtime_error(
                    "solve_coupled_momentum_continuity: invalid momentum diagonal at cell " +
                    std::to_string(r) + " component " + std::to_string(component));

            bool has_diagonal = false;
            for (std::uint32_t k = ro[r]; k < ro[r + 1]; ++k) {
                const std::size_t column = component * nc + co[k];
                A.push_back(row_base + r, column, va[k]);
                has_diagonal = has_diagonal || (column == row_base + r);
            }
            // The scalar equation keeps its authoritative diagonal separately
            // from the CSR entries. If a transport row reaches the coupled
            // assembler without an explicit diagonal entry, restore that
            // coefficient from the authoritative metadata instead of creating
            // an empty coupled momentum row.
            if (!has_diagonal)
                A.push_back(row_base + r, row_base + r, diagonal);
            b(row_base + r) =
                eq.rhs(r) + grad_p_old.component_data(component)[r] * geometry.cell_volumes[r];
        }
    };
    add_momentum_block(ex, 0);
    add_momentum_block(ey, 1);
    add_momentum_block(ez, 2);

    // Move the current explicit pressure term from b back into the linear
    // system. The FV Gauss gradient is linear in the pressure field, so its
    // coefficients can be assembled face-by-face.
    for (std::size_t c = 0; c < nc; ++c) {
        const Offset off = mesh.cells().offsets_data()[c];
        const Offset count = mesh.cells().offsets_data()[c + 1] - off;
        for (Offset k = 0; k < count; ++k) {
            const std::size_t f = mesh.cells().faces_data()[off + k];
            const bool owner = mesh.ownership().owner(f) == c;
            const double sign = owner ? 1.0 : -1.0;
            const Vec3 Sf = geometry.face_area_vectors[f] * sign;
            const auto nr = mesh.ownership().neighbour(f);

            if (nr >= 0) {
                // Same owner/neighbour rule as the continuity block below:
                // if c is the neighbour, nr points back to c itself.
                const int other_cell = owner
                    ? nr
                    : static_cast<int>(mesh.ownership().owner(f));
                if (other_cell < 0 ||
                    static_cast<std::size_t>(other_cell) >= nc ||
                    static_cast<std::size_t>(other_cell) == c)
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid face topology on face " +
                        std::to_string(f) + " cell " + std::to_string(c));
                const std::size_t ncell = static_cast<std::size_t>(other_cell);
                // Match gauss_gradient_with_boundary() exactly: the
                // cell-centred pressure gradient uses a distance-weighted
                // face value, not an unconditional arithmetic average.
                // On a skewed mesh dc and dn need not be equal. Using 0.5 here
                // therefore assembled a different pressure-gradient operator
                // from the one used by the segregated momentum equation and
                // from the independent final momentum audit.
                const double dc =
                    (geometry.face_centres[f] - geometry.cell_centres[c]).mag();
                const double dn =
                    (geometry.cell_centres[ncell] - geometry.face_centres[f]).mag();
                if (!(dc > 0.0) || !(dn > 0.0) || !std::isfinite(dc) ||
                    !std::isfinite(dn))
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid face-centre distances on face " +
                        std::to_string(f));
                const double weight_c = dn / (dc + dn);
                const double weight_n = dc / (dc + dn);
                A.push_back(c, 3*nc + c, Sf.x * weight_c);
                A.push_back(c, 3*nc + ncell, Sf.x * weight_n);
                A.push_back(nc + c, 3*nc + c, Sf.y * weight_c);
                A.push_back(nc + c, 3*nc + ncell, Sf.y * weight_n);
                A.push_back(2*nc + c, 3*nc + c, Sf.z * weight_c);
                A.push_back(2*nc + c, 3*nc + ncell, Sf.z * weight_n);
            } else {
                const std::size_t patch = geometry.face_patch[f];
                const ScalarBoundaryCondition* bc = nullptr;
                ScalarBoundaryCondition local_bc;
                if (patch < mesh.boundary().n_patches()) {
                    const auto& name = mesh.boundary().patch(patch).name;
                    const auto it = pressure_bcs.find(name);
                    if (it != pressure_bcs.end()) {
                        local_bc = it->second;
                        bc = &local_bc;
                    }
                }
                if (bc && bc->type == ScalarBoundaryType::FIXED_VALUE) {
                    // Known boundary pressure is moved to the momentum RHS.
                    b(c) -= Sf.x * bc->value;
                    b(nc + c) -= Sf.y * bc->value;
                    b(2*nc + c) -= Sf.z * bc->value;
                } else {
                    // The transport-layer Gauss gradient treats zero-gradient
                    // pressure as owner-cell extrapolation. Retain that exact
                    // operator here so the coupled and segregated residuals
                    // remain comparable.
                    A.push_back(c, 3*nc + c, Sf.x);
                    A.push_back(nc + c, 3*nc + c, Sf.y);
                    A.push_back(2*nc + c, 3*nc + c, Sf.z);
                }
            }
        }
    }

    // Continuity is assembled from the same face velocity interpolation used
    // by the pressure-based solver. The pressure Schur diagonal below is the
    // exact diagonal of the same momentum-weighted (Rhie-Chow) face-flux
    // derivative and is passed explicitly to the coupled preconditioner.
    std::vector<double> coupled_pressure_schur_diagonal(nc, 0.0);
    for (std::size_t c = 0; c < nc; ++c) {
        const std::size_t row = nv + c;
        const Offset off = mesh.cells().offsets_data()[c];
        const Offset count = mesh.cells().offsets_data()[c + 1] - off;
        for (Offset k = 0; k < count; ++k) {
            const std::size_t f = mesh.cells().faces_data()[off + k];
            const bool owner = mesh.ownership().owner(f) == c;
            const double sign = owner ? 1.0 : -1.0;
            const Vec3 Sf = geometry.face_area_vectors[f] * sign;
            const auto nr = mesh.ownership().neighbour(f);

            if (nr >= 0) {
                // neighbour(f) is the opposite cell only on the owner side.
                // When c is the neighbour, neighbour(f) == c; using it directly
                // creates a self-coupled face (zero centre distance) and corrupts
                // the coupled pressure Schur coefficient.
                const int other_cell = owner
                    ? nr
                    : static_cast<int>(mesh.ownership().owner(f));
                if (other_cell < 0 ||
                    static_cast<std::size_t>(other_cell) >= nc ||
                    static_cast<std::size_t>(other_cell) == c)
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid face topology on face " +
                        std::to_string(f) + " cell " + std::to_string(c));
                const std::size_t ncell = static_cast<std::size_t>(other_cell);
                const double half_rho = 0.5 * rho;
                A.push_back(row, c, half_rho * Sf.x);
                A.push_back(row, ncell, half_rho * Sf.x);
                A.push_back(row, nc + c, half_rho * Sf.y);
                A.push_back(row, nc + ncell, half_rho * Sf.y);
                A.push_back(row, 2*nc + c, half_rho * Sf.z);
                A.push_back(row, 2*nc + ncell, half_rho * Sf.z);

                const Vec3 rawSf = geometry.face_area_vectors[f];
                const double area = rawSf.mag();
                // Keep the geometric direction aligned with the stored
                // owner-to-neighbour face orientation even when this row is
                // the neighbour cell; the matrix coefficient must be positive
                // and identical in both continuity rows.
                const std::size_t face_owner = mesh.ownership().owner(f);
                const std::size_t face_neighbour = static_cast<std::size_t>(nr);
                const Vec3 dvec = geometry.cell_centres[face_neighbour] -
                                  geometry.cell_centres[face_owner];
                const double d = dvec.mag();
                if (!(area > 0.0) || !std::isfinite(area) ||
                    !(d > 0.0) || !std::isfinite(d))
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid face geometry on face " +
                        std::to_string(f));
                const Vec3 e{dvec.x / d, dvec.y / d, dvec.z / d};
                const double aox = ex.diagonal[c], anx = ex.diagonal[ncell];
                const double aoy = ey.diagonal[c], any = ey.diagonal[ncell];
                const double aoz = ez.diagonal[c], anz = ez.diagonal[ncell];
                if (!(aox > 0.0) || !(anx > 0.0) ||
                    !(aoy > 0.0) || !(any > 0.0) ||
                    !(aoz > 0.0) || !(anz > 0.0) ||
                    !std::isfinite(aox) || !std::isfinite(anx) ||
                    !std::isfinite(aoy) || !std::isfinite(any) ||
                    !std::isfinite(aoz) || !std::isfinite(anz))
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid momentum diagonal on face " +
                        std::to_string(f) + " cells " + std::to_string(c) + "/" +
                        std::to_string(ncell));
                // The coupled Schur complement uses the same dAU = V/A_P
                // pressure response as the segregated Rhie-Chow operator.
                const double rAUx_f = 0.5 * (
                    geometry.cell_volumes[c] / aox +
                    geometry.cell_volumes[ncell] / anx);
                const double rAUy_f = 0.5 * (
                    geometry.cell_volumes[c] / aoy +
                    geometry.cell_volumes[ncell] / any);
                const double rAUz_f = 0.5 * (
                    geometry.cell_volumes[c] / aoz +
                    geometry.cell_volumes[ncell] / anz);
                const double rfn = rAUx_f*e.x*e.x + rAUy_f*e.y*e.y + rAUz_f*e.z*e.z;
                if (!(rfn > 0.0) || !std::isfinite(rfn))
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid face rAU on face " +
                        std::to_string(f));
                const double D = rho * pressure_face_orthogonal_coefficient(
                    rawSf, dvec, {rAUx_f, rAUy_f, rAUz_f});
                if (!(d > 0.0) || !(area > 0.0) ||
                    !std::isfinite(d) || !std::isfinite(area) ||
                    !(D > 0.0) || !std::isfinite(D)) {
                    if (diagnostics.coupled_matrix_summary) {
                        std::cerr << "COUPLED_FACE_DEBUG face=" << f
                                  << " cells=" << c << "/" << ncell
                                  << " area=" << area
                                  << " distance=" << d
                                  << " a=(" << aox << "," << anx
                                  << ";" << aoy << "," << any
                                  << ";" << aoz << "," << anz << ")"
                                  << " rAU=(" << rAUx_f << "," << rAUy_f
                                  << "," << rAUz_f << ")"
                                  << " rfn=" << rfn
                                  << " D=" << D << "\n";
                    }
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid pressure coefficient on face " +
                        std::to_string(f) + " cells " + std::to_string(c) + "/" +
                        std::to_string(ncell));
                }
                // phi_p = D (p_P - p_N). Only this orthogonal derivative is
                // implicit in the coupled matrix. The non-orthogonal remainder
                // is a deferred correction evaluated from p_old below, exactly
                // as in the segregated Rhie-Chow operator.
                A.push_back(row, nv + c, D);
                A.push_back(row, nv + ncell, -D);
                coupled_pressure_schur_diagonal[c] += D;

                const Vec3 gp{
                    0.5*(grad_p_old.component_data(0)[c] +
                         grad_p_old.component_data(0)[ncell]),
                    0.5*(grad_p_old.component_data(1)[c] +
                         grad_p_old.component_data(1)[ncell]),
                    0.5*(grad_p_old.component_data(2)[c] +
                         grad_p_old.component_data(2)[ncell])};
                // Reuse the same deferred tensor-response operator as SIMPLE
                // and the pressure-flux reconstruction.
                const double nonorth_flux = rho * pressure_face_nonorthogonal_flux(
                    rawSf, dvec, {rAUx_f, rAUy_f, rAUz_f}, gp);
                // The continuity equation is div(phi)=0 and the pressure
                // correction is written as phi = phi_U - D grad(p).
                // Therefore the explicit non-orthogonal pressure flux is
                // subtracted from the current continuity residual.
                // nonorth_flux is evaluated with the raw owner->neighbour
                // face orientation, so the neighbour-side visit must flip
                // its sign: like the implicit pressure part, which enters
                // every row as D*(p_self - p_other), the deferred flux
                // must enter this cell's outward divergence with the
                // face sign.
                b(row) += sign * nonorth_flux;
            } else {
                const std::size_t patch = geometry.face_patch[f];
                if (patch < mesh.boundary().n_patches() &&
                    mesh.boundary().patch(patch).type == PatchType::EMPTY)
                    continue;
                const VelocityBoundaryCondition* vbc = nullptr;
                VelocityBoundaryCondition local_vbc;
                const ScalarBoundaryCondition* pbc = nullptr;
                ScalarBoundaryCondition local_pbc;
                if (patch < mesh.boundary().n_patches()) {
                    const auto& name = mesh.boundary().patch(patch).name;
                    const auto vit = velocity_bcs.find(name);
                    if (vit != velocity_bcs.end()) {
                        local_vbc = vit->second;
                        vbc = &local_vbc;
                    }
                    const auto pit = pressure_bcs.find(name);
                    if (pit != pressure_bcs.end()) {
                        local_pbc = pit->second;
                        pbc = &local_pbc;
                    }
                }

                if (vbc && vbc->type == VelocityBoundaryCondition::Type::FIXED_VALUE) {
                    b(row) -= rho * vbc->value.dot(Sf);
                } else {
                    A.push_back(row, c, rho * Sf.x);
                    A.push_back(row, nc + c, rho * Sf.y);
                    A.push_back(row, 2*nc + c, rho * Sf.z);
                }

                if (pbc && pbc->type == ScalarBoundaryType::FIXED_VALUE) {
                    const Vec3 rawSf = geometry.face_area_vectors[f];
                    const Vec3 dvec = geometry.face_centres[f] -
                                      geometry.cell_centres[c];
                    const double a_x = ex.diagonal[c];
                    const double a_y = ey.diagonal[c];
                    const double a_z = ez.diagonal[c];
                    if (!(a_x > 0.0) || !(a_y > 0.0) || !(a_z > 0.0) ||
                        !std::isfinite(a_x) || !std::isfinite(a_y) || !std::isfinite(a_z))
                        throw std::runtime_error(
                            "solve_coupled_momentum_continuity: invalid boundary momentum diagonal on face " +
                            std::to_string(f) + " cell " + std::to_string(c));
                    const double rAUx = geometry.cell_volumes[c] / a_x;
                    const double rAUy = geometry.cell_volumes[c] / a_y;
                    const double rAUz = geometry.cell_volumes[c] / a_z;
                    const double D = rho * pressure_face_orthogonal_coefficient(
                        rawSf, dvec, {rAUx, rAUy, rAUz});
                    if (!(D > 0.0) || !std::isfinite(D))
                        throw std::runtime_error(
                            "solve_coupled_momentum_continuity: invalid boundary pressure coefficient on face " +
                            std::to_string(f));
                    A.push_back(row, nv + c, D);
                    coupled_pressure_schur_diagonal[c] += D;
                    b(row) += D * pbc->value;
                    const Vec3 gp_boundary{
                        grad_p_old.component_data(0)[c],
                        grad_p_old.component_data(1)[c],
                        grad_p_old.component_data(2)[c]};
                    const double deferred_boundary = pressure_face_nonorthogonal_flux(
                        rawSf, dvec, {rAUx, rAUy, rAUz}, gp_boundary);
                    // Match the boundary's deferred pressure response to the
                    // internal-face treatment; otherwise COUPLED drops the
                    // tangential pressure flux at fixed-pressure boundaries.
                    b(row) += rho * deferred_boundary;
                }
            }
        }
    }

    // A pure-Neumann pressure field has one gauge null mode. Eliminate the
    // reference pressure column and replace the reference continuity row by
    // p_ref = reference_value. This keeps the coupled system deterministic.
    bool has_fixed_pressure = false;
    for (const auto& [name, bc] : pressure_bcs) {
        (void)name;
        has_fixed_pressure = has_fixed_pressure ||
            bc.type == ScalarBoundaryType::FIXED_VALUE;
    }

    if (!has_fixed_pressure) {
        // Finalize the assembled COO matrix before reading its CSR arrays.
        // Without this explicit finalize(), row_offsets_/values_/columns_ still
        // describe an empty CSR shell, so gauge elimination would silently
        // discard every assembled momentum and continuity coefficient.
        A.finalize();

        // Rebuild the matrix while eliminating the known gauge column.
        SparseMatrix reduced(n, n);
        const auto* ro = A.row_offsets_data();
        const auto* co = A.columns_data();
        const auto* va = A.values_data();
        for (std::size_t r = 0; r < n; ++r) {
            if (r == nv + reference_cell) {
                reduced.push_back(r, nv + reference_cell, 1.0);
                b(r) = reference_value;
                continue;
            }
            for (std::uint32_t k = ro[r]; k < ro[r + 1]; ++k) {
                const std::size_t col = co[k];
                if (col == nv + reference_cell) {
                    b(r) -= va[k] * reference_value;
                } else {
                    reduced.push_back(r, col, va[k]);
                }
            }
        }
        reduced.finalize();
        A = std::move(reduced);
    } else {
        A.finalize();
    }

    Vector x(n, 0.0);
    for (std::size_t c = 0; c < nc; ++c) {
        x(c) = U_old.component_data(0)[c];
        x(nc + c) = U_old.component_data(1)[c];
        x(2*nc + c) = U_old.component_data(2)[c];
        x(nv + c) = p_old(c);
    }

    // PR #832 fix: zero the gauge row of the initial guess before GMRES.
    // The gauge row is an identity row with rhs = reference_value (default 0.0).
    // If x0(gauge) = p_old(reference_cell) is non-zero, the right-preconditioned
    // GMRES may struggle to correct it because the preconditioner preserves
    // the gauge component. Setting x0(gauge) = 0.0 gives the solver the correct
    // starting point and avoids the stagnation observed in PR #830.
    x(nv + reference_cell) = controls.pressure_reference_value;

    // Reject non-finite coefficients before entering GMRES. A bad local
    // coefficient must be diagnosed at assembly time rather than surfacing as
    // an opaque Krylov failure.
    {
        const auto* values = A.values_data();
        for (std::size_t k = 0; k < A.nnz(); ++k) {
            if (!std::isfinite(values[k]))
                throw std::runtime_error(
                    "solve_coupled_momentum_continuity: non-finite matrix coefficient at nnz " +
                    std::to_string(k));
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (!std::isfinite(b(i)))
                throw std::runtime_error(
                    "solve_coupled_momentum_continuity: non-finite RHS at row " +
                    std::to_string(i));
        }
    }

    // Evidence-only: preserve the raw (pre-pin) Rhie-Chow diagonal for the
    // coupled_assembly_callback when one is installed. The pinned form
    // stored in coupled_pressure_schur_diagonal is diagnostic-only and is
    // never consumed by production preconditioning.
    std::vector<double> pressure_schur_diagonal_raw_copy;
    if (controls.coupled_assembly_callback)
        pressure_schur_diagonal_raw_copy = coupled_pressure_schur_diagonal;

    // Optional algebraic microscope. It is deliberately computed from
    // the exact post-gauge matrix passed to GMRES, so it exposes the actual
    // M/G/D/C blocks rather than reconstructed proxy operators.
    if (diagnostics.coupled_matrix_summary) {
        double m2 = 0.0, g2 = 0.0, d2 = 0.0, c2 = 0.0;
        std::size_t zero_rows = 0;
        std::vector<std::size_t> row_nnz(n, 0);
        std::vector<std::size_t> col_nnz(n, 0);
        const auto* ro = A.row_offsets_data();
        const auto* co = A.columns_data();
        const auto* va = A.values_data();

        for (std::size_t r = 0; r < n; ++r) {
            row_nnz[r] = ro[r + 1] - ro[r];
            if (row_nnz[r] == 0) ++zero_rows;
            for (std::size_t k = ro[r]; k < ro[r + 1]; ++k) {
                const std::size_t col = co[k];
                const double value = va[k];
                col_nnz[col]++;
                if (r < nv && col < nv)
                    m2 += value * value;
                else if (r < nv && col >= nv)
                    g2 += value * value;
                else if (r >= nv && col < nv)
                    d2 += value * value;
                else
                    c2 += value * value;
            }
        }
        std::size_t zero_cols = 0;
        for (const auto count : col_nnz)
            if (count == 0) ++zero_cols;

        double schur_diag_min = std::numeric_limits<double>::infinity();
        double schur_diag_max = 0.0;
        double schur_mismatch_max = 0.0;
        std::size_t schur_near_zero = 0;
        for (std::size_t c = 0; c < nc; ++c) {
            const std::size_t pr = nv + c;
            double effective = 0.0;
            for (std::size_t k = ro[pr]; k < ro[pr + 1]; ++k) {
                const std::size_t u = co[k];
                if (u >= nv) continue;
                const double d = va[k];
                double inv_m = 0.0;
                bool found = false;
                for (std::size_t mk = ro[u]; mk < ro[u + 1]; ++mk) {
                    if (co[mk] == u) {
                        const double diag = va[mk];
                        if (std::isfinite(diag) && std::abs(diag) > 1e-30) {
                            inv_m = 1.0 / diag;
                            found = true;
                        }
                        break;
                    }
                }
                if (!found) continue;
                for (std::size_t gk = ro[u]; gk < ro[u + 1]; ++gk) {
                    const std::size_t q = co[gk];
                    if (q == pr)
                        effective -= d * inv_m * va[gk];
                }
            }
            for (std::size_t k = ro[pr]; k < ro[pr + 1]; ++k)
                if (co[k] == pr) effective += va[k];

            schur_diag_min = std::min(schur_diag_min, effective);
            schur_diag_max = std::max(schur_diag_max, std::abs(effective));
            if (std::abs(effective) <=
                1e-12 * std::max(1.0, schur_diag_max))
                ++schur_near_zero;
            if (c < coupled_pressure_schur_diagonal.size() &&
                std::isfinite(coupled_pressure_schur_diagonal[c]))
                schur_mismatch_max = std::max(
                    schur_mismatch_max,
                    std::abs(effective - coupled_pressure_schur_diagonal[c]));
        }

        std::cerr << "COUPLED_MATRIX_SUMMARY n=" << n
                  << " nnz=" << A.nnz()
                  << " M_l2=" << std::sqrt(m2)
                  << " G_l2=" << std::sqrt(g2)
                  << " D_l2=" << std::sqrt(d2)
                  << " C_l2=" << std::sqrt(c2)
                  << " zero_rows=" << zero_rows
                  << " zero_cols=" << zero_cols
                  << " schur_diag_min=" << schur_diag_min
                  << " schur_diag_absmax=" << schur_diag_max
                  << " schur_near_zero=" << schur_near_zero
                  << " rhie_chow_schur_vs_algebraic_diag_delta=" << schur_mismatch_max
                  << " gauge_row=" << (nv + reference_cell)
                  << "\n";
    }

    // The coupled matrix is a saddle-point system. A cell-local 4x4
    // inverse is generally singular because continuity rows have no
    // pressure-pressure diagonal. Select the N8 linear plan before constructing
    // its preconditioner so the production path uses exactly the resolved policy.
    coupled_pressure_schur_diagonal[reference_cell] = 1.0;
    // The production Couette system is only 512 unknowns. A restart of
    // 128 can stagnate on the nonsymmetric saddle-point spectrum even when
    // the assembled system is nonsingular. Use one full Krylov space for
    // small/medium coupled systems so a restart boundary cannot manufacture
    // an artificial MAX_ITER failure.
    // Feed the dispatcher the characteristics measured from the assembled
    // production CSR matrix. The convenience overload only knows the problem
    // kind and equation count; it cannot see coefficient scaling, symmetry,
    // diagonal dominance or the actual sparsity pattern.
    const auto matrix_characteristics =
        cfdx::core::measure_matrix_characteristics(A, true);
    const auto solver_plan = select_linear_solver(
        LinearProblemKind::CoupledPressureVelocity,
        matrix_characteristics,
        solver_request);
    if (solver_plan.krylov != KrylovModel::GMRES &&
        solver_plan.krylov != KrylovModel::FGMRES)
        throw std::invalid_argument("coupled solver currently requires GMRES or FGMRES");
    if (solver_plan.preconditioner != PreconditionerModel::CoupledBlockSchur &&
        solver_plan.preconditioner != PreconditionerModel::PCD &&
        solver_plan.preconditioner != PreconditionerModel::MGR)
        throw std::invalid_argument(
            "coupled solver requires coupled_block_schur, pcd, or mgr");

    // N8 owns the linear-solver/preconditioner policy. The coupled FV
    // assembler remains N9-owned, but an explicit BlockSchur request must use
    // the N8 AMG-backed Schur implementation rather than the legacy scalar
    // preconditioner. MGR remains a separate N8-qualified path until its
    // coupled dispatch metadata is integrated.
    //
    // Selection reads the *resolved* plan, not the request, so an automatically
    // resolved policy takes the same branch an explicit one would.
    const bool use_n8_block_schur =
        solver_plan.preconditioner == PreconditionerModel::CoupledBlockSchur;
    // The Schur family is the production method selector. The outer
    // preconditioner family selects the surrounding block infrastructure;
    // PCD therefore remains reachable even when the outer plan is the generic
    // coupled-block Schur container.
    const bool use_pcd =
        schur_model == CoupledSchurModel::PCD;
    const bool use_lsc_bfbt =
        schur_model == CoupledSchurModel::LSC ||
        schur_model == CoupledSchurModel::BFBT;
    const bool use_simple_schur =
        schur_model == CoupledSchurModel::SIMPLE ||
        schur_model == CoupledSchurModel::SIMPLEC;
    std::unique_ptr<Preconditioner> coupled_preconditioner;
    // Evidence-only alias of the production PCD pressure operators, kept
    // alive for the coupled_assembly_callback (the branch-local pcd_ops
    // shared_ptr goes out of scope below). Assigned only when a callback is
    // installed and the plan resolves to PCD; null otherwise.
    std::shared_ptr<PcdPressureOperators> evidence_pcd_ops;
    const auto extract_coupled_block = [&](std::size_t row_block,
                                           std::size_t col_block,
                                           std::size_t row_size,
                                           std::size_t col_size) {
        SparseMatrix block(row_size, col_size);
        const std::size_t velocity_size = 3 * nc;
        const std::size_t row_offset = row_block == 0 ? 0 : velocity_size;
        const std::size_t col_offset = col_block == 0 ? 0 : velocity_size;
        for (std::size_t r = 0; r < row_size; ++r) {
            const std::size_t gr = row_offset + r;
            for (std::size_t k = A.row_offsets_data()[gr];
                 k < A.row_offsets_data()[gr + 1]; ++k) {
                const std::size_t gc = A.columns_data()[k];
                if (gc >= col_offset && gc < col_offset + col_size)
                    block.push_back(r, gc - col_offset, A.values_data()[k]);
            }
        }
        block.finalize();
        return block;
    };
    if (use_simple_schur) {
        if (!use_n8_block_schur)
            throw std::invalid_argument(
                "SIMPLE/SIMPLEC coupled Schur models require the coupled_block_schur linear plan");

        const auto mode =
            schur_model == CoupledSchurModel::SIMPLE
                ? SimplerSchurMode::SIMPLE
                : SimplerSchurMode::SIMPLEC;
        // The coupled matrix is deliberately unrelaxed, but the SIMPLEC
        // consistent denominator is defined on the under-relaxed momentum
        // row A_P/alpha - sum(a_PN), exactly as the segregated SIMPLEC
        // assembles it from relax_momentum_equation() output. On the
        // unrelaxed steady row the plain row sum A_P - sum(a_PN) is
        // numerically zero for interior cells (only boundary/source terms
        // survive), which made SimplerSchurApproximation::setup() reject
        // the matrix and abort the N9 campaign at the first COUPLED
        // simplec case. Feed the live momentum under-relaxation to the
        // SIMPLEC Schur approximation; SIMPLE keeps its plain diag(Auu)
        // denominator unchanged.
        auto simpler = std::make_unique<SimplerSchurApproximation>(
            mode, std::nullopt,
            schur_model == CoupledSchurModel::SIMPLEC
                ? momentum_under_relaxation
                : 1.0);
        CoupledBlockSchurOptions options;
        options.factorization = CoupledSchurFactorization::Full;
        options.velocity_approximation = CoupledSchurVelocityApproximation::Block;
        options.schur_approximation =
            schur_model == CoupledSchurModel::SIMPLE
                ? CoupledSchurApproximationModel::SIMPLE
                : CoupledSchurApproximationModel::SIMPLEC;

        auto schur = std::make_unique<CoupledBlockSchurAMGPreconditioner>(
            nc, options, has_fixed_pressure
                ? std::numeric_limits<std::size_t>::max() : reference_cell);
        schur->set_simpler_schur(std::move(simpler));
        if (!schur->setup(A)) {
            throw std::runtime_error(
                std::string("N9 SIMPLE/SIMPLEC Schur setup failed: ") +
                schur->last_error());
        }
        coupled_preconditioner = std::move(schur);
    } else if (use_lsc_bfbt) {
        if (!use_n8_block_schur)
            throw std::invalid_argument(
                "LSC/BFBt coupled Schur models require the coupled_block_schur linear plan");

        const SparseMatrix Auu = extract_coupled_block(0, 0, 3 * nc, 3 * nc);
        const SparseMatrix G = extract_coupled_block(0, 1, 3 * nc, nc);
        const SparseMatrix D = extract_coupled_block(1, 0, nc, 3 * nc);
        const SparseMatrix C = extract_coupled_block(1, 1, nc, nc);

        std::vector<double> q_diagonal(3 * nc, 1.0);
        if (schur_model == CoupledSchurModel::LSC) {
            for (std::size_t i = 0; i < 3 * nc; ++i) {
                double diagonal = 0.0;
                for (std::size_t k = Auu.row_offsets_data()[i];
                     k < Auu.row_offsets_data()[i + 1]; ++k) {
                    if (Auu.columns_data()[k] == i)
                        diagonal += Auu.values_data()[k];
                }
                if (!(diagonal > 0.0) || !std::isfinite(diagonal))
                    throw std::runtime_error(
                        "LSC Schur setup requires positive finite momentum diagonal");
                q_diagonal[i] = diagonal;
            }
        }

        const BlockOperator lsc_blocks(Auu, G, D, C);
        SparseMatrix pressure_operator =
            LscBfbtSchurApproximation::assemble_stabilized_pressure_operator(
                lsc_blocks,
                schur_model == CoupledSchurModel::LSC
                    ? LscBfbtSchurApproximation::Mode::LSC
                    : LscBfbtSchurApproximation::Mode::BFBT,
                q_diagonal);

        // Gauge the nested LSC/BFBt pressure operator ONLY when no fixed
        // pressure boundary already removes the constant null mode (same
        // sentinel policy as the reduced coupled system and the outer
        // preconditioner), and mirror that reduced system exactly:
        // reference row AND reference column eliminated with a unit
        // diagonal. The nested solve keeps the raw right-hand side, so
        // the gauge component passes through to the outer Krylov
        // iteration (pinning it freezes the coupled gauge row residual;
        // see apply_nested_schur_pressure_gauge).
        const std::size_t lsc_gauge_cell = has_fixed_pressure
            ? std::numeric_limits<std::size_t>::max() : reference_cell;
        apply_nested_schur_pressure_gauge(pressure_operator, lsc_gauge_cell);

        auto pressure_preconditioner =
            std::shared_ptr<cfdx::core::Preconditioner>(
                cfdx::core::make_scalar_preconditioner(
                    cfdx::core::PreconditionerModel::NativeAMG).release());
        if (!pressure_preconditioner ||
            !pressure_preconditioner->setup(pressure_operator)) {
            throw std::runtime_error(
                "N8 LSC/BFBt pressure AMG setup failed");
        }

        auto solve_pressure = [pressure_operator,
                               pressure_preconditioner,
                               max_iterations, tolerance](
            const Vector& rhs, Vector& x) mutable {
            const int restart = static_cast<int>(std::min<std::size_t>(
                64, pressure_operator.n_rows()));
            const auto result = solve_gmres(
                pressure_operator, rhs, x, restart,
                max_iterations, tolerance,
                pressure_preconditioner.get());
            return nested_schur_pressure_solve_usable(
                result, x, pressure_operator.n_rows());
        };

        auto algebraic_schur = std::make_unique<LscBfbtSchurApproximation>(
            schur_model == CoupledSchurModel::LSC
                ? LscBfbtSchurApproximation::Mode::LSC
                : LscBfbtSchurApproximation::Mode::BFBT,
            std::move(solve_pressure),
            std::move(q_diagonal));

        CoupledBlockSchurOptions options;
        options.factorization = CoupledSchurFactorization::Full;
        options.velocity_approximation = CoupledSchurVelocityApproximation::Block;
        options.schur_approximation =
            schur_model == CoupledSchurModel::LSC
                ? CoupledSchurApproximationModel::LSC
                : CoupledSchurApproximationModel::BFBT;

        auto schur = std::make_unique<CoupledBlockSchurAMGPreconditioner>(
            nc, options, has_fixed_pressure
                ? std::numeric_limits<std::size_t>::max() : reference_cell);
        schur->set_algebraic_schur(std::move(algebraic_schur));
        if (!schur->setup(A)) {
            throw std::runtime_error(
                std::string("N9 LSC/BFBt Schur setup failed: ") +
                schur->last_error());
        }
        coupled_preconditioner = std::move(schur);
    } else if (use_pcd) {
        auto pcd_ops = std::make_shared<PcdPressureOperators>(
            assemble_pcd_pressure_operators(
                mesh, geometry, U_old, rho, kinematic_viscosity, velocity_bcs, pressure_bcs,
                reference_cell, ConvectionScheme::UPWIND/* exact matrix for every scheme, see assemble_pcd_pressure_operators */));

        const auto prepare_pcd_solver =
            [](const cfdx::core::SparseMatrix& matrix,
               const cfdx::core::LinearSolverRequest& request) {
                if (request.krylov == cfdx::core::KrylovModel::Auto ||
                    request.preconditioner == cfdx::core::PreconditionerModel::Auto)
                    throw std::invalid_argument(
                        "PCD nested solvers require explicit Krylov and preconditioner models");
                auto preconditioner =
                    cfdx::core::make_scalar_preconditioner(request.preconditioner);
                if (!preconditioner || !preconditioner->setup(matrix))
                    throw std::runtime_error(
                        "PCD nested pressure preconditioner setup failed");
                return preconditioner;
            };

        pcd_ops->laplacian_preconditioner =
            prepare_pcd_solver(pcd_ops->laplacian, controls.pcd_laplacian_linear_solver);
        pcd_ops->convection_diffusion_preconditioner =
            prepare_pcd_solver(
                pcd_ops->convection_diffusion,
                controls.pcd_convection_diffusion_linear_solver);

        const auto solve_kp = [pcd_ops, controls](
            const Vector& rhs, Vector& x) {
            const auto& matrix = pcd_ops->laplacian;
            auto* preconditioner = pcd_ops->laplacian_preconditioner.get();
            const int restart = std::min<int>(
                controls.pcd_laplacian_linear_solver.gmres_restart,
                static_cast<int>(matrix.n_rows()));
            cfdx::core::SolverResult result;
            if (controls.pcd_laplacian_linear_solver.krylov ==
                cfdx::core::KrylovModel::CG) {
                result = solve_cg_prepared(
                    matrix, rhs, x, *preconditioner,
                    controls.pcd_laplacian_max_iterations,
                    controls.pcd_laplacian_tolerance);
            } else if (controls.pcd_laplacian_linear_solver.krylov ==
                       cfdx::core::KrylovModel::GMRES) {
                result = solve_gmres(
                    matrix, rhs, x, restart,
                    controls.pcd_laplacian_max_iterations,
                    controls.pcd_laplacian_tolerance,
                    preconditioner, {}, {}, false);
            } else if (controls.pcd_laplacian_linear_solver.krylov ==
                       cfdx::core::KrylovModel::FGMRES) {
                result = solve_fgmres(
                    matrix, rhs, x, restart,
                    controls.pcd_laplacian_max_iterations,
                    controls.pcd_laplacian_tolerance,
                    preconditioner, {}, {}, false);
            } else {
                throw std::invalid_argument(
                    "unsupported PCD Laplacian Krylov model");
            }
            return nested_schur_pressure_solve_usable(
                result, x, matrix.n_rows());
        };

        const auto solve_fp = [pcd_ops, controls](
            const Vector& rhs, Vector& x) {
            const auto& matrix = pcd_ops->convection_diffusion;
            auto* preconditioner = pcd_ops->convection_diffusion_preconditioner.get();
            const int restart = std::min<int>(
                controls.pcd_convection_diffusion_linear_solver.gmres_restart,
                static_cast<int>(matrix.n_rows()));
            cfdx::core::SolverResult result;
            if (controls.pcd_convection_diffusion_linear_solver.krylov ==
                cfdx::core::KrylovModel::GMRES) {
                result = solve_gmres(
                    matrix, rhs, x, restart,
                    controls.pcd_convection_diffusion_max_iterations,
                    controls.pcd_convection_diffusion_tolerance,
                    preconditioner, {}, {}, false);
            } else if (controls.pcd_convection_diffusion_linear_solver.krylov ==
                       cfdx::core::KrylovModel::FGMRES) {
                result = solve_fgmres(
                    matrix, rhs, x, restart,
                    controls.pcd_convection_diffusion_max_iterations,
                    controls.pcd_convection_diffusion_tolerance,
                    preconditioner, {}, {}, false);
            } else {
                throw std::invalid_argument(
                    "PCD convection-diffusion solve requires GMRES or FGMRES");
            }
            if (result.status != SolverStatus::CONVERGED &&
                controls.diagnostics.coupled_matrix_summary) {
                std::cerr << "COUPLED_NESTED_PCD_BEST_EFFORT operator=Fp"
                          << " status=" << static_cast<int>(result.status)
                          << " iterations=" << result.iterations
                          << " residual=" << result.residual
                          << " relative=" << result.residual_relative
                          << "\n";
            }
            return nested_schur_pressure_solve_usable(
                result, x, matrix.n_rows());
        };

        auto pcd = std::make_unique<PcdSchurApproximation>(
            pcd_ops->mass, pcd_ops->laplacian, pcd_ops->convection_diffusion,
            std::move(solve_kp), std::move(solve_fp),
            std::nullopt,
            has_fixed_pressure ? std::optional<std::size_t>()
                               : std::optional<std::size_t>(reference_cell));
        // Evidence-only: expose the production operators to the
        // coupled_assembly_callback when one is installed (kept alive via
        // evidence_pcd_ops past this branch scope).
        if (controls.coupled_assembly_callback)
            evidence_pcd_ops = pcd_ops;

        CoupledBlockSchurOptions options;
        options.factorization = CoupledSchurFactorization::Full;
        options.velocity_approximation = CoupledSchurVelocityApproximation::Block;
        switch (schur_model) {
            case CoupledSchurModel::BlockLocal:
                options.schur_approximation = CoupledSchurApproximationModel::BlockLocal;
                break;
            case CoupledSchurModel::PCD:
                options.schur_approximation = CoupledSchurApproximationModel::PCD;
                break;
            case CoupledSchurModel::LSC:
            case CoupledSchurModel::BFBT:
            case CoupledSchurModel::SIMPLE:
            case CoupledSchurModel::SIMPLEC:
                throw std::invalid_argument(
                    "selected algebraic Schur model is handled by its dedicated production path");
            default:
                throw std::invalid_argument("unsupported coupled Schur model");
        }

        auto schur = std::make_unique<CoupledBlockSchurAMGPreconditioner>(
            nc, options, has_fixed_pressure
                ? std::numeric_limits<std::size_t>::max() : reference_cell);
        if (schur_model == CoupledSchurModel::PCD)
            schur->set_pcd_schur(std::move(pcd));
        // Set the preconditioner up here for the same reason the BlockSchur
        // branch below does. solve_fgmres also calls setup() on its
        // preconditioner before the first apply, so this is not what makes PCD
        // reachable; performing it explicitly is what makes a setup failure
        // report its own cause. Without it, a rejected PCD setup surfaces only
        // as SolverStatus::NOT_APPLICABLE from the Krylov driver and the
        // PcdSchurApproximation last_error() is lost. reset() preserves
        // pcd_schur_, so the driver's second setup is clean.
        if (!schur->setup(A)) {
            throw std::runtime_error(
                std::string("N8 coupled PCD Schur setup failed: ") +
                schur->last_error());
        }
        coupled_preconditioner = std::move(schur);
    } else if (use_n8_block_schur) {
        auto schur_amg = std::make_unique<CoupledBlockSchurAMGPreconditioner>(
            nc, CoupledBlockSchurOptions{}, has_fixed_pressure
                ? std::numeric_limits<std::size_t>::max() : reference_cell);
        if (!schur_amg->setup(A)) {
            throw std::runtime_error(
                std::string("N8 coupled BlockSchur AMG setup failed: ") +
                schur_amg->last_error());
        }
        coupled_preconditioner = std::move(schur_amg);
    } else {
        std::vector<std::size_t> fine_variables(nv);
        std::vector<std::size_t> coarse_variables(np);
        std::iota(fine_variables.begin(), fine_variables.end(), 0);
        std::iota(coarse_variables.begin(), coarse_variables.end(), nv);
        coupled_preconditioner =
            std::make_unique<NativeMGRPreconditioner>(
                std::move(fine_variables), std::move(coarse_variables));
    }

    const bool automatic_coupled = solver_request.krylov == KrylovModel::Auto &&
        solver_request.preconditioner == PreconditionerModel::Auto;
    if (resolved_linear_plan != nullptr)
        *resolved_linear_plan = solver_plan;
    const int gmres_restart = automatic_coupled
        ? static_cast<int>(A.n_rows() <= 4096
                               ? std::min<std::size_t>(A.n_rows(),
                                                       std::max<std::size_t>(1024, max_iterations))
                               : std::min<std::size_t>(A.n_rows(), 512))
        : std::min<int>(solver_request.gmres_restart, static_cast<int>(A.n_rows()));
    // The coupled acceptance problem is small enough for a full Krylov space.
    // Do not let the global adaptive-restart policy silently clamp the restart
    // back to its generic [10,40] range: that was the direct cause of the
    // observed residual plateau around 1.33e-6. The restart must also be the
    // FULL space for small/medium coupled systems (<= 1024 unknowns): the
    // 12x16 channel system has 768 unknowns, so a 512 cap truncated the
    // Krylov space below the system dimension and manufactured a plateau
    // near 1e-7 on the Poiseuille Picard system, converting a solvable
    // 1e-12 gate into an artificial MAX_ITER failure. Medium coupled systems
    // (1024 < n <= 4096) get a restart of max(1024, max_iterations), capped
    // by the system dimension: the 16x32 cavity Picard system has 4096
    // unknowns, and a 512 cap truncated the Krylov space mid-budget
    // (1000 iterations = 512 + 488), stalling the true relative residual at
    // 3.66e-12 against the 1e-12 gate with three identical relaxation
    // retries. A 1024 restart is never truncated within the iteration
    // budget, so convergence is decided by the budget, not by a restart
    // boundary. Systems above 4096 unknowns keep the 512 cap for memory.
    // Keep the restart fixed for this block solve; adaptive restart remains
    // available to generic GMRES callers.
    cfdx::core::KrylovControls coupled_gmres_controls;
    coupled_gmres_controls.restart_min = gmres_restart;
    coupled_gmres_controls.restart_max = gmres_restart;
    coupled_gmres_controls.adaptive_restart = false;
    if (diagnostics.coupled_matrix_summary) {
        std::cerr << "COUPLED_PRECONDITIONER name=" << coupled_preconditioner->name()
                  << " restart=" << gmres_restart
                  << " adaptive_restart=0\n";
    }

    // Evidence-only capture hook, at the exact Krylov call point (null by
    // default; see IncompressibleSolverControls::coupled_assembly_callback).
    // Everything the solver consumes is now final: assembly, gauge
    // treatment, finite checks and preconditioner setup are done, and no
    // mutation happens between this hook and the solve calls below. The
    // snapshot therefore aliases the very A, b and x objects handed to the
    // Krylov driver - the captured system IS the solved system, by
    // construction rather than by convention. The consumer must still copy
    // anything it wants to keep: the references die with this scope, and x
    // becomes the solution below.
    if (controls.coupled_assembly_callback) {
        const CoupledAssemblySnapshot snapshot{
            &A, &b, &x,
            &coupled_pressure_schur_diagonal,
            &pressure_schur_diagonal_raw_copy,
            evidence_pcd_ops ? &evidence_pcd_ops->mass : nullptr,
            evidence_pcd_ops ? &evidence_pcd_ops->laplacian : nullptr,
            evidence_pcd_ops ? &evidence_pcd_ops->convection_diffusion : nullptr,
            &U_old, reference_cell, has_fixed_pressure};
        controls.coupled_assembly_callback(snapshot);
    }

    SolverResult result;
    if (solver_plan.krylov == KrylovModel::FGMRES) {
        result = solve_fgmres(
            A, b, x, gmres_restart, max_iterations, tolerance,
            coupled_preconditioner.get(), coupled_gmres_controls);
    } else {
        result = solve_gmres(
            A, b, x, gmres_restart, max_iterations, tolerance,
            coupled_preconditioner.get(), coupled_gmres_controls);
    }

    // Keep coupled-solver failures visible. Do not replace a failed Schur
    // setup by identity-preconditioned GMRES: that would hide the defect in
    // the production preconditioner and make the numerical gate non-diagnostic.
    if (result.status != SolverStatus::CONVERGED)
        return result;

    // Independently recompute ||b-Ax|| after the linear solve.  The Krylov
    // result already carries a true residual, but keeping an explicit
    // post-solve matrix check here makes the coupled algebraic gate independent
    // from the nonlinear FV diagnostics and exposes any stale/inconsistent
    // solver bookkeeping.
    double coupled_matrix_residual_sq = 0.0;
    double coupled_rhs_norm_sq = 0.0;
    for (std::size_t row = 0; row < A.n_rows(); ++row) {
        double ax = 0.0;
        const auto rb = A.row_offsets_data()[row];
        const auto re = A.row_offsets_data()[row + 1];
        for (std::size_t k = rb; k < re; ++k)
            ax += A.values_data()[k] * x(A.columns_data()[k]);
        const double r = b(row) - ax;
        coupled_matrix_residual_sq += r * r;
        coupled_rhs_norm_sq += b(row) * b(row);
    }
    const double coupled_matrix_residual = std::sqrt(coupled_matrix_residual_sq);
    const double coupled_rhs_norm = std::sqrt(coupled_rhs_norm_sq);
    const double coupled_matrix_relative =
        coupled_matrix_residual / (coupled_rhs_norm > 0.0 ? coupled_rhs_norm : 1.0);
    if (!std::isfinite(coupled_matrix_residual) ||
        !std::isfinite(coupled_matrix_relative)) {
        result.status = SolverStatus::DIVERGED;
        result.residual = std::numeric_limits<double>::infinity();
        result.residual_relative = std::numeric_limits<double>::infinity();
        return result;
    }
    // The independently recomputed residual is the authoritative algebraic
    // gate. Never accept a Krylov-reported convergence if b-Ax is still above
    // the requested coupled linear tolerance.
    result.residual = coupled_matrix_residual;
    result.residual_relative = coupled_matrix_relative;
    if (coupled_matrix_relative > tolerance) {
        result.status = SolverStatus::DIVERGED;
        return result;
    }
    if (diagnostics.coupled_matrix_summary) {
        std::cerr << "COUPLED_LINEAR_TRUE_RESIDUAL norm="
                  << coupled_matrix_residual
                  << " relative=" << coupled_matrix_relative
                  << " gmres_reported=" << result.residual
                  << "\n";
    }

    for (std::size_t c = 0; c < nc; ++c) {
        U.component_data(0)[c] = x(c);
        U.component_data(1)[c] = x(nc + c);
        U.component_data(2)[c] = x(2*nc + c);
        p(c) = x(nv + c);
    }
    return result;
}

inline IncompressibleSolveResult solve_steady_incompressible(
    const cfdx::core::Mesh& mesh,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    const VelocityBoundaryConditions& velocity_bcs,
    const ScalarBoundaryConditions& pressure_bcs,
    const IncompressibleSolverControls& controls = {},
    const std::string& restart_path = {},
    cfdx::io::DatRestartFields restart_fields = {})
{
    using namespace cfdx::core;
    if (U.dimension() != 3 || U.size() != mesh.n_cells() ||
        p.dimension() != 1 || p.size() != mesh.n_cells())
        throw std::invalid_argument("solve_steady_incompressible: invalid fields");

    // Validate all controls and derived material properties before initialization
    // can mutate the caller-provided solution fields.
    validate_incompressible_controls(controls, mesh.n_cells());
    // Validate the caller-owned SST bundle before initialization can mutate
    // the carrier fields. Derived SST quantities must be refreshed only after
    // initialization/restart has established the actual velocity field.
    if (controls.sst_runtime)
        validate_sst_production_state(*controls.sst_runtime->state, mesh.n_cells());

    apply_initialization(
        controls.initialization, mesh, U, p, restart_fields, restart_path);
    if (controls.initialized_state_callback)
        controls.initialized_state_callback(mesh, U, p);

    // Authoritative dynamic viscosity used by predictor assembly, nonlinear
    // acceptance residuals, and the independent momentum conservation audit.
    std::vector<double> mu_eff(mesh.n_cells());
    if (controls.sst_runtime) {
        std::map<std::string, std::array<double, 3>> fixed_velocity_values;
        for (const auto& [patch, bc] : velocity_bcs) {
            if (bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                fixed_velocity_values.emplace(
                    patch, std::array<double, 3>{bc.value.x, bc.value.y, bc.value.z});
        }
        // Initialization may replace U (uniform, profile, or DAT restart).
        // Refresh wall distance, strain and nu_t against that final field so
        // the very first momentum predictor cannot use stale SST viscosity.
        refresh_sst_production_derived_state_from_mesh(
            *controls.sst_runtime->state, mesh, controls.sst_runtime->transport,
            U, fixed_velocity_values, controls.sst_runtime->wall_patch_names);
    }
    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
        const double nu_t = controls.sst_runtime
            ? controls.sst_runtime->state->nu_t(cell)
            : (controls.cell_turbulent_viscosity.empty()
                ? controls.turbulent_viscosity
                : controls.cell_turbulent_viscosity[cell]);
        mu_eff[cell] = controls.density * (controls.kinematic_viscosity + nu_t);
        if (!(mu_eff[cell] >= 0.0) || !std::isfinite(mu_eff[cell]))
            throw std::invalid_argument(
                "effective cellwise dynamic viscosity must be finite and non-negative");
    }

    const FvGeometry geometry = build_fv_geometry(mesh);
    cfdx::core::ConvergenceMonitorControls monitor_controls;
    monitor_controls.residual.absolute_tolerance = controls.convergence.absolute_tolerance;
    monitor_controls.residual.relative_tolerance = controls.convergence.relative_tolerance;
    monitor_controls.conservation_tolerance = controls.convergence.continuity_tolerance;
    monitor_controls.minimum_iterations = std::min<std::size_t>(2, controls.convergence.max_iterations);
    monitor_controls.max_iterations = controls.convergence.max_iterations;
    // Conservative defaults: stagnation is diagnostic, not a shortcut around
    // the configured residual/conservation gates. Divergence is detected only
    // after repeated growth and never by relaxing a tolerance.
    monitor_controls.stagnation_window = 10;
    monitor_controls.stagnation_relative_improvement =
        controls.convergence.stagnation_relative_improvement;
    monitor_controls.divergence_window = 3;
    monitor_controls.divergence_growth_factor = 10.0;
    monitor_controls.qoi_gates = controls.qoi_gates;
    cfdx::core::ConvergenceMonitor convergence_monitor(monitor_controls);

    bool has_fixed_pressure_boundary = false;
    for (const auto& [name, bc] : pressure_bcs) {
        (void)name;
        if (bc.type == ScalarBoundaryType::FIXED_VALUE) {
            has_fixed_pressure_boundary = true;
            break;
        }
    }

    const std::size_t configured_pcorr =
        (controls.algorithm == PressureVelocityAlgorithm::PISO ||
         controls.algorithm == PressureVelocityAlgorithm::PIMPLE)
            ? static_cast<std::size_t>(controls.coupling.n_pressure_correctors)
            : (controls.algorithm == PressureVelocityAlgorithm::FRACTIONAL_STEP
                ? static_cast<std::size_t>(controls.coupling.n_fractional_steps)
                : 1u);
    const std::size_t minimum_outer_correctors =
        controls.algorithm == PressureVelocityAlgorithm::PIMPLE
            ? static_cast<std::size_t>(controls.coupling.n_outer_correctors)
            : 1u;

    IncompressibleSolveResult result;
    result.requested_coupled_schur_model = controls.coupling.schur_model;
    result.resolved_coupled_schur_model = controls.coupling.schur_model;
    result.coupled_schur_model_resolved = true;
    LinearSolverPlan pressure_plan;
    std::unique_ptr<Preconditioner> pressure_preconditioner;
    std::unique_ptr<ReusableCgContext> pressure_context;
    std::unique_ptr<NullSpaceProjector> pressure_null_space;
    const bool use_projected_pressure_null_space =
        !has_fixed_pressure_boundary &&
        controls.pressure_linear_solver.null_space ==
            NullSpaceModel::Constant;
    const bool reduce_pressure_gauge =
        !has_fixed_pressure_boundary &&
        !use_projected_pressure_null_space &&
        mesh.n_cells() > 1;
    if (controls.algorithm != PressureVelocityAlgorithm::COUPLED) {
        pressure_plan = select_linear_solver(
            LinearProblemKind::PressurePoisson,
            reduce_pressure_gauge ? mesh.n_cells() - 1 : mesh.n_cells(),
            controls.pressure_linear_solver);
        result.pressure_linear_plan = pressure_plan;
        result.pressure_linear_plan_resolved = true;
        if (pressure_plan.krylov == KrylovModel::CG) {
            pressure_preconditioner =
                make_scalar_preconditioner(
                    pressure_plan.preconditioner,
                    pressure_plan.null_space == NullSpaceModel::Constant);
            pressure_context = std::make_unique<ReusableCgContext>(
                *pressure_preconditioner);
            if (pressure_plan.null_space == NullSpaceModel::Constant) {
                pressure_null_space = std::make_unique<NullSpaceProjector>(
                    NullSpaceProjector::constant(mesh.n_cells()));
            }
        }
    }

    auto build_hbya = [&](const ScalarEquation& eq,
                          const Vector& field,
                          const Field<double, Location::CELL>& gradp,
                          std::size_t component,
                          const std::vector<double>& rAU) {
        std::vector<double> hbya(mesh.n_cells(), 0.0);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            double h = eq.rhs(c) + gradp.component_data(component)[c] * geometry.cell_volumes[c];
            const auto begin = eq.matrix.row_offsets_data()[c];
            const auto end = eq.matrix.row_offsets_data()[c + 1];
            for (std::uint32_t k = begin; k < end; ++k) {
                const auto col = eq.matrix.columns_data()[k];
                if (col != c)
                    h -= eq.matrix.values_data()[k] * field(col);
            }
            // HbyA is always the plain momentum-inverse combination
            // rAU*H, for SIMPLEC included. The SIMPLEC consistent
            // operator (rAtU) belongs to the pressure equation and its
            // face fluxes only: the over-relaxation terms added to
            // phiHbyA and to the corrected flux are sized so that
            // (dAtU - dAU) cancels the rAtU response of the pressure
            // operator, which collapses the SIMPLEC fixed point back
            // onto the plain Rhie-Chow state built with rAU. Folding
            // (rAtU - rAU)*grad(p) into the cell HbyA instead leaves
            // the converged state on a displaced continuity operator
            // carrying the spurious source
            // div(rho*(dAtU - dAU)*snGrad(p)*|Sf|) of O(Laplacian(p))
            // scale: invisible on Couette (constant p) and Poiseuille
            // (linear p) but far outside the algorithm-equivalence
            // gate on the lid-driven cavity.
            hbya[c] = rAU[c] * h;
            if (!std::isfinite(hbya[c]))
                throw std::runtime_error("solve_steady_incompressible: non-finite HbyA");
        }
        return hbya;
    };

    auto make_hbya_field = [&](const std::array<std::vector<double>,3>& hbya) {
        Field<double, Location::CELL> field(
            mesh.n_cells(), "HbyA", "m/s", 3);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            field.component_data(0)[c] = hbya[0][c];
            field.component_data(1)[c] = hbya[1][c];
            field.component_data(2)[c] = hbya[2][c];
        }
        return field;
    };

    // Face pressure projection coefficient.
    //
    // Momentum diagonals are integrated finite-volume coefficients A_P [N s/m].
    // rAU = 1/A_P therefore has units of m/(N s), while the velocity response
    // to a pressure gradient is dAU = V*rAU.  The same dAU operator must be
    // used by cell velocity reconstruction, Rhie-Chow face fluxes and the
    // pressure-correction matrix.
    auto directional_face_coefficient =
        [&](const std::array<std::vector<double>,3>& coeff,
            std::size_t o, std::size_t n,
            const Vec3& Sf, const Vec3& dvec) {
            const double cx = 0.5 * (
                geometry.cell_volumes[o] * coeff[0][o] +
                geometry.cell_volumes[n] * coeff[0][n]);
            const double cy = 0.5 * (
                geometry.cell_volumes[o] * coeff[1][o] +
                geometry.cell_volumes[n] * coeff[1][n]);
            const double cz = 0.5 * (
                geometry.cell_volumes[o] * coeff[2][o] +
                geometry.cell_volumes[n] * coeff[2][n]);
            return pressure_face_orthogonal_coefficient(Sf, dvec, {cx, cy, cz});
        };

    Field<double, Location::FACE> mass_flux =
        make_mass_flux(mesh, geometry, U, controls.density, velocity_bcs);

    const bool freeze_state = controls.diagnostics.freeze_state_probe;
    Field<double, Location::FACE> phi_used_in_momentum;
    std::array<std::vector<double>, 3> frozen_rAU;
    std::array<std::vector<double>, 3> frozen_hbyA;
    Field<double, Location::CELL> frozen_grad_p;
    bool frozen_state_valid = false;

    double effective_alpha_u = controls.coupling.alpha_u;
    double effective_alpha_p = controls.coupling.alpha_p;
    std::size_t adaptive_u_improvement_streak = 0;
    std::size_t adaptive_u_degradation_streak = 0;
    std::size_t adaptive_p_improvement_streak = 0;
    std::size_t adaptive_p_degradation_streak = 0;
    std::size_t adaptive_u_recovery_cooldown = 0;
    std::size_t adaptive_p_recovery_cooldown = 0;

    for (std::size_t iter = 1; iter <= controls.convergence.max_iterations; ++iter) {
        bool stop_after_iteration = false;
        NonlinearRetryController retry_controller(controls.nonlinear_retry);
        NonlinearStateRollback transaction(U, p);
        const auto mass_flux_snapshot = mass_flux;
        bool iteration_completed = false;
        while (!iteration_completed) {
            // A retry must restore not only physical fields but also any
            // convergence/history state tentatively updated before the outer
            // transaction commits.
            stop_after_iteration = false;
            transaction.begin();
            const auto convergence_monitor_snapshot = convergence_monitor;
            const auto result_history_size = result.history.size();
            const auto probe_samples_size = result.probe_samples.size();
            const auto previous_result_converged = result.converged;
            const auto previous_result_status = result.convergence_status;
            const auto previous_result_reason = result.convergence_reason;
            const auto previous_result_iterations = result.iterations;
            const auto sst_attempt_snapshot = controls.sst_runtime
                ? std::optional<SSTProductionState>(
                    snapshot_sst_production_state(*controls.sst_runtime->state))
                : std::nullopt;
            try {
        double linear_tolerance = controls.linear_tolerance;
        if (controls.acceleration.adaptive_linear_tolerance) {
            const auto& forcing = controls.acceleration.linear_forcing;
            if (result.history.empty()) {
                linear_tolerance = forcing.eta_max;
            } else {
                const auto nonlinear_residual = std::max(
                    result.history.back().momentum_residual,
                    result.history.back().continuity_normalized);
                const auto previous_residual = result.history.size() > 1
                    ? std::max(result.history[result.history.size()-2].momentum_residual,
                               result.history[result.history.size()-2].continuity_normalized)
                    : std::max(nonlinear_residual, forcing.residual_floor);
                linear_tolerance = nonlinear_forcing_tolerance(
                    nonlinear_residual,
                    std::max(previous_residual, forcing.residual_floor), forcing);
            }
        }

        if (controls.adaptive_relaxation.enabled) {
            // Adapt the momentum and pressure relaxation independently. SIMPLE
            // couples the two equations, but their residuals measure different
            // numerical mechanisms: momentum relaxation is controlled from the
            // momentum-equation residual, while pressure relaxation is
            // controlled from the normalized continuity imbalance. A single
            // aggregate metric must not force both factors in the same direction.
            const auto window = controls.adaptive_relaxation.adaptation_window;
            if (window > 0 && result.history.size() > window &&
                (result.history.size() - 1) % window == 0) {
                const auto& current = result.history.back();
                const auto& reference = result.history[result.history.size() - 1 - window];

                if (adaptive_u_recovery_cooldown > 0) {
                    --adaptive_u_recovery_cooldown;
                    adaptive_u_improvement_streak = 0;
                    adaptive_u_degradation_streak = 0;
                } else if (std::isfinite(current.momentum_equation_residual_relative) &&
                           std::isfinite(reference.momentum_equation_residual_relative) &&
                           reference.momentum_equation_residual_relative > 0.0) {
                    effective_alpha_u = adapt_relaxation_factor_windowed(
                        effective_alpha_u,
                        reference.momentum_equation_residual_relative,
                        current.momentum_equation_residual_relative,
                        controls.adaptive_relaxation.min_alpha_u,
                        controls.adaptive_relaxation.max_alpha_u,
                        controls.adaptive_relaxation,
                        adaptive_u_improvement_streak,
                        adaptive_u_degradation_streak);
                } else {
                    adaptive_u_improvement_streak = 0;
                    adaptive_u_degradation_streak = 0;
                }

                if (adaptive_p_recovery_cooldown > 0) {
                    --adaptive_p_recovery_cooldown;
                    adaptive_p_improvement_streak = 0;
                    adaptive_p_degradation_streak = 0;
                } else if (std::isfinite(current.continuity_normalized) &&
                           std::isfinite(reference.continuity_normalized) &&
                           reference.continuity_normalized > 0.0) {
                    effective_alpha_p = adapt_relaxation_factor_windowed(
                        effective_alpha_p,
                        reference.continuity_normalized,
                        current.continuity_normalized,
                        controls.adaptive_relaxation.min_alpha_p,
                        controls.adaptive_relaxation.max_alpha_p,
                        controls.adaptive_relaxation,
                        adaptive_p_improvement_streak,
                        adaptive_p_degradation_streak);
                } else {
                    adaptive_p_improvement_streak = 0;
                    adaptive_p_degradation_streak = 0;
                }
            }
        }

        std::size_t pressure_correctors = configured_pcorr;
        const bool supports_multiple_pressure_corrections =
            controls.algorithm == PressureVelocityAlgorithm::PISO ||
            controls.algorithm == PressureVelocityAlgorithm::PIMPLE ||
            controls.algorithm == PressureVelocityAlgorithm::FRACTIONAL_STEP;
        if (supports_multiple_pressure_corrections &&
            controls.acceleration.adaptive_pressure_correctors) {
            pressure_correctors = result.history.empty()
                ? controls.acceleration.pressure_correctors.max_correctors
                : choose_pressure_correctors(
                      result.history.back().continuity_normalized,
                      controls.acceleration.pressure_correctors);
        }
        const auto U_old = U;
        const auto p_old = p;

        // Exact conservative flux consumed by the momentum assembly.
        if (freeze_state)
            phi_used_in_momentum = mass_flux;

        auto grad_p = gauss_gradient_with_boundary(p, mesh, geometry, pressure_bcs);

        Field<double, Location::CELL> body_x(mesh.n_cells(), "body_x", "N/m3", 1);
        Field<double, Location::CELL> body_y(mesh.n_cells(), "body_y", "N/m3", 1);
        Field<double, Location::CELL> body_z(mesh.n_cells(), "body_z", "N/m3", 1);
        const bool cellwise_body_force = !controls.cell_body_force.empty();
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            const auto& cell_force = cellwise_body_force
                ? controls.cell_body_force[c]
                : controls.body_force;
            body_x(c) = controls.density * cell_force.x;
            body_y(c) = controls.density * cell_force.y;
            body_z(c) = controls.density * cell_force.z;
        }

        ScalarBoundaryConditions ubc_x, ubc_y, ubc_z;
        Field<double, Location::CELL> u_x_field(mesh.n_cells(), "Ux_reconstruction", "m/s", 1);
        Field<double, Location::CELL> u_y_field(mesh.n_cells(), "Uy_reconstruction", "m/s", 1);
        Field<double, Location::CELL> u_z_field(mesh.n_cells(), "Uz_reconstruction", "m/s", 1);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            u_x_field(c) = U.component_data(0)[c];
            u_y_field(c) = U.component_data(1)[c];
            u_z_field(c) = U.component_data(2)[c];
        }
        for (const auto& [name, bc] : velocity_bcs) {
            ubc_x[name] = {bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE
                               ? ScalarBoundaryType::FIXED_VALUE : ScalarBoundaryType::ZERO_GRADIENT,
                           bc.value.x, 0.0};
            ubc_y[name] = {bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE
                               ? ScalarBoundaryType::FIXED_VALUE : ScalarBoundaryType::ZERO_GRADIENT,
                           bc.value.y, 0.0};
            ubc_z[name] = {bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE
                               ? ScalarBoundaryType::FIXED_VALUE : ScalarBoundaryType::ZERO_GRADIENT,
                           bc.value.z, 0.0};
        }

        // Momentum predictor uses the conservative corrected flux from the
        // previous outer iteration. This preserves the SIMPLE phi state.

        auto ex = assemble_momentum_component(
            mesh, geometry, mass_flux, grad_p, body_x, mu_eff, ubc_x, 0,
            controls.use_bounded_convection, controls.convection_scheme, &u_x_field);
        auto ey = assemble_momentum_component(
            mesh, geometry, mass_flux, grad_p, body_y, mu_eff, ubc_y, 1,
            controls.use_bounded_convection, controls.convection_scheme, &u_y_field);
        auto ez = assemble_momentum_component(
            mesh, geometry, mass_flux, grad_p, body_z, mu_eff, ubc_z, 2,
            controls.use_bounded_convection, controls.convection_scheme, &u_z_field);

        if (controls.transient.enabled) {
            const auto& u_n = *controls.transient.previous;
            const auto& u_nm1 = controls.transient.history_valid
                ? *controls.transient.previous_previous : u_n;
            const bool use_bdf2 = controls.transient.scheme == cfdx::core::TimeScheme::BDF2 &&
                                  controls.transient.history_valid;
            add_dual_time_mass_term(ex, u_n, u_nm1, geometry, controls.density,
                                    controls.transient.dt, controls.transient.dt_previous,
                                    use_bdf2, 0);
            add_dual_time_mass_term(ey, u_n, u_nm1, geometry, controls.density,
                                    controls.transient.dt, controls.transient.dt_previous,
                                    use_bdf2, 1);
            add_dual_time_mass_term(ez, u_n, u_nm1, geometry, controls.density,
                                    controls.transient.dt, controls.transient.dt_previous,
                                    use_bdf2, 2);
        }

        // Preserve the pre-solve velocity for equation relaxation. HbyA is
        // reconstructed after the momentum solve from the solved neighbour
        // values, which is the algebraic H/A split of the assembled equation.
        // The pre-solve state is used only for the under-relaxation source term.
        Vector ux(mesh.n_cells()), uy(mesh.n_cells()), uz(mesh.n_cells());
        Vector ux_old(mesh.n_cells()), uy_old(mesh.n_cells()), uz_old(mesh.n_cells());
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            ux(c) = U.component_data(0)[c];
            uy(c) = U.component_data(1)[c];
            uz(c) = U.component_data(2)[c];
            ux_old(c) = ux(c);
            uy_old(c) = uy(c);
            uz_old(c) = uz(c);
        }

        cfdx::core::SolverResult rx{}, ry{}, rz{};
        double pressure_residual = std::numeric_limits<double>::infinity();
        std::size_t pressure_iterations = 0;
        std::array<double, 3> simplec_signed_offdiag_sum_audit{0.0, 0.0, 0.0};
        std::array<double, 3> simplec_abs_offdiag_sum_audit{0.0, 0.0, 0.0};
        std::array<double, 3> simplec_denominator_signed_audit{0.0, 0.0, 0.0};
        std::array<double, 3> simplec_denominator_abs_audit{0.0, 0.0, 0.0};
        std::array<double, 3> simplec_denominator_gap_audit{0.0, 0.0, 0.0};
        std::array<std::size_t, 3> simplec_denominator_worst_cell_audit{0, 0, 0};
        // Rhie-Chow operator identity audit (common scope): each algorithm
        // branch evaluates it right after its authoritative flux rebuild.
        // Evidence channel only; no acceptance criterion consumes it.
        bool rhie_chow_identity_evaluated = false;
        double rhie_chow_identity_audit_linf = 0.0;
        if (controls.algorithm == PressureVelocityAlgorithm::COUPLED) {
            const auto coupled_result = solve_coupled_momentum_continuity(
                mesh, geometry, ex, ey, ez, U_old, p_old,
                velocity_bcs, pressure_bcs, controls.density,
                controls.kinematic_viscosity,
                controls.pressure_reference_cell,
                controls.pressure_reference_value,
                controls.coupling.coupled_max_iterations,
                controls.coupling.coupled_linear_tolerance,
                controls.diagnostics,
                controls.coupling.schur_model,
                controls.coupled_linear_solver,
                controls,
                retry_controller.alpha_u(effective_alpha_u),
                U, p, &result.coupled_linear_plan);
            result.coupled_linear_plan_resolved = true;
            if (coupled_result.status != cfdx::core::SolverStatus::CONVERGED) {
                // Report the stalled linear state in scientific notation.
                // std::to_string truncates to six fixed decimals, so a
                // MAX_ITER stall printed as "residual=0.000000" hid the
                // plateau level and made the failure impossible to calibrate
                // against the coupled linear gate or the nested PCD budget.
                const auto coupled_linear_failure_value = [](double value) {
                    std::ostringstream out;
                    out << std::scientific << std::setprecision(6) << value;
                    return out.str();
                };
                throw NonlinearRetryableFailure(
                    "solve_steady_incompressible: coupled momentum-continuity solve did not converge "
                    "(status=" + std::to_string(static_cast<int>(coupled_result.status)) +
                    ", iterations=" + std::to_string(coupled_result.iterations) +
                    ", residual=" + coupled_linear_failure_value(coupled_result.residual) +
                    ", relative=" + coupled_linear_failure_value(coupled_result.residual_relative) + ")");
            }

            // The coupled linear system is a Picard/Newton linearization of the
            // current nonlinear iterate. Keep the nonlinear field update explicit
            // and bounded; there is no hidden scalar fallback.
            for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
                for (std::size_t d = 0; d < 3; ++d)
                    U.component_data(d)[c] =
                        U_old.component_data(d)[c] +
                        retry_controller.alpha_u(effective_alpha_u) *
                        (U.component_data(d)[c] - U_old.component_data(d)[c]);
                p(c) = p_old(c) +
                    retry_controller.alpha_p(effective_alpha_p) * (p(c) - p_old(c));
            }
            apply_pressure_gauge_policy(
                p, geometry, controls.pressure_gauge_policy,
                controls.pressure_reference_cell, controls.pressure_reference_value,
                has_fixed_pressure_boundary);

            // The coupled solve has already enforced the discrete continuity
            // equation. Rebuild the authoritative flux from the same RC
            // interpolation used by the block matrix for diagnostics and for
            // the next nonlinear iteration.
            std::array<std::vector<double>,3> coupled_rAU;
            coupled_rAU[0].resize(mesh.n_cells());
            coupled_rAU[1].resize(mesh.n_cells());
            coupled_rAU[2].resize(mesh.n_cells());
            for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
                if (!(ex.diagonal[c] > 0.0) || !std::isfinite(ex.diagonal[c]) ||
                    !(ey.diagonal[c] > 0.0) || !std::isfinite(ey.diagonal[c]) ||
                    !(ez.diagonal[c] > 0.0) || !std::isfinite(ez.diagonal[c]))
                    throw std::runtime_error(
                        "solve_steady_incompressible: invalid coupled momentum diagonal after solve");
                // make_rhie_chow_mass_flux() and rhie_chow_pressure_flux_internal()
                // take the INVERSE momentum diagonal rAU = 1/A_P and apply the
                // cell volume internally (cx = 0.5*(V*rAU) = dAU). Passing V/A_P
                // here double-counted the volume and shrank the Rhie-Chow
                // pressure correction of the authoritative coupled flux by a
                // factor V, silently removing the checkerboard suppression.
                // On the constant/linear-pressure channel cases the correction
                // vanishes and the defect was invisible (Couette/Poiseuille
                // passed); on the lid-driven cavity the assembled continuity
                // rows and the authoritative flux disagreed by div(interp(U)),
                // a permanent continuity floor (~8.5e-3 against the 1e-7 gate)
                // that no iteration budget or linear tolerance could remove:
                // the velocity froze at the Picard fixed point and the
                // stagnation detector terminated the case.
                coupled_rAU[0][c] = 1.0 / ex.diagonal[c];
                coupled_rAU[1][c] = 1.0 / ey.diagonal[c];
                coupled_rAU[2][c] = 1.0 / ez.diagonal[c];
            }
            mass_flux = make_rhie_chow_mass_flux(
                mesh, geometry, U, p, coupled_rAU,
                controls.density, velocity_bcs, pressure_bcs);

            // Operator identity audit: re-derive the authoritative flux from
            // the same state through the same production call and compare
            // against the stored field. Zero unless the stored flux was
            // mutated after the rebuild.
            {
                const auto rederived_flux = make_rhie_chow_mass_flux(
                    mesh, geometry, U, p, coupled_rAU,
                    controls.density, velocity_bcs, pressure_bcs);
                for (std::size_t f = 0; f < mesh.n_faces(); ++f)
                    rhie_chow_identity_audit_linf = std::max(
                        rhie_chow_identity_audit_linf,
                        std::abs(mass_flux(f) - rederived_flux(f)));
                rhie_chow_identity_evaluated = true;
            }

            // A coupled iteration has no segregated momentum/pressure Krylov
            // sub-solves. The block residual is the authoritative linear metric.
            // The independent nonlinear momentum residual below remains the
            // final physical acceptance metric.
        } else {
        relax_momentum_equation(ex, ux_old, retry_controller.alpha_u(effective_alpha_u));
        relax_momentum_equation(ey, uy_old, retry_controller.alpha_u(effective_alpha_u));
        relax_momentum_equation(ez, uz_old, retry_controller.alpha_u(effective_alpha_u));

        // Equation relaxation is already part of the matrix. The Krylov solve
        // therefore returns the solution of the relaxed equation directly;
        // applying alpha_u again here would double-relax the predictor.
        ScalarSolveControls momentum_solve_controls;
        momentum_solve_controls.max_iterations = controls.linear_max_iterations;
        momentum_solve_controls.tolerance = linear_tolerance;
        momentum_solve_controls.relaxation = 1.0;
        momentum_solve_controls.problem = LinearProblemKind::Momentum;
        momentum_solve_controls.linear_solver = controls.momentum_linear_solver;
        rx = solve_scalar_equation(ex, ux, momentum_solve_controls);
        ry = solve_scalar_equation(ey, uy, momentum_solve_controls);
        rz = solve_scalar_equation(ez, uz, momentum_solve_controls);

        auto require_linear_convergence = [](const char* component, const auto& solve) {
            if (solve.status != cfdx::core::SolverStatus::CONVERGED) {
                std::ostringstream message;
                message << "solve_steady_incompressible: " << component
                        << " momentum solve did not converge (status="
                        << static_cast<int>(solve.status)
                        << ", iterations=" << solve.iterations
                        << ", residual=" << std::scientific << std::setprecision(3)
                        << solve.residual
                        << ", relative=" << solve.residual_relative << ")";
                throw NonlinearRetryableFailure(message.str());
            }
        };
        require_linear_convergence("Ux", rx);
        require_linear_convergence("Uy", ry);
        require_linear_convergence("Uz", rz);

        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            U.component_data(0)[c] = ux(c);
            U.component_data(1)[c] = uy(c);
            U.component_data(2)[c] = uz(c);
        }

        std::array<std::vector<double>,3> rAU, rAtU, hbya;
        for (std::size_t d = 0; d < 3; ++d) {
            rAU[d].resize(mesh.n_cells());
            rAtU[d].resize(mesh.n_cells());
            for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
                const ScalarEquation* eqs[3] = {&ex, &ey, &ez};
                const double a = eqs[d]->diagonal[c];
                if (!(a > 0.0) || !std::isfinite(a))
                    throw std::runtime_error(
                        "solve_steady_incompressible: invalid momentum diagonal");
                rAU[d][c] = 1.0 / a;
                rAtU[d][c] = rAU[d][c];
                if (controls.algorithm == PressureVelocityAlgorithm::SIMPLEC) {
                    double signed_offdiag_sum = 0.0;
                    double abs_offdiag_sum = 0.0;
                    const auto begin = eqs[d]->matrix.row_offsets_data()[c];
                    const auto end = eqs[d]->matrix.row_offsets_data()[c + 1];
                    for (std::uint32_t k = begin; k < end; ++k) {
                        const auto col = eqs[d]->matrix.columns_data()[k];
                        // Boundary contributions are assembled into the owner
                        // diagonal/source; off-diagonal entries are therefore
                        // genuine internal neighbour coefficients.
                        if (col != c) {
                            const double value = eqs[d]->matrix.values_data()[k];
                            signed_offdiag_sum += value;
                            abs_offdiag_sum += std::abs(value);
                        }
                    }
                    const double diagonal = 1.0 / rAU[d][c];
                    const double denom =
                        simplec_consistent_diagonal(diagonal, signed_offdiag_sum);
                    // Independent audit counterpart. With the expected non-positive
                    // stored off-diagonals this equals the signed-contract
                    // denominator; a nonzero gap flags a matrix-sign mismatch.
                    const double abs_denom =
                        diagonal - abs_offdiag_sum;
                    simplec_signed_offdiag_sum_audit[d] =
                        std::max(simplec_signed_offdiag_sum_audit[d],
                                 std::abs(signed_offdiag_sum));
                    simplec_abs_offdiag_sum_audit[d] =
                        std::max(simplec_abs_offdiag_sum_audit[d], abs_offdiag_sum);
                    simplec_denominator_signed_audit[d] =
                        std::max(simplec_denominator_signed_audit[d], std::abs(denom));
                    simplec_denominator_abs_audit[d] =
                        std::max(simplec_denominator_abs_audit[d], std::abs(abs_denom));
                    const double gap = std::abs(denom - abs_denom);
                    if (gap > simplec_denominator_gap_audit[d])
                        simplec_denominator_worst_cell_audit[d] = c;
                    simplec_denominator_gap_audit[d] =
                        std::max(simplec_denominator_gap_audit[d], gap);
                    rAtU[d][c] = 1.0 / denom;
                }
            }
        }

        hbya[0] = build_hbya(ex, ux, grad_p, 0, rAU[0]);
        hbya[1] = build_hbya(ey, uy, grad_p, 1, rAU[1]);
        hbya[2] = build_hbya(ez, uz, grad_p, 2, rAU[2]);
        auto HbyA = make_hbya_field(hbya);

        if (freeze_state) {
            frozen_rAU = rAU;
            frozen_hbyA = hbya;
            frozen_grad_p = grad_p;
            frozen_state_valid = true;
        }

        // rAU is the inverse integrated momentum diagonal (1/A_P).
        // Pressure-gradient terms in the integrated momentum equation therefore
        // carry the cell volume explicitly when reconstructing cell velocity.
        // phiHbyA is the authoritative predictor flux. The Rhie-Chow term
        // uses the same directional inverse momentum coefficient as the
        // pressure equation, avoiding a second, incompatible face operator.
        auto phiHbyA = make_mass_flux(
            mesh, geometry, HbyA, controls.density, velocity_bcs);

        if (controls.algorithm == PressureVelocityAlgorithm::SIMPLEC) {
            for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                const auto nr = mesh.ownership().neighbour(f);
                const std::size_t o = mesh.ownership().owner(f);
                if (nr < 0) {
                    const std::size_t patch = geometry.face_patch[f];
                    if (patch >= mesh.boundary().n_patches() ||
                        mesh.boundary().patch(patch).type == PatchType::EMPTY)
                        continue;
                    const auto pbc = pressure_bcs.find(
                        mesh.boundary().patch(patch).name);
                    if (pbc == pressure_bcs.end() ||
                        pbc->second.type != ScalarBoundaryType::FIXED_VALUE)
                        continue;
                    const Vec3 Sf = geometry.face_area_vectors[f];
                    const Vec3 dvec = geometry.face_centres[f] -
                                      geometry.cell_centres[o];
                    const Vec3 delta_dAU{
                        geometry.cell_volumes[o] * (rAtU[0][o] - rAU[0][o]),
                        geometry.cell_volumes[o] * (rAtU[1][o] - rAU[1][o]),
                        geometry.cell_volumes[o] * (rAtU[2][o] - rAU[2][o])};
                    const Vec3 grad_p_cell{
                        grad_p.component_data(0)[o],
                        grad_p.component_data(1)[o],
                        grad_p.component_data(2)[o]};
                    const double correction_coeff = pressure_face_orthogonal_coefficient(
                        Sf, dvec, delta_dAU);
                    const double correction_deferred = pressure_face_nonorthogonal_flux(
                        Sf, dvec, delta_dAU, grad_p_cell);
                    phiHbyA(f) += controls.density * (
                        correction_coeff * (pbc->second.value - p(o)) +
                        correction_deferred);
                    continue;
                }
                const std::size_t n = static_cast<std::size_t>(nr);
                const auto Sf = geometry.face_area_vectors[f];
                const Vec3 dvec = geometry.cell_centres[n] - geometry.cell_centres[o];
                const double dx = 0.5*(geometry.cell_volumes[o]*rAtU[0][o] +
                                       geometry.cell_volumes[n]*rAtU[0][n]) -
                                  0.5*(geometry.cell_volumes[o]*rAU[0][o] +
                                       geometry.cell_volumes[n]*rAU[0][n]);
                const double dy = 0.5*(geometry.cell_volumes[o]*rAtU[1][o] +
                                       geometry.cell_volumes[n]*rAtU[1][n]) -
                                  0.5*(geometry.cell_volumes[o]*rAU[1][o] +
                                       geometry.cell_volumes[n]*rAU[1][n]);
                const double dz = 0.5*(geometry.cell_volumes[o]*rAtU[2][o] +
                                       geometry.cell_volumes[n]*rAtU[2][n]) -
                                  0.5*(geometry.cell_volumes[o]*rAU[2][o] +
                                       geometry.cell_volumes[n]*rAU[2][n]);
                const Vec3 delta_dAU{dx, dy, dz};
                const double correction_coeff = pressure_face_orthogonal_coefficient(
                    Sf, dvec, delta_dAU);
                const Vec3 grad_p_face{
                    0.5 * (grad_p.component_data(0)[o] + grad_p.component_data(0)[n]),
                    0.5 * (grad_p.component_data(1)[o] + grad_p.component_data(1)[n]),
                    0.5 * (grad_p.component_data(2)[o] + grad_p.component_data(2)[n])};
                // SIMPLEC changes the full face response, not only its
                // centre-line coefficient. Include the matching deferred
                // difference or skewed-face continuity retains an O(grad p)
                // mismatch between rAtU and rAU.
                const double correction_deferred = pressure_face_nonorthogonal_flux(
                    Sf, dvec, delta_dAU, grad_p_face);
                phiHbyA(f) += controls.density * (
                    correction_coeff * (p(n)-p(o)) + correction_deferred);
            }
        }

        for (std::size_t corr = 0; corr < pressure_correctors; ++corr) {
            const std::size_t nc = mesh.n_cells();
            SparseMatrix A(nc, nc);
            Vector b(nc, 0.0);
            std::vector<std::map<std::size_t, double>> rows(nc);
            std::vector<double> diag(nc, 0.0);
            std::vector<double> continuity(nc, 0.0);
            const auto current_grad_p =
                gauss_gradient_with_boundary(p, mesh, geometry, pressure_bcs);

            const auto* cell_faces = mesh.cells().faces_data();
            const auto* cell_offsets = mesh.cells().offsets_data();
            for (std::size_t c = 0; c < nc; ++c) {
                const Offset off = cell_offsets[c];
                const Offset count = cell_offsets[c + 1] - off;
                for (Offset k = 0; k < count; ++k) {
                    const std::size_t f = cell_faces[off+k];
                    continuity[c] += mesh.ownership().owner(f) == c
                        ? phiHbyA(f) : -phiHbyA(f);
                    if (mesh.ownership().neighbour(f) >= 0) {
                        const double phi_nonorth = rhie_chow_pressure_flux_internal(
                            mesh, geometry, f, p, current_grad_p,
                            controls.algorithm == PressureVelocityAlgorithm::SIMPLEC ? rAtU : rAU,
                            controls.density);
                        continuity[c] += mesh.ownership().owner(f) == c
                            ? -phi_nonorth : phi_nonorth;
                    }
                }
            }

            for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                const auto nr = mesh.ownership().neighbour(f);
                if (nr < 0) continue;
                const std::size_t o = mesh.ownership().owner(f);
                const std::size_t n = static_cast<std::size_t>(nr);
                const Vec3 Sf = geometry.face_area_vectors[f];
                const Vec3 dvec = geometry.cell_centres[n] - geometry.cell_centres[o];
                const auto& pressure_coeff =
                    controls.algorithm == PressureVelocityAlgorithm::SIMPLEC ? rAtU : rAU;
                const double coeff = controls.density * directional_face_coefficient(
                    pressure_coeff, o, n, Sf, dvec);
                diag[o] += coeff;
                diag[n] += coeff;
                rows[o][n] -= coeff;
                rows[n][o] -= coeff;
            }

            for (std::size_t c = 0; c < nc; ++c) {
                rows[c][c] += diag[c];
                b(c) = -continuity[c];
            }

            for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                if (mesh.ownership().neighbour(f) >= 0) continue;
                const std::size_t patch = geometry.face_patch[f];
                if (patch >= mesh.boundary().n_patches()) continue;
                const auto& name = mesh.boundary().patch(patch).name;
                const auto it = pressure_bcs.find(name);
                if (it == pressure_bcs.end() ||
                    it->second.type != ScalarBoundaryType::FIXED_VALUE) continue;
                const std::size_t o = mesh.ownership().owner(f);
                const Vec3 dvec = geometry.face_centres[f] - geometry.cell_centres[o];
                const Vec3 Sf = geometry.face_area_vectors[f];
                const auto& pressure_coeff =
                    controls.algorithm == PressureVelocityAlgorithm::SIMPLEC ? rAtU : rAU;
                const Vec3 dAU_boundary{
                    geometry.cell_volumes[o] * pressure_coeff[0][o],
                    geometry.cell_volumes[o] * pressure_coeff[1][o],
                    geometry.cell_volumes[o] * pressure_coeff[2][o]};
                const double boundary_coeff = pressure_face_orthogonal_coefficient(
                    Sf, dvec, dAU_boundary);
                rows[o][o] += controls.density * boundary_coeff;
                // Explicit pressure flux at fixed-pressure boundaries (#894):
                // the half-Dirichlet diagonal term above is only consistent
                // when the boundary mass flux is reconstructed as
                //   phi_b = phiHbyA_b - rho*boundary_coeff*(p_bc - p_P).
                // boundary_coeff uses the same centre-line projection as the
                // boundary flux reconstruction below. Without the matching
                // source term the discretised equation forces the boundary
                // flux to zero and continuity stalls on a nonzero floor.
                const Vec3 gp_boundary{
                    current_grad_p.component_data(0)[o],
                    current_grad_p.component_data(1)[o],
                    current_grad_p.component_data(2)[o]};
                const double deferred_boundary = pressure_face_nonorthogonal_flux(
                    Sf, dvec, dAU_boundary, gp_boundary);
                b(o) += controls.density * (
                    boundary_coeff * (it->second.value - p(o)) +
                    deferred_boundary);
            }

            Vector p_corr(nc, 0.0);
            cfdx::core::SolverResult rp;
            const auto solve_pressure_system =
                [&](const SparseMatrix& matrix,
                    const Vector& rhs,
                    Vector& correction) {
                    if (pressure_context) {
                        return pressure_context->solve(
                            matrix, rhs, correction,
                            controls.linear_max_iterations,
                            linear_tolerance, {}, {},
                            pressure_null_space.get());
                    }
                    return solve_linear_system(
                        matrix, rhs, correction,
                        LinearProblemKind::PressurePoisson,
                        controls.pressure_linear_solver,
                        controls.linear_max_iterations,
                        linear_tolerance).result;
                };
            if (has_fixed_pressure_boundary ||
                use_projected_pressure_null_space) {
                // With a prescribed pressure boundary the pressure-correction
                // operator is SPD. For pure Neumann pressure, an explicitly
                // requested constant null-space model instead keeps the
                // singular operator and projects every CG vector.
                for (std::size_t row = 0; row < nc; ++row)
                    for (const auto& [col, value] : rows[row])
                        A.push_back(row, col, value);
                A.finalize();
                rp = solve_pressure_system(A, b, p_corr);
            } else if (nc == 1) {
                // Eliminating the only pressure unknown would produce an empty
                // linear system. Keep the trivial gauge equation instead.
                A.push_back(0, 0, 1.0);
                A.finalize();
                b(0) = 0.0;
                rp = solve_pressure_system(A, b, p_corr);
            } else {
                // Pure-Neumann pressure has one gauge null mode. Replacing one
                // row by p_ref=0 breaks symmetry. Eliminate the reference DOF
                // instead: the remaining principal submatrix is SPD for a
                // connected cavity and remains eligible for the pressure
                // solver policy and its true-residual convergence checks.
                const std::size_t ref = controls.pressure_reference_cell;
                SparseMatrix reduced(nc - 1, nc - 1);
                Vector reduced_b(nc - 1, 0.0);
                const auto reduced_index = [ref](std::size_t cell) {
                    return cell < ref ? cell : cell - 1;
                };
                for (std::size_t row = 0; row < nc; ++row) {
                    if (row == ref) continue;
                    const std::size_t rr = reduced_index(row);
                    reduced_b(rr) = b(row);
                    for (const auto& [col, value] : rows[row]) {
                        if (col == ref) continue;
                        reduced.push_back(rr, reduced_index(col), value);
                    }
                }
                reduced.finalize();

                Vector reduced_corr(nc - 1, 0.0);
                rp = solve_pressure_system(
                    reduced, reduced_b, reduced_corr);
                if (rp.status == cfdx::core::SolverStatus::CONVERGED) {
                    for (std::size_t cell = 0; cell < nc; ++cell)
                        if (cell != ref)
                            p_corr(cell) = reduced_corr(reduced_index(cell));
                    p_corr(ref) = 0.0;
                }
            }
            pressure_residual = rp.residual_relative;
            pressure_iterations = rp.iterations;
            if (rp.status != cfdx::core::SolverStatus::CONVERGED) {
                std::ostringstream message;
                message << "solve_steady_incompressible: pressure-correction solve did not converge "
                        << "(status=" << static_cast<int>(rp.status)
                        << ", iterations=" << rp.iterations
                        << ", residual=" << std::scientific
                        << std::setprecision(std::numeric_limits<double>::max_digits10)
                        << rp.residual
                        << ", relative=" << rp.residual_relative << ")";
                throw NonlinearRetryableFailure(message.str());
            }

            // Relax the physical pressure exactly once. The reference is a
            // gauge: enforce it by a uniform shift, never by overwriting one
            // cell and creating an artificial local pressure jump.
            const double alpha_p = retry_controller.alpha_p(effective_alpha_p);
            for (std::size_t c = 0; c < nc; ++c)
                p(c) += alpha_p * p_corr(c);
            apply_pressure_gauge_policy(
                p, geometry, controls.pressure_gauge_policy,
                controls.pressure_reference_cell, controls.pressure_reference_value,
                has_fixed_pressure_boundary);

            auto corrected_grad_p =
                gauss_gradient_with_boundary(p, mesh, geometry, pressure_bcs);

            // Rebuild the corrected velocity from the same rAU operator
            // that built the cell HbyA, so the converged velocity
            // satisfies the momentum equation. The SIMPLEC consistent
            // operator (rAtU) belongs to the pressure equation and its
            // face fluxes only. This is the momentum correction; no
            // least-squares flux fitting is permitted.
            for (std::size_t c = 0; c < nc; ++c) {
                U.component_data(0)[c] =
                    HbyA.component_data(0)[c] - rAU[0][c] * corrected_grad_p.component_data(0)[c] * geometry.cell_volumes[c];
                U.component_data(1)[c] =
                    HbyA.component_data(1)[c] - rAU[1][c] * corrected_grad_p.component_data(1)[c] * geometry.cell_volumes[c];
                U.component_data(2)[c] =
                    HbyA.component_data(2)[c] - rAU[2][c] * corrected_grad_p.component_data(2)[c] * geometry.cell_volumes[c];
            }

            // The conservative face flux is rebuilt from the same
            // momentum-weighted interpolation used by the pressure equation.
            // This is the single authoritative phi state for continuity.
            mass_flux = make_rhie_chow_mass_flux(
                mesh, geometry, HbyA, p,
                controls.algorithm == PressureVelocityAlgorithm::SIMPLEC ? rAtU : rAU,
                controls.density, velocity_bcs, pressure_bcs);

            if (controls.algorithm == PressureVelocityAlgorithm::SIMPLEC) {
                // The SIMPLEC pressure equation balances the over-relaxed
                // predictor flux phiHbyA + rho*(dAtU - dAU)*(p_N - p_P),
                // not the bare interpolated predictor. The corrected flux
                // must therefore carry the same consistent-response term,
                // evaluated at the corrected pressure; otherwise it differs
                // from the operator the pressure equation just balanced by
                // exactly div(rho*(dAtU - dAU)*dp) and a continuity floor of
                // that size survives on every case with a non-constant
                // pressure field (e.g. the lid-driven cavity) even though
                // the momentum residual is zero. With the term included the
                // face operator collapses to the plain Rhie-Chow flux built
                // with rAU: (dAtU - dAU) cancels the rAtU response, leaving
                // the conservative flux of the reconstructed cell velocity.
                for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                    const auto nr = mesh.ownership().neighbour(f);
                    if (nr < 0) continue;
                    const std::size_t o = mesh.ownership().owner(f);
                    const std::size_t n = static_cast<std::size_t>(nr);
                    const auto Sf = geometry.face_area_vectors[f];
                    const Vec3 dvec = geometry.cell_centres[n] - geometry.cell_centres[o];
                    const double dx = 0.5*(geometry.cell_volumes[o]*rAtU[0][o] +
                                           geometry.cell_volumes[n]*rAtU[0][n]) -
                                      0.5*(geometry.cell_volumes[o]*rAU[0][o] +
                                           geometry.cell_volumes[n]*rAU[0][n]);
                    const double dy = 0.5*(geometry.cell_volumes[o]*rAtU[1][o] +
                                           geometry.cell_volumes[n]*rAtU[1][n]) -
                                      0.5*(geometry.cell_volumes[o]*rAU[1][o] +
                                           geometry.cell_volumes[n]*rAU[1][n]);
                    const double dz = 0.5*(geometry.cell_volumes[o]*rAtU[2][o] +
                                           geometry.cell_volumes[n]*rAtU[2][n]) -
                                      0.5*(geometry.cell_volumes[o]*rAU[2][o] +
                                           geometry.cell_volumes[n]*rAU[2][n]);
                    const double correction_coeff = pressure_face_orthogonal_coefficient(
                        Sf, dvec, {dx, dy, dz});
                    mass_flux(f) += controls.density * correction_coeff * (p(n)-p(o));
                }
            }

            // Issa corrector form: every pressure corrector is driven by the
            // ORIGINAL predictor flux together with the current pressure,
            // div(phiHbyA - RC(p)); the assembly at the top of this loop
            // already subtracts the Rhie-Chow pressure flux of the current
            // p. Overwriting phiHbyA with the corrected flux here would
            // subtract RC(p) twice in the next corrector's continuity source
            // (div(phiHbyA) - 2*div(RC(p))), which displaces the fixed point
            // and leaves a permanent continuity floor on every
            // multi-corrector algorithm (PISO, PIMPLE, FRACTIONAL_STEP).
            // The momentum matrix stays frozen inside the inner loop; the
            // relaxed pressure update is what carries the first correction
            // into the next one.
        }

        // Operator identity audit: re-derive the authoritative flux from the
        // CURRENT HbyA/p/rAtU(rAU) state through the same operator primitives
        // the pressure equation balanced (including the SIMPLEC
        // consistent-response term) and compare with the stored mass_flux.
        // Expected zero up to round-off; a nonzero value flags a stale or
        // mutated flux, or an omitted/duplicated SIMPLEC response term.
        {
            auto rederived_flux = make_rhie_chow_mass_flux(
                mesh, geometry, HbyA, p,
                controls.algorithm == PressureVelocityAlgorithm::SIMPLEC
                    ? rAtU : rAU,
                controls.density, velocity_bcs, pressure_bcs);
            if (controls.algorithm == PressureVelocityAlgorithm::SIMPLEC) {
                for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                    const auto nr = mesh.ownership().neighbour(f);
                    if (nr < 0) continue;
                    const std::size_t o = mesh.ownership().owner(f);
                    const std::size_t n = static_cast<std::size_t>(nr);
                    const auto Sf = geometry.face_area_vectors[f];
                    const Vec3 dvec =
                        geometry.cell_centres[n] - geometry.cell_centres[o];
                    const double dx = 0.5*(geometry.cell_volumes[o]*rAtU[0][o] +
                                          geometry.cell_volumes[n]*rAtU[0][n]) -
                                      0.5*(geometry.cell_volumes[o]*rAU[0][o] +
                                           geometry.cell_volumes[n]*rAU[0][n]);
                    const double dy = 0.5*(geometry.cell_volumes[o]*rAtU[1][o] +
                                          geometry.cell_volumes[n]*rAtU[1][n]) -
                                      0.5*(geometry.cell_volumes[o]*rAU[1][o] +
                                           geometry.cell_volumes[n]*rAU[1][n]);
                    const double dz = 0.5*(geometry.cell_volumes[o]*rAtU[2][o] +
                                          geometry.cell_volumes[n]*rAtU[2][n]) -
                                      0.5*(geometry.cell_volumes[o]*rAU[2][o] +
                                           geometry.cell_volumes[n]*rAU[2][n]);
                    const double correction_coeff =
                        pressure_face_orthogonal_coefficient(
                            Sf, dvec, {dx, dy, dz});
                    rederived_flux(f) +=
                        controls.density * correction_coeff * (p(n)-p(o));
                }
            }
            for (std::size_t f = 0; f < mesh.n_faces(); ++f)
                rhie_chow_identity_audit_linf = std::max(
                    rhie_chow_identity_audit_linf,
                    std::abs(mass_flux(f) - rederived_flux(f)));
            rhie_chow_identity_evaluated = true;
        }

        }
        
        // Independent final momentum residual: evaluate the actual corrected
        // state against freshly assembled equations. Do not use a U reconstructed
        // from the conservative flux as the residual state.
        Field<double, Location::CELL> final_ux_field(mesh.n_cells(), "Ux_final", "m/s", 1);
        Field<double, Location::CELL> final_uy_field(mesh.n_cells(), "Uy_final", "m/s", 1);
        Field<double, Location::CELL> final_uz_field(mesh.n_cells(), "Uz_final", "m/s", 1);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            final_ux_field(c) = U.component_data(0)[c];
            final_uy_field(c) = U.component_data(1)[c];
            final_uz_field(c) = U.component_data(2)[c];
        }
        auto final_grad_p = gauss_gradient_with_boundary(p, mesh, geometry, pressure_bcs);
        auto final_ex = assemble_momentum_component(
            mesh, geometry, mass_flux, final_grad_p, body_x, mu_eff, ubc_x, 0,
            controls.use_bounded_convection, controls.convection_scheme, &final_ux_field);
        auto final_ey = assemble_momentum_component(
            mesh, geometry, mass_flux, final_grad_p, body_y, mu_eff, ubc_y, 1,
            controls.use_bounded_convection, controls.convection_scheme, &final_uy_field);
        auto final_ez = assemble_momentum_component(
            mesh, geometry, mass_flux, final_grad_p, body_z, mu_eff, ubc_z, 2,
            controls.use_bounded_convection, controls.convection_scheme, &final_uz_field);

        // The final nonlinear acceptance residual must represent the same
        // production equation that was solved. For transient runs the
        // physical-time derivative is part of that equation; omitting it here
        // would measure only the steady spatial operator and reject every
        // legitimate transient state whose temporal acceleration balances the
        // spatial momentum residual. Reuse the exact production mass-term
        // assembly, including variable-step BDF2 coefficients and the accepted
        // history, so the acceptance metric is algebraically consistent with
        // the coupled solve above.
        if (controls.transient.enabled) {
            const auto& u_n = *controls.transient.previous;
            const auto& u_nm1 = controls.transient.history_valid
                ? *controls.transient.previous_previous : u_n;
            const bool use_bdf2 =
                controls.transient.scheme == cfdx::core::TimeScheme::BDF2 &&
                controls.transient.history_valid;
            add_dual_time_mass_term(
                final_ex, u_n, u_nm1, geometry, controls.density,
                controls.transient.dt, controls.transient.dt_previous,
                use_bdf2, 0);
            add_dual_time_mass_term(
                final_ey, u_n, u_nm1, geometry, controls.density,
                controls.transient.dt, controls.transient.dt_previous,
                use_bdf2, 1);
            add_dual_time_mass_term(
                final_ez, u_n, u_nm1, geometry, controls.density,
                controls.transient.dt, controls.transient.dt_previous,
                use_bdf2, 2);
        }

        Vector final_ux(mesh.n_cells()), final_uy(mesh.n_cells()), final_uz(mesh.n_cells());
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            final_ux(c) = U.component_data(0)[c];
            final_uy(c) = U.component_data(1)[c];
            final_uz(c) = U.component_data(2)[c];
        }
        const double final_momentum_residual = std::max({
            scalar_equation_residual_inf(final_ex, final_ux),
            scalar_equation_residual_inf(final_ey, final_uy),
            scalar_equation_residual_inf(final_ez, final_uz)
        });

        // Recompute the cell-wise residuals explicitly so a large independent
        // residual can be localized to a component and to an interior/boundary
        // cell. This is intentionally redundant with scalar_equation_residual_inf:
        // the purpose is forensic diagnostics, not a second numerical path.
        const auto cell_has_boundary_face = [&](std::size_t cell) {
            const Offset off = mesh.cells().offsets_data()[cell];
            const Offset count = mesh.cells().offsets_data()[cell + 1] - off;
            for (Offset k = 0; k < count; ++k) {
                const std::size_t face = mesh.cells().faces_data()[off + k];
                if (mesh.ownership().neighbour(face) < 0) {
                    const std::size_t patch = geometry.face_patch[face];
                    if (patch < mesh.boundary().n_patches() &&
                        mesh.boundary().patch(patch).type == PatchType::EMPTY)
                        continue;
                    return true;
                }
            }
            return false;
        };
        struct ResidualDiagnostic {
            double global = 0.0;
            double interior = 0.0;
            double boundary = 0.0;
            std::size_t max_cell = 0;
            double no_pressure = 0.0;
            double pressure_contribution = 0.0;
            std::string patch;
        };
        const auto residual_diagnostics =
            [&](const ScalarEquation& equation,
                const Vector& solution,
                std::size_t component) {
                ResidualDiagnostic diagnostic;
                for (std::size_t row = 0; row < mesh.n_cells(); ++row) {
                    double ax = 0.0;
                    const auto begin = equation.matrix.row_offsets_data()[row];
                    const auto end = equation.matrix.row_offsets_data()[row + 1];
                    for (std::uint32_t k = begin; k < end; ++k)
                        ax += equation.matrix.values_data()[k] *
                              solution(equation.matrix.columns_data()[k]);
                    const double full_signed = ax - equation.rhs(row);
                    const double residual = std::abs(full_signed);
                    if (residual > diagnostic.global) {
                        diagnostic.global = residual;
                        diagnostic.max_cell = row;
                    }
                    if (cell_has_boundary_face(row))
                        diagnostic.boundary = std::max(diagnostic.boundary, residual);
                    else
                        diagnostic.interior = std::max(diagnostic.interior, residual);

                    // final_ex is assembled with source = body - grad(p).
                    // Remove the pressure source only for diagnosis, never for
                    // the acceptance metric, to distinguish pressure/transport
                    // inconsistency from boundary transport inconsistency.
                    const double pressure = final_grad_p.component_data(component)[row] *
                                            geometry.cell_volumes[row];
                    const double transport_signed = full_signed - pressure;
                    diagnostic.no_pressure =
                        std::max(diagnostic.no_pressure, std::abs(transport_signed));
                    diagnostic.pressure_contribution =
                        std::max(diagnostic.pressure_contribution, std::abs(pressure));
                }
                const std::size_t row = diagnostic.max_cell;
                const Offset off = mesh.cells().offsets_data()[row];
                const Offset count = mesh.cells().offsets_data()[row + 1] - off;
                for (Offset k = 0; k < count; ++k) {
                    const std::size_t face = mesh.cells().faces_data()[off + k];
                    if (mesh.ownership().neighbour(face) < 0) {
                        const std::size_t patch = geometry.face_patch[face];
                        if (patch >= mesh.boundary().n_patches())
                            continue;
                        // EMPTY faces are geometric extrusion faces, not
                        // physical boundary patches. Do not let the first
                        // front/back face mask the physical patch adjacent to
                        // the worst momentum-residual cell.
                        if (mesh.boundary().patch(patch).type == PatchType::EMPTY)
                            continue;
                        diagnostic.patch = mesh.boundary().patch(patch).name;
                        break;
                    }
                }
                return diagnostic;
            };
        const auto rx_diag = residual_diagnostics(final_ex, final_ux, 0);
        const auto ry_diag = residual_diagnostics(final_ey, final_uy, 1);
        const auto rz_diag = residual_diagnostics(final_ez, final_uz, 2);

        double pressure_gradient_linf = 0.0;
        double pressure_gradient_l2_sum = 0.0;
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            const double gx = final_grad_p.component_data(0)[c];
            const double gy = final_grad_p.component_data(1)[c];
            const double gz = final_grad_p.component_data(2)[c];
            const double gp2 = gx*gx + gy*gy + gz*gz;
            pressure_gradient_linf = std::max(pressure_gradient_linf, std::sqrt(gp2));
            pressure_gradient_l2_sum += gp2;
        }
        const double pressure_gradient_l2 =
            std::sqrt(pressure_gradient_l2_sum /
                      std::max<std::size_t>(mesh.n_cells(), 1));

        double l1 = 0.0, linf = 0.0;
        const auto* cell_faces = mesh.cells().faces_data();
        const auto* cell_offsets = mesh.cells().offsets_data();
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            double div = 0.0;
            const Offset off = cell_offsets[c];
            const Offset count = cell_offsets[c+1] - off;
            for (Offset k = 0; k < count; ++k) {
                const std::size_t f = cell_faces[off+k];
                div += mesh.ownership().owner(f) == c ? mass_flux(f) : -mass_flux(f);
            }
            l1 += std::abs(div);
            linf = std::max(linf, std::abs(div));
        }

        constexpr double tiny = 1e-300;
        double velocity_change_inf = 0.0, pressure_change_inf = 0.0;
        double velocity_scale = 0.0;
        double pressure_min = std::numeric_limits<double>::infinity();
        double pressure_max = -std::numeric_limits<double>::infinity();
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            for (std::size_t d = 0; d < 3; ++d) {
                velocity_change_inf = std::max(
                    velocity_change_inf,
                    std::abs(U.component_data(d)[c]-U_old.component_data(d)[c]));
                velocity_scale = std::max(
                    velocity_scale, std::abs(U.component_data(d)[c]));
            }
            const double pressure_delta = has_fixed_pressure_boundary
                ? std::abs(p(c) - p_old(c))
                : std::abs((p(c) - p(controls.pressure_reference_cell)) -
                            (p_old(c) - p_old(controls.pressure_reference_cell)));
            pressure_change_inf = std::max(pressure_change_inf, pressure_delta);
            pressure_min = std::min(pressure_min, p(c));
            pressure_max = std::max(pressure_max, p(c));
        }
        velocity_change_inf /= std::max(velocity_scale, tiny);
        const double pressure_range = pressure_max - pressure_min;
        // Use a physical pressure scale when the converged solution is nearly uniform.
        // Otherwise round-off-level pressure changes are divided by a vanishing range
        // and can prevent convergence forever on constant-pressure cases.
        const double pressure_scale = std::max(
            {pressure_range, controls.density * velocity_scale * velocity_scale, 1.0e-30});
        pressure_change_inf /= pressure_scale;

        const std::array<const ScalarEquation*, 3> momentum_equations{
            &final_ex, &final_ey, &final_ez};
        std::array<double, 3> momentum_rhs_scale_components{0.0, 0.0, 0.0};
        for (std::size_t d = 0; d < momentum_equations.size(); ++d) {
            for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
                momentum_rhs_scale_components[d] = std::max(
                    momentum_rhs_scale_components[d],
                    std::abs(momentum_equations[d]->rhs(c)));
            }
            momentum_rhs_scale_components[d] =
                std::max(momentum_rhs_scale_components[d], tiny);
        }
        const double momentum_rhs_scale = std::max({
            momentum_rhs_scale_components[0],
            momentum_rhs_scale_components[1],
            momentum_rhs_scale_components[2]});

        IncompressibleIteration h;
        h.iteration = iter;
        h.momentum_residual = final_momentum_residual / momentum_rhs_scale;
        h.pressure_residual = pressure_residual;
        h.continuity_l1 = l1;
        h.continuity_linf = linf;
        h.momentum_equation_residual = final_momentum_residual;
        h.momentum_rhs_scale = momentum_rhs_scale;
        h.momentum_rhs_scale_components = momentum_rhs_scale_components;
        h.momentum_equation_residual_rhs_relative =
            final_momentum_residual / momentum_rhs_scale;
        h.momentum_equation_residual_rhs_relative_components = {
            rx_diag.global / momentum_rhs_scale_components[0],
            ry_diag.global / momentum_rhs_scale_components[1],
            rz_diag.global / momentum_rhs_scale_components[2]};
        double flux_reference = 0.0;
        for (std::size_t f = 0; f < mesh.n_faces(); ++f)
            flux_reference += std::abs(mass_flux(f));
        h.continuity_normalized = l1 / std::max(flux_reference, tiny);
        if (controls.algorithm == PressureVelocityAlgorithm::COUPLED)
            h.pressure_residual = h.continuity_normalized;
        if (iter == 1)
            result.reference_momentum_residual =
                std::max(final_momentum_residual, tiny);
        h.momentum_equation_residual_reference =
            std::max(result.reference_momentum_residual, tiny);
        h.momentum_equation_residual_relative =
            final_momentum_residual / h.momentum_equation_residual_reference;
        h.velocity_change_inf = velocity_change_inf;
        h.pressure_change_inf = pressure_change_inf;
        // Report the relaxation factors actually used on this attempt, including
        // any retry-controller scaling. These are the values needed to replay
        // an accepted iteration deterministically from its history record.
        h.effective_alpha_u = retry_controller.alpha_u(effective_alpha_u);
        h.effective_alpha_p = retry_controller.alpha_p(effective_alpha_p);
        h.simplec_signed_offdiag_sum_linf = simplec_signed_offdiag_sum_audit;
        h.simplec_abs_offdiag_sum_linf = simplec_abs_offdiag_sum_audit;
        h.simplec_denominator_signed_linf = simplec_denominator_signed_audit;
        h.simplec_denominator_abs_linf = simplec_denominator_abs_audit;
        h.simplec_denominator_gap_linf = simplec_denominator_gap_audit;
        h.simplec_denominator_worst_cell = simplec_denominator_worst_cell_audit;
        h.momentum_linear_iterations =
            controls.algorithm == PressureVelocityAlgorithm::COUPLED
                ? 0 : std::max({rx.iterations, ry.iterations, rz.iterations});
        h.pressure_linear_iterations =
            controls.algorithm == PressureVelocityAlgorithm::COUPLED ? 0 : pressure_iterations;
        h.pressure_correctors_used =
            controls.algorithm == PressureVelocityAlgorithm::COUPLED ? 0 : pressure_correctors;
        h.linear_tolerance_used = controls.algorithm == PressureVelocityAlgorithm::COUPLED
            ? controls.coupling.coupled_linear_tolerance : linear_tolerance;
        h.corrected_flux_continuity_linf = linf;
        // Independent flux/velocity audit (restores the measurement channel
        // that commit 7c7dadf5 made tautological with
        // reconstructed_flux = mass_flux). Two reconstruction paths are
        // compared:
        //   channel RC (authoritative): mass_flux, the Rhie-Chow flux built
        //     from HbyA/p plus, for SIMPLEC, the consistent-response term;
        //   channel U (independent): the flux rebuilt from the corrected
        //     cell velocity with the plain linear interpolation of
        //     make_mass_flux, which shares no operator with Rhie-Chow.
        // The mismatch is split internal / fixed-pressure boundary faces:
        // on fixed-pressure boundary faces the two channels differ
        // analytically by the Rhie-Chow boundary correction plus the
        // predictor/corrector velocity response, so a persistent value there
        // is an expected artifact of the boundary contract, not a
        // conservation error. A persistent internal-face mismatch at
        // convergence is an evidence channel about the interpolation
        // identity. Neither channel is an acceptance metric: no criterion
        // consumes them (see the N9/S9 test suite for the sensitivity
        // proofs).
        auto reconstructed_flux = make_mass_flux(
            mesh, geometry, U, controls.density, velocity_bcs);
        double reconstructed_linf = 0.0;
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            double div_reconstructed = 0.0;
            const Offset off = cell_offsets[c];
            const Offset count = cell_offsets[c+1] - off;
            for (Offset k = 0; k < count; ++k) {
                const std::size_t f = cell_faces[off+k];
                div_reconstructed +=
                    mesh.ownership().owner(f) == c ? reconstructed_flux(f) : -reconstructed_flux(f);
            }
            reconstructed_linf = std::max(reconstructed_linf, std::abs(div_reconstructed));
        }
        double flux_mismatch_linf = 0.0;
        double flux_mismatch_boundary_linf = 0.0;
        double flux_mismatch_scale = 0.0;
        double max_velocity_abs = 0.0;
        double max_face_area = 0.0;
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            for (std::size_t d = 0; d < 3; ++d)
                max_velocity_abs = std::max(
                    max_velocity_abs, std::abs(U.component_data(d)[c]));
        for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
            const double mismatch =
                std::abs(mass_flux(f) - reconstructed_flux(f));
            flux_mismatch_scale = std::max(flux_mismatch_scale, std::abs(mass_flux(f)));
            max_face_area = std::max(
                max_face_area, geometry.face_area_vectors[f].mag());
            if (mesh.ownership().neighbour(f) >= 0) {
                flux_mismatch_linf = std::max(flux_mismatch_linf, mismatch);
                continue;
            }
            const std::size_t patch = geometry.face_patch[f];
            if (patch >= mesh.boundary().n_patches())
                continue;
            if (mesh.boundary().patch(patch).type == PatchType::EMPTY)
                continue;
            const auto& patch_name = mesh.boundary().patch(patch).name;
            const auto pbc = pressure_bcs.find(patch_name);
            if (pbc != pressure_bcs.end() &&
                pbc->second.type == ScalarBoundaryType::FIXED_VALUE)
                flux_mismatch_boundary_linf =
                    std::max(flux_mismatch_boundary_linf, mismatch);
        }
        // Normalized internal mismatch with an explicit physical floor:
        // E_phi = ||phi_RC - phi_U||_inf / max(||phi_RC||_inf,
        // rho * max|U| * max|Sf|). Never divides by a vanishing scale.
        const double flux_mismatch_denominator = std::max(
            flux_mismatch_scale,
            controls.density * max_velocity_abs * max_face_area);
        const double flux_mismatch_normalized =
            flux_mismatch_denominator > 0.0
                ? flux_mismatch_linf / flux_mismatch_denominator
                : 0.0;
        h.reconstructed_velocity_continuity_linf = reconstructed_linf;
        h.flux_velocity_mismatch_linf = flux_mismatch_linf;
        h.flux_velocity_mismatch_boundary_linf = flux_mismatch_boundary_linf;
        h.flux_velocity_mismatch_normalized = flux_mismatch_normalized;
        h.rhie_chow_identity_linf =
            rhie_chow_identity_evaluated ? rhie_chow_identity_audit_linf : 0.0;
        h.momentum_equation_residual_components =
            std::array<double, 3>{rx_diag.global, ry_diag.global, rz_diag.global};
        h.momentum_equation_residual_internal =
            std::max({rx_diag.interior, ry_diag.interior, rz_diag.interior});
        h.momentum_equation_residual_boundary =
            std::max({rx_diag.boundary, ry_diag.boundary, rz_diag.boundary});
        h.pressure_gradient_linf = pressure_gradient_linf;
        h.pressure_gradient_l2 = pressure_gradient_l2;
        const ResidualDiagnostic* worst_diag = &rx_diag;
        if (ry_diag.global > worst_diag->global) worst_diag = &ry_diag;
        if (rz_diag.global > worst_diag->global) worst_diag = &rz_diag;
        h.momentum_residual_cell = worst_diag->max_cell;
        h.momentum_residual_no_pressure = worst_diag->no_pressure;
        h.momentum_pressure_contribution = worst_diag->pressure_contribution;
        h.momentum_residual_patch = worst_diag->patch;

        // Optional forensic microscope for one cell. This is deliberately
        // disabled unless CFDX_DEBUG_CELL is set, so production/CI output is
        // unchanged. Use CFDX_DEBUG_CELL=<index> for a fixed cell or
        // CFDX_DEBUG_CELL=auto to inspect the independently reassembled
        // worst-residual cell. The dump is emitted once per algorithm at
        // convergence or at the configured iteration limit.
        const bool debug_cell_auto = controls.diagnostics.debug_cell_auto;
        std::size_t debug_cell = 0;
        const bool debug_cell_enabled =
            debug_cell_auto ||
            (controls.diagnostics.debug_cell >= 0 &&
             static_cast<std::size_t>(controls.diagnostics.debug_cell) < mesh.n_cells());
        if (debug_cell_auto)
            debug_cell = h.momentum_residual_cell;
        else if (debug_cell_enabled)
            debug_cell = static_cast<std::size_t>(controls.diagnostics.debug_cell);

        // The first-iteration-relative equation residual is diagnostic only:
        // if the initial residual is already at round-off, its denominator is
        // ill-conditioned and the ratio can grow while the actual equation
        // residual remains negligible. Use the RHS-normalized residual for the
        // nonlinear convergence monitor and keep the reference-relative value
        // available in history for diagnostics/adaptive-relaxation analysis.
        // The first-iteration normalization is useful only when its reference
        // is resolved relative to the momentum RHS. If the reference is at
        // round-off, keep it for diagnostics but do not let its ratio dominate
        // the nonlinear monitor.
        const double momentum_reference_roundoff =
            100.0 * std::numeric_limits<double>::epsilon() *
            std::max(h.momentum_rhs_scale, tiny);
        const bool momentum_reference_ill_conditioned =
            h.momentum_equation_residual_reference <= momentum_reference_roundoff;
        double nonlinear_metric = momentum_reference_ill_conditioned
            ? std::max({
                  h.momentum_residual,
                  h.pressure_residual,
                  h.velocity_change_inf,
                  h.pressure_change_inf})
            : std::max({
                  h.momentum_residual,
                  h.momentum_equation_residual_relative,
                  h.pressure_residual,
                  h.velocity_change_inf,
                  h.pressure_change_inf});
        h.nonlinear_convergence_metric = nonlinear_metric;

        if (controls.diagnostics.iteration_trace &&
            (iter == 1 ||
             iter % std::max<std::size_t>(1, controls.diagnostics.iteration_trace_frequency) == 0)) {
            // Diagnostic-only trace: expose the authoritative nonlinear history
            // before the convergence gate. No criterion or tolerance is changed.
            std::cerr << "INCOMPRESSIBLE_ITER iter=" << iter
                      << " momentum_residual_rhs_relative=" << h.momentum_equation_residual_rhs_relative
                      << " momentum_residual_reference_relative=" << h.momentum_equation_residual_relative
                      << " momentum_residual_raw=" << h.momentum_equation_residual
                      << " momentum_rhs_scale=" << h.momentum_rhs_scale
                      << " momentum_residual_reference=" << h.momentum_equation_residual_reference
                      << " momentum_rhs_scale_components=("
                      << h.momentum_rhs_scale_components[0] << ","
                      << h.momentum_rhs_scale_components[1] << ","
                      << h.momentum_rhs_scale_components[2] << ")"
                      << " momentum_rhs_relative_components=("
                      << h.momentum_equation_residual_rhs_relative_components[0] << ","
                      << h.momentum_equation_residual_rhs_relative_components[1] << ","
                      << h.momentum_equation_residual_rhs_relative_components[2] << ")"
                      << " continuity_normalized=" << h.continuity_normalized
                      << " velocity_change_inf=" << h.velocity_change_inf
                      << " pressure_change_inf=" << h.pressure_change_inf
                      << " alpha_u=" << h.effective_alpha_u
                      << " alpha_p=" << h.effective_alpha_p
                      << " pressure_linear_iterations=" << pressure_iterations
                      << " pressure_relative_residual=" << pressure_residual
                      << " momentum_linear_iterations="
                      << ((controls.algorithm == PressureVelocityAlgorithm::COUPLED)
                              ? 0 : std::max({rx.iterations, ry.iterations, rz.iterations}))
                      << " flux_velocity_mismatch_linf=" << h.flux_velocity_mismatch_linf
                      << " worst_momentum_cell=" << h.momentum_residual_cell
                      << " worst_momentum_patch=" << h.momentum_residual_patch
                      << "\n";
        }
        // Acceptance stays per-criterion. The monitor divides its synthetic
        // maximum by the maximum seen at the first iteration, so a metric whose
        // absolute scale sits far below that first maximum would be exempted
        // from the gate. nonlinear_metric is therefore a diagnostic summary and
        // never a replacement for the individual physical criteria.
        const bool momentum_residual_gate =
            h.momentum_equation_residual_rhs_relative <= controls.convergence.relative_tolerance;
        const bool momentum_equation_residual_gate =
            momentum_reference_ill_conditioned ||
            h.momentum_equation_residual_relative <= controls.convergence.relative_tolerance;
        const bool pressure_residual_gate =
            h.pressure_residual <= controls.convergence.relative_tolerance;
        const bool velocity_change_gate =
            h.velocity_change_inf <= controls.convergence.relative_tolerance;
        const bool pressure_change_gate =
            h.pressure_change_inf <= controls.convergence.relative_tolerance;
        const bool per_criterion_gates_ok =
            momentum_residual_gate && momentum_equation_residual_gate &&
            pressure_residual_gate && velocity_change_gate && pressure_change_gate;

        // Reject severe adaptive overshoots before they reach the convergence
        // monitor. A retry is the same nonlinear iteration, so convergence
        // history must only observe the accepted state.
        if (controls.adaptive_relaxation.enabled && !result.history.empty()) {
            const auto& previous = result.history.back();
            const bool severe_u =
                std::isfinite(previous.momentum_equation_residual_relative) &&
                std::isfinite(h.momentum_equation_residual_relative) &&
                previous.momentum_equation_residual_relative > 0.0 &&
                adaptive_relaxation_severe_degradation(
                    previous.momentum_equation_residual_relative,
                    h.momentum_equation_residual_relative,
                    previous.effective_alpha_u,
                    h.effective_alpha_u,
                    controls.adaptive_relaxation);
            const bool severe_p =
                std::isfinite(previous.continuity_normalized) &&
                std::isfinite(h.continuity_normalized) &&
                previous.continuity_normalized > 0.0 &&
                adaptive_relaxation_severe_degradation(
                    previous.continuity_normalized,
                    h.continuity_normalized,
                    previous.effective_alpha_p,
                    h.effective_alpha_p,
                    controls.adaptive_relaxation);

            if (severe_u || severe_p) {
                if (!retry_controller.can_retry())
                    throw NonlinearRetryableFailure(
                        "solve_steady_incompressible: adaptive relaxation severe overshoot "
                        "retry limit exhausted");
                if (severe_u) {
                    effective_alpha_u = std::max(
                        controls.adaptive_relaxation.min_alpha_u,
                        effective_alpha_u - controls.adaptive_relaxation.decrease_step);
                    adaptive_u_improvement_streak = 0;
                    adaptive_u_degradation_streak = 0;
                    adaptive_u_recovery_cooldown =
                        controls.adaptive_relaxation.recovery_cooldown_windows;
                }
                if (severe_p) {
                    effective_alpha_p = std::max(
                        controls.adaptive_relaxation.min_alpha_p,
                        effective_alpha_p - controls.adaptive_relaxation.decrease_step);
                    adaptive_p_improvement_streak = 0;
                    adaptive_p_degradation_streak = 0;
                    adaptive_p_recovery_cooldown =
                        controls.adaptive_relaxation.recovery_cooldown_windows;
                }
                transaction.reject();
                mass_flux = mass_flux_snapshot;
                if (controls.sst_runtime && sst_attempt_snapshot)
                    restore_sst_production_state(*controls.sst_runtime->state, *sst_attempt_snapshot);
                frozen_state_valid = false;
                // The severe guard has already reduced the responsible effective
                // relaxation factor. Count the retry without applying the generic
                // retry-controller shrink a second time.
                retry_controller.retry_without_relaxation();
                continue;
            }
        }
        std::map<std::string, double> qoi_values;
        if (!controls.qoi_gates.empty()) {
            for (const auto& gate : controls.qoi_gates)
                qoi_values[gate.name] =
                    sample_incompressible_probe(
                        *std::find_if(
                            controls.probes.begin(), controls.probes.end(),
                            [&](const IncompressiblePointProbe& probe) {
                                return probe.name == gate.name;
                            }),
                        mesh, geometry, U, p);
        }
        SSTProductionStepResult sst_step;
        double sst_nu_t_relative_change = 0.0;
        if (controls.sst_runtime) {
            std::map<std::string, std::array<double, 3>> fixed_velocity_values;
            for (const auto& [patch, bc] : velocity_bcs) {
                if (bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                    fixed_velocity_values.emplace(
                        patch, std::array<double, 3>{bc.value.x, bc.value.y, bc.value.z});
            }
            sst_step = advance_sst_production_iterate_transactionally(
                mesh, U, p, mass_flux, *controls.sst_runtime->state,
                controls.sst_runtime->transport, fixed_velocity_values,
                controls.sst_runtime->wall_patch_names,
                controls.sst_runtime->k_boundary_conditions,
                controls.sst_runtime->omega_boundary_conditions,
                controls.sst_runtime->max_transport_iterations,
                controls.sst_runtime->transport_tolerance,
                &controls.sst_runtime->omega_face_values);
            if (!sst_step.committed)
                throw NonlinearRetryableFailure(
                    "solve_steady_incompressible: SST update rejected: " + sst_step.reason);
            // Turbulence convergence is an additional gate; it never relaxes
            // the existing carrier residual, continuity, or QoI criteria.
            if (!std::isfinite(sst_step.transport.k_residual) ||
                !std::isfinite(sst_step.transport.second_residual) ||
                sst_step.transport.nonfinite_bound_events != 0)
                throw NonlinearRetryableFailure(
                    "solve_steady_incompressible: SST residual/boundedness gate failed");
            // The transport solver's converged flag is based on its linear
            // solver status, bound audit, and field-update criterion. Its raw
            // equation residuals have dimensional scales and must not be mixed
            // directly with the normalized carrier convergence metric.
            if (!sst_attempt_snapshot)
                throw std::logic_error("SST production attempt is missing its state snapshot");
            sst_nu_t_relative_change = sst_relative_eddy_viscosity_change(
                *sst_attempt_snapshot, *controls.sst_runtime->state,
                controls.kinematic_viscosity);
        }
        const auto convergence_report = convergence_monitor.update({
            iter,
            nonlinear_metric,
            h.continuity_normalized,
            std::move(qoi_values)});
        const bool sst_residual_gate = !controls.sst_runtime ||
            (sst_step.committed && sst_step.transport.converged &&
             std::isfinite(sst_step.transport.k_residual) &&
             std::isfinite(sst_step.transport.second_residual) &&
             sst_step.transport.nonfinite_bound_events == 0 &&
             std::isfinite(sst_nu_t_relative_change) &&
             sst_nu_t_relative_change <= controls.convergence.relative_tolerance);
        const bool converged_now =
            iter >= minimum_outer_correctors &&
            convergence_report.status == cfdx::core::ConvergenceStatus::CONVERGED &&
            per_criterion_gates_ok && sst_residual_gate;
        if (debug_cell_enabled && (converged_now || iter == controls.convergence.max_iterations)) {
            std::cerr << "\n=== CFDX MOMENTUM MICROSCOPE cell=" << debug_cell
                      << " cell_source=" << (debug_cell_auto ? "auto_worst" : "explicit")
                      << " algorithm=" << static_cast<int>(controls.algorithm)
                      << " iteration=" << iter << " ===\n";
            const auto dump_equation = [&](const char* name,
                                           const ScalarEquation& eq,
                                           const Vector& solution) {
                double ax = 0.0;
                const auto begin = eq.matrix.row_offsets_data()[debug_cell];
                const auto end = eq.matrix.row_offsets_data()[debug_cell + 1];
                std::cerr << name << " matrix row: rhs=" << eq.rhs(debug_cell)
                          << " diagonal=" << eq.diagonal[debug_cell]
                          << " nnz=" << (end - begin) << "\n";
                for (std::uint32_t k = begin; k < end; ++k) {
                    const auto col = eq.matrix.columns_data()[k];
                    const auto value = eq.matrix.values_data()[k];
                    ax += value * solution(col);
                    std::cerr << "  col=" << col << " a=" << value
                              << " x=" << solution(col)
                              << " ax=" << value * solution(col) << "\n";
                }
                const double residual = ax - eq.rhs(debug_cell);
                std::cerr << "  exact_row_ax=" << ax
                          << " exact_row_residual=" << residual
                          << " abs=" << std::abs(residual) << "\n";
            };
            dump_equation("Ux", final_ex, final_ux);
            dump_equation("Uy", final_ey, final_uy);
            dump_equation("Uz", final_ez, final_uz);

            // Reconstruct the same momentum predictor used by the pressure-
            // velocity coupling from the final, independently assembled
            // momentum rows. This is diagnostic only; it does not alter the
            // converged state. dAU = V/A_P is the exact pressure-response
            // coefficient used by the segregated Rhie-Chow operator.
            std::array<std::vector<double>, 3> diag_rAU;
            std::array<std::vector<double>, 3> hbyA;
            const ScalarEquation* final_eqs[3] = {&final_ex, &final_ey, &final_ez};
            const Vector* final_u[3] = {&final_ux, &final_uy, &final_uz};
            for (std::size_t d = 0; d < 3; ++d) {
                diag_rAU[d].resize(mesh.n_cells());
                hbyA[d].resize(mesh.n_cells());
                for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
                    const double ap = final_eqs[d]->diagonal[c];
                    if (!(ap > 0.0) || !std::isfinite(ap))
                        throw std::runtime_error(
                            "CFDX momentum microscope: invalid final momentum diagonal");
                    diag_rAU[d][c] = 1.0 / ap;
                    double h_explicit = final_eqs[d]->rhs(c) +
                        final_grad_p.component_data(d)[c] * geometry.cell_volumes[c];
                    const auto begin = final_eqs[d]->matrix.row_offsets_data()[c];
                    const auto end = final_eqs[d]->matrix.row_offsets_data()[c + 1];
                    for (std::uint32_t k = begin; k < end; ++k) {
                        const auto col = final_eqs[d]->matrix.columns_data()[k];
                        if (col != c)
                            h_explicit -=
                                final_eqs[d]->matrix.values_data()[k] * (*final_u[d])(col);
                    }
                    hbyA[d][c] = diag_rAU[d][c] * h_explicit;
                }
            }

            const auto reconstructed_flux =
                make_mass_flux(mesh, geometry, U, controls.density, velocity_bcs);
            double local_div_authoritative = 0.0;
            double local_div_reconstructed = 0.0;
            double max_delta_phi = 0.0;
            double max_delta_grad_p = 0.0;
            double max_delta_u_face = 0.0;
            double max_momentum_pressure_flux = 0.0;
            double max_correction_pressure_flux = 0.0;
            double max_rhie_chow_actual = 0.0;
            double max_rhie_chow_expected = 0.0;
            double max_rhie_chow_mismatch = 0.0;
            if (freeze_state && frozen_state_valid) {
                double phi_delta_linf = 0.0;
                double hbyA_delta_linf = 0.0;
                std::size_t hbyA_delta_cell = 0;
                for (std::size_t f = 0; f < mesh.n_faces(); ++f)
                    phi_delta_linf = std::max(phi_delta_linf,
                        std::abs(mass_flux(f) - phi_used_in_momentum(f)));
                for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
                    double delta = 0.0;
                    for (std::size_t d = 0; d < 3; ++d)
                        delta = std::max(delta, std::abs(hbyA[d][c] - frozen_hbyA[d][c]));
                    if (delta > hbyA_delta_linf) {
                        hbyA_delta_linf = delta;
                        hbyA_delta_cell = c;
                    }
                }
                const std::size_t comparison_cell =
                    debug_cell == 33 && mesh.n_cells() > 70 ? 70 : 33;
                auto dump_state_cell = [&](std::size_t cell, const char* label) {
                    if (cell >= mesh.n_cells()) return;
                    double phi_delta_l1 = 0.0;
                    double phi_used_l1 = 0.0;
                    double phi_reassembled_l1 = 0.0;
                    double hbyA_delta = 0.0;
                    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                        const auto owner = mesh.ownership().owner(f);
                        const auto neighbour = mesh.ownership().neighbour(f);
                        if (owner != cell && (neighbour < 0 || static_cast<std::size_t>(neighbour) != cell)) continue;
                        const double sign = owner == cell ? 1.0 : -1.0;
                        const double used = sign * phi_used_in_momentum(f);
                        const double current = sign * mass_flux(f);
                        phi_used_l1 += std::abs(used);
                        phi_reassembled_l1 += std::abs(current);
                        phi_delta_l1 += std::abs(used-current);
                    }
                    for (std::size_t d = 0; d < 3; ++d)
                        hbyA_delta = std::max(hbyA_delta,
                            std::abs(hbyA[d][cell] - frozen_hbyA[d][cell]));
                    std::cerr << "FROZEN_STATE_CELL label=" << label
                              << " cell=" << cell
                              << " phi_used_in_momentum_assembly=" << phi_used_l1
                              << " phi_used_in_reassembly=" << phi_reassembled_l1
                              << " phi_delta_l1=" << phi_delta_l1
                              << " HbyA_delta_linf=" << hbyA_delta
                              << " frozen_rAU=(" << frozen_rAU[0][cell] << ","
                              << frozen_rAU[1][cell] << "," << frozen_rAU[2][cell] << ")"
                              << " final_rAU=(" << diag_rAU[0][cell] << ","
                              << diag_rAU[1][cell] << "," << diag_rAU[2][cell] << ")\n";
                };
                std::cerr << "FROZEN_STATE_COMPARE phi_global_linf=" << phi_delta_linf
                          << " HbyA_global_linf=" << hbyA_delta_linf
                          << " HbyA_max_cell=" << hbyA_delta_cell << "\n";
                dump_state_cell(debug_cell, "debug");
                if (comparison_cell != debug_cell)
                    dump_state_cell(comparison_cell, "comparison");
            }
            const Offset off = mesh.cells().offsets_data()[debug_cell];
            const Offset count = mesh.cells().offsets_data()[debug_cell + 1] - off;
            std::cerr << "faces=" << count << "\n";
            for (Offset k = 0; k < count; ++k) {
                const std::size_t face = mesh.cells().faces_data()[off + k];
                const auto owner = mesh.ownership().owner(face);
                const auto neighbour = mesh.ownership().neighbour(face);
                const bool owner_side = owner == debug_cell;
                const double sign = owner_side ? 1.0 : -1.0;
                const double phi_auth = sign * mass_flux(face);
                const double phi_reconstructed = sign * reconstructed_flux(face);
                local_div_authoritative += phi_auth;
                local_div_reconstructed += phi_reconstructed;
                const Vec3 Sf = geometry.face_area_vectors[face] * sign;
                std::size_t patch = geometry.face_patch[face];
                std::string patch_name = "INTERNAL";
                int patch_type = -1;
                if (neighbour < 0 && patch < mesh.boundary().n_patches()) {
                    patch_name = mesh.boundary().patch(patch).name;
                    patch_type = static_cast<int>(mesh.boundary().patch(patch).type);
                }

                Vec3 u_face_auth{0.0, 0.0, 0.0};
                Vec3 u_face_reconstructed{0.0, 0.0, 0.0};
                Vec3 grad_face_momentum{0.0, 0.0, 0.0};
                Vec3 grad_face_correction{0.0, 0.0, 0.0};
                Vec3 hbyA_face{0.0, 0.0, 0.0};
                double dAU_n = 0.0;
                double pressure_flux_correction = 0.0;
                double rhie_chow_expected = 0.0;
                double rhie_chow_actual = 0.0;

                if (neighbour >= 0) {
                    const std::size_t ncell = static_cast<std::size_t>(neighbour);
                    u_face_auth = Vec3{
                        0.5 * (final_ux(debug_cell) + final_ux(ncell)),
                        0.5 * (final_uy(debug_cell) + final_uy(ncell)),
                        0.5 * (final_uz(debug_cell) + final_uz(ncell))};
                    hbyA_face = Vec3{
                        0.5 * (hbyA[0][debug_cell] + hbyA[0][ncell]),
                        0.5 * (hbyA[1][debug_cell] + hbyA[1][ncell]),
                        0.5 * (hbyA[2][debug_cell] + hbyA[2][ncell])};
                    // Build the local face frame from the debug-cell side.
                    // For a neighbour-side debug cell, the global owner-oriented
                    // face vector and the centre-to-centre vector have opposite
                    // orientations. Mixing those frames was the source of the
                    // NaN/incorrect pressure-flux microscope records and could
                    // make a valid internal face look non-commuting.
                    const std::size_t other_cell = owner_side
                        ? ncell : owner;
                    const Vec3 dvec = geometry.cell_centres[other_cell] -
                                      geometry.cell_centres[debug_cell];
                    const double d = dvec.mag();
                    if (!(d > 0.0) || !std::isfinite(d))
                        throw std::runtime_error(
                            "CFDX momentum microscope: degenerate internal face distance");
                    const Vec3 e{dvec.x/d, dvec.y/d, dvec.z/d};
                    const double cx = 0.5 * (
                        geometry.cell_volumes[debug_cell] * diag_rAU[0][debug_cell] +
                        geometry.cell_volumes[other_cell] * diag_rAU[0][other_cell]);
                    const double cy = 0.5 * (
                        geometry.cell_volumes[debug_cell] * diag_rAU[1][debug_cell] +
                        geometry.cell_volumes[other_cell] * diag_rAU[1][other_cell]);
                    const double cz = 0.5 * (
                        geometry.cell_volumes[debug_cell] * diag_rAU[2][debug_cell] +
                        geometry.cell_volumes[other_cell] * diag_rAU[2][other_cell]);
                    dAU_n = cx*e.x*e.x + cy*e.y*e.y + cz*e.z*e.z;
                    const Vec3 gp{
                        0.5 * (final_grad_p.component_data(0)[debug_cell] +
                               final_grad_p.component_data(0)[other_cell]),
                        0.5 * (final_grad_p.component_data(1)[debug_cell] +
                               final_grad_p.component_data(1)[other_cell]),
                        0.5 * (final_grad_p.component_data(2)[debug_cell] +
                               final_grad_p.component_data(2)[other_cell])};
                    // Compare the authoritative face-flux correction against
                    // the exact Rhie-Chow operator used by the pressure equation.
                    // This is deliberately evaluated from HbyA, not the final
                    // velocity, because phi = HbyA_f.Sf - RC(p) is the
                    // conservative pressure-velocity coupling operator.
                    // Use the exact inverse momentum diagonals reconstructed
                    // above from the final momentum rows.  The microscope is
                    // diagnostic only and must not introduce a separate
                    // SIMPLEC/SIMPLE coefficient path.
                    rhie_chow_expected = rhie_chow_pressure_flux_internal(
                        mesh, geometry, face, p, final_grad_p,
                        diag_rAU, controls.density);
                    const double hbyA_flux = controls.density * hbyA_face.dot(Sf);
                    rhie_chow_actual = hbyA_flux - phi_auth;
                    max_rhie_chow_actual = std::max(
                        max_rhie_chow_actual, std::abs(rhie_chow_actual));
                    max_rhie_chow_expected = std::max(
                        max_rhie_chow_expected, std::abs(rhie_chow_expected));
                    max_rhie_chow_mismatch = std::max(
                        max_rhie_chow_mismatch,
                        std::abs(rhie_chow_actual - rhie_chow_expected));
                    const Vec3 localSf = Sf;
                    const double normal_dp = (p(other_cell) - p(debug_cell)) / d;
                    const Vec3 dAU_face{cx, cy, cz};
                    grad_face_correction =
                        e * normal_dp + gp - e * gp.dot(e);
                    grad_face_momentum = gp;
                    pressure_flux_correction = controls.density * (
                        pressure_face_orthogonal_coefficient(
                            localSf, dvec, dAU_face) *
                            (p(other_cell) - p(debug_cell)) +
                        pressure_face_nonorthogonal_flux(
                            localSf, dvec, dAU_face, gp));
                    u_face_reconstructed = hbyA_face -
                        Vec3{cx * gp.x, cy * gp.y, cz * gp.z};
                } else {
                    if (patch < mesh.boundary().n_patches()) {
                        const auto& name = mesh.boundary().patch(patch).name;
                        const auto pit = pressure_bcs.find(name);
                        const auto uit = velocity_bcs.find(name);
                        if (uit != velocity_bcs.end() &&
                            uit->second.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                            u_face_auth = uit->second.value;
                        else
                            u_face_auth = Vec3{
                                final_ux(debug_cell), final_uy(debug_cell), final_uz(debug_cell)};
                        if (pit != pressure_bcs.end() &&
                            pit->second.type == ScalarBoundaryType::FIXED_VALUE) {
                            const Vec3 Sf_raw = geometry.face_area_vectors[face];
                            const Vec3 dvec = geometry.face_centres[face] -
                                              geometry.cell_centres[debug_cell];
                            const double d = dvec.mag();
                            if (!(d > 0.0) || !std::isfinite(d))
                                throw std::runtime_error(
                                    "CFDX momentum microscope: degenerate boundary face distance");
                            const Vec3 e{dvec.x/d, dvec.y/d, dvec.z/d};
                            const Vec3 dAU_boundary{
                                geometry.cell_volumes[debug_cell] * diag_rAU[0][debug_cell],
                                geometry.cell_volumes[debug_cell] * diag_rAU[1][debug_cell],
                                geometry.cell_volumes[debug_cell] * diag_rAU[2][debug_cell]};
                            dAU_n = dAU_boundary.x*e.x*e.x +
                                    dAU_boundary.y*e.y*e.y +
                                    dAU_boundary.z*e.z*e.z;
                            const double normal_dp =
                                (pit->second.value - p(debug_cell)) / d;
                            const Vec3 gp{
                                final_grad_p.component_data(0)[debug_cell],
                                final_grad_p.component_data(1)[debug_cell],
                                final_grad_p.component_data(2)[debug_cell]};
                            pressure_flux_correction = controls.density * (
                                pressure_face_orthogonal_coefficient(
                                    Sf_raw, dvec, dAU_boundary) *
                                    (pit->second.value - p(debug_cell)) +
                                pressure_face_nonorthogonal_flux(
                                    Sf_raw, dvec, dAU_boundary, gp));
                            grad_face_correction = e * normal_dp + gp - e * gp.dot(e);
                        } else {
                            grad_face_correction = Vec3{
                                final_grad_p.component_data(0)[debug_cell],
                                final_grad_p.component_data(1)[debug_cell],
                                final_grad_p.component_data(2)[debug_cell]};
                        }
                        grad_face_momentum = Vec3{
                            final_grad_p.component_data(0)[debug_cell],
                            final_grad_p.component_data(1)[debug_cell],
                            final_grad_p.component_data(2)[debug_cell]};
                        hbyA_face = Vec3{
                            hbyA[0][debug_cell],
                            hbyA[1][debug_cell],
                            hbyA[2][debug_cell]};
                        const Vec3 gp = grad_face_momentum;
                        u_face_reconstructed = hbyA_face -
                            Vec3{dAU_n*gp.x, dAU_n*gp.y, dAU_n*gp.z};
                    }
                }

                const double momentum_pressure_flux = grad_face_momentum.dot(Sf);
                max_delta_phi = std::max(
                    max_delta_phi, std::abs(phi_auth - phi_reconstructed));
                max_delta_grad_p = std::max(
                    max_delta_grad_p,
                    (grad_face_momentum - grad_face_correction).mag());
                max_delta_u_face = std::max(
                    max_delta_u_face, (u_face_auth - u_face_reconstructed).mag());
                max_momentum_pressure_flux = std::max(
                    max_momentum_pressure_flux, std::abs(momentum_pressure_flux));
                max_correction_pressure_flux = std::max(
                    max_correction_pressure_flux, std::abs(pressure_flux_correction));

                std::cerr << "  face=" << face
                          << " owner=" << owner
                          << " neighbour=" << neighbour
                          << " patch=" << patch_name
                          << " patch_type=" << patch_type
                          << " Sf=(" << Sf.x << "," << Sf.y << "," << Sf.z << ")"
                          << " phi_authoritative=" << phi_auth
                          << " phi_reconstructed_from_Uface=" << phi_reconstructed
                          << " delta_phi=" << (phi_auth - phi_reconstructed)
                          << " U_face_authoritative=(" << u_face_auth.x << ","
                          << u_face_auth.y << "," << u_face_auth.z << ")"
                          << " U_face_reconstructed_HbyA_dAUgradp=("
                          << u_face_reconstructed.x << "," << u_face_reconstructed.y
                          << "," << u_face_reconstructed.z << ")"
                          << " delta_U_face=" << (u_face_auth-u_face_reconstructed).mag()
                          << " dAU=" << dAU_n
                          << " gradp_face_momentum=(" << grad_face_momentum.x << ","
                          << grad_face_momentum.y << "," << grad_face_momentum.z << ")"
                          << " gradp_face_correction=(" << grad_face_correction.x << ","
                          << grad_face_correction.y << "," << grad_face_correction.z << ")"
                          << " delta_gradp="
                          << (grad_face_momentum-grad_face_correction).mag()
                          << " momentum_pressure_flux=" << momentum_pressure_flux
                          << " correction_pressure_flux=" << pressure_flux_correction
                          << " rhie_chow_expected=" << rhie_chow_expected
                          << " rhie_chow_actual=" << rhie_chow_actual
                          << " delta_rhie_chow="
                          << (rhie_chow_actual-rhie_chow_expected)
                          << "\n";
            }

            std::cerr << "local_mass_flux_divergence=" << local_div_authoritative
                      << " local_reconstructed_flux_divergence=" << local_div_reconstructed
                      << " reconstructed_velocity_continuity_linf="
                      << h.reconstructed_velocity_continuity_linf
                      << " delta_phi_linf=" << max_delta_phi
                      << " delta_u_face_linf=" << max_delta_u_face
                      << " delta_gradp_linf=" << max_delta_grad_p
                      << " momentum_pressure_flux_linf=" << max_momentum_pressure_flux
                      << " correction_pressure_flux_linf=" << max_correction_pressure_flux
                      << " rhie_chow_expected_linf=" << max_rhie_chow_expected
                      << " rhie_chow_actual_linf=" << max_rhie_chow_actual
                      << " delta_rhie_chow_linf=" << max_rhie_chow_mismatch
                      << "\n";
            std::cerr << "diagnostic_global=(" << h.momentum_equation_residual_components[0]
                      << "," << h.momentum_equation_residual_components[1]
                      << "," << h.momentum_equation_residual_components[2] << ")"
                      << " internal=" << h.momentum_equation_residual_internal
                      << " boundary=" << h.momentum_equation_residual_boundary
                      << " worst_cell=" << h.momentum_residual_cell
                      << " worst_patch=" << h.momentum_residual_patch << "\n";
            std::cerr << "no_pressure=" << h.momentum_residual_no_pressure
                      << " pressure_contribution=" << h.momentum_pressure_contribution
                      << " gradp_linf=" << h.pressure_gradient_linf
                      << "\n";
            std::cerr << "=== END CFDX MOMENTUM MICROSCOPE ===\n";
        }
        {
            const auto conservation =
                cfdx::core::audit_face_flux_conservation(mesh, mass_flux);
            const double scale =
                std::max(conservation.global_abs_boundary_flux,
                         std::max(conservation.l1_local_imbalance, 1.0));
            h.mass_boundary_flux = conservation.global_boundary_flux;
            h.mass_global_cell_balance = conservation.global_cell_balance;
            h.mass_local_l1 = conservation.l1_local_imbalance;
            h.mass_local_linf = conservation.max_local_imbalance;
            h.mass_local_l2 = conservation.l2_local_imbalance;
            h.mass_normalized_imbalance =
                conservation.l1_local_imbalance / scale;
            h.mass_worst_cell = conservation.worst_cell;
            h.mass_nonfinite_faces = conservation.nonfinite_faces;
            for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
                for (std::size_t d = 0; d < 3; ++d) {
                    if (!std::isfinite(U.component_data(d)[c]))
                        ++h.boundedness_nonfinite_velocity;
                }
            }

            // Independent post-solve momentum conservation reconstruction.
            // The flux is rebuilt from the accepted final field, independently
            // of the assembled matrix residual. Internal faces therefore have
            // one owner-oriented flux and are audited separately for balance.
            for (std::size_t d = 0; d < 3; ++d) {
                cfdx::core::Field<double, cfdx::core::Location::CELL> component(
                    mesh.n_cells(), "U_reconstructed", "m/s", 1);
                cfdx::core::Field<double, cfdx::core::Location::CELL> source(
                    mesh.n_cells(), "momentum_source_reconstructed", "N/m3", 1);
                ScalarBoundaryConditions scalar_bcs;
                for (const auto& [name, bc] : velocity_bcs) {
                    if (bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                        scalar_bcs[name] = {ScalarBoundaryType::FIXED_VALUE, bc.value.x, 0.0};
                    else
                        scalar_bcs[name] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
                    if (d == 1 && bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                        scalar_bcs[name].value = bc.value.y;
                    if (d == 2 && bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                        scalar_bcs[name].value = bc.value.z;
                }
                for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
                    component(cell) = U.component_data(d)[cell];
                    // Same channel selection as the momentum assembly
                    // (cellwise replaces the uniform vector when
                    // supplied), preserving this audit's own convention
                    // for the source term.
                    const bool cellwise_body_force =
                        !controls.cell_body_force.empty();
                    const double body_component =
                        cellwise_body_force
                            ? (d == 0
                                   ? controls.cell_body_force[cell].x
                                   : (d == 1
                                          ? controls.cell_body_force[cell].y
                                          : controls.cell_body_force[cell].z))
                            : (d == 0
                                   ? controls.body_force.x
                                   : (d == 1
                                          ? controls.body_force.y
                                          : controls.body_force.z));
                    source(cell) = body_component - final_grad_p.component_data(d)[cell];
                }
                const auto flux = reconstruct_scalar_transport_flux(
                    mesh, geometry, mass_flux, component, 0.0, scalar_bcs,
                    controls.use_bounded_convection, controls.convection_scheme,
                    nullptr, &mu_eff);
                const auto balance = cfdx::core::audit_integrated_balance(
                    mesh, flux, source, geometry.cell_volumes);
                h.momentum_conservation_residual[d] = balance.residual;
                h.momentum_conservation_normalized[d] = balance.normalized_residual;
                h.momentum_conservation_worst_cell[d] = balance.worst_cell;
            }
        }
        result.history.push_back(h);

        // Samples are retained regardless of probe_callback. A probe that was
        // requested and sampled but not retained left no artifact behind, so
        // there was nothing to export or validate after the run.
        for (const auto& probe : controls.probes) {
            if (probe.name.empty())
                throw std::invalid_argument("solve_steady_incompressible: probe name must not be empty");
            const double value = sample_incompressible_probe(probe, mesh, geometry, U, p);
            if (!std::isfinite(value))
                throw std::runtime_error("solve_steady_incompressible: probe value is not finite");
            result.probe_samples.push_back(
                IncompressibleProbeSample{probe.name, iter, 0.0, value});
            if (controls.probe_callback)
                controls.probe_callback(result.probe_samples.back());
        }

        if (controls.iteration_output_callback &&
            !controls.iteration_output_callback(iter, 0.0, mesh, U, p)) {
            result.iterations = iter;
            stop_after_iteration = true;
        }

        result.convergence_status = convergence_report.status;
        result.convergence_reason = convergence_report.reason;
        if ((convergence_report.status == cfdx::core::ConvergenceStatus::CONVERGED ||
             convergence_report.status == cfdx::core::ConvergenceStatus::MAX_ITERATIONS) &&
            !sst_residual_gate) {
            // The carrier monitor does not know about the additional SST
            // closure/transport gate. Do not publish a converged status while
            // the accepted momentum residual still corresponds to an earlier
            // effective-viscosity field.
            result.convergence_status =
                iter == controls.convergence.max_iterations
                    ? cfdx::core::ConvergenceStatus::MAX_ITERATIONS
                    : cfdx::core::ConvergenceStatus::CONTINUE;
            result.convergence_reason =
                "carrier criteria satisfied but SST transport/closure convergence gate is not satisfied";
        }
        if (converged_now) {
            result.converged = true;
            result.iterations = iter;
            stop_after_iteration = true;
        }
        if (convergence_report.status == cfdx::core::ConvergenceStatus::DIVERGED ||
            convergence_report.status == cfdx::core::ConvergenceStatus::STAGNATED) {
            result.iterations = iter;
            stop_after_iteration = true;
        }
        result.iterations = iter;
            transaction.commit();
            if (controls.sst_runtime) {
                for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell)
                    mu_eff[cell] = controls.density *
                        (controls.kinematic_viscosity +
                         controls.sst_runtime->state->nu_t(cell));
            }
            // Only publish a face flux after the nonlinear transaction has
            // accepted the corresponding U/p state. mass_flux is mutable during
            // a tentative iteration and is not part of NonlinearStateRollback;
            // publishing it after a rejected attempt would therefore expose a
            // flux belonging to a state that was rolled back.
            result.authoritative_mass_flux = mass_flux;
            iteration_completed = true;
            } catch (const NonlinearRetryableFailure& error) {
                transaction.reject();
                mass_flux = mass_flux_snapshot;
                if (controls.sst_runtime && sst_attempt_snapshot)
                    restore_sst_production_state(*controls.sst_runtime->state, *sst_attempt_snapshot);
                if (controls.sst_runtime) {
                    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell)
                        mu_eff[cell] = controls.density *
                            (controls.kinematic_viscosity +
                             controls.sst_runtime->state->nu_t(cell));
                }
                // The monitor and exported iteration evidence are tentative
                // until the carrier/SST transaction commits. Replaying a retry
                // with the same iteration index must not see stale samples or
                // a previously published converged flag.
                convergence_monitor = convergence_monitor_snapshot;
                result.history.resize(result_history_size);
                result.probe_samples.resize(probe_samples_size);
                result.converged = previous_result_converged;
                result.convergence_status = previous_result_status;
                result.convergence_reason = previous_result_reason;
                result.iterations = previous_result_iterations;
                stop_after_iteration = false;
                frozen_state_valid = false;
                if (!retry_controller.can_retry())
                    throw;
                retry_controller.reject();
                std::cerr << "CFDX nonlinear iteration " << iter
                          << " rejected; retry=" << retry_controller.retries()
                          << " alpha_u=" << retry_controller.alpha_u(effective_alpha_u)
                          << " alpha_p=" << retry_controller.alpha_p(effective_alpha_p)
                          << " reason=" << error.what() << "\n";
            } catch (...) {
                transaction.reject();
                mass_flux = mass_flux_snapshot;
                if (controls.sst_runtime && sst_attempt_snapshot)
                    restore_sst_production_state(*controls.sst_runtime->state, *sst_attempt_snapshot);
                throw;
            }
        }
        if (stop_after_iteration)
            break;
    }

    // authoritative_mass_flux is updated only after an accepted nonlinear
    // transaction above. In particular, a rejected retry cannot overwrite the
    // evidence for the last accepted U/p state. Do not reconstruct it from U.
    if (pressure_context)
        result.pressure_linear_context = pressure_context->stats();
    return result;
}

} // namespace cfdx::physics


```

### `src/cfdx/physics/wall_distance.h`

SHA256 `355b00d3c75ad6388571ccf497683e47eb0ed77e79d477ddc56f6d584882d4ad`

```cpp
#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::physics {

struct WallDistanceVec3 {
    double x{0.0}, y{0.0}, z{0.0};
    WallDistanceVec3 operator+(const WallDistanceVec3& b) const { return {x+b.x,y+b.y,z+b.z}; }
    WallDistanceVec3 operator-(const WallDistanceVec3& b) const { return {x-b.x,y-b.y,z-b.z}; }
    WallDistanceVec3 operator*(double s) const { return {x*s,y*s,z*s}; }
};

inline double wd_dot(const WallDistanceVec3& a,const WallDistanceVec3& b) {
    return a.x*b.x+a.y*b.y+a.z*b.z;
}
inline double wd_norm2(const WallDistanceVec3& a) { return wd_dot(a,a); }
inline double wd_norm(const WallDistanceVec3& a) { return std::sqrt(wd_norm2(a)); }

struct WallTriangle { std::array<std::size_t,3> v{}; };

struct WallSurface {
    std::vector<WallDistanceVec3> points;
    std::vector<WallTriangle> triangles;
};

enum class WallDistanceConvergenceStatus {
    DIRECT,
    CONVERGED,
    MAX_ITER,
    RESIDUAL_TOO_HIGH,
    RESIDUAL_STAGNATION,
    NOT_REPORTED,
    INVALID_FIELD
};

struct WallDistanceResult {
    std::vector<double> distance;
    std::vector<unsigned char> valid;
    std::string method;
    std::size_t iterations{0};
    double residual_inf{0.0};
    std::size_t auxiliary_iterations{0};
    double auxiliary_residual_inf{0.0};
    std::size_t poisson_wall_ray_hits{0};
    std::size_t poisson_wall_ray_misses{0};
    std::size_t poisson_wall_fallbacks{0};
    std::size_t poisson_wall_bad_alignment{0};
    double poisson_wall_min_alignment{1.0};
    double poisson_phi_min{0.0};
    double poisson_phi_max{0.0};
    double poisson_grad_min{0.0};
    double poisson_grad_max{0.0};
    double poisson_distance_l2_error{0.0};
    double poisson_distance_linf_error{0.0};
    // F-W1 (N11 backlog, #480): number of nearest_normal() calls that hit a
    // degenerate (zero-area) nearest triangle and returned the {1,0,0}
    // fallback instead of a real surface normal. Filled by the methods that
    // query the wall normal (directional mesh wave, Poisson cut-face audit);
    // zero means every queried normal was a genuine geometric normal.
    std::size_t nearest_normal_fallbacks{0};
    // F-W2 (N11 backlog, #480): number of pseudo-time updates in the explicit
    // transport steppers (Hamilton-Jacobi, hybrid Poisson+HJ) whose candidate
    // distance was negative and was clamped to zero. The clamp itself is
    // unchanged; this counter only makes it attributable in the result.
    std::size_t positivity_clamp_activations{0};
    bool converged{false};
    std::string stopping_reason{"not_reported"};
    WallDistanceConvergenceStatus convergence_status{WallDistanceConvergenceStatus::NOT_REPORTED};
};

enum class WallDistanceMethod {
    EXACT_GEOMETRIC,
    SEARCH_BASED,
    MESH_WAVE,
    DIRECTIONAL_MESH_WAVE,
    POISSON,
    EIKONAL,
    HAMILTON_JACOBI,
    ADVECTION_DIFFUSION,
    HYBRID_POISSON_EIKONAL
};

// Recommended default wall-distance method.
//
// EIKONAL (fast sweeping of |grad d| = 1) is the accuracy/cost sweet spot
// documented in the literature (Zhao, Math. Comp. 2005; Sethian 1996) and
// verified by our own complex-geometry benchmark:
//
//   * exact_geometric / search_based : 0 % L2 error, but O(N_cell * N_face)
//     for exact and requires a global BVH sweep for search.
//   * eikonal                        : ~7 % L2 error in 3 fast-sweep passes.
//   * hamilton_jacobi                : ~13 % L2 error, ~1000+ iterations
//     (O(n^2) pseudo-time). Its H-J correction only helps on smooth
//     single-body shapes; on multi-body wing-body-tail it slightly worsens
//     accuracy compared to plain eikonal because the Poisson-integrated
//     "seed" front is contaminated by contributions from every wall.
//   * advection_diffusion / hybrid   : 10-33 % L2 error, thousand+
//     iterations. Same multi-body limitation as HJ; the hybrid variant
//     compounds it because it starts from a Poisson-reconstructed field
//     and the reconstruction is intrinsically 15-30 % accurate in the
//     literature (Tucker 2003, Fares-Schroder 2002) for this geometry
//     class.
//
// The Poisson family is retained because it produces a smooth d(x) with
// analytic dphi/dn, which is what Spalart-Allmaras and k-omega SST source
// terms actually need -- OpenFOAM's PoissonPatchDist has the same
// rationale. Do not use it when the value of d itself matters.
inline constexpr WallDistanceMethod default_wall_distance_method =
    WallDistanceMethod::EIKONAL;

inline const char* wall_distance_method_name(WallDistanceMethod m) {
    switch(m) {
        case WallDistanceMethod::EXACT_GEOMETRIC: return "exact_geometric";
        case WallDistanceMethod::SEARCH_BASED: return "search_based";
        case WallDistanceMethod::MESH_WAVE: return "mesh_wave";
        case WallDistanceMethod::DIRECTIONAL_MESH_WAVE: return "directional_mesh_wave";
        case WallDistanceMethod::POISSON: return "poisson";
        case WallDistanceMethod::EIKONAL: return "eikonal";
        case WallDistanceMethod::HAMILTON_JACOBI: return "hamilton_jacobi";
        case WallDistanceMethod::ADVECTION_DIFFUSION: return "advection_diffusion";
        case WallDistanceMethod::HYBRID_POISSON_EIKONAL: return "hybrid_poisson_eikonal";
    }
    return "unknown";
}

// Iteration budget for the iterative methods, indexed by a linear resolution n.
//
// A single constant budget across all resolutions and methods silently turns a
// convergence gate into a resolution gate. Each solver has its own cost law:
//
//   * Poisson is a Krylov solve; its iteration count grows like O(n). Measured
//     46 / 98 / 207 / 461 at n=16 / 32 / 64 / 128, a clean factor 2.2 per
//     refinement.
//   * Hamilton-Jacobi, advection-diffusion and the hybrid all take an explicit
//     pseudo-time step bounded by dt = O(h^2), so the number of steps to reach
//     steady state grows like O(n^2). Measured with the old fixed budget of
//     500: advection-diffusion needed 68 steps at n=16 and 224 at n=32
//     (factor 3.3), then hit the ceiling at n=64 with residual 4.3e-5.
//
// Budgets below are those cost laws, normalized at n=16 with roughly a factor 4
// of margin over the counts measured there (46 for Poisson, 208 for the slowest
// pseudo-time method). Methods that are direct or purely graph-based ignore the
// budget entirely.
inline std::size_t wall_distance_iteration_budget(WallDistanceMethod method,
                                                  std::size_t n) {
    const double ratio=static_cast<double>(n)/16.0;
    switch(method) {
    case WallDistanceMethod::POISSON:
        return static_cast<std::size_t>(200.0*ratio);
    case WallDistanceMethod::HAMILTON_JACOBI:
    case WallDistanceMethod::ADVECTION_DIFFUSION:
    case WallDistanceMethod::HYBRID_POISSON_EIKONAL:
        return static_cast<std::size_t>(1000.0*ratio*ratio);
    default:
        return 500;
    }
}

inline double point_triangle_distance2(const WallDistanceVec3& p,
                                       const WallDistanceVec3& a,
                                       const WallDistanceVec3& b,
                                       const WallDistanceVec3& c) {
    const WallDistanceVec3 ab=b-a, ac=c-a, ap=p-a;
    const double d1=wd_dot(ab,ap), d2=wd_dot(ac,ap);
    if(d1<=0.0 && d2<=0.0) return wd_norm2(ap);
    const WallDistanceVec3 bp=p-b;
    const double d3=wd_dot(ab,bp), d4=wd_dot(ac,bp);
    if(d3>=0.0 && d4<=d3) return wd_norm2(bp);
    const double vc=d1*d4-d3*d2;
    if(vc<=0.0 && d1>=0.0 && d3<=0.0) {
        const double v=d1/(d1-d3);
        return wd_norm2(p-(a+ab*v));
    }
    const WallDistanceVec3 cp=p-c;
    const double d5=wd_dot(ab,cp), d6=wd_dot(ac,cp);
    if(d6>=0.0 && d5<=d6) return wd_norm2(cp);
    const double vb=d5*d2-d1*d6;
    if(vb<=0.0 && d2>=0.0 && d6<=0.0) {
        const double w=d2/(d2-d6);
        return wd_norm2(p-(a+ac*w));
    }
    const double va=d3*d6-d5*d4;
    if(va<=0.0 && (d4-d3)>=0.0 && (d5-d6)>=0.0) {
        const double w=(d4-d3)/((d4-d3)+(d5-d6));
        return wd_norm2(p-(b+(c-b)*w));
    }
    const double denom=1.0/(va+vb+vc);
    const double v=vb*denom, w=vc*denom;
    return wd_norm2(p-(a+ab*v+ac*w));
}

struct WallDistanceAabb {
    WallDistanceVec3 lo{std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::infinity()};
    WallDistanceVec3 hi{-std::numeric_limits<double>::infinity(),
                         -std::numeric_limits<double>::infinity(),
                         -std::numeric_limits<double>::infinity()};
};

inline WallDistanceAabb wall_distance_expand(WallDistanceAabb box,
                                              const WallDistanceVec3& p) {
    box.lo.x=std::min(box.lo.x,p.x); box.lo.y=std::min(box.lo.y,p.y); box.lo.z=std::min(box.lo.z,p.z);
    box.hi.x=std::max(box.hi.x,p.x); box.hi.y=std::max(box.hi.y,p.y); box.hi.z=std::max(box.hi.z,p.z);
    return box;
}

inline double wall_distance_aabb_distance2(const WallDistanceAabb& box,
                                            const WallDistanceVec3& p) {
    const auto axis=[](double x,double lo,double hi) {
        return x<lo ? lo-x : (x>hi ? x-hi : 0.0);
    };
    const double dx=axis(p.x,box.lo.x,box.hi.x);
    const double dy=axis(p.y,box.lo.y,box.hi.y);
    const double dz=axis(p.z,box.lo.z,box.hi.z);
    return dx*dx+dy*dy+dz*dz;
}

class WallDistanceBvh {
public:
    explicit WallDistanceBvh(const WallSurface& surface) : surface_(surface) {
        if(surface_.points.empty() || surface_.triangles.empty())
            throw std::invalid_argument("wall distance BVH surface is empty");
        indices_.resize(surface_.triangles.size());
        for(std::size_t i=0;i<indices_.size();++i) indices_[i]=i;
        nodes_.reserve(indices_.size()*2);
        build(0,indices_.size());
    }

    double nearest_distance(const WallDistanceVec3& p) const {
        double best2=std::numeric_limits<double>::infinity();
        nearest(0,p,best2);
        return std::sqrt(best2);
    }

    double nearest_distance_within(const WallDistanceVec3& p,double threshold) const {
        if(threshold<0.0) return std::numeric_limits<double>::infinity();
        double best2=threshold*threshold;
        bool found=false;
        nearest_within(0,p,best2,found);
        return found ? std::sqrt(best2) : std::numeric_limits<double>::infinity();
    }

    // Distance from p to the first wall intersection along a grid-face ray.
    // This is the correct geometric quantity for a cut Dirichlet face:
    // the boundary lies on the segment joining the fluid node to the solid
    // neighbour.  nearest_distance()/nearest_normal() are not sufficient
    // near edges/corners because the closest surface point need not lie on
    // that face-normal ray.
    double ray_distance(const WallDistanceVec3& p,
                        const WallDistanceVec3& direction,
                        double max_distance) const {
        const double dn=wd_norm(direction);
        if(!(dn>1e-30) || !(max_distance>0.0)) return std::numeric_limits<double>::infinity();
        const WallDistanceVec3 d=direction*(1.0/dn);
        const auto cross=[](const WallDistanceVec3& a,const WallDistanceVec3& b) {
            return WallDistanceVec3{
                a.y*b.z-a.z*b.y,
                a.z*b.x-a.x*b.z,
                a.x*b.y-a.y*b.x
            };
        };
        double best=std::numeric_limits<double>::infinity();
        const double eps=1e-12*std::max(1.0,max_distance);

        for(const auto& tri:surface_.triangles) {
            const auto& v0=surface_.points[tri.v[0]];
            const auto& v1=surface_.points[tri.v[1]];
            const auto& v2=surface_.points[tri.v[2]];
            const auto e1=v1-v0;
            const auto e2=v2-v0;
            const auto h=cross(d,e2);
            const double a=wd_dot(e1,h);
            if(std::abs(a)<1e-14) continue;
            const double inv_a=1.0/a;
            const auto s=p-v0;
            const double u=inv_a*wd_dot(s,h);
            if(u < -1e-12 || u > 1.0+1e-12) continue;
            const auto q=cross(s,e1);
            const double v=inv_a*wd_dot(d,q);
            if(v < -1e-12 || u+v > 1.0+1e-12) continue;
            const double t=inv_a*wd_dot(e2,q);
            if(t>=-eps && t<=max_distance+eps) best=std::min(best,std::max(0.0,t));
        }
        return best;
    }

    std::size_t triangle_count() const { return surface_.triangles.size(); }
    std::size_t node_count() const { return nodes_.size(); }

    // F-W1 (N11 backlog, #480): count of nearest_normal() calls whose nearest
    // triangle has a degenerate (zero-area) normal, so the {1,0,0} fallback
    // was returned instead of a real surface normal. The fallback value and
    // direction are unchanged; this counter only makes the silent fallback
    // attributable by the caller. Per BVH instance, monotonic, not reset.
    //
    // Concurrency contract (audit of #865): all queries of a constructed
    // BVH are read-only, and this counter is the only mutable state. It is
    // a std::atomic incremented with relaxed ordering, so concurrent
    // const nearest_normal() calls on the same instance are well defined:
    // the counter equals the exact number of fallbacks after any storm of
    // calls, and the returned fallback normal is unaffected. Relaxed
    // ordering is sufficient because the counter is a diagnostic tally,
    // never a synchronization point. Note: making the counter atomic
    // deletes the implicit copy constructor/assignment of WallDistanceBvh;
    // no current caller copies an instance (construction from a
    // WallSurface plus const-reference use only).
    std::size_t degenerate_normal_fallbacks() const {
        return degenerate_normal_fallbacks_.load(std::memory_order_relaxed);
    }

    WallDistanceVec3 nearest_normal(const WallDistanceVec3& p) const {
        double best2=std::numeric_limits<double>::infinity();
        std::size_t best=0;
        nearest_triangle(0,p,best2,best);
        const auto& t=surface_.triangles[best];
        const auto ab=surface_.points[t.v[1]]-surface_.points[t.v[0]];
        const auto ac=surface_.points[t.v[2]]-surface_.points[t.v[0]];
        const WallDistanceVec3 n{ab.y*ac.z-ab.z*ac.y,
                                 ab.z*ac.x-ab.x*ac.z,
                                 ab.x*ac.y-ab.y*ac.x};
        const double nn=wd_norm(n);
        if(nn>1e-30) return n*(1.0/nn);
        // Degenerate triangle (collinear vertices): the cross product is the
        // zero vector and no geometric normal exists. The historical
        // behaviour is preserved exactly (unit +x fallback); the counter
        // above is what makes the fallback observable (F-W1). Relaxed
        // fetch_add keeps the tally exact under concurrent const calls
        // (see the concurrency contract on degenerate_normal_fallbacks()).
        degenerate_normal_fallbacks_.fetch_add(1,std::memory_order_relaxed);
        return WallDistanceVec3{1.0,0.0,0.0};
    }

private:
    struct Node {
        WallDistanceAabb box{};
        std::size_t begin{0}, end{0};
        std::size_t left{std::numeric_limits<std::size_t>::max()};
        std::size_t right{std::numeric_limits<std::size_t>::max()};
        bool leaf{false};
    };

    const WallSurface& surface_;
    std::vector<std::size_t> indices_;
    std::vector<Node> nodes_;
    // F-W1: mutated through const nearest_normal() so the degenerate-normal
    // fallback stays observable without changing the const query contract.
    // Atomic (relaxed) so the const-query API is also well defined under
    // concurrent calls on one instance; see the contract note on
    // degenerate_normal_fallbacks().
    mutable std::atomic<std::size_t> degenerate_normal_fallbacks_{0};
    static constexpr std::size_t leaf_size_=4;

    WallDistanceAabb triangle_box(std::size_t ti) const {
        const auto& t=surface_.triangles[ti];
        WallDistanceAabb box;
        box=wall_distance_expand(box,surface_.points[t.v[0]]);
        box=wall_distance_expand(box,surface_.points[t.v[1]]);
        box=wall_distance_expand(box,surface_.points[t.v[2]]);
        return box;
    }

    WallDistanceVec3 triangle_centroid(std::size_t ti) const {
        const auto& t=surface_.triangles[ti];
        return (surface_.points[t.v[0]]+surface_.points[t.v[1]]+surface_.points[t.v[2]])*(1.0/3.0);
    }

    std::size_t build(std::size_t begin,std::size_t end) {
        const std::size_t node_id=nodes_.size();
        nodes_.push_back({});
        Node& node=nodes_.back();
        node.begin=begin; node.end=end;
        for(std::size_t i=begin;i<end;++i) {
            const auto box=triangle_box(indices_[i]);
            node.box=wall_distance_expand(node.box,box.lo);
            node.box=wall_distance_expand(node.box,box.hi);
        }
        const std::size_t count=end-begin;
        if(count<=leaf_size_) {
            node.leaf=true;
            return node_id;
        }
        const auto extent=node.box.hi-node.box.lo;
        int axis=0;
        if(extent.y>extent.x && extent.y>=extent.z) axis=1;
        else if(extent.z>extent.x && extent.z>extent.y) axis=2;
        const std::size_t mid=begin+count/2;
        auto coord=[axis](const WallDistanceVec3& p) {
            return axis==0?p.x:(axis==1?p.y:p.z);
        };
        std::nth_element(indices_.begin()+static_cast<std::ptrdiff_t>(begin),
                         indices_.begin()+static_cast<std::ptrdiff_t>(mid),
                         indices_.begin()+static_cast<std::ptrdiff_t>(end),
                         [&](std::size_t a,std::size_t b) {
                             return coord(triangle_centroid(a))<coord(triangle_centroid(b));
                         });
        const std::size_t left=build(begin,mid);
        const std::size_t right=build(mid,end);
        nodes_[node_id].left=left;
        nodes_[node_id].right=right;
        return node_id;
    }

    void nearest_triangle(std::size_t node_id,const WallDistanceVec3& p,
                          double& best2,std::size_t& best) const {
        const Node& node=nodes_[node_id];
        const double box_d2=wall_distance_aabb_distance2(node.box,p);
        if(box_d2>=best2) return;
        if(node.leaf) {
            for(std::size_t i=node.begin;i<node.end;++i) {
                const std::size_t ti=indices_[i];
                const auto& t=surface_.triangles[ti];
                const double d2=point_triangle_distance2(
                    p,surface_.points[t.v[0]],surface_.points[t.v[1]],surface_.points[t.v[2]]);
                if(d2<best2) { best2=d2; best=ti; }
            }
            return;
        }
        const double dl=wall_distance_aabb_distance2(nodes_[node.left].box,p);
        const double dr=wall_distance_aabb_distance2(nodes_[node.right].box,p);
        if(dl<dr) { nearest_triangle(node.left,p,best2,best); nearest_triangle(node.right,p,best2,best); }
        else { nearest_triangle(node.right,p,best2,best); nearest_triangle(node.left,p,best2,best); }
    }

    void nearest_within(std::size_t node_id,const WallDistanceVec3& p,
                         double& best2,bool& found) const {
        const Node& node=nodes_[node_id];
        const double box_d2=wall_distance_aabb_distance2(node.box,p);
        if(box_d2>best2) return;
        if(node.leaf) {
            for(std::size_t i=node.begin;i<node.end;++i) {
                const auto& t=surface_.triangles[indices_[i]];
                const double d2=point_triangle_distance2(
                    p,surface_.points[t.v[0]],surface_.points[t.v[1]],surface_.points[t.v[2]]);
                if(d2<=best2) {
                    best2=d2;
                    found=true;
                }
            }
            return;
        }
        const Node& left=nodes_[node.left];
        const Node& right=nodes_[node.right];
        const double dl=wall_distance_aabb_distance2(left.box,p);
        const double dr=wall_distance_aabb_distance2(right.box,p);
        if(dl<dr) {
            nearest_within(node.left,p,best2,found);
            nearest_within(node.right,p,best2,found);
        } else {
            nearest_within(node.right,p,best2,found);
            nearest_within(node.left,p,best2,found);
        }
    }

    void nearest(std::size_t node_id,const WallDistanceVec3& p,double& best2) const {
        const Node& node=nodes_[node_id];
        if(wall_distance_aabb_distance2(node.box,p)>=best2) return;
        if(node.leaf) {
            for(std::size_t i=node.begin;i<node.end;++i) {
                const auto& t=surface_.triangles[indices_[i]];
                best2=std::min(best2,point_triangle_distance2(
                    p,surface_.points[t.v[0]],surface_.points[t.v[1]],surface_.points[t.v[2]]));
            }
            return;
        }
        const Node& left=nodes_[node.left];
        const Node& right=nodes_[node.right];
        const double dl=wall_distance_aabb_distance2(left.box,p);
        const double dr=wall_distance_aabb_distance2(right.box,p);
        if(dl<dr) {
            nearest(node.left,p,best2);
            nearest(node.right,p,best2);
        } else {
            nearest(node.right,p,best2);
            nearest(node.left,p,best2);
        }
    }
};

inline double exact_point_distance(const WallSurface& s,const WallDistanceVec3& p) {
    if(s.points.empty() || s.triangles.empty())
        throw std::invalid_argument("wall distance surface is empty");
    double d2=std::numeric_limits<double>::infinity();
    for(const auto& t:s.triangles) {
        d2=std::min(d2,point_triangle_distance2(p,s.points[t.v[0]],s.points[t.v[1]],s.points[t.v[2]]));
    }
    return std::sqrt(d2);
}

struct WallDistanceGrid {
    std::size_t nx{0}, ny{0}, nz{0};
    WallDistanceVec3 origin{};
    WallDistanceVec3 spacing{1.0,1.0,1.0};
    std::vector<WallDistanceVec3> points;
    std::vector<unsigned char> solid;
    std::vector<std::vector<std::size_t>> neighbours;

    std::size_t index(std::size_t i,std::size_t j,std::size_t k) const {
        return (k*ny+j)*nx+i;
    }
};

inline WallDistanceGrid make_wall_distance_grid(std::size_t nx,std::size_t ny,std::size_t nz,
                                                 const WallDistanceVec3& origin,
                                                 const WallDistanceVec3& spacing,
                                                 const std::function<bool(const WallDistanceVec3&)>& inside) {
    if(nx<3 || ny<3 || nz<3) throw std::invalid_argument("wall distance grid too small");
    WallDistanceGrid g; g.nx=nx; g.ny=ny; g.nz=nz; g.origin=origin; g.spacing=spacing;
    const std::size_t n=nx*ny*nz;
    g.points.resize(n); g.solid.resize(n,0); g.neighbours.resize(n);
    for(std::size_t k=0;k<nz;++k) for(std::size_t j=0;j<ny;++j) for(std::size_t i=0;i<nx;++i) {
        const std::size_t id=g.index(i,j,k);
        g.points[id]={origin.x+spacing.x*static_cast<double>(i),
                      origin.y+spacing.y*static_cast<double>(j),
                      origin.z+spacing.z*static_cast<double>(k)};
        g.solid[id]=inside(g.points[id])?1:0;
        constexpr std::array<std::array<int,3>,6> dirs={{{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}}};
        for(const auto& d:dirs) {
            const long ii=static_cast<long>(i)+d[0], jj=static_cast<long>(j)+d[1], kk=static_cast<long>(k)+d[2];
            if(ii>=0 && jj>=0 && kk>=0 && ii<static_cast<long>(nx) && jj<static_cast<long>(ny) && kk<static_cast<long>(nz))
                g.neighbours[id].push_back(g.index(static_cast<std::size_t>(ii),static_cast<std::size_t>(jj),static_cast<std::size_t>(kk)));
        }
    }
    return g;
}

inline std::vector<double> exact_reference(const WallSurface& s,const WallDistanceGrid& g) {
    std::vector<double> d(g.points.size(),std::numeric_limits<double>::infinity());
    for(std::size_t i=0;i<g.points.size();++i) if(!g.solid[i]) d[i]=exact_point_distance(s,g.points[i]);
    return d;
}

inline std::vector<double> search_based_reference(const WallSurface& s,const WallDistanceGrid& g,double threshold) {
    // Exact surface search accelerated by a BVH. The threshold argument is
    // retained for API compatibility, but it must never switch the method to
    // a vertex-distance approximation: the wall-distance contract is exact
    // point-to-triangle distance for every fluid query point.
    (void)threshold;
    const WallDistanceBvh bvh(s);
    std::vector<double> d(g.points.size(),std::numeric_limits<double>::infinity());
    for(std::size_t i=0;i<g.points.size();++i)
        if(!g.solid[i]) d[i]=bvh.nearest_distance(g.points[i]);
    return d;
}

inline std::vector<std::size_t> wall_seed_nodes(const WallDistanceBvh& bvh,const WallDistanceGrid& g,double seed_width) {
    std::vector<std::size_t> seeds;
    for(std::size_t i=0;i<g.points.size();++i)
        if(!g.solid[i] && std::isfinite(bvh.nearest_distance_within(g.points[i],seed_width)))
            seeds.push_back(i);
    return seeds;
}

inline std::vector<std::size_t> wall_seed_nodes(const WallSurface& s,const WallDistanceGrid& g,double seed_width) {
    const WallDistanceBvh bvh(s);
    return wall_seed_nodes(bvh,g,seed_width);
}

inline WallDistanceResult graph_wave(const WallSurface& s,const WallDistanceGrid& g,bool directional) {
    const double h=std::min({g.spacing.x,g.spacing.y,g.spacing.z});
    const WallDistanceBvh bvh(s);
    const auto seeds=wall_seed_nodes(bvh,g,1.6*h);
    std::vector<double> d(g.points.size(),std::numeric_limits<double>::infinity());
    using Item=std::pair<double,std::size_t>;
    std::priority_queue<Item,std::vector<Item>,std::greater<Item>> q;
    for(const auto id:seeds) { d[id]=bvh.nearest_distance(g.points[id]); q.push({d[id],id}); }
    while(!q.empty()) {
        const auto [du,u]=q.top(); q.pop();
        if(du!=d[u]) continue;
        for(const auto v:g.neighbours[u]) {
            if(g.solid[v]) continue;
            const auto dv=g.points[v]-g.points[u];
            const double r=wd_norm(dv);
            double w=r;
            if(directional) {
                // Penalize propagation tangential to the local wall normal.
                // This is a real directional metric; the previous implementation
                // used dot(dv,-dv), which is identically -|dv|^2.
                const auto n=bvh.nearest_normal(g.points[u]);
                const double align=std::abs(wd_dot(dv,n))/std::max(r,1e-30);
                w*=1.0+0.15*(1.0-align);
            }
            if(du+w<d[v]) { d[v]=du+w; q.push({d[v],v}); }
        }
    }
    WallDistanceResult r;
    r.distance=std::move(d);
    r.valid.assign(g.points.size(),1);
    r.method=directional?"directional_mesh_wave":"mesh_wave";
    r.iterations=0;
    // F-W1: attribute the degenerate-normal fallbacks of the directional
    // metric (the plain mesh wave never queries the surface normal).
    r.nearest_normal_fallbacks=bvh.degenerate_normal_fallbacks();
    return r;
}

inline double poisson_wall_offset(const WallDistanceBvh& bvh,
                                  const WallDistanceVec3& p,
                                  const WallDistanceVec3& q,
                                  double h);

// Per-axis 1D Shortley-Weller stencil (Shortley & Weller 1938).
//
// For a fluid node P with a neighbour on each side of the axis, each
// neighbour is one of:
//   - absent (outside the grid): homogeneous Neumann, drop that side
//     entirely (matches the pre-existing outer-boundary treatment);
//   - fluid: distance = grid spacing h, value = f[q];
//   - solid (wall): distance = delta (BVH ray, clamped to
//     wall_offset_relative_floor*h), value = phi_wall = 0.
//
// When both sides are present, the true 2nd-order irregular-boundary
// stencil is
//     f''_P = 2/(h_-*(h_-+h_+)) * f_-
//           + 2/(h_+*(h_-+h_+)) * f_+
//           - 2/(h_-*h_+)       * f_P
// which is exact for a quadratic and reduces to (f_- + f_+ - 2 f_P)/h^2
// when h_- = h_+ = h.
//
// The previous per-face stencil (-f_P/(h*delta) on a solid face plus
// (f_q - f_P)/h^2 on the opposite fluid face) is inconsistent when
// delta != h: its relative error on the coefficient is (1 - delta/h)/2,
// independent of h. This produced ~97% L2 error on the complex wing-body
// benchmark. Tracked in #466.
//
// When only one side is present (grid edge), fall back to the one-sided
// per-face contribution so the outer-boundary treatment is unchanged.
inline double poisson_wall_offset(const WallDistanceBvh& bvh,
                                      const WallDistanceVec3& p,
                                      const WallDistanceVec3& q,
                                      double h);

namespace poisson_stencil_detail {

struct AxisNeighbour {
    bool present{false};
    double distance{0.0};  // h for fluid, delta for solid, ignored if !present
    double value{0.0};     // f[q] for fluid, 0 for solid (wall Dirichlet)
};


inline AxisNeighbour classify_neighbour(const WallDistanceBvh& bvh,
                                        const std::vector<double>& f,
                                        const WallDistanceGrid& g,
                                        const WallDistanceVec3& p,
                                        std::size_t q, bool exists, double h) {
    AxisNeighbour out;
    if (!exists) return out;
    if (g.solid[q]) {
        out.present = true;
        out.distance = poisson_wall_offset(bvh, p, g.points[q], h);
        out.value = 0.0;
    } else if (std::isfinite(f[q])) {
        out.present = true;
        out.distance = h;
        out.value = f[q];
    }
    return out;
}

// Contribution of one axis to the Laplacian at f[id].
inline double axis_laplacian(const AxisNeighbour& m,
                             const AxisNeighbour& p,
                             double h, double f_id) {
    if (m.present && p.present) {
        // Shortley-Weller 2nd-order irregular stencil.
        const double sum = m.distance + p.distance;
        return 2.0 * (m.value / (m.distance * sum)
                      + p.value / (p.distance * sum)
                      - f_id  / (m.distance * p.distance));
    }
    if (m.present) return (m.value - f_id) / (h * m.distance);
    if (p.present) return (p.value - f_id) / (h * p.distance);
    return 0.0;
}

// Contribution of one axis to the diagonal of A = -Laplacian.
inline double axis_diagonal(const AxisNeighbour& m,
                            const AxisNeighbour& p,
                            double h) {
    if (m.present && p.present) return 2.0 / (m.distance * p.distance);
    if (m.present) return 1.0 / (h * m.distance);
    if (p.present) return 1.0 / (h * p.distance);
    return 0.0;
}

// 2nd-order gradient at f_P for asymmetric neighbour distances.
//
// Derived by eliminating f''(P) from a two-point Taylor expansion at
// distances h_- and h_+ (each either the grid spacing h for a fluid
// neighbour or the wall offset delta for a solid neighbour):
//
//   f'(P) = [(h_+^2 - h_-^2)*f_P - h_+^2*f_- + h_-^2*f_+]
//         / [h_- * h_+ * (h_- + h_+)]
//
// Reduces to the standard centred difference (f_+ - f_-)/(2h) when
// h_- = h_+ = h, and to the correct O(h^2) one-sided cut-face formula
// when one side is a wall at distance delta. This is the gradient
// analogue of axis_laplacian above and must be used in the distance
// reconstruction d = sqrt(|grad phi|^2 + 2 phi) - |grad phi|, because
// that formula amplifies gradient error 1:1 into the reconstructed
// distance near the wall.
//
// The previous per-face gradient used -f_P/delta one-sided at cut
// cells, which is O(1) accurate as delta shrinks and produced ~713%
// near-wall error on the complex wing-body benchmark. Tracked in #466.
inline double axis_gradient(const AxisNeighbour& m,
                            const AxisNeighbour& p,
                            double h, double f_id) {
    if (m.present && p.present) {
        const double hm = m.distance;
        const double hp = p.distance;
        const double sum = hm + hp;
        return ((hp*hp - hm*hm) * f_id
                - hp*hp * m.value
                + hm*hm * p.value)
             / (hm * hp * sum);
    }
    if (m.present) return (f_id - m.value) / m.distance;
    if (p.present) return (p.value - f_id) / p.distance;
    (void)h;
    return 0.0;
}

} // namespace poisson_stencil_detail

// Immutable per-grid cache for cut-face offsets used by iterative transport
// stencils. Geometry, solid mask and spacing are fixed during a solve, so
// these values need not be recomputed by BVH ray queries on every iteration.
class WallDistanceCutFaceOffsetCache {
public:
    WallDistanceCutFaceOffsetCache(const WallDistanceBvh& bvh,
                                   const WallDistanceGrid& g)
        : node_offset_index_(g.points.size(), invalid_index),
          spacing_{{g.spacing.x,g.spacing.y,g.spacing.z}} {
        for(std::size_t id=0; id<g.points.size(); ++id) {
            if(g.solid[id]) continue;
            const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny);
            const std::size_t j=rem/g.nx, i=rem%g.nx;
            const auto fill=[&](std::size_t slot,std::size_t q,bool exists,double h) {
                if(!exists || !g.solid[q]) return;
                if(node_offset_index_[id]==invalid_index) {
                    node_offset_index_[id]=offsets_.size();
                    offsets_.push_back({{g.spacing.x,g.spacing.x,g.spacing.y,
                                         g.spacing.y,g.spacing.z,g.spacing.z}});
                }
                offsets_[node_offset_index_[id]][slot]=
                    poisson_wall_offset(bvh,g.points[id],g.points[q],h);
            };
            fill(0,i>0?g.index(i-1,j,k):0,i>0,g.spacing.x);
            fill(1,i+1<g.nx?g.index(i+1,j,k):0,i+1<g.nx,g.spacing.x);
            fill(2,j>0?g.index(i,j-1,k):0,j>0,g.spacing.y);
            fill(3,j+1<g.ny?g.index(i,j+1,k):0,j+1<g.ny,g.spacing.y);
            fill(4,k>0?g.index(i,j,k-1):0,k>0,g.spacing.z);
            fill(5,k+1<g.nz?g.index(i,j,k+1):0,k+1<g.nz,g.spacing.z);
        }
    }

    double delta(std::size_t id,std::size_t face) const {
        const std::size_t index=node_offset_index_[id];
        return index==invalid_index ? spacing_[face/2] : offsets_[index][face];
    }

    std::size_t cut_node_count() const { return offsets_.size(); }

    std::size_t storage_bytes() const {
        return node_offset_index_.capacity()*sizeof(std::size_t)
             + offsets_.capacity()*sizeof(std::array<double,6>);
    }

private:
    static constexpr std::size_t invalid_index=std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> node_offset_index_;
    std::vector<std::array<double,6>> offsets_;
    std::array<double,3> spacing_;
};





inline double eikonal_update(const std::array<double,3>& a,
                                  const std::array<double,3>& h,
                                  double rhs) {
    if(rhs<=0.0) return 0.0;
    std::array<std::pair<double,double>,3> v{};
    std::size_t n=0;
    for(std::size_t q=0;q<3;++q) if(std::isfinite(a[q]) && h[q]>0.0) v[n++]={a[q],h[q]};
    if(n==0) return std::numeric_limits<double>::infinity();
    // n<=3 by construction.  Use a bounded insertion sort instead of
    // std::sort on a fixed-size array; GCC otherwise emits a spurious
    // -Warray-bounds warning after aggressive inlining of the 3-element case.
    for(std::size_t q=1;q<n;++q) {
        auto key=v[q];
        std::size_t j=q;
        while(j>0 && key.first<v[j-1].first) {
            v[j]=v[j-1];
            --j;
        }
        v[j]=key;
    }
    double A=0.0,B=0.0,C=-rhs*rhs;
    double x=v[0].first+v[0].second*rhs;
    for(std::size_t m=0;m<n;++m) {
        const double ai=v[m].first, hi=v[m].second;
        A+=1.0/(hi*hi); B+=-2.0*ai/(hi*hi); C+=ai*ai/(hi*hi);
        const double disc=std::max(0.0,B*B-4.0*A*C);
        x=(-B+std::sqrt(disc))/(2.0*A);
        if(m+1==n || x<=v[m+1].first) return x;
    }
    return x;
}

inline std::vector<double> eikonal_fast_sweep(const WallDistanceBvh& bvh,const WallDistanceGrid& g,
                                               std::size_t max_iter,std::size_t* used_iter=nullptr) {
    const double h=std::min({g.spacing.x,g.spacing.y,g.spacing.z});
    const auto seeds=wall_seed_nodes(bvh,g,1.6*h);
    std::vector<unsigned char> fixed(g.points.size(),0);
    std::vector<double> d(g.points.size(),std::numeric_limits<double>::infinity());
    for(auto id:seeds) { d[id]=bvh.nearest_distance(g.points[id]); fixed[id]=1; }
    const std::array<int,2> signs={-1,1};
    std::size_t it_used=max_iter;
    for(std::size_t it=0;it<max_iter;++it) {
        double max_change=0.0;
        for(int sx:signs) for(int sy:signs) for(int sz:signs)
            for(std::size_t kk=0;kk<g.nz;++kk) for(std::size_t jj=0;jj<g.ny;++jj) for(std::size_t ii=0;ii<g.nx;++ii) {
                const std::size_t i=sx>0?ii:g.nx-1-ii, j=sy>0?jj:g.ny-1-jj, k=sz>0?kk:g.nz-1-kk;
                const std::size_t id=g.index(i,j,k);
                if(g.solid[id] || fixed[id]) continue;
                std::array<double,3> a={std::numeric_limits<double>::infinity(),
                                        std::numeric_limits<double>::infinity(),
                                        std::numeric_limits<double>::infinity()};
                if(i>0) a[0]=std::min(a[0],d[g.index(i-1,j,k)]);
                if(i+1<g.nx) a[0]=std::min(a[0],d[g.index(i+1,j,k)]);
                if(j>0) a[1]=std::min(a[1],d[g.index(i,j-1,k)]);
                if(j+1<g.ny) a[1]=std::min(a[1],d[g.index(i,j+1,k)]);
                if(k>0) a[2]=std::min(a[2],d[g.index(i,j,k-1)]);
                if(k+1<g.nz) a[2]=std::min(a[2],d[g.index(i,j,k+1)]);
                const double cand=eikonal_update(a,{g.spacing.x,g.spacing.y,g.spacing.z},1.0);
                if(std::isfinite(cand) && (!std::isfinite(d[id]) || cand<d[id])) {
                    const double old=d[id]; d[id]=cand;
                    if(std::isfinite(old)) max_change=std::max(max_change,std::abs(cand-old));
                    else max_change=std::max(max_change,cand);
                }
            }
        if(max_change<1e-10*h) { it_used=it+1; break; }
    }
    if(used_iter) *used_iter=it_used;
    return d;
}

inline std::vector<double> eikonal_fast_sweep(const WallSurface& s,const WallDistanceGrid& g,
                                               std::size_t max_iter,std::size_t* used_iter=nullptr) {
    const WallDistanceBvh bvh(s);
    return eikonal_fast_sweep(bvh,g,max_iter,used_iter);
}

inline double laplacian_at(const std::vector<double>& f,const WallDistanceGrid& g,
                           std::size_t id) {
    // Finite-volume-compatible seven-point operator on the active fluid graph.
    //
    // Every fluid-fluid or fluid-solid face contributes exactly once with
    // weight 1/h^2.  A missing outer neighbour is a zero-flux Neumann face and
    // therefore contributes nothing.  This is deliberately different from
    // the old mirrored-ghost stencil (2*(f_nb-f_i)): that stencil makes the
    // matrix non-symmetric at an outer boundary (boundary row -2/h^2 versus
    // interior row -1/h^2), so it is not a valid Euclidean-inner-product PCG
    // operator.
    const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
    double l=0.0;
    auto add_face=[&](std::size_t q,bool exists,double h) {
        if(!exists) return; // homogeneous Neumann: zero flux
        const double w=1.0/(h*h);
        const double fq=g.solid[q] ? 0.0 : f[q];
        if(std::isfinite(fq)) l+=(fq-f[id])*w;
    };
    const std::size_t xm=i>0?g.index(i-1,j,k):0, xp=i+1<g.nx?g.index(i+1,j,k):0;
    const std::size_t ym=j>0?g.index(i,j-1,k):0, yp=j+1<g.ny?g.index(i,j+1,k):0;
    const std::size_t zm=k>0?g.index(i,j,k-1):0, zp=k+1<g.nz?g.index(i,j,k+1):0;
    add_face(xm,i>0,g.spacing.x); add_face(xp,i+1<g.nx,g.spacing.x);
    add_face(ym,j>0,g.spacing.y); add_face(yp,j+1<g.ny,g.spacing.y);
    add_face(zm,k>0,g.spacing.z); add_face(zp,k+1<g.nz,g.spacing.z);
    return l;
}

// BVH-aware laplacian used by the Hamilton-Jacobi and hybrid pseudo-time
// integrators. Wall neighbours contribute a Dirichlet zero-value flux at
// distance delta rather than at grid distance h. Shares the same
// Shortley-Weller irregular stencil as poisson_laplacian_at for
// consistency between the elliptic solve and the pseudo-time relaxation.
inline double laplacian_at(const WallDistanceBvh& bvh,
                           const std::vector<double>& f,
                           const WallDistanceGrid& g,
                           std::size_t id) {
    using poisson_stencil_detail::classify_neighbour;
    using poisson_stencil_detail::axis_laplacian;
    const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
    const auto& p=g.points[id];
    const std::size_t xm=i>0?g.index(i-1,j,k):0, xp=i+1<g.nx?g.index(i+1,j,k):0;
    const std::size_t ym=j>0?g.index(i,j-1,k):0, yp=j+1<g.ny?g.index(i,j+1,k):0;
    const std::size_t zm=k>0?g.index(i,j,k-1):0, zp=k+1<g.nz?g.index(i,j,k+1):0;
    return axis_laplacian(classify_neighbour(bvh,f,g,p,xm,i>0,g.spacing.x),
                          classify_neighbour(bvh,f,g,p,xp,i+1<g.nx,g.spacing.x),
                          g.spacing.x,f[id])
         + axis_laplacian(classify_neighbour(bvh,f,g,p,ym,j>0,g.spacing.y),
                          classify_neighbour(bvh,f,g,p,yp,j+1<g.ny,g.spacing.y),
                          g.spacing.y,f[id])
         + axis_laplacian(classify_neighbour(bvh,f,g,p,zm,k>0,g.spacing.z),
                          classify_neighbour(bvh,f,g,p,zp,k+1<g.nz,g.spacing.z),
                          g.spacing.z,f[id]);
}

inline double godunov_gradient_at(const WallDistanceBvh& bvh,
                                  const std::vector<double>& d,const WallDistanceGrid& g,
                                  std::size_t id,
                                  const WallDistanceCutFaceOffsetCache* offsets=nullptr) {
    const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
    auto component=[&](int axis)->double {
        double dm=0.0, dp=0.0;
        auto face=[&](std::size_t q,bool minus_side,double hh,std::size_t face_index) {
            if(g.solid[q]) {
                const double delta=offsets ? offsets->delta(id,face_index)
                                           : poisson_wall_offset(bvh,g.points[id],g.points[q],hh);
                if(delta>0.0 && std::isfinite(delta)) {
                    if(minus_side) dm=d[id]/delta;
                    else dp=-d[id]/delta;
                }
            } else if(std::isfinite(d[q])) {
                if(minus_side) dm=(d[id]-d[q])/hh;
                else dp=(d[q]-d[id])/hh;
            }
        };
        if(axis==0) {
            if(i>0) face(g.index(i-1,j,k),true,g.spacing.x,0);
            if(i+1<g.nx) face(g.index(i+1,j,k),false,g.spacing.x,1);
        } else if(axis==1) {
            if(j>0) face(g.index(i,j-1,k),true,g.spacing.y,2);
            if(j+1<g.ny) face(g.index(i,j+1,k),false,g.spacing.y,3);
        } else {
            if(k>0) face(g.index(i,j,k-1),true,g.spacing.z,4);
            if(k+1<g.nz) face(g.index(i,j,k+1),false,g.spacing.z,5);
        }
        // The upwind Eikonal Hamiltonian chooses one direction per axis.
        // Adding both squared slopes double-counts a medial ridge and
        // disagrees with eikonal_update(), which uses the nearest neighbour.
        return std::max({dm,-dp,0.0});
    };
    const double gx=component(0), gy=component(1), gz=component(2);
    return std::sqrt(gx*gx+gy*gy+gz*gz);
}

inline double godunov_gradient_at(const std::vector<double>& d,const WallDistanceGrid& g,
                                  std::size_t id) {
    const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
    auto component=[&](int axis)->double {
        double dm=0.0, dp=0.0;
        if(axis==0) {
            if(i>0) {
                const auto m=g.index(i-1,j,k);
                if(g.solid[m]) dm=d[id]/g.spacing.x;
                else if(std::isfinite(d[m])) dm=(d[id]-d[m])/g.spacing.x;
            }
            if(i+1<g.nx) {
                const auto p=g.index(i+1,j,k);
                if(g.solid[p]) dp=-d[id]/g.spacing.x;
                else if(std::isfinite(d[p])) dp=(d[p]-d[id])/g.spacing.x;
            }
        } else if(axis==1) {
            if(j>0) {
                const auto m=g.index(i,j-1,k);
                if(g.solid[m]) dm=d[id]/g.spacing.y;
                else if(std::isfinite(d[m])) dm=(d[id]-d[m])/g.spacing.y;
            }
            if(j+1<g.ny) {
                const auto p=g.index(i,j+1,k);
                if(g.solid[p]) dp=-d[id]/g.spacing.y;
                else if(std::isfinite(d[p])) dp=(d[p]-d[id])/g.spacing.y;
            }
        } else {
            if(k>0) {
                const auto m=g.index(i,j,k-1);
                if(g.solid[m]) dm=d[id]/g.spacing.z;
                else if(std::isfinite(d[m])) dm=(d[id]-d[m])/g.spacing.z;
            }
            if(k+1<g.nz) {
                const auto p=g.index(i,j,k+1);
                if(g.solid[p]) dp=-d[id]/g.spacing.z;
                else if(std::isfinite(d[p])) dp=(d[p]-d[id])/g.spacing.z;
            }
        }
        // The upwind Eikonal Hamiltonian chooses one direction per axis.
        // Adding both squared slopes double-counts a medial ridge and
        // disagrees with eikonal_update(), which uses the nearest neighbour.
        return std::max({dm,-dp,0.0});
    };
    const double gx=component(0), gy=component(1), gz=component(2);
    return std::sqrt(gx*gx+gy*gy+gz*gz);
}


inline double laplacian_at(const WallDistanceBvh& bvh,
                           const std::vector<double>& f,
                           const WallDistanceGrid& g,
                           std::size_t id,
                           const WallDistanceCutFaceOffsetCache& offsets) {
    using poisson_stencil_detail::AxisNeighbour;
    using poisson_stencil_detail::axis_laplacian;
    const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
    const auto side=[&](std::size_t q,bool exists,double h,std::size_t face) {
        AxisNeighbour out;
        if(!exists) return out;
        if(g.solid[q]) {
            out.present=true;
            out.distance=offsets.delta(id,face);
            out.value=0.0;
        } else if(std::isfinite(f[q])) {
            out.present=true;
            out.distance=h;
            out.value=f[q];
        }
        return out;
    };
    (void)bvh;
    const auto xm=i>0?g.index(i-1,j,k):0, xp=i+1<g.nx?g.index(i+1,j,k):0;
    const auto ym=j>0?g.index(i,j-1,k):0, yp=j+1<g.ny?g.index(i,j+1,k):0;
    const auto zm=k>0?g.index(i,j,k-1):0, zp=k+1<g.nz?g.index(i,j,k+1):0;
    return axis_laplacian(side(xm,i>0,g.spacing.x,0),side(xp,i+1<g.nx,g.spacing.x,1),g.spacing.x,f[id])
         + axis_laplacian(side(ym,j>0,g.spacing.y,2),side(yp,j+1<g.ny,g.spacing.y,3),g.spacing.y,f[id])
         + axis_laplacian(side(zm,k>0,g.spacing.z,4),side(zp,k+1<g.nz,g.spacing.z,5),g.spacing.z,f[id]);
}

inline std::vector<double> hamilton_jacobi_distance(const WallSurface& s,const WallDistanceGrid& g,
                                                     std::size_t max_iter,double epsilon=0.25,
                                                     double relaxation=0.25,
                                                     std::size_t* used_iter=nullptr,
                                                     double* residual_out=nullptr,
                                                     std::size_t* clamp_activations=nullptr) {
    // Modified Hamilton-Jacobi equation:
    //   |grad d| = 1 + Gamma(d) laplacian(d), Gamma(d)=epsilon*d.
    // Solve it by pseudo-time relaxation with an upwind/Godunov gradient
    // using a CFL-scale explicit step for the advection/diffusion update.
    // and a central Laplacian. The local pseudo-time step is diffusion-limited.
    const double h=std::min({g.spacing.x,g.spacing.y,g.spacing.z});
    const WallDistanceBvh bvh(s);
    const WallDistanceCutFaceOffsetCache offsets(bvh,g);
    const auto seeds=wall_seed_nodes(bvh,g,1.6*h);
    std::vector<unsigned char> fixed(g.points.size(),0);
    // eikonal_fast_sweep uses the same seed criterion and already assigns exact
    // BVH distances to those nodes. Do not repeat the nearest-triangle queries.
    std::vector<double> d=eikonal_fast_sweep(bvh,g,std::max<std::size_t>(8,max_iter));
    for(auto id:seeds) fixed[id]=1;
    // Preserve the original ascending in-place update order, but avoid scanning
    // solid, fixed-seed and unreachable nodes on every pseudo-time iteration.
    std::vector<std::size_t> active_ids;
    active_ids.reserve(d.size());
    for(std::size_t id=0;id<d.size();++id)
        if(!g.solid[id] && !fixed[id] && std::isfinite(d[id])) active_ids.push_back(id);

    const double inverse_spacing_sum=1.0/g.spacing.x+1.0/g.spacing.y+1.0/g.spacing.z;
    const double dt_adv_constant=0.9/inverse_spacing_sum;
    const double relaxation_factor=std::clamp(relaxation,0.1,1.0);
    const std::size_t steps=std::max<std::size_t>(20,max_iter);
    std::size_t used=steps;
    std::size_t clamps=0;
    double final_residual=std::numeric_limits<double>::infinity();
    const double inv_h2=1.0/(g.spacing.x*g.spacing.x)
                      +1.0/(g.spacing.y*g.spacing.y)
                      +1.0/(g.spacing.z*g.spacing.z);
    for(std::size_t it=0;it<steps;++it) {
        double max_change=0.0;
        for(const std::size_t id:active_ids) {
            const double grad=godunov_gradient_at(bvh,d,g,id,&offsets);
            const double gamma=epsilon*std::max(d[id],0.0);
            const double lap=laplacian_at(bvh,d,g,id,offsets);
            const double residual=grad-1.0-gamma*lap;
            // Use the full explicit CFL limits of the upwind advection and
            // central diffusion terms. The previous 0.25 safety factor was
            // unnecessarily restrictive and made H-J convergence ~4x slower.
            // 3-D CFL: dt*h^{-1} must account for all active
            // upwind directions.  0.9*h is only a 1-D bound and can violate
            // the monotone explicit stability limit when several components
            // of the Godunov gradient are active.
            const double dt_diff=gamma>0.0 ? 0.9/(2.0*gamma*inv_h2)
                                           : std::numeric_limits<double>::infinity();
            const double dt=relaxation_factor*std::min(dt_adv_constant,dt_diff);
            const double candidate=d[id]-dt*residual;
            if(candidate<0.0) ++clamps;  // F-W2: observable positivity clamp
            const double nd=std::max(0.0,candidate);
            max_change=std::max(max_change,std::abs(nd-d[id]));
            d[id]=nd;
        }
        // Pure Poisson mode deliberately does not overwrite the PDE solution with BVH seeds.
        if((it&3u)==3u || max_change<1e-10*h) {
            final_residual=0.0;
            for(const std::size_t id:active_ids) {
                final_residual=std::max(final_residual,std::abs(
                    godunov_gradient_at(bvh,d,g,id,&offsets)-1.0-
                    epsilon*std::max(d[id],0.0)*laplacian_at(bvh,d,g,id,offsets)));
            }
            if(final_residual<1e-8 && max_change<1e-9*h) { used=it+1; break; }
        }
    }
    if(used_iter) *used_iter=used;
    if(residual_out) *residual_out=final_residual;
    if(clamp_activations) *clamp_activations=clamps;
    return d;
}

inline double poisson_laplacian_at(const WallDistanceBvh& bvh,
                                   const std::vector<double>& f,
                                   const WallDistanceGrid& g,
                                   std::size_t id);

inline double poisson_residual_inf(const WallDistanceBvh& bvh,
                                      const std::vector<double>& phi,
                                      const WallDistanceGrid& g,
                                      const std::vector<unsigned char>& fixed) {
    // Public diagnostic: use exactly the same cut-face operator as PCG.
    double rmax=0.0;
    for(std::size_t id=0;id<phi.size();++id) {
        if(g.solid[id] || fixed[id] || !std::isfinite(phi[id])) continue;
        rmax=std::max(rmax,std::abs(-poisson_laplacian_at(bvh,phi,g,id)-1.0));
    }
    return rmax;
}

inline double poisson_wall_offset(const WallDistanceBvh& bvh,
                                      const WallDistanceVec3& p,
                                      const WallDistanceVec3& q,
                                      double h);
inline double poisson_laplacian_at(const WallDistanceBvh& bvh,
                                   const std::vector<double>& f,
                                   const WallDistanceGrid& g,
                                   std::size_t id);
inline double poisson_diagonal(const WallDistanceBvh& bvh,
                               const WallDistanceGrid& g,
                               std::size_t id);

struct WallDistancePoissonAudit {
    std::size_t fluid_nodes{0};
    std::size_t fluid_fluid_faces{0};
    std::size_t solid_faces{0};
    std::size_t outer_faces{0};
    double min_diagonal{std::numeric_limits<double>::infinity()};
    double max_diagonal{0.0};
    double min_diagonal_dominance{std::numeric_limits<double>::infinity()};
    double symmetry_error{0.0};
};

struct WallDistancePoissonOffsetAudit {
    double min_delta{std::numeric_limits<double>::infinity()};
    double max_delta{0.0};
    double min_delta_over_h{std::numeric_limits<double>::infinity()};
    double max_delta_over_h{0.0};
    double max_wall_coefficient{0.0};
    std::size_t max_wall_coefficient_cell{0};
    double max_wall_coefficient_delta{0.0};
    double max_wall_coefficient_h{0.0};
    WallDistanceVec3 max_wall_coefficient_point{};
    std::size_t degenerate_count{0};
    std::size_t solid_face_count{0};
    std::size_t ray_hit_count{0};
    std::size_t ray_miss_count{0};
    std::size_t fallback_count{0};
    std::size_t bad_alignment_count{0};
    double min_alignment{1.0};
};

struct WallDistancePoissonWallOffsetDiagnostic {
    double delta{0.0};
    bool ray_hit{false};
    bool fallback{false};
    bool bad_alignment{false};
    double alignment{1.0};
};

inline WallDistancePoissonWallOffsetDiagnostic poisson_wall_offset_diagnostic(
    const WallDistanceBvh& bvh,const WallDistanceVec3& p,
    const WallDistanceVec3& q,double h);

inline WallDistancePoissonOffsetAudit audit_poisson_wall_offsets(
    const WallDistanceBvh& bvh, const WallDistanceGrid& g) {
    WallDistancePoissonOffsetAudit a;
    for(std::size_t id=0; id<g.points.size(); ++id) {
        if(g.solid[id]) continue;
        const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
        auto inspect=[&](std::size_t q, bool exists, double h) {
            if(!exists || !g.solid[q]) return;
            ++a.solid_face_count;
            const auto wd=poisson_wall_offset_diagnostic(bvh,g.points[id],g.points[q],h);
            const double delta=wd.delta;
            if(wd.ray_hit) ++a.ray_hit_count;
            else ++a.ray_miss_count;
            if(wd.fallback) ++a.fallback_count;
            if(wd.bad_alignment) ++a.bad_alignment_count;
            a.min_alignment=std::min(a.min_alignment,wd.alignment);
            if(!(std::isfinite(delta) && delta>0.0)) {
                ++a.degenerate_count;
                return;
            }
            const double ratio=delta/h;
            const double coeff=1.0/(h*delta);
            a.min_delta=std::min(a.min_delta,delta);
            a.max_delta=std::max(a.max_delta,delta);
            a.min_delta_over_h=std::min(a.min_delta_over_h,ratio);
            a.max_delta_over_h=std::max(a.max_delta_over_h,ratio);
            if(coeff>a.max_wall_coefficient) {
                a.max_wall_coefficient=coeff;
                a.max_wall_coefficient_cell=id;
                a.max_wall_coefficient_delta=delta;
                a.max_wall_coefficient_h=h;
                a.max_wall_coefficient_point=g.points[id];
            }
        };
        inspect(i>0?g.index(i-1,j,k):0,i>0,g.spacing.x);
        inspect(i+1<g.nx?g.index(i+1,j,k):0,i+1<g.nx,g.spacing.x);
        inspect(j>0?g.index(i,j-1,k):0,j>0,g.spacing.y);
        inspect(j+1<g.ny?g.index(i,j+1,k):0,j+1<g.ny,g.spacing.y);
        inspect(k>0?g.index(i,j,k-1):0,k>0,g.spacing.z);
        inspect(k+1<g.nz?g.index(i,j,k+1):0,k+1<g.nz,g.spacing.z);
    }
    if(!std::isfinite(a.min_delta)) a.min_delta=0.0;
    if(!std::isfinite(a.min_delta_over_h)) a.min_delta_over_h=0.0;
    return a;
}

inline WallDistancePoissonAudit audit_poisson_operator(const WallDistanceBvh& bvh,
                                                   const WallDistanceGrid& g) {
    WallDistancePoissonAudit a;
    const std::size_t n=g.points.size();

    auto offdiag=[&](std::size_t row,std::size_t col)->double {
        if(g.solid[row] || g.solid[col]) return 0.0;
        const std::size_t rk=row/(g.nx*g.ny), rr=row%(g.nx*g.ny), rj=rr/g.nx, ri=rr%g.nx;
        const std::size_t ck=col/(g.nx*g.ny), cr=col%(g.nx*g.ny), cj=cr/g.nx, ci=cr%g.nx;
        if(rj==cj && rk==ck && (ri+1==ci || ci+1==ri))
            return -1.0/(g.spacing.x*g.spacing.x);
        if(ri==ci && rk==ck && (rj+1==cj || cj+1==rj))
            return -1.0/(g.spacing.y*g.spacing.y);
        if(ri==ci && rj==cj && (rk+1==ck || ck+1==rk))
            return -1.0/(g.spacing.z*g.spacing.z);
        return 0.0;
    };

    for(std::size_t id=0;id<n;++id) {
        if(g.solid[id]) continue;
        ++a.fluid_nodes;
        const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
        const double diag=poisson_diagonal(bvh,g,id);
        a.min_diagonal=std::min(a.min_diagonal,diag);
        a.max_diagonal=std::max(a.max_diagonal,diag);

        double abs_offdiag_sum=0.0;
        auto inspect=[&](std::size_t q,bool exists) {
            if(!exists) {
                ++a.outer_faces;
                return;
            }
            if(g.solid[q]) {
                ++a.solid_faces;
                return;
            }
            if(q>id) ++a.fluid_fluid_faces;
            abs_offdiag_sum += std::abs(offdiag(id,q));
            a.symmetry_error=std::max(
                a.symmetry_error,
                std::abs(offdiag(id,q)-offdiag(q,id)));
        };

        inspect(i>0?g.index(i-1,j,k):0,i>0);
        inspect(i+1<g.nx?g.index(i+1,j,k):0,i+1<g.nx);
        inspect(j>0?g.index(i,j-1,k):0,j>0);
        inspect(j+1<g.ny?g.index(i,j+1,k):0,j+1<g.ny);
        inspect(k>0?g.index(i,j,k-1):0,k>0);
        inspect(k+1<g.nz?g.index(i,j,k+1):0,k+1<g.nz);

        // For the actual operator A=-L:
        //   A_ii = sum(fluid-face weights) + sum(solid-face weights)
        //   A_ij = -fluid-face weight.
        // Therefore the diagonal dominance margin is exactly the total
        // Dirichlet contribution. This uses poisson_diagonal(), including
        // the BVH-based cut-face/wall-offset coefficient actually used by PCG.
        a.min_diagonal_dominance=std::min(
            a.min_diagonal_dominance, diag-abs_offdiag_sum);
    }

    if(!std::isfinite(a.min_diagonal)) a.min_diagonal=0.0;
    if(!std::isfinite(a.min_diagonal_dominance)) a.min_diagonal_dominance=0.0;
    return a;
}

// Minimum wall offset as a fraction of the local grid spacing.
//
// The cut-face Dirichlet contribution scales as 1/(h*delta). If delta is
// allowed to approach zero (wall arbitrarily close to a grid node), the
// diagonal coefficient of the Poisson operator diverges: on the complex
// benchmark we observed max_wall_coefficient = 3.3e12 against typical
// interior coefficients of ~23, producing a condition number of ~10^11
// and a PCG that converges to the wrong solution (residual small, phi
// pinned to zero over the corrupted row).
//
// Clamping delta from below at ~5% of h keeps the diagonal bounded by
// 20/h^2 while introducing at most a 5%*h positioning error at the wall,
// which is below the discretisation error of the one-sided ghost-cell
// closure itself. This is the standard treatment in cut-cell Poisson
// solvers (Shortley-Weller with epsilon-clamp).
inline constexpr double wall_offset_relative_floor = 0.05;

inline WallDistancePoissonWallOffsetDiagnostic poisson_wall_offset_diagnostic(
    const WallDistanceBvh& bvh,const WallDistanceVec3& p,
    const WallDistanceVec3& q,double h) {
    WallDistancePoissonWallOffsetDiagnostic out;
    const double floor=wall_offset_relative_floor*h;
    const auto segment=q-p;
    const double length=wd_norm(segment);
    const double ray=bvh.ray_distance(p,segment,length);
    if(std::isfinite(ray) && ray>0.0) {
        out.delta=std::clamp(ray,floor,h);
        out.ray_hit=true;
        return out;
    }

    out.fallback=true;
    const double d=bvh.nearest_distance(p);
    if(!std::isfinite(d) || d<=0.0) {
        out.delta=h;
        out.bad_alignment=true;
        out.alignment=0.0;
        return out;
    }

    const auto dir=segment*(1.0/std::max(length,1e-30));
    const auto n=bvh.nearest_normal(p);
    out.alignment=std::abs(wd_dot(n,dir));
    if(out.alignment<0.25) {
        out.bad_alignment=true;
        out.delta=std::clamp(d,floor,h);
    } else {
        out.delta=std::clamp(d/out.alignment,floor,h);
    }
    return out;
}

inline double poisson_wall_offset(const WallDistanceBvh& bvh,
                                      const WallDistanceVec3& p,
                                      const WallDistanceVec3& q,
                                      double h) {
    return poisson_wall_offset_diagnostic(bvh,p,q,h).delta;
}


inline double poisson_laplacian_at(const WallDistanceBvh& bvh,
                                   const std::vector<double>& f,
                                   const WallDistanceGrid& g,
                                   std::size_t id) {
    using poisson_stencil_detail::classify_neighbour;
    using poisson_stencil_detail::axis_laplacian;
    const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
    const auto& p=g.points[id];
    const std::size_t xm=i>0?g.index(i-1,j,k):0, xp=i+1<g.nx?g.index(i+1,j,k):0;
    const std::size_t ym=j>0?g.index(i,j-1,k):0, yp=j+1<g.ny?g.index(i,j+1,k):0;
    const std::size_t zm=k>0?g.index(i,j,k-1):0, zp=k+1<g.nz?g.index(i,j,k+1):0;
    const auto xn=classify_neighbour(bvh,f,g,p,xm,i>0,g.spacing.x);
    const auto xP=classify_neighbour(bvh,f,g,p,xp,i+1<g.nx,g.spacing.x);
    const auto yn=classify_neighbour(bvh,f,g,p,ym,j>0,g.spacing.y);
    const auto yP=classify_neighbour(bvh,f,g,p,yp,j+1<g.ny,g.spacing.y);
    const auto zn=classify_neighbour(bvh,f,g,p,zm,k>0,g.spacing.z);
    const auto zP=classify_neighbour(bvh,f,g,p,zp,k+1<g.nz,g.spacing.z);
    return axis_laplacian(xn,xP,g.spacing.x,f[id])
         + axis_laplacian(yn,yP,g.spacing.y,f[id])
         + axis_laplacian(zn,zP,g.spacing.z,f[id]);
}


inline double poisson_laplacian_at(const WallDistanceBvh& bvh,
                                   const std::vector<double>& f,
                                   const WallDistanceGrid& g,
                                   std::size_t id,
                                   const WallDistanceCutFaceOffsetCache& offsets) {
    return laplacian_at(bvh,f,g,id,offsets);
}

inline double poisson_diagonal(const WallDistanceBvh& bvh,
                               const WallDistanceGrid& g,
                               std::size_t id) {
    using poisson_stencil_detail::AxisNeighbour;
    using poisson_stencil_detail::axis_diagonal;
    const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
    const auto& p=g.points[id];
    // The diagonal only depends on presence and distance; values are
    // irrelevant. Inline a lightweight classifier so we never touch f[q].
    auto side=[&](std::size_t q, bool exists, double h) {
        AxisNeighbour out;
        if(!exists) return out;
        out.present=true;
        out.distance=g.solid[q]
            ? poisson_wall_offset(bvh,p,g.points[q],h)
            : h;
        return out;
    };
    const std::size_t xm=i>0?g.index(i-1,j,k):0, xp=i+1<g.nx?g.index(i+1,j,k):0;
    const std::size_t ym=j>0?g.index(i,j-1,k):0, yp=j+1<g.ny?g.index(i,j+1,k):0;
    const std::size_t zm=k>0?g.index(i,j,k-1):0, zp=k+1<g.nz?g.index(i,j,k+1):0;
    return axis_diagonal(side(xm,i>0,g.spacing.x),
                         side(xp,i+1<g.nx,g.spacing.x),g.spacing.x)
         + axis_diagonal(side(ym,j>0,g.spacing.y),
                         side(yp,j+1<g.ny,g.spacing.y),g.spacing.y)
         + axis_diagonal(side(zm,k>0,g.spacing.z),
                         side(zp,k+1<g.nz,g.spacing.z),g.spacing.z);
}


inline double poisson_diagonal(const WallDistanceBvh& bvh,
                               const WallDistanceGrid& g,
                               std::size_t id,
                               const WallDistanceCutFaceOffsetCache& offsets) {
    using poisson_stencil_detail::AxisNeighbour;
    using poisson_stencil_detail::axis_diagonal;
    const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
    const auto side=[&](bool exists,std::size_t face) {
        AxisNeighbour out;
        if(!exists) return out;
        out.present=true;
        out.distance=offsets.delta(id,face);
        return out;
    };
    (void)bvh;
    return axis_diagonal(side(i>0,0),side(i+1<g.nx,1),g.spacing.x)
         + axis_diagonal(side(j>0,2),side(j+1<g.ny,3),g.spacing.y)
         + axis_diagonal(side(k>0,4),side(k+1<g.nz,5),g.spacing.z);
}


inline std::vector<double> poisson_potential(const WallDistanceBvh& bvh,const WallDistanceGrid& g,
                                               std::size_t max_iter,double,
                                               std::size_t* used_iter=nullptr,
                                               double* residual_out=nullptr,
                                               const WallDistanceCutFaceOffsetCache* offsets=nullptr) {
    // Solve the actual elliptic boundary-value problem
    //
    //     -laplacian(phi) = 1,    phi = 0 on the wall,
    //
    // with homogeneous Neumann treatment at the outer computational boundary.
    //
    // Solid neighbours are Dirichlet values phi_b=0 in the discrete operator.
    // Outer boundaries are homogeneous Neumann zero-flux faces.
    //
    // The discrete operator uses the 2nd-order Shortley-Weller irregular
    // stencil at cut cells (see poisson_laplacian_at). That stencil is
    // consistent to O(h^2) even on a wall arbitrarily close to a grid node,
    // but it is NOT symmetric when the two axis neighbours of a fluid node
    // see different distances (typical of cut cells). The previous PCG
    // solver here required an SPD operator; with the S-W stencil we use
    // BiCGStab, which handles non-symmetric operators. The system stays
    // positive-definite (a similarity transform of an SPD operator), so
    // BiCGStab remains an appropriate outer Krylov method. Diagonal Jacobi
    // preconditioning uses poisson_diagonal, which for the S-W stencil
    // gives the true 2/(hm*hp) diagonal entry.
    const std::size_t n=g.points.size();
    std::vector<double> phi(n,0.0);

    auto is_fluid=[](const WallDistanceGrid& grid,std::size_t id) {
        return !grid.solid[id];
    };

    auto diagonal=[&](std::size_t id) {
        const double d=offsets ? poisson_diagonal(bvh,g,id,*offsets)
                               : poisson_diagonal(bvh,g,id);
        return d>0.0 ? d : 1.0;
    };

    // A * x, evaluated only on fluid degrees of freedom.
    auto apply=[&](const std::vector<double>& x,std::vector<double>& y) {
        std::fill(y.begin(),y.end(),0.0);
        for(std::size_t id=0;id<n;++id) {
            if(!is_fluid(g,id)) continue;
            y[id] = offsets ? -poisson_laplacian_at(bvh,x,g,id,*offsets)
                            : -poisson_laplacian_at(bvh,x,g,id);
        }
    };

    // Diagonal Jacobi preconditioner: M^-1 z = y  <=>  z = D^-1 y.
    auto precondition=[&](const std::vector<double>& in,std::vector<double>& out) {
        std::fill(out.begin(),out.end(),0.0);
        for(std::size_t id=0;id<n;++id)
            if(is_fluid(g,id)) out[id]=in[id]/diagonal(id);
    };

    // BiCGStab state (Saad, Iterative Methods for Sparse Linear Systems,
    // Algorithm 7.7).
    std::vector<double> r(n,0.0), rhat(n,0.0), p(n,0.0), v(n,0.0);
    std::vector<double> s(n,0.0), t(n,0.0);
    std::vector<double> phat(n,0.0), shat(n,0.0);
    std::vector<double> Ax(n,0.0);

    // r0 = b - A x0 = 1 - A phi  (phi starts at zero, so r0 = 1 on fluid DOFs)
    apply(phi,Ax);
    double residual_inf=0.0;
    std::size_t fluid_count=0;
    for(std::size_t id=0;id<n;++id) {
        if(!is_fluid(g,id)) continue;
        r[id]=1.0-Ax[id];
        rhat[id]=r[id];
        residual_inf=std::max(residual_inf,std::abs(r[id]));
        ++fluid_count;
    }

    if(fluid_count==0) {
        if(used_iter) *used_iter=0;
        if(residual_out) *residual_out=0.0;
        return phi;
    }

    constexpr double abs_tol=1e-10;
    constexpr double breakdown_tol=1e-30;
    std::size_t used=0;

    if(residual_inf<=abs_tol) {
        if(used_iter) *used_iter=0;
        if(residual_out) *residual_out=residual_inf;
        return phi;
    }

    double rho_old=1.0, alpha=1.0, omega=1.0;

    for(std::size_t it=0;it<max_iter;++it) {
        double rho=0.0;
        for(std::size_t id=0;id<n;++id)
            if(is_fluid(g,id)) rho += rhat[id]*r[id];

        if(!std::isfinite(rho) || std::abs(rho)<=breakdown_tol) {
            used=it; break;
        }

        if(it==0) {
            for(std::size_t id=0;id<n;++id) if(is_fluid(g,id)) p[id]=r[id];
        } else {
            if(std::abs(rho_old)<=breakdown_tol || std::abs(omega)<=breakdown_tol) {
                used=it; break;
            }
            const double beta=(rho/rho_old)*(alpha/omega);
            for(std::size_t id=0;id<n;++id)
                if(is_fluid(g,id)) p[id]=r[id]+beta*(p[id]-omega*v[id]);
        }

        precondition(p,phat);
        apply(phat,v);

        double rhat_v=0.0;
        for(std::size_t id=0;id<n;++id)
            if(is_fluid(g,id)) rhat_v += rhat[id]*v[id];
        if(std::abs(rhat_v)<=breakdown_tol || !std::isfinite(rhat_v)) {
            used=it; break;
        }
        alpha=rho/rhat_v;

        // s = r - alpha v
        for(std::size_t id=0;id<n;++id)
            if(is_fluid(g,id)) s[id]=r[id]-alpha*v[id];

        double s_inf=0.0;
        for(std::size_t id=0;id<n;++id)
            if(is_fluid(g,id)) s_inf=std::max(s_inf,std::abs(s[id]));
        if(s_inf<=abs_tol) {
            // x <- x + alpha * phat, residual has been driven below tol
            for(std::size_t id=0;id<n;++id)
                if(is_fluid(g,id)) phi[id]+=alpha*phat[id];
            residual_inf=s_inf;
            used=it+1;
            break;
        }

        precondition(s,shat);
        apply(shat,t);

        double tt=0.0, ts=0.0;
        for(std::size_t id=0;id<n;++id) {
            if(!is_fluid(g,id)) continue;
            tt += t[id]*t[id];
            ts += t[id]*s[id];
        }
        if(!(tt>breakdown_tol) || !std::isfinite(tt)) {
            used=it; break;
        }
        omega=ts/tt;

        residual_inf=0.0;
        double max_change=0.0;
        for(std::size_t id=0;id<n;++id) {
            if(!is_fluid(g,id)) continue;
            const double delta=alpha*phat[id]+omega*shat[id];
            phi[id]+=delta;
            r[id]=s[id]-omega*t[id];
            residual_inf=std::max(residual_inf,std::abs(r[id]));
            max_change=std::max(max_change,std::abs(delta));
        }

        used=it+1;
        if(residual_inf<=abs_tol) break;
        if(std::abs(omega)<=breakdown_tol) break;

        // Stagnation is a numerical failure, not convergence.
        if(max_change==0.0 && residual_inf>abs_tol) break;

        rho_old=rho;
    }

    if(used_iter) *used_iter=used;
    if(residual_out) *residual_out=residual_inf;

    return phi;
}

struct WallDistancePoissonReconstructionAudit {
    double phi_min{std::numeric_limits<double>::infinity()};
    double phi_max{-std::numeric_limits<double>::infinity()};
    double grad_min{std::numeric_limits<double>::infinity()};
    double grad_max{0.0};
    double distance_l2_error{0.0};
    double distance_linf_error{0.0};
};


inline double poisson_laplacian_at(const WallDistanceBvh& bvh,
                                   const std::vector<double>& f,
                                   const WallDistanceGrid& g,
                                   std::size_t id);

// P0/P1 numerical-audit contract: the gradient used by the Poisson
// reconstruction is kept in one place and uses the same wall offset as the
// Poisson face flux.  This does NOT make the Poisson-to-distance model exact;
// it only removes an implementation inconsistency between the elliptic solve
// and the post-processing operator.
inline double poisson_reconstruction_gradient_component(
    const WallDistanceBvh& bvh,const std::vector<double>& phi,
    const WallDistanceGrid& g,std::size_t id,int axis) {
    using poisson_stencil_detail::classify_neighbour;
    using poisson_stencil_detail::axis_gradient;
    const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
    const auto& p=g.points[id];
    if(axis==0) {
        const std::size_t xm=i>0?g.index(i-1,j,k):0, xp=i+1<g.nx?g.index(i+1,j,k):0;
        return axis_gradient(classify_neighbour(bvh,phi,g,p,xm,i>0,g.spacing.x),
                             classify_neighbour(bvh,phi,g,p,xp,i+1<g.nx,g.spacing.x),
                             g.spacing.x,phi[id]);
    }
    if(axis==1) {
        const std::size_t ym=j>0?g.index(i,j-1,k):0, yp=j+1<g.ny?g.index(i,j+1,k):0;
        return axis_gradient(classify_neighbour(bvh,phi,g,p,ym,j>0,g.spacing.y),
                             classify_neighbour(bvh,phi,g,p,yp,j+1<g.ny,g.spacing.y),
                             g.spacing.y,phi[id]);
    }
    const std::size_t zm=k>0?g.index(i,j,k-1):0, zp=k+1<g.nz?g.index(i,j,k+1):0;
    return axis_gradient(classify_neighbour(bvh,phi,g,p,zm,k>0,g.spacing.z),
                         classify_neighbour(bvh,phi,g,p,zp,k+1<g.nz,g.spacing.z),
                         g.spacing.z,phi[id]);
}

inline double poisson_reconstructed_distance(const WallDistanceBvh& bvh,
                                              const WallDistanceGrid& g,
                                              const std::vector<double>& phi,
                                              std::size_t id) {
    const double gx=poisson_reconstruction_gradient_component(bvh,phi,g,id,0);
    const double gy=poisson_reconstruction_gradient_component(bvh,phi,g,id,1);
    const double gz=poisson_reconstruction_gradient_component(bvh,phi,g,id,2);
    const double grad=std::sqrt(gx*gx+gy*gy+gz*gz);
    return std::max(0.0,std::sqrt(std::max(0.0,grad*grad+2.0*phi[id]))-grad);
}

struct WallDistancePoissonManufacturedAudit {
    double operator_residual_inf{0.0};
    double gradient_error_inf{0.0};
    double reconstruction_l2_relative{0.0};
    double reconstruction_linf_relative{0.0};
    double eikonal_residual_inf{0.0};
    std::size_t samples{0};
};

// This audit deliberately separates three error sources:
//   (1) A_h(phi_exact)-b: consistency of the discrete Poisson operator;
//   (2) grad_h(phi_exact)-grad(phi_exact): reconstruction-gradient error;
//   (3) reconstructed d_h-d_exact: the Poisson-to-distance model itself.
// A small linear-solver residual alone cannot distinguish these effects.
inline WallDistancePoissonManufacturedAudit audit_poisson_manufactured(
    const WallDistanceBvh& bvh,const WallDistanceGrid& g,
    const std::vector<double>& phi_exact,
    const std::function<double(const WallDistanceVec3&)>& exact_distance,
    const std::function<WallDistanceVec3(const WallDistanceVec3&)>& exact_gradient) {
    if(phi_exact.size()!=g.points.size())
        throw std::invalid_argument("manufactured Poisson field size mismatch");
    WallDistancePoissonManufacturedAudit a;
    double e2=0.0, ref2=0.0, emax=0.0, refmax=0.0;
    for(std::size_t id=0;id<g.points.size();++id) {
        if(g.solid[id] || !std::isfinite(phi_exact[id])) continue;
        const auto x=g.points[id];
        const auto ge=exact_gradient(x);
        const double gx=poisson_reconstruction_gradient_component(bvh,phi_exact,g,id,0);
        const double gy=poisson_reconstruction_gradient_component(bvh,phi_exact,g,id,1);
        const double gz=poisson_reconstruction_gradient_component(bvh,phi_exact,g,id,2);
        const double ge_err=std::max({std::abs(gx-ge.x),std::abs(gy-ge.y),std::abs(gz-ge.z)});
        a.gradient_error_inf=std::max(a.gradient_error_inf,ge_err);
        const double op=std::abs(-poisson_laplacian_at(bvh,phi_exact,g,id)-1.0);
        a.operator_residual_inf=std::max(a.operator_residual_inf,op);
        const double d=poisson_reconstructed_distance(bvh,g,phi_exact,id);
        const double de=exact_distance(x);
        const double e=std::abs(d-de);
        e2+=e*e; ref2+=de*de; emax=std::max(emax,e); refmax=std::max(refmax,std::abs(de));
        a.eikonal_residual_inf=std::max(a.eikonal_residual_inf,std::abs(
            std::sqrt(gx*gx+gy*gy+gz*gz)-1.0));
        ++a.samples;
    }
    a.reconstruction_l2_relative=std::sqrt(e2/std::max(1e-30,ref2));
    a.reconstruction_linf_relative=emax/std::max(1e-30,refmax);
    return a;
}

inline std::vector<double> poisson_distance(const WallDistanceBvh& bvh,const WallDistanceGrid& g,
                                             std::size_t max_iter,double smooth,
                                             std::size_t* used_iter=nullptr,
                                             double* residual_out=nullptr,
                                             WallDistancePoissonReconstructionAudit* audit_out=nullptr) {
    std::vector<double> phi=poisson_potential(bvh,g,max_iter,smooth,used_iter,residual_out);
    std::vector<double> d(g.points.size(),std::numeric_limits<double>::infinity());
    for(std::size_t id=0;id<g.points.size();++id) if(!g.solid[id])
        d[id]=poisson_reconstructed_distance(bvh,g,phi,id);

    if(audit_out) {
        WallDistancePoissonReconstructionAudit a;
        double e2=0.0, ref2=0.0, emax=0.0, refmax=0.0;
        for(std::size_t id=0;id<g.points.size();++id) {
            if(g.solid[id] || !std::isfinite(phi[id]) || !std::isfinite(d[id])) continue;
            const double gx=poisson_reconstruction_gradient_component(bvh,phi,g,id,0);
            const double gy=poisson_reconstruction_gradient_component(bvh,phi,g,id,1);
            const double gz=poisson_reconstruction_gradient_component(bvh,phi,g,id,2);
            const double grad=std::sqrt(gx*gx+gy*gy+gz*gz);
            a.phi_min=std::min(a.phi_min,phi[id]);
            a.phi_max=std::max(a.phi_max,phi[id]);
            a.grad_min=std::min(a.grad_min,grad);
            a.grad_max=std::max(a.grad_max,grad);
            const double ref=bvh.nearest_distance(g.points[id]), e=d[id]-ref;
            e2+=e*e; ref2+=ref*ref;
            emax=std::max(emax,std::abs(e)); refmax=std::max(refmax,ref);
        }
        if(!std::isfinite(a.phi_min)) a.phi_min=0.0;
        if(!std::isfinite(a.grad_min)) a.grad_min=0.0;
        a.distance_l2_error=std::sqrt(e2/std::max(1e-30,ref2));
        a.distance_linf_error=emax/std::max(1e-30,refmax);
        *audit_out=a;
    }
    return d;
}

inline std::vector<double> poisson_distance(const WallSurface& s,const WallDistanceGrid& g,
                                             std::size_t max_iter,double smooth,
                                             std::size_t* used_iter=nullptr,
                                             double* residual_out=nullptr,
                                             WallDistancePoissonReconstructionAudit* audit_out=nullptr) {
    const WallDistanceBvh bvh(s);
    return poisson_distance(bvh,g,max_iter,smooth,used_iter,residual_out,audit_out);
}

inline double wall_distance_pde_residual_inf(WallDistanceMethod method,
                                                        const std::vector<double>& d,
                                                        const WallDistanceGrid& g,
                                                        const std::vector<unsigned char>* fixed=nullptr,
                                                        double advection_diffusion_gamma=0.05,
                                                        const WallDistanceBvh* bvh=nullptr);

inline std::vector<double> advection_diffusion_distance(const WallSurface& s,const WallDistanceGrid& g,
                                                        std::size_t max_iter,double gamma=0.05,
                                                        std::size_t* used_iter=nullptr,
                                                        double* residual_out=nullptr) {
    // NASA/Tucker transport form: U·∇d = 1 + Gamma ∇²d, with U derived
    // from the current distance field. Advection is first-order upwind and
    // diffusion is second-order central, matching the published formulation.
    const double h=std::min({g.spacing.x,g.spacing.y,g.spacing.z});
    const WallDistanceBvh bvh(s);
    const auto seeds=wall_seed_nodes(bvh,g,1.6*h);
    std::vector<unsigned char> fixed(g.points.size(),0);
    std::vector<double> d=eikonal_fast_sweep(bvh,g,std::max<std::size_t>(8,max_iter/2));
    for(std::size_t id=0;id<d.size();++id) if(g.solid[id]) d[id]=0.0;
    for(auto id:seeds) { d[id]=bvh.nearest_distance(g.points[id]); fixed[id]=1; }
    std::size_t used=max_iter;
    double final_residual=std::numeric_limits<double>::infinity();
    for(std::size_t it=0;it<max_iter;++it) {
        double max_change=0.0;
        for(std::size_t id=0;id<d.size();++id) {
            if(g.solid[id]||fixed[id]||!std::isfinite(d[id])) continue;
            const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
            const auto one=[&](int axis)->double {
                double gm=0.0,gp=0.0;
                if(axis==0) { if(i>0) gm=(d[g.index(i,j,k)]-d[g.index(i-1,j,k)])/g.spacing.x; if(i+1<g.nx) gp=(d[g.index(i+1,j,k)]-d[g.index(i,j,k)])/g.spacing.x; }
                if(axis==1) { if(j>0) gm=(d[g.index(i,j,k)]-d[g.index(i,j-1,k)])/g.spacing.y; if(j+1<g.ny) gp=(d[g.index(i,j+1,k)]-d[g.index(i,j,k)])/g.spacing.y; }
                if(axis==2) { if(k>0) gm=(d[g.index(i,j,k)]-d[g.index(i,j,k-1)])/g.spacing.z; if(k+1<g.nz) gp=(d[g.index(i,j,k+1)]-d[g.index(i,j,k)])/g.spacing.z; }
                return std::max(gm,0.0)+std::min(gp,0.0);
            };
            const double gx=one(0),gy=one(1),gz=one(2),gn=std::sqrt(gx*gx+gy*gy+gz*gz);
            const double ux=gn>1e-14?gx/gn:0.0, uy=gn>1e-14?gy/gn:0.0, uz=gn>1e-14?gz/gn:0.0;
            double diag=0.0,rhs=1.0;
            auto add_axis=[&](std::size_t m,bool has_m,std::size_t p,bool has_p,
                                   double hh,double u) {
                auto face_delta=[&](std::size_t q)->double {
                    return g.solid[q] ? poisson_wall_offset(bvh,g.points[id],g.points[q],hh) : hh;
                };
                if(has_m) {
                    const double delta=face_delta(m);
                    const double w=gamma/(hh*delta);
                    const double dm=g.solid[m]?0.0:d[m];
                    diag+=w; rhs+=w*dm;
                    if(u>0.0) {
                        const double a=u/delta;
                        diag+=a; if(!g.solid[m]) rhs+=a*d[m];
                    }
                }
                if(has_p) {
                    const double delta=face_delta(p);
                    const double w=gamma/(hh*delta);
                    const double dp=g.solid[p]?0.0:d[p];
                    diag+=w; rhs+=w*dp;
                    if(u<0.0) {
                        const double a=-u/delta;
                        diag+=a; if(!g.solid[p]) rhs+=a*d[p];
                    }
                }
            };
            const std::size_t xm=i>0?g.index(i-1,j,k):0, xp=i+1<g.nx?g.index(i+1,j,k):0;
            const std::size_t ym=j>0?g.index(i,j-1,k):0, yp=j+1<g.ny?g.index(i,j+1,k):0;
            const std::size_t zm=k>0?g.index(i,j,k-1):0, zp=k+1<g.nz?g.index(i,j,k+1):0;
            add_axis(xm,i>0,xp,i+1<g.nx,g.spacing.x,ux);
            add_axis(ym,j>0,yp,j+1<g.ny,g.spacing.y,uy);
            add_axis(zm,k>0,zp,k+1<g.nz,g.spacing.z,uz);
            if(diag>0.0) { const double nd=rhs/diag; max_change=std::max(max_change,std::abs(nd-d[id])); d[id]=nd; }
        }
        for(auto id:seeds) d[id]=bvh.nearest_distance(g.points[id]);
        if((it&3u)==3u || max_change<1e-10*h) {
            final_residual=wall_distance_pde_residual_inf(
                WallDistanceMethod::ADVECTION_DIFFUSION,d,g,&fixed,gamma,&bvh);
            if(final_residual<1e-6 && max_change<1e-8*h) {
                used=it+1;
                break;
            }
        }
    }
    if(used_iter) *used_iter=used;
    if(residual_out) *residual_out=final_residual;
    return d;
}

inline std::vector<double> hybrid_poisson_hamilton_jacobi_distance(
    const WallSurface& s,const WallDistanceGrid& g,std::size_t max_iter,
    double epsilon=0.25,double relaxation=0.9,
    std::size_t* used_iter=nullptr,double* residual_out=nullptr,
    double poisson_weight=1.0,
    std::size_t* poisson_iter_out=nullptr,double* poisson_residual_out=nullptr,
    std::size_t* clamp_activations=nullptr) {
    // Tucker 2011: Poisson supplies an auxiliary front-propagation direction
    // (effectively a wall normal); H-J then propagates distance with that
    // direction, optionally blended with the evolving Eikonal direction.
    const double h=std::min({g.spacing.x,g.spacing.y,g.spacing.z});
    const WallDistanceBvh bvh(s);
    const WallDistanceCutFaceOffsetCache offsets(bvh,g);
    const auto seeds=wall_seed_nodes(bvh,g,1.6*h);
    std::vector<unsigned char> fixed(g.points.size(),0);
    for(auto id:seeds) fixed[id]=1;

    std::size_t poisson_iter=0; double poisson_residual=0.0;
    const std::size_t poisson_max_iter=std::max<std::size_t>(200,4*max_iter);
    const auto phi=poisson_potential(bvh,g,poisson_max_iter,1.5,
                                     &poisson_iter,&poisson_residual,&offsets);
    if(poisson_iter_out) *poisson_iter_out=poisson_iter;
    if(poisson_residual_out) *poisson_residual_out=poisson_residual;

    std::vector<double> d(g.points.size(),std::numeric_limits<double>::infinity());
    auto grad_comp=[&](const std::vector<double>& f,std::size_t id,int a)->double {
        const std::size_t k=id/(g.nx*g.ny),rem=id%(g.nx*g.ny),j=rem/g.nx,i=rem%g.nx;
        if(a==0) {
            const bool hm=i>0, hp=i+1<g.nx;
            const auto im=hm?g.index(i-1,j,k):0, ip=hp?g.index(i+1,j,k):0;
            const bool fm=hm && !g.solid[im], fp=hp && !g.solid[ip];
            if(fm&&fp) return (f[ip]-f[im])/(2*g.spacing.x);
            if(fp) return (f[ip]-f[id])/g.spacing.x;
            if(fm) return (f[id]-f[im])/g.spacing.x;
        } else if(a==1) {
            const bool hm=j>0, hp=j+1<g.ny;
            const auto im=hm?g.index(i,j-1,k):0, ip=hp?g.index(i,j+1,k):0;
            const bool fm=hm && !g.solid[im], fp=hp && !g.solid[ip];
            if(fm&&fp) return (f[ip]-f[im])/(2*g.spacing.y);
            if(fp) return (f[ip]-f[id])/g.spacing.y;
            if(fm) return (f[id]-f[im])/g.spacing.y;
        } else {
            const bool hm=k>0, hp=k+1<g.nz;
            const auto im=hm?g.index(i,j,k-1):0, ip=hp?g.index(i,j,k+1):0;
            const bool fm=hm && !g.solid[im], fp=hp && !g.solid[ip];
            if(fm&&fp) return (f[ip]-f[im])/(2*g.spacing.z);
            if(fp) return (f[ip]-f[id])/g.spacing.z;
            if(fm) return (f[id]-f[im])/g.spacing.z;
        }
        return 0.0;
    };
    for(std::size_t id=0;id<d.size();++id) if(!g.solid[id]) {
        const double gx=grad_comp(phi,id,0),gy=grad_comp(phi,id,1),gz=grad_comp(phi,id,2);
        const double gn=std::sqrt(gx*gx+gy*gy+gz*gz);
        d[id]=gn>1e-14 ? std::max(0.0,std::sqrt(std::max(0.0,gn*gn+2*phi[id]))-gn) : 0.0;
    }
    for(auto id:seeds) d[id]=bvh.nearest_distance(g.points[id]);
    // The fixed set and finite-distance mask are invariant during transport.
    // Cache active IDs in ascending order to preserve the in-place Gauss-Seidel sweep.
    std::vector<std::size_t> active_ids;
    active_ids.reserve(d.size());
    for(std::size_t id=0;id<d.size();++id)
        if(!g.solid[id] && !fixed[id] && std::isfinite(d[id])) active_ids.push_back(id);

    const double inv_h2=1.0/(g.spacing.x*g.spacing.x)
                       +1.0/(g.spacing.y*g.spacing.y)
                       +1.0/(g.spacing.z*g.spacing.z);
    const double relaxation_factor=std::clamp(relaxation,0.1,1.0);
    const double alpha=std::clamp(poisson_weight,0.0,1.0);
    const std::size_t steps=std::max<std::size_t>(20,max_iter);
    std::size_t used=steps; double final_residual=std::numeric_limits<double>::infinity();
    std::size_t clamps=0;

    // The Poisson field is immutable during transport. Cache its normalized
    // gradient once instead of rebuilding three finite-difference components
    // and a norm for every active node on every pseudo-time iteration.
    // Keep one entry per active ID to bound storage to the transport domain.
    std::vector<WallDistanceVec3> phi_direction;
    phi_direction.reserve(active_ids.size());
    for(const std::size_t id:active_ids) {
        const double px=grad_comp(phi,id,0),py=grad_comp(phi,id,1),pz=grad_comp(phi,id,2);
        const double pn=std::sqrt(px*px+py*py+pz*pz);
        phi_direction.push_back({pn>1e-14?px/pn:0.0,
                                 pn>1e-14?py/pn:0.0,
                                 pn>1e-14?pz/pn:0.0});
    }
    auto velocity=[&](std::size_t active_index,std::size_t id,
                      double& ux,double& uy,double& uz) {
        const auto& base=phi_direction[active_index];
        ux=base.x; uy=base.y; uz=base.z;
        if(alpha<1.0) {
            const double dx=grad_comp(d,id,0),dy=grad_comp(d,id,1),dz=grad_comp(d,id,2);
            const double dn=std::sqrt(dx*dx+dy*dy+dz*dz);
            if(dn>1e-14) {
                ux=alpha*ux+(1-alpha)*dx/dn; uy=alpha*uy+(1-alpha)*dy/dn; uz=alpha*uz+(1-alpha)*dz/dn;
                const double n=std::sqrt(ux*ux+uy*uy+uz*uz);
                if(n>1e-14){ux/=n;uy/=n;uz/=n;}
            }
        }
    };

    for(std::size_t it=0;it<steps;++it) {
        double max_change=0.0;
        for(std::size_t active_index=0;active_index<active_ids.size();++active_index) {
            const std::size_t id=active_ids[active_index];
            const std::size_t k=id/(g.nx*g.ny),rem=id%(g.nx*g.ny),j=rem/g.nx,i=rem%g.nx;
            double ux,uy,uz; velocity(active_index,id,ux,uy,uz);
            const auto up=[&](double u,int a)->double {
                auto value=[&](std::size_t q)->double {
                    return g.solid[q] ? 0.0 : d[q];
                };
                if(a==0) {
                    if(u>=0) {
                        if(i==0) return 0.0;
                        const auto q=g.index(i-1,j,k);
                        const double delta=offsets.delta(id,0);
                        return u*(d[id]-value(q))/delta;
                    }
                    if(i+1==g.nx) return 0.0;
                    const auto q=g.index(i+1,j,k);
                    const double delta=offsets.delta(id,1);
                    return u*(value(q)-d[id])/delta;
                }
                if(a==1) {
                    if(u>=0) {
                        if(j==0) return 0.0;
                        const auto q=g.index(i,j-1,k);
                        const double delta=offsets.delta(id,2);
                        return u*(d[id]-value(q))/delta;
                    }
                    if(j+1==g.ny) return 0.0;
                    const auto q=g.index(i,j+1,k);
                    const double delta=offsets.delta(id,3);
                    return u*(value(q)-d[id])/delta;
                }
                if(u>=0) {
                    if(k==0) return 0.0;
                    const auto q=g.index(i,j,k-1);
                    const double delta=offsets.delta(id,4);
                    return u*(d[id]-value(q))/delta;
                }
                if(k+1==g.nz) return 0.0;
                const auto q=g.index(i,j,k+1);
                const double delta=offsets.delta(id,5);
                return u*(value(q)-d[id])/delta;
            };
            const double adv=up(ux,0)+up(uy,1)+up(uz,2);
            const double gamma=epsilon*std::max(d[id],0.0);
            const double residual=adv-1.0-gamma*laplacian_at(bvh,d,g,id,offsets);
            // The hybrid transport step is also a 3-D explicit upwind
            // update. Use its actual velocity components for the CFL bound.
            const double cfl_rate=std::abs(ux)/g.spacing.x+
                                  std::abs(uy)/g.spacing.y+
                                  std::abs(uz)/g.spacing.z;
            const double dt_adv=cfl_rate>1e-14 ? 0.9/cfl_rate :
                                                std::numeric_limits<double>::infinity();
            const double dt_diff=gamma>0.0?0.9/(2*gamma*inv_h2):std::numeric_limits<double>::infinity();
            const double dt=relaxation_factor*std::min(dt_adv,dt_diff);
            const double candidate=d[id]-dt*residual;
            if(candidate<0.0) ++clamps;  // F-W2: observable positivity clamp
            const double nd=std::max(0.0,candidate);
            max_change=std::max(max_change,std::abs(nd-d[id])); d[id]=nd;
        }
        // Seed IDs are fixed and excluded from active_ids, so their exact
        // distances cannot change during transport. Re-querying the BVH for
        // every seed on every iteration was redundant work.
        if((it&3u)==3u||max_change<1e-10*h) {
            // Evaluate exactly the same upwind transport operator used by the
            // pseudo-time update. The previous diagnostic used a central
            // gradient here, so the reported residual could disagree with
            // the nonlinear operator actually being iterated.
            final_residual=0.0;
            for(std::size_t active_index=0;active_index<active_ids.size();++active_index) {
            const std::size_t id=active_ids[active_index];
                double ux,uy,uz; velocity(active_index,id,ux,uy,uz);
                const std::size_t k=id/(g.nx*g.ny),rem=id%(g.nx*g.ny),j=rem/g.nx,i=rem%g.nx;
                auto upwind=[&](double u,int axis)->double {
                    auto value=[&](std::size_t q)->double { return g.solid[q] ? 0.0 : d[q]; };
                    if(axis==0) {
                        if(u>=0.0) {
                            if(i==0) return 0.0;
                            const auto q=g.index(i-1,j,k);
                            return u*(d[id]-value(q))/offsets.delta(id,0);
                        }
                        if(i+1==g.nx) return 0.0;
                        const auto q=g.index(i+1,j,k);
                        return u*(value(q)-d[id])/offsets.delta(id,1);
                    }
                    if(axis==1) {
                        if(u>=0.0) {
                            if(j==0) return 0.0;
                            const auto q=g.index(i,j-1,k);
                            return u*(d[id]-value(q))/offsets.delta(id,2);
                        }
                        if(j+1==g.ny) return 0.0;
                        const auto q=g.index(i,j+1,k);
                        return u*(value(q)-d[id])/offsets.delta(id,3);
                    }
                    if(u>=0.0) {
                        if(k==0) return 0.0;
                        const auto q=g.index(i,j,k-1);
                        return u*(d[id]-value(q))/offsets.delta(id,4);
                    }
                    if(k+1==g.nz) return 0.0;
                    const auto q=g.index(i,j,k+1);
                    return u*(value(q)-d[id])/offsets.delta(id,5);
                };
                const double adv=upwind(ux,0)+upwind(uy,1)+upwind(uz,2);
                final_residual=std::max(final_residual,
                    std::abs(adv-1.0-epsilon*std::max(d[id],0.0)*laplacian_at(bvh,d,g,id,offsets)));
            }
            if(final_residual<1e-8&&max_change<1e-9*h){used=it+1;break;}
        }
    }
    if(used_iter) *used_iter=used;
    if(residual_out) *residual_out=final_residual;
    if(clamp_activations) *clamp_activations=clamps;
    return d;
}

inline double wall_distance_pde_residual_inf(WallDistanceMethod method,
                                                        const std::vector<double>& d,
                                                        const WallDistanceGrid& g,
                                                        const std::vector<unsigned char>* fixed,
                                                        double advection_diffusion_gamma,
                                                        const WallDistanceBvh* bvh) {
    if(d.size()!=g.points.size() || g.solid.size()!=g.points.size() ||
       (fixed && fixed->size()!=g.points.size()))
        throw std::invalid_argument("wall_distance_pde_residual_inf: field/grid size mismatch");
    double rmax=0.0;
    std::size_t fluid_samples=0;
    const auto accumulate=[&](double residual) {
        rmax=std::isfinite(residual)
            ? std::max(rmax,std::abs(residual))
            : std::numeric_limits<double>::infinity();
    };
    for(std::size_t id=0;id<d.size();++id) {
        if(g.solid[id]) continue;
        ++fluid_samples;
        // A missing or corrupt value is not a zero PDE residual. Check fixed
        // seeds too: otherwise an invalid boundary field can certify a solve.
        if(!std::isfinite(d[id]) || d[id]<0.0)
            return std::numeric_limits<double>::infinity();
        if(fixed && (*fixed)[id]) continue;
        if(method==WallDistanceMethod::EIKONAL) {
            const double grad=bvh ? godunov_gradient_at(*bvh,d,g,id)
                                  : godunov_gradient_at(d,g,id);
            accumulate(grad-1.0);
            continue;
        }
        if(method==WallDistanceMethod::HAMILTON_JACOBI ||
           method==WallDistanceMethod::HYBRID_POISSON_EIKONAL) {
            const double gamma=0.25*std::max(d[id],0.0);
            const double grad=bvh ? godunov_gradient_at(*bvh,d,g,id)
                                  : godunov_gradient_at(d,g,id);
            const double lap=bvh ? laplacian_at(*bvh,d,g,id)
                                 : laplacian_at(d,g,id);
            accumulate(grad-1.0-gamma*lap);
            continue;
        }
        if(method==WallDistanceMethod::ADVECTION_DIFFUSION) {
            const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
            auto value=[&](std::size_t q)->double { return g.solid[q] ? 0.0 : d[q]; };
            auto one=[&](int axis)->double {
                double gm=0.0,gp=0.0;
                if(axis==0) {
                    if(i>0) gm=(d[id]-value(g.index(i-1,j,k)))/g.spacing.x;
                    if(i+1<g.nx) gp=(value(g.index(i+1,j,k))-d[id])/g.spacing.x;
                } else if(axis==1) {
                    if(j>0) gm=(d[id]-value(g.index(i,j-1,k)))/g.spacing.y;
                    if(j+1<g.ny) gp=(value(g.index(i,j+1,k))-d[id])/g.spacing.y;
                } else {
                    if(k>0) gm=(d[id]-value(g.index(i,j,k-1)))/g.spacing.z;
                    if(k+1<g.nz) gp=(value(g.index(i,j,k+1))-d[id])/g.spacing.z;
                }
                return std::max(gm,0.0)+std::min(gp,0.0);
            };
            const double gx=one(0), gy=one(1), gz=one(2);
            const double gn=std::sqrt(gx*gx+gy*gy+gz*gz);
            const double ux=gn>1e-14 ? gx/gn : 0.0;
            const double uy=gn>1e-14 ? gy/gn : 0.0;
            const double uz=gn>1e-14 ? gz/gn : 0.0;
            // Certify the same discrete equation as the in-place update in
            // advection_diffusion_distance(): normalized signed upwind U,
            // first-order face transport, and gamma/(h*delta) diffusion.
            // Godunov |grad d| plus a Shortley-Weller Laplacian is a different
            // operator and can stay O(1) after this solve reaches roundoff.
            double lhs=0.0;
            const auto face=[&](std::size_t q,bool exists,bool minus_side,double h,double u) {
                if(!exists) return;
                const double delta=(bvh && g.solid[q])
                    ? poisson_wall_offset(*bvh,g.points[id],g.points[q],h) : h;
                const double difference=d[id]-value(q);
                lhs+=advection_diffusion_gamma*difference/(h*delta);
                if((minus_side && u>0.0) || (!minus_side && u<0.0))
                    lhs+=std::abs(u)*difference/delta;
            };
            face(i>0?g.index(i-1,j,k):0,i>0,true,g.spacing.x,ux);
            face(i+1<g.nx?g.index(i+1,j,k):0,i+1<g.nx,false,g.spacing.x,ux);
            face(j>0?g.index(i,j-1,k):0,j>0,true,g.spacing.y,uy);
            face(j+1<g.ny?g.index(i,j+1,k):0,j+1<g.ny,false,g.spacing.y,uy);
            face(k>0?g.index(i,j,k-1):0,k>0,true,g.spacing.z,uz);
            face(k+1<g.nz?g.index(i,j,k+1):0,k+1<g.nz,false,g.spacing.z,uz);
            accumulate(lhs-1.0);
        }
    }
    return fluid_samples>0 ? rmax : std::numeric_limits<double>::infinity();
}

// Default arguments deliberately not applied: callers must state the method
// they want. The overload below (compute_wall_distance(surface, grid, iter))
// exposes the recommended default at the call site so the choice remains
// visible in reviews.
inline WallDistanceResult compute_wall_distance(WallDistanceMethod method,const WallSurface& s,
                                                 const WallDistanceGrid& g,
                                                 std::size_t iterations=500) {
    WallDistanceResult r;
    r.method=wall_distance_method_name(method);
    r.valid.assign(g.points.size(),1);
    switch(method) {
        case WallDistanceMethod::EXACT_GEOMETRIC:
            r.distance=exact_reference(s,g);
            r.converged=true;
            r.stopping_reason="direct";
            r.convergence_status=WallDistanceConvergenceStatus::DIRECT;
            break;
        case WallDistanceMethod::SEARCH_BASED:
            r.distance=search_based_reference(s,g,4.0*std::min({g.spacing.x,g.spacing.y,g.spacing.z}));
            r.converged=true;
            r.stopping_reason="direct";
            break;
        case WallDistanceMethod::MESH_WAVE:
            r=graph_wave(s,g,false);
            r.converged=true;
            r.stopping_reason="graph_complete";
            r.convergence_status=WallDistanceConvergenceStatus::DIRECT;
            break;
        case WallDistanceMethod::DIRECTIONAL_MESH_WAVE:
            r=graph_wave(s,g,true);
            r.converged=true;
            r.stopping_reason="graph_complete";
            break;
        case WallDistanceMethod::POISSON: {
            WallDistancePoissonReconstructionAudit pa;
            r.distance=poisson_distance(s,g,iterations,1.5,&r.iterations,&r.residual_inf,&pa);
            const WallDistanceBvh pbvh(s);
            const auto wa=audit_poisson_wall_offsets(pbvh,g);
            r.poisson_wall_ray_hits=wa.ray_hit_count;
            r.poisson_wall_ray_misses=wa.ray_miss_count;
            r.poisson_wall_fallbacks=wa.fallback_count;
            r.poisson_wall_bad_alignment=wa.bad_alignment_count;
            // F-W1: degenerate-normal fallbacks observed on the audited cut
            // faces. nearest_normal() is only reached on the ray-miss
            // fallback path of poisson_wall_offset_diagnostic(), so this
            // counts the audited offsets whose fallback normal was degenerate.
            r.nearest_normal_fallbacks=pbvh.degenerate_normal_fallbacks();
            r.poisson_wall_min_alignment=wa.min_alignment;
            r.poisson_phi_min=pa.phi_min;
            r.poisson_phi_max=pa.phi_max;
            r.poisson_grad_min=pa.grad_min;
            r.poisson_grad_max=pa.grad_max;
            r.poisson_distance_l2_error=pa.distance_l2_error;
            r.poisson_distance_linf_error=pa.distance_linf_error;
            r.converged = std::isfinite(r.residual_inf) && r.residual_inf < 1e-6;
            r.convergence_status = r.converged ? WallDistanceConvergenceStatus::CONVERGED :
                                   (r.iterations >= iterations ? WallDistanceConvergenceStatus::MAX_ITER : WallDistanceConvergenceStatus::RESIDUAL_TOO_HIGH);
            r.stopping_reason = r.converged ? "converged" :
                                (r.iterations >= iterations ? "max_iter" : "residual_too_high");
            break;
        }
        case WallDistanceMethod::EIKONAL: {
            const WallDistanceBvh bvh(s);
            r.distance=eikonal_fast_sweep(s,g,iterations,&r.iterations);
            std::vector<unsigned char> fixed(g.points.size(),0);
            for(const auto id:wall_seed_nodes(bvh,g,1.6*std::min({g.spacing.x,g.spacing.y,g.spacing.z}))) fixed[id]=1;
            r.residual_inf=wall_distance_pde_residual_inf(WallDistanceMethod::EIKONAL,r.distance,g,&fixed,0.05,&bvh);
            r.converged = std::isfinite(r.residual_inf) && r.residual_inf < 1e-6;
            r.convergence_status = r.converged ? WallDistanceConvergenceStatus::CONVERGED :
                                   (r.iterations >= iterations ? WallDistanceConvergenceStatus::MAX_ITER : WallDistanceConvergenceStatus::RESIDUAL_TOO_HIGH);
            r.stopping_reason = r.converged ? "converged" :
                                (r.iterations >= iterations ? "max_iter" : "residual_too_high");
            break;
        }
        case WallDistanceMethod::HAMILTON_JACOBI: {
            r.distance=hamilton_jacobi_distance(s,g,iterations,0.25,0.7,&r.iterations,&r.residual_inf,
                &r.positivity_clamp_activations);
            r.converged = std::isfinite(r.residual_inf) && r.residual_inf < 1e-6;
            r.convergence_status = r.converged ? WallDistanceConvergenceStatus::CONVERGED :
                                   (r.iterations >= iterations ? WallDistanceConvergenceStatus::MAX_ITER : WallDistanceConvergenceStatus::RESIDUAL_TOO_HIGH);
            r.stopping_reason = r.converged ? "converged" :
                                (r.iterations >= iterations ? "max_iter" : "residual_too_high");
            break;
        }
        case WallDistanceMethod::ADVECTION_DIFFUSION: {
            r.distance=advection_diffusion_distance(s,g,iterations,0.05,&r.iterations,&r.residual_inf);
            r.converged = std::isfinite(r.residual_inf) && r.residual_inf < 1e-6;
            r.convergence_status = r.converged ? WallDistanceConvergenceStatus::CONVERGED :
                                   (r.iterations >= iterations ? WallDistanceConvergenceStatus::MAX_ITER : WallDistanceConvergenceStatus::RESIDUAL_TOO_HIGH);
            r.stopping_reason = r.converged ? "converged" :
                                (r.iterations >= iterations ? "max_iter" : "residual_too_high");
            break;
        }
        case WallDistanceMethod::HYBRID_POISSON_EIKONAL: {
            std::size_t poisson_iter=0; double poisson_residual=0.0;
            r.distance=hybrid_poisson_hamilton_jacobi_distance(s,g,iterations,0.25,0.9,
                &r.iterations,&r.residual_inf,1.0,&poisson_iter,&poisson_residual,
                &r.positivity_clamp_activations);
            r.auxiliary_iterations=poisson_iter;
            r.auxiliary_residual_inf=poisson_residual;
            r.converged = std::isfinite(r.residual_inf) && r.residual_inf < 1e-6;
            r.convergence_status = r.converged ? WallDistanceConvergenceStatus::CONVERGED :
                                   (r.iterations >= iterations ? WallDistanceConvergenceStatus::MAX_ITER : WallDistanceConvergenceStatus::RESIDUAL_TOO_HIGH);
            r.stopping_reason = r.converged ? "converged" :
                                (r.iterations >= iterations ? "max_iter" : "residual_too_high");
            break;
        }
    }
    std::size_t fluid_samples=0, invalid_fluid_samples=0;
    for(std::size_t i=0;i<r.distance.size();++i) {
        if(g.solid[i]) { r.valid[i]=0; continue; }
        ++fluid_samples;
        if(!std::isfinite(r.distance[i]) || r.distance[i]<0.0) {
            r.valid[i]=0;
            ++invalid_fluid_samples;
        }
    }
    if(fluid_samples==0 || invalid_fluid_samples!=0) {
        r.converged=false;
        r.residual_inf=std::numeric_limits<double>::infinity();
        r.stopping_reason=fluid_samples==0 ? "no_fluid_samples" : "invalid_distance_field";
        r.convergence_status=WallDistanceConvergenceStatus::INVALID_FIELD;
    }
    return r;
}

// Convenience overload using the recommended default method (eikonal).
// Prefer this at the top-level call site unless the caller has an explicit
// reason to pick another method (e.g. requiring the smooth derivative of
// Poisson for SA/SST source terms).
inline WallDistanceResult compute_wall_distance(const WallSurface& s,
                                                 const WallDistanceGrid& g,
                                                 std::size_t iterations=500) {
    return compute_wall_distance(default_wall_distance_method, s, g, iterations);
}

struct WallDistanceBenchmarkMetrics {
    double l2_relative{0.0};
    double linf_relative{0.0};
    double near_wall_l2_relative{0.0};
    double monotonicity_violations{0.0};
    std::size_t samples{0};
    std::size_t invalid_samples{0};
    std::size_t near_wall_samples{0};
};

inline WallDistanceBenchmarkMetrics compare_wall_distance(const WallDistanceGrid& g,
                                                          const std::vector<double>& ref,
                                                          const std::vector<double>& value,
                                                          double near_wall) {
    if(ref.size()!=g.points.size() || value.size()!=g.points.size() ||
       g.solid.size()!=g.points.size() || g.neighbours.size()!=g.points.size())
        throw std::invalid_argument("compare_wall_distance: field/grid size mismatch");
    if(!std::isfinite(near_wall) || near_wall<0.0)
        throw std::invalid_argument("compare_wall_distance: invalid near-wall width");
    double sum2=0.0, ref2=0.0, maxe=0.0, max_ref=0.0, near_sum2=0.0, near_ref2=0.0;
    std::size_t n=0, near_n=0, violations=0, invalid=0;
    for(std::size_t i=0;i<ref.size();++i) {
        if(g.solid[i]) continue;
        if(!std::isfinite(ref[i]) || !std::isfinite(value[i]) || ref[i]<0.0 || value[i]<0.0) {
            ++invalid;
            continue;
        }
        const double e=std::abs(value[i]-ref[i]); sum2+=e*e; ref2+=ref[i]*ref[i]; maxe=std::max(maxe,e); max_ref=std::max(max_ref,ref[i]); ++n;
        if(ref[i]<=near_wall) { near_sum2+=e*e; near_ref2+=ref[i]*ref[i]; ++near_n; }
        for(auto nb:g.neighbours[i]) if(!g.solid[nb] && ref[nb]>ref[i] && value[nb]+1e-12<value[i]) ++violations;
    }
    const double unavailable=std::numeric_limits<double>::infinity();
    // Never turn an empty set or a partially invalid field into an accuracy
    // success. Counts are retained so the caller can diagnose the failure.
    if(n==0 || invalid!=0 || !std::isfinite(sum2) || !std::isfinite(ref2) ||
       !std::isfinite(near_sum2) || !std::isfinite(near_ref2))
        return {unavailable,unavailable,unavailable,static_cast<double>(violations),n,invalid,near_n};
    return {std::sqrt(sum2/std::max(1e-30,ref2)),
            maxe/std::max(1e-30,max_ref),
            near_n>0 ? std::sqrt(near_sum2/std::max(1e-30,near_ref2)) : unavailable,
            static_cast<double>(violations),n,invalid,near_n};
}

} // namespace cfdx::physics



```

### `tests/numerical/test_temporal.cpp`

SHA256 `7e2bf643a3c12a459b33c214af107602be6ee5b1a199088b830524bac52b993e`

```cpp
// M0.13-T02/T03 — Tests for Temporal Discretization (Local Time Stepping, Adaptive Time Step)

#include "cfdx/core/numerics/temporal.h"
#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/boundary.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include "common/test_harness.h"
#include <cmath>

using namespace cfdx::core;
using namespace cfdx::testing;

// Four closed hexahedra: dx=1/4, dy=dz=1/10, V=1/400 each.
// The former coplanar fixture had zero volume and only exercised the 1e30
// fallback in the CFL/local-step routines.
Mesh make_test_mesh() {
    constexpr std::size_t cells=4;
    Mesh mesh;
    mesh.points().resize(4*(cells+1));
    for(std::size_t i=0;i<=cells;++i) {
        const auto b=4*i;
        const double x=.25*static_cast<double>(i);
        mesh.points().set(b,  x,0.0,0.0);
        mesh.points().set(b+1,x,.1,0.0);
        mesh.points().set(b+2,x,.1,.1);
        mesh.points().set(b+3,x,0.0,.1);
        if(i==0) mesh.faces().push_face({b,b+3,b+2,b+1});
        else mesh.faces().push_face({b,b+1,b+2,b+3});
    }
    for(std::size_t i=0;i<cells;++i) {
        const auto b=4*i;
        const auto side=mesh.n_faces();
        mesh.faces().push_face({b,b+4,b+7,b+3});
        mesh.faces().push_face({b+1,b+2,b+6,b+5});
        mesh.faces().push_face({b,b+1,b+5,b+4});
        mesh.faces().push_face({b+3,b+7,b+6,b+2});
        mesh.cells().push_cell({i,i+1,side,side+1,side+2,side+3});
    }
    mesh.ownership().resize(mesh.n_faces());
    for(std::size_t i=0;i<=cells;++i) {
        mesh.ownership().set_owner(i,i==0 ? 0 : i-1);
        mesh.ownership().set_neighbour(i,i>0 && i<cells
            ? static_cast<CellIndex>(i) : FaceOwnership::BOUNDARY);
    }
    for(std::size_t i=0;i<cells;++i)
        for(std::size_t side=0;side<4;++side) {
            const auto face=cells+1+4*i+side;
            mesh.ownership().set_owner(face,i);
            mesh.ownership().set_neighbour(face,FaceOwnership::BOUNDARY);
        }
    BoundaryPatches patches;
    patches.add_patch("inlet",PatchType::INLET);
    patches.patch(0).face_ids={0};
    patches.add_patch("outlet",PatchType::OUTLET);
    patches.patch(1).face_ids={cells};
    patches.add_patch("walls",PatchType::WALL);
    for(std::size_t face=cells+1;face<mesh.n_faces();++face)
        patches.patch(2).face_ids.push_back(face);
    mesh.set_boundary(patches);
    return mesh;
}

// Simple RHS for testing: dφ/dt = -φ (decay)
void decay_rhs(const Field<double, Location::CELL>& phi, Field<double, Location::CELL>& rhs) {
    const std::size_t n = phi.size();
    const std::size_t dim = phi.dimension();
    for (std::size_t c = 0; c < n; ++c) {
        for (std::size_t d = 0; d < dim; ++d) {
            rhs.component_data(d)[c] = -phi.component_data(d)[c];
        }
    }
}

int main() {
    run_case("euler_explicit_decay", [&]() {
        Field<double, Location::CELL> phi(4, "phi", "1", 1);
        phi.fill(1.0);

        double dt = 0.1;
        auto phi_new = advance_time(phi, dt, decay_rhs, TimeScheme::EULER_EXPLICIT);

        // Exact: φ(t) = φ0 * exp(-t)
        // Euler: φ^{n+1} = φ^n * (1 - dt)
        EXPECT_NEAR(phi_new(0), 1.0 * (1.0 - dt), 1e-12);
    });

    run_case("euler_implicit_decay", []() {
        Field<double,Location::CELL> phi(4,"phi","1",2);
        for(std::size_t cell=0;cell<phi.size();++cell) {
            phi(cell,0)=1.0+.1*static_cast<double>(cell);
            phi(cell,1)=.5+.2*static_cast<double>(cell);
        }
        constexpr double dt=.1;
        const auto result=advance_time(phi,dt,decay_rhs,TimeScheme::EULER_IMPLICIT);
        for(std::size_t cell=0;cell<phi.size();++cell)
            for(std::size_t component=0;component<phi.dimension();++component)
                EXPECT_NEAR(result(cell,component),phi(cell,component)/(1.0+dt),1e-12);
    });

    run_case("crank_nicolson_decay", [&]() {
        Field<double, Location::CELL> phi(4, "phi", "1", 1);
        phi.fill(1.0);

        double dt = 0.1;
        auto phi_new = advance_time(phi, dt, decay_rhs, TimeScheme::CRANK_NICOLSON);

        // Solve (1+dt/2) phi^{n+1} = (1-dt/2) phi^n.
        // The former oracle was Heun's explicit predictor-corrector, not CN.
        const double expected=(1.0-.5*dt)/(1.0+.5*dt);
        for(std::size_t cell=0;cell<phi.size();++cell) {
            EXPECT_NEAR(phi_new(cell),expected,1e-12);
            EXPECT_NEAR(phi_new(cell)-phi(cell)+.5*dt*(phi(cell)+phi_new(cell)),0.0,1e-12);
        }
    });

    run_case("bdf2_requires_history", [&]() {
        Field<double, Location::CELL> phi(4, "phi", "1", 1);
        phi.fill(1.0);

        TimeIntegrationContext ctx(4, 1, "phi");
        ctx.initialize(phi);

        double dt = 0.1;
        // Bootstrap solves (1+dt) phi1 = phi0 by implicit Euler.
        EXPECT_FALSE(ctx.has_prev);
        const auto phi1=advance_time(phi,dt,decay_rhs,TimeScheme::BDF2,&ctx);
        EXPECT_TRUE(ctx.has_prev);
        EXPECT_NEAR(ctx.dt_prev,dt,1e-15);
        for(std::size_t cell=0;cell<phi.size();++cell) {
            EXPECT_NEAR(phi1(cell),1.0/(1.0+dt),1e-12);
            EXPECT_NEAR(ctx.phi_prev(cell),phi(cell),1e-15);
            EXPECT_NEAR(ctx.phi_curr(cell),phi1(cell),1e-15);
        }

        // BDF2: (3+2dt) phi2 = 4 phi1 - phi0.
        const auto phi2=advance_time(phi1,dt,decay_rhs,TimeScheme::BDF2,&ctx);
        for(std::size_t cell=0;cell<phi.size();++cell) {
            const double expected=(4.0*phi1(cell)-phi(cell))/(3.0+2.0*dt);
            EXPECT_NEAR(phi2(cell),expected,1e-12);
            EXPECT_NEAR((3.0+2.0*dt)*phi2(cell)-4.0*phi1(cell)+phi(cell),0.0,1e-12);
            EXPECT_NEAR(ctx.phi_prev(cell),phi1(cell),1e-15);
            EXPECT_NEAR(ctx.phi_curr(cell),phi2(cell),1e-15);
        }
    });

    run_case("bdf2_variable_step_implicit_decay", []() {
        Field<double,Location::CELL> phi(4,"phi","1",1);
        phi.fill(1.0);
        TimeIntegrationContext ctx(4,1,"phi");
        ctx.initialize(phi);
        const auto first=advance_time(phi,.1,decay_rhs,TimeScheme::BDF2,&ctx);
        // dt_new / dt_old = 2: 5*phi2 - 9*phi1 + 4*phi0 = .6*RHS(phi2).
        // This closed-form oracle does not call bdf2_coefficients().
        const auto second=advance_time(first,.2,decay_rhs,TimeScheme::BDF2,&ctx);
        for(std::size_t cell=0;cell<phi.size();++cell)
            EXPECT_NEAR(second(cell),(9.0*first(cell)-4.0*phi(cell))/5.6,1e-12);
        EXPECT_NEAR(ctx.dt_prev,.2,1e-15);
    });

    run_case("compute_cfl_time_step", [&]() {
        Mesh mesh = make_test_mesh();
        Field<double, Location::CELL> velocity(4, "U", "m/s", 3);
        // Set uniform velocity of 10 m/s in x direction
        for (std::size_t c = 0; c < 4; ++c) {
            velocity(c, 0) = 10.0;  // ux
            velocity(c, 1) = 0.0;   // uy
            velocity(c, 2) = 0.0;   // uz
        }

        double dt = compute_cfl_time_step(mesh, velocity, 0.5, 0.0);

        // The existing API uses h=V^(1/3), with V=.25*.1*.1.
        const double expected=.5*std::cbrt(.0025)/10.0;
        EXPECT_NEAR(dt,expected,1e-12);
        EXPECT_NEAR(compute_cfl_time_step(mesh,velocity,.5,2.0),
                    .5*std::cbrt(.0025)/12.0,1e-12);
    });

    run_case("adaptive_time_stepper", [&]() {
        AdaptiveTimeStepper stepper;
        stepper.dt = 0.1;
        stepper.dt_min = 1e-6;
        stepper.dt_max = 1.0;
        stepper.cfl_target = 0.5;

        // If CFL achieved is higher than target, reduce dt
        stepper.update(0.8);  // CFL = 0.8 > 0.5
        EXPECT_TRUE(stepper.dt < 0.1);

        // If CFL achieved is lower than target, increase dt
        stepper.dt = 0.1;
        stepper.update(0.2);  // CFL = 0.2 < 0.5
        EXPECT_TRUE(stepper.dt > 0.1);

        // Test reject
        stepper.dt = 0.1;
        stepper.reject();
        EXPECT_NEAR(stepper.dt, 0.05, 1e-12);  // shrink_factor = 0.5
    });

    run_case("local_time_stepping", [&]() {
        Mesh mesh = make_test_mesh();
        Field<double, Location::CELL> velocity(4, "U", "m/s", 3);
        // Varying velocities
        velocity(0, 0) = 10.0; velocity(0, 1) = 0.0; velocity(0, 2) = 0.0;
        velocity(1, 0) = 20.0; velocity(1, 1) = 0.0; velocity(1, 2) = 0.0;
        velocity(2, 0) = 5.0;  velocity(2, 1) = 0.0; velocity(2, 2) = 0.0;
        velocity(3, 0) = 1.0;  velocity(3, 1) = 0.0; velocity(3, 2) = 0.0;

        auto dt_local = compute_local_time_steps(mesh, velocity, 0.5);

        // Equal positive volumes make the expected local dt independent of
        // the geometry routine under test. A 1e30 fallback must fail.
        for(std::size_t c=0;c<4;++c)
            EXPECT_NEAR(dt_local[c],.5*std::cbrt(.0025)/velocity(c,0),1e-12);

        // Apply local time stepping
        Field<double, Location::CELL> phi(4, "phi", "1", 1);
        phi.fill(1.0);
        Field<double, Location::CELL> phi_new(4, "phi_new", "1", 1);

        apply_local_time_stepping(phi, dt_local, decay_rhs, phi_new);

        // Each cell updated with its own dt: φ_new = φ + dt_local * (-φ)
        for (std::size_t c = 0; c < 4; ++c) {
            double expected = 1.0 * (1.0 - dt_local[c]);
            EXPECT_NEAR(phi_new(c), expected, 1e-12);
        }
    });

    run_case("time_scheme_enum", [&]() {
        EXPECT_TRUE(std::string(to_string(TimeScheme::EULER_EXPLICIT)) == "euler_explicit");
        EXPECT_TRUE(std::string(to_string(TimeScheme::EULER_IMPLICIT)) == "euler_implicit");
        EXPECT_TRUE(std::string(to_string(TimeScheme::CRANK_NICOLSON)) == "crank_nicolson");
        EXPECT_TRUE(std::string(to_string(TimeScheme::BDF2)) == "bdf2");

        EXPECT_TRUE(time_scheme_from_string("euler_explicit") == TimeScheme::EULER_EXPLICIT);
        EXPECT_TRUE(time_scheme_from_string("euler_implicit") == TimeScheme::EULER_IMPLICIT);
        EXPECT_TRUE(time_scheme_from_string("crank_nicolson") == TimeScheme::CRANK_NICOLSON);
        EXPECT_TRUE(time_scheme_from_string("bdf2") == TimeScheme::BDF2);
    });

    return run_all();
}
```

### `tests/python/test_su2_naca0012_tutorial.py`

SHA256 `80ccf72e1b662db47da3f46d304b985671dd2d99ac59a021c088f9b0d0d49f47`

```python
"""The public NACA reader command cannot claim a CFDX flow validation."""
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("naca_reader_cli", ROOT / "validation/run_su2_naca0012.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_real_fixture_reports_reader_and_production_separately(tmp_path):
    report = module.qualify_reader(tmp_path)
    assert report["reader_passed"]
    assert all(report["reader_checks"].values())
    assert not report["canonical_conversion_passed"]
    blocking = [item for item in report["canonical_conversion_findings"]["findings"]
                if item["severity"] == "unsupported_blocking"]
    assert any(item["feature"] == "canonical_selection" and "JST" in item["detail"]
               for item in blocking)
    assert report["boundary_initial_condition_equivalence_passed"] == all(report["physical_metadata_checks"].values())
    assert not report["production_supported"]
    assert report["physical_validation"] == "NOT_RUN"
    assert not list(tmp_path.glob("*.cfdx.h5"))


def test_default_command_cannot_exit_success_for_an_unrun_flow(tmp_path):
    assert module.main(["--output-dir", str(tmp_path)]) == 2
    assert module.main(["--output-dir", str(tmp_path), "--reader-only"]) == 0

```

### `tests/python/test_tmr_naca0012_import.py`

SHA256 `f0a8ca8c73136eec1c312f3722360cf25f001acf720aac69beaf250190a2e715`

```python
"""Executable topology regressions for the NASA/TMR NACA C-grid importer."""
from pathlib import Path
import gzip
import sys

import h5py
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from import_tmr_naca0012 import convert, read_plot3d_2d, main
from generate_2d_validation_meshes import write_mesh


def miniature_c_grid():
    # Two paired wake edges on each side, a four-edge airfoil, one radial
    # layer. This is independent of NASA's ni/nj ratios and is deliberately
    # small enough for topology counts to be an exact independent oracle.
    x = np.array([[3, 2, 1, .5, 0, .5, 1, 2, 3],
                  [3, 2, 1, 0, -1, 0, 1, 2, 3]], dtype=float)
    y = np.array([[0, 0, 0, -.1, 0, .1, 0, 0, 0],
                  [-1, -1, -1, -1, 0, 1, 1, 1, 1]], dtype=float)
    return x, y


def plot3d_text(x, y):
    nj, ni = x.shape
    return f"1\n{ni} {nj}\n" + "\n".join(str(v) for v in np.r_[x.ravel(), y.ravel()])


@pytest.mark.parametrize("compressed", [False, True])
def test_reads_official_single_block_header_and_gzip(tmp_path, compressed):
    x, y = miniature_c_grid()
    path = tmp_path / ("n0012_9-2.p2dfmt.gz" if compressed else "n0012_9-2.p2dfmt")
    text = plot3d_text(x, y)
    if compressed:
        with gzip.open(path, "wt", encoding="ascii") as f:
            f.write(text)
    else:
        path.write_text(text, encoding="ascii")
    actual_x, actual_y = read_plot3d_2d(path)
    np.testing.assert_array_equal(actual_x, x)
    np.testing.assert_array_equal(actual_y, y)


@pytest.mark.parametrize("mutation,diagnostic", [
    (lambda s: s.replace("1\n", "2\n", 1), "single-block"),
    (lambda s: s + "\n123", "numeric tokens"),
    (lambda s: s.replace("3.0", "nan", 1), "finite"),
])
def test_rejects_unsupported_or_malformed_plot3d(tmp_path, mutation, diagnostic):
    x, y = miniature_c_grid()
    path = tmp_path / "n0012_9-2.p2dfmt"
    path.write_text(mutation(plot3d_text(x, y)), encoding="ascii")
    with pytest.raises(ValueError, match=diagnostic):
        read_plot3d_2d(path)


def test_wake_is_internal_and_extrusion_is_empty(tmp_path):
    points, cells, boundaries = convert(*miniature_c_grid())
    assert len(points) == 30  # 18 plane points minus 3 duplicate pairs, twice.
    assert len(cells) == 8
    assert {name: len(faces) for name, faces in boundaries.items()} == {
        "airfoil": 4, "farfield": 10, "front": 8, "back": 8,
    }
    path = tmp_path / "miniature.h5"
    write_mesh(path, points, cells, boundaries)
    with h5py.File(path, "r") as f:
        neighbours = f["neighbour"][:]
        assert np.count_nonzero(neighbours >= 0) == 9  # 7 adjacent + 2 joined wake.
        assert np.count_nonzero(neighbours == -1) == 30
        patch_ids = f["patch_face_ids"][:]
        assert len(np.unique(patch_ids)) == len(patch_ids)
        assert set(patch_ids) == set(np.flatnonzero(neighbours == -1))
        patch_types = {v.split(":")[0]: int(v.split(":")[-1])
                       for v in f.attrs["boundary_patches"].split(";")}
        assert patch_types == {"airfoil": 0, "farfield": 7, "front": 6, "back": 6}
        points = f["points"][:]
        face_vertices, face_offsets = f["face_vertices"][:], f["face_offsets"][:]
        centres = np.array([points[face_vertices[face_offsets[i]:face_offsets[i+1]]].mean(axis=0)
                            for i in range(len(neighbours))])
        on_wake = (centres[:, 0] > 1) & (np.abs(centres[:, 1]) < 1e-14) & (centres[:, 2] == .5)
        assert np.count_nonzero(on_wake) == 2
        assert (neighbours[on_wake] >= 0).all()
        # Signed polygon area need not share a winding convention; its
        # nonzero magnitude is independently checked for every input cell.
        for cell in cells:
            q = points[np.asarray(cell[:4]), :2]
            area = .5 * np.sum(q[:, 0] * np.roll(q[:, 1], -1) - q[:, 1] * np.roll(q[:, 0], -1))
            assert abs(area) > 0


def test_rejects_nonmatching_wake():
    x, y = miniature_c_grid()
    y[0, -1] = .01
    with pytest.raises(ValueError, match="coincident"):
        convert(x, y)


def test_command_records_provenance_without_physical_promotion(tmp_path, monkeypatch):
    source = tmp_path / "n0012_9-2.p2dfmt"
    output = tmp_path / "grid.h5"
    source.write_text(plot3d_text(*miniature_c_grid()), encoding="ascii")
    monkeypatch.setattr(sys, "argv", ["import_tmr_naca0012", "--input", str(source), "--output", str(output)])
    assert main() == 0
    with h5py.File(output, "r") as f:
        assert f.attrs["validation_wake_pairs"] == 2
        assert f.attrs["validation_wall_faces"] == 4
        assert f.attrs["validation_status"] == "imported-provenance-only"
        assert len(f.attrs["validation_source_sha256"]) == 64


def test_partial_wake_corruption_cannot_become_solid_wall():
    x, y = miniature_c_grid()
    y[0, -3] = .01  # break only the trailing-edge pair, after two good pairs.
    with pytest.raises(ValueError, match="trailing edge"):
        convert(x, y)

```

### `tests/unit/test_cg_solver.cpp`

SHA256 `41f1da4d6b9eae26ed3b1890be9807e5d571a404773c6ae244ea693b0b1783e7`

```cpp
// M0.8-T04 — Tests for CG solver

#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"
#include <limits>
#include <cmath>

using namespace cfdx::core;
using namespace cfdx::testing;

// Crée une matrice symétrique définie positive (diagonale dominante).
SparseMatrix make_poisson_1d(std::size_t n) {
    SparseMatrix A(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        A.push_back(i, i, 2.0);
        if (i > 0) A.push_back(i, i - 1, -1.0);
        if (i + 1 < n) A.push_back(i, i + 1, -1.0);
    }
    A.finalize();
    return A;
}

namespace {
enum class FailurePoint { Setup, InitialApply, IterationApply, ReplacementDot };

class FailingPreconditioner final : public Preconditioner {
public:
    explicit FailingPreconditioner(FailurePoint point) : point_(point) {}
    bool setup(const SparseMatrix&) override { return point_ != FailurePoint::Setup; }
    bool apply(const Vector& r, Vector& z) const override {
        ++applications;
        if ((point_ == FailurePoint::InitialApply && applications == 1) ||
            (point_ == FailurePoint::IterationApply && applications == 2))
            return false;
        const double sign =
            point_ == FailurePoint::ReplacementDot && applications >= 2 ? -1.0 : 1.0;
        for (std::size_t i = 0; i < r.size(); ++i) z(i) = sign * r(i);
        return true;
    }
    const char* name() const override { return "Failure injection"; }
    mutable std::size_t applications = 0;
private:
    FailurePoint point_;
};

void check_preconditioner_failure_residual(FailurePoint point) {
    SparseMatrix A(2, 2);
    A.push_back(0, 0, 1.0);
    A.push_back(1, 1, 2.0);
    A.finalize();
    Vector b(2, 1.0), x(2, 0.0);
    const bool before_iteration =
        point == FailurePoint::Setup || point == FailurePoint::InitialApply;
    if (before_iteration) {
        x(0) = 0.25;
        x(1) = -0.5;
    }
    FailingPreconditioner pc(point);
    const auto result = solve_cg(A, b, x, pc, 10, 1e-12);
    EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE);
    EXPECT_TRUE(result.iterations == (before_iteration ? 0u : 1u));
    if (point == FailurePoint::ReplacementDot) {
        // The third application follows residual replacement and rejects
        // r^T M^-1 r < 0, as in the captured NACA pressure failure.
        EXPECT_TRUE(pc.applications == 3);
        EXPECT_TRUE(result.residual_replacements == 1);
    }
    if (!before_iteration) {
        EXPECT_NEAR(x(0), 2.0 / 3.0, 1e-15);
        EXPECT_NEAR(x(1), 2.0 / 3.0, 1e-15);
    }
    // Independent two-equation oracle, evaluated at the returned iterate.
    const double expected = std::hypot(1.0 - x(0), 1.0 - 2.0 * x(1));
    EXPECT_TRUE(expected > 0.0);
    EXPECT_NEAR(result.residual, expected, 1e-14);
    EXPECT_NEAR(result.residual_relative, expected / std::sqrt(2.0), 1e-14);
}
} // namespace

int main() {
    run_case("cg_failed_setup_reports_current_true_residual", [] {
        check_preconditioner_failure_residual(FailurePoint::Setup);
    });
    run_case("cg_failed_initial_apply_reports_current_true_residual", [] {
        check_preconditioner_failure_residual(FailurePoint::InitialApply);
    });
    run_case("cg_failed_iteration_apply_reports_updated_true_residual", [] {
        check_preconditioner_failure_residual(FailurePoint::IterationApply);
    });
    run_case("cg_failed_replacement_reports_updated_true_residual", [] {
        check_preconditioner_failure_residual(FailurePoint::ReplacementDot);
    });
    run_case("cg_failed_zero_rhs_reports_nonzero_true_residual", [] {
        SparseMatrix A(1, 1);
        A.push_back(0, 0, 2.0);
        A.finalize();
        Vector b(1, 0.0), x(1, 3.0);
        FailingPreconditioner pc(FailurePoint::Setup);
        const auto result = solve_cg(A, b, x, pc, 10, 1e-12);
        EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE);
        EXPECT_TRUE(result.iterations == 0);
        EXPECT_NEAR(x(0), 3.0, 0.0);
        EXPECT_NEAR(result.residual, 6.0, 1e-14);
        EXPECT_TRUE(std::isinf(result.residual_relative));
        EXPECT_TRUE(result.residual_relative > 0.0);
    });
    run_case("cg_nonpositive_curvature_reports_true_residual", [] {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0); A.push_back(0, 1, 2.0);
        A.push_back(1, 0, 2.0); A.push_back(1, 1, 1.0);
        A.finalize();
        Vector b(2), x(2, 0.0);
        b(0) = 1.0; b(1) = -1.0;
        const auto result = solve_cg(A, b, x, 10, 1e-12);
        EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE);
        EXPECT_TRUE(result.iterations == 1);
        EXPECT_NEAR(x(0), 0.0, 0.0);
        EXPECT_NEAR(x(1), 0.0, 0.0);
        EXPECT_NEAR(result.residual, std::sqrt(2.0), 1e-14);
        EXPECT_NEAR(result.residual_relative, 1.0, 1e-14);
    });

    run_case("cg_solve_diagonal", []() {
        // A = 2I, b = [2, 2, 2], solution attendue x = [1, 1, 1].
        SparseMatrix A(3, 3);
        A.push_back(0, 0, 2.0);
        A.push_back(1, 1, 2.0);
        A.push_back(2, 2, 2.0);
        A.finalize();

        Vector b(3);
        b(0) = 2.0; b(1) = 2.0; b(2) = 2.0;

        Vector x(3, 0.0);
        auto result = solve_cg(A, b, x);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(result.iterations > 0);
        EXPECT_TRUE(result.residual_relative < 1e-10);
        EXPECT_NEAR(x(0), 1.0, 1e-10);
        EXPECT_NEAR(x(1), 1.0, 1e-10);
        EXPECT_NEAR(x(2), 1.0, 1e-10);
    });

    run_case("cg_solve_poisson_1d", []() {
        // Système 1D : -u'' = f avec u(0)=u(n-1)=0.
        // A = tridiag(-1, 2, -1), b = [0, 1, 1, 1, 0].
        const std::size_t n = 5;
        auto A = make_poisson_1d(n);

        Vector b(n);
        b(0) = 0.0; b(1) = 1.0; b(2) = 1.0; b(3) = 1.0; b(4) = 0.0;

        Vector x(n, 0.0);
        auto result = solve_cg(A, b, x);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(result.residual_relative < 1e-10);
        // Vérification : A x ≈ b.
        auto Ax = A.matvec(x);
        for (std::size_t i = 0; i < n; ++i) {
            EXPECT_NEAR(Ax[i], b(i), 1e-10);
        }
    });

    run_case("cg_already_converged", []() {
        // x = 0, b = 0 → déjà convergé.
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 1, 1.0);
        A.finalize();

        Vector b(2, 0.0);
        Vector x(2, 0.0);
        auto result = solve_cg(A, b, x);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(result.iterations == 0);
    });

    run_case("cg_max_iter_reached", []() {
        // Système mal conditionné avec très peu d'itérations.
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1e10);
        A.push_back(1, 1, 1e-10);
        A.finalize();

        Vector b(2);
        b(0) = 1.0; b(1) = 1.0;
        Vector x(2, 0.0);
        auto result = solve_cg(A, b, x, 1, 1e-20);
        // Soit convergé, soit max_iter atteint.
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED ||
                    result.status == SolverStatus::MAX_ITER_REACHED);
    });

    run_case("cg_not_applicable_non_square", []() {
        SparseMatrix A(2, 3);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 1, 1.0);
        A.finalize();

        Vector b(2);
        Vector x(3, 0.0);
        auto result = solve_cg(A, b, x);
        EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE);
    });

    run_case("cg_rejects_nonfinite_system", []() {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 1, 1.0);
        A.finalize();
        A.values_data()[1] = std::numeric_limits<double>::quiet_NaN();
        Vector b(2, 1.0);
        Vector x(2, 0.0);
        const auto result = solve_cg(A, b, x);
        EXPECT_TRUE(result.status == SolverStatus::DIVERGED);
    });

    run_case("cg_dimension_mismatch", []() {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 1, 1.0);
        A.finalize();

        Vector b(3);  // mauvaise taille
        Vector x(2, 0.0);
        auto result = solve_cg(A, b, x);
        EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE);
    });

    run_case("cg_exact_spd_3x3_reference_solution", []() {
        SparseMatrix A(3, 3);
        A.push_back(0, 0, 4.0); A.push_back(0, 1, 1.0);
        A.push_back(1, 0, 1.0); A.push_back(1, 1, 3.0); A.push_back(1, 2, 1.0);
        A.push_back(2, 1, 1.0); A.push_back(2, 2, 2.0);
        A.finalize();

        Vector x_exact(3);
        x_exact(0) = 1.0; x_exact(1) = -2.0; x_exact(2) = 3.0;
        const auto b_values = A.matvec(x_exact);
        Vector b(3);
        for (std::size_t i = 0; i < 3; ++i) b(i) = b_values[i];
        Vector x(3, 0.0);

        const auto result = solve_cg(A, b, x, 20, 1e-12);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(result.iterations <= 3);
        EXPECT_TRUE(result.residual_relative < 1e-12);
        const auto Ax = A.matvec(x);
        double true_residual_sq = 0.0;
        double b_sq = 0.0;
        for (std::size_t i = 0; i < 3; ++i) {
            const double r = Ax[i] - b(i);
            true_residual_sq += r * r;
            b_sq += b(i) * b(i);
            EXPECT_TRUE(std::isfinite(x(i)));
            EXPECT_NEAR(x(i), x_exact(i), 1e-11);
        }
        EXPECT_TRUE(std::sqrt(true_residual_sq / b_sq) < 1e-10);
    });

    run_case("cg_reports_non_spd_breakdown", []() {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0);
        A.push_back(1, 1, -1.0);
        A.finalize();

        Vector b(2);
        b(0) = 1.0; b(1) = 1.0;
        Vector x(2, 0.0);

        const auto result = solve_cg(A, b, x, 10, 1e-12);
        EXPECT_TRUE(result.status == SolverStatus::NOT_APPLICABLE ||
                    result.status == SolverStatus::DIVERGED);
    });

    return run_all();
}

```

### `tests/unit/test_core_temporal.cpp`

SHA256 `25b8eeb82a6bea1b99c85e0a2203cacce41d3983120362d55f70d60db595359c`

```cpp
#include "cfdx/core/numerics/temporal.h"
#include "common/test_harness.h"
#include <cmath>
#include <stdexcept>
#include <limits>

using namespace cfdx::core;

int main() {
    using cfdx::testing::run_case;

    run_case("implicit_euler", [] {
        ScalarCellField phi(1, "phi", "1", 1);
        phi(0) = 1.0;
        RhsFunction rhs = [](const ScalarCellField& in, ScalarCellField& out) {
            out(0) = -2.0 * in(0);
        };
        const auto result = advance_time(phi, 0.1, rhs, TimeScheme::EULER_IMPLICIT);
        EXPECT_NEAR(result(0), 1.0 / 1.2, 1e-10);
    });

    run_case("crank_nicolson", [] {
        ScalarCellField phi(1, "phi", "1", 1);
        phi(0) = 1.0;
        RhsFunction rhs = [](const ScalarCellField& in, ScalarCellField& out) {
            out(0) = -2.0 * in(0);
        };
        const auto result = advance_time(phi, 0.1, rhs, TimeScheme::CRANK_NICOLSON);
        EXPECT_NEAR(result(0), 0.9 / 1.1, 1e-10);
    });

    run_case("context_shape_mismatch", [] {
        ScalarCellField phi(2, "phi", "1", 1);
        TimeIntegrationContext ctx(1, 1, "phi");
        ctx.initialize(ScalarCellField(1, "initial", "1", 1));
        RhsFunction rhs = [](const ScalarCellField&, ScalarCellField& out) { out(0) = 0.0; };
        EXPECT_THROW(advance_time(phi, 0.1, rhs, TimeScheme::BDF2, &ctx), std::invalid_argument);
    });

    run_case("nonfinite_rhs_rejected", [] {
        ScalarCellField phi(1, "phi", "1", 1);
        phi(0) = 1.0;
        RhsFunction rhs = [](const ScalarCellField&, ScalarCellField& out) {
            out(0) = std::numeric_limits<double>::quiet_NaN();
        };
        EXPECT_THROW(advance_time(phi, 0.1, rhs, TimeScheme::EULER_EXPLICIT), std::runtime_error);
    });

    run_case("invalid_dt", [] {
        ScalarCellField phi(1, "phi", "1", 1);
        RhsFunction rhs = [](const ScalarCellField&, ScalarCellField& out) {
            out(0) = 0.0;
        };
        EXPECT_THROW(
            advance_time(phi, 0.0, rhs, TimeScheme::EULER_IMPLICIT),
            std::invalid_argument);
    });

    run_case("bdf2_bootstrap_and_second_step", [] {
        ScalarCellField phi(1, "phi", "1", 1);
        phi(0) = 1.0;
        TimeIntegrationContext ctx(1, 1, "phi");
        ctx.initialize(phi);
        RhsFunction rhs = [](const ScalarCellField& in, ScalarCellField& out) {
            out(0) = -2.0 * in(0);
        };

        const auto first = advance_time(phi, 0.1, rhs, TimeScheme::BDF2, &ctx);
        EXPECT_NEAR(first(0), 1.0 / 1.2, 1e-10);

        const auto second = advance_time(first, 0.1, rhs, TimeScheme::BDF2, &ctx);
        const double expected = (4.0 * first(0) - 1.0) / 3.4;
        EXPECT_NEAR(second(0), expected, 1e-10);
        EXPECT_TRUE(std::isfinite(second(0)));
        EXPECT_TRUE(second(0) > 0.0);
    });

    for (const auto scheme : {TimeScheme::EULER_IMPLICIT,
                              TimeScheme::CRANK_NICOLSON, TimeScheme::BDF2}) {
        run_case(std::string("intermediate_nonfinite_rhs_rejected_") + to_string(scheme), [scheme] {
            ScalarCellField phi(1, "phi", "1", 1);
            phi(0) = 1.0;
            TimeIntegrationContext ctx(1, 1, "phi");
            ctx.initialize(phi);
            const RhsFunction rhs = [](const ScalarCellField& in, ScalarCellField& out) {
                out(0) = -std::sqrt(in(0));
            };
            const double dt = scheme == TimeScheme::CRANK_NICOLSON ? 4.0 : 2.0;
            EXPECT_THROW(advance_time(phi, dt, rhs, scheme, &ctx), std::runtime_error);
            EXPECT_NEAR(ctx.phi_curr(0), 1.0, 0.0);
            EXPECT_FALSE(ctx.has_prev);
        });
    }

    run_case("nonfinite_initial_field_rejected_with_finite_rhs", [] {
        ScalarCellField phi(1, "phi", "1", 1);
        phi(0) = std::numeric_limits<double>::quiet_NaN();
        const RhsFunction rhs = [](const ScalarCellField&, ScalarCellField& out) { out(0) = 0.0; };
        EXPECT_THROW(advance_time(phi, 0.1, rhs, TimeScheme::EULER_EXPLICIT), std::runtime_error);
    });

    run_case("explicit_overflow_rejected_before_history_commit", [] {
        ScalarCellField phi(1, "phi", "1", 1);
        phi(0) = std::numeric_limits<double>::max();
        TimeIntegrationContext ctx(1, 1, "phi");
        ctx.initialize(phi);
        const RhsFunction rhs = [](const ScalarCellField&, ScalarCellField& out) {
            out(0) = std::numeric_limits<double>::max();
        };
        EXPECT_THROW(advance_time(phi, 1.0, rhs, TimeScheme::EULER_EXPLICIT, &ctx), std::runtime_error);
        EXPECT_NEAR(ctx.phi_curr(0), phi(0), 0.0);
        EXPECT_FALSE(ctx.has_prev);
    });

    run_case("bdf2_history_path_rejects_nonfinite_rhs_without_committing", [] {
        ScalarCellField phi(1, "phi", "1", 1);
        phi(0) = 1.0;
        TimeIntegrationContext ctx(1, 1, "phi");
        ctx.initialize(phi);
        ctx.shift(phi);
        ctx.dt_prev = 2.0;
        const RhsFunction rhs = [](const ScalarCellField& in, ScalarCellField& out) {
            out(0) = -std::sqrt(in(0));
        };
        EXPECT_THROW(advance_time(phi, 2.0, rhs, TimeScheme::BDF2, &ctx), std::runtime_error);
        EXPECT_NEAR(ctx.phi_curr(0), 1.0, 0.0);
        EXPECT_NEAR(ctx.phi_prev(0), 1.0, 0.0);
        EXPECT_NEAR(ctx.dt_prev, 2.0, 0.0);
        EXPECT_TRUE(ctx.has_prev);
    });

    return cfdx::testing::run_all();
}

```

### `tests/unit/test_matrix_characteristics.cpp`

SHA256 `b885b0ea28a8066237db26f811d4776a9d89c9ef0f969226283ca95cbc1be1b4`

```cpp
#include "cfdx/core/linalg/matrix_diagnostics.h"
#include "cfdx/core/linalg/linear_solver_models.h"
#include "common/test_harness.h"

#include <cmath>

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("measured_characteristics_feed_coupled_dispatch", [] {
        SparseMatrix matrix(3, 3);
        matrix.push_back(0, 0, 4.0);
        matrix.push_back(0, 1, -1.0);
        matrix.push_back(1, 0, -1.0);
        matrix.push_back(1, 1, 4.0);
        matrix.push_back(1, 2, -1.0);
        matrix.push_back(2, 1, -1.0);
        matrix.push_back(2, 2, 4.0);
        matrix.finalize();

        const auto characteristics = measure_matrix_characteristics(matrix, true);
        EXPECT_TRUE(characteristics.square);
        EXPECT_TRUE(characteristics.saddle_point);
        EXPECT_TRUE(characteristics.numerically_symmetric);
        EXPECT_TRUE(characteristics.diagonally_dominant);
        EXPECT_TRUE(std::abs(characteristics.average_nnz_per_row - 7.0 / 3.0) < 1e-14);

        LinearSolverRequest request;
        const auto plan = select_linear_solver(
            LinearProblemKind::CoupledPressureVelocity,
            characteristics,
            request);
        EXPECT_TRUE(plan.krylov == KrylovModel::FGMRES);
        EXPECT_TRUE(plan.preconditioner == PreconditionerModel::CoupledBlockSchur);
    });

    run_case("measured_characteristics_detect_nonsymmetry_and_scaling", [] {
        SparseMatrix matrix(2, 2);
        matrix.push_back(0, 0, 1.0e9);
        matrix.push_back(0, 1, 1.0);
        matrix.push_back(1, 0, 0.5);
        matrix.push_back(1, 1, 1.0e-2);
        matrix.finalize();

        const auto characteristics = measure_matrix_characteristics(matrix);
        EXPECT_TRUE(characteristics.square);
        EXPECT_TRUE(!characteristics.numerically_symmetric);
        EXPECT_TRUE(characteristics.strongly_scaled);
        EXPECT_TRUE(characteristics.coefficient_range > 1.0e10);
    });

    return run_all();
}

```

### `tests/unit/test_mesh.cpp`

SHA256 `7033031a1b6e6aad5c279309a04ad4e2b106d56cc409a54d84ff488b53d50033`

```cpp
// M0.3-T01 — Mesh validator tests (topology + geometry + quality + conservation)

#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/geometry/mesh_validator.h"
#include "common/test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

// Construit un cube unité [0,1]^3 — 1 cellule, 6 faces, 8 sommets.
Mesh make_unit_cube() {
    Mesh m;

    m.points().resize(8);
    m.points().set(0, 0.0, 0.0, 0.0);
    m.points().set(1, 1.0, 0.0, 0.0);
    m.points().set(2, 1.0, 1.0, 0.0);
    m.points().set(3, 0.0, 1.0, 0.0);
    m.points().set(4, 0.0, 0.0, 1.0);
    m.points().set(5, 1.0, 0.0, 1.0);
    m.points().set(6, 1.0, 1.0, 1.0);
    m.points().set(7, 0.0, 1.0, 1.0);

    // 6 faces (CCW vues de l'extérieur)
    m.faces().push_face({0, 3, 2, 1});  // 0 : bas
    m.faces().push_face({4, 5, 6, 7});  // 1 : haut
    m.faces().push_face({0, 1, 5, 4});  // 2 : avant
    m.faces().push_face({3, 7, 6, 2});  // 3 : arrière
    m.faces().push_face({0, 4, 7, 3});  // 4 : gauche
    m.faces().push_face({1, 2, 6, 5});  // 5 : droite

    m.ownership().resize(6);
    for (std::size_t i = 0; i < 6; ++i) {
        m.ownership().set_owner(i, 0);
        m.ownership().set_neighbour(i, FaceOwnership::BOUNDARY);
    }

    m.cells().push_cell({0, 1, 2, 3, 4, 5});

    BoundaryPatches bp;
    const char* names[6] = {"bottom", "top", "front", "back", "left", "right"};
    const PatchType types[6] = {PatchType::WALL, PatchType::WALL, PatchType::INLET,
                                 PatchType::OUTLET, PatchType::SYMMETRY, PatchType::WALL};
    for (int i = 0; i < 6; ++i) {
        Patch p;
        p.name = names[i];
        p.type = types[i];
        p.face_ids = {static_cast<std::uint32_t>(i)};
        bp.add_patch(p);
    }
    m.set_boundary(bp);

    return m;
}

int main() {
    run_case("default_empty", []() {
        Mesh m;
        EXPECT_TRUE(m.n_points() == 0);
        EXPECT_TRUE(m.n_faces() == 0);
        EXPECT_TRUE(m.n_cells() == 0);
        EXPECT_TRUE(m.topo_validate().ok);
    });

    run_case("unit_cube_stats", []() {
        Mesh m = make_unit_cube();
        EXPECT_TRUE(m.n_points() == 8);
        EXPECT_TRUE(m.n_faces() == 6);
        EXPECT_TRUE(m.n_cells() == 1);
        EXPECT_TRUE(m.stats().n_patches == 6);
        EXPECT_TRUE(m.stats().n_boundary_faces == 6);
        EXPECT_TRUE(m.stats().n_internal_faces == 0);
    });

    run_case("unit_cube_topo_validate_ok", []() {
        Mesh m = make_unit_cube();
        auto result = m.topo_validate();
        if (!result.ok) {
            for (const auto& e : result.errors) {
                std::fprintf(stderr, "  error: %s\n", e.c_str());
            }
        }
        EXPECT_TRUE(result.ok);
    });

    run_case("invalid_points_nan", []() {
        Mesh m = make_unit_cube();
        m.points().set(0, std::nan(""), 0.0, 0.0);
        EXPECT_FALSE(m.topo_validate().ok);
    });

    run_case("invalid_owner_out_of_range", []() {
        Mesh m = make_unit_cube();
        m.ownership().set_owner(0, 99);
        EXPECT_FALSE(m.topo_validate().ok);
    });

    run_case("invalid_cell_face_refs", []() {
        Mesh m = make_unit_cube();
        m.cells().clear();
        m.cells().push_cell({0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
        EXPECT_FALSE(m.topo_validate().ok);
    });

    run_case("boundary_patch_consistent", []() {
        Mesh m = make_unit_cube();
        EXPECT_TRUE(m.boundary().is_consistent(m.n_faces()));
    });

    run_case("boundary_patch_overlap", []() {
        Mesh m = make_unit_cube();
        Patch p;
        p.name = "dup";
        p.type = PatchType::WALL;
        p.face_ids = {0};
        m.boundary().add_patch(p);
        EXPECT_FALSE(m.boundary().is_consistent(m.n_faces()));
    });

    run_case("boundary_mesh_membership_complete", []() {
        const Mesh m = make_unit_cube();
        EXPECT_TRUE(m.boundary().is_consistent_with_mesh(m));
    });

    run_case("boundary_mesh_membership_orphan", []() {
        Mesh m = make_unit_cube();
        m.boundary().remove_patch(0);
        EXPECT_FALSE(m.boundary().is_consistent_with_mesh(m));
    });

    run_case("boundary_mesh_membership_duplicate", []() {
        Mesh m = make_unit_cube();
        m.boundary().patch(0).face_ids.push_back(0);
        EXPECT_FALSE(m.boundary().is_consistent_with_mesh(m));
    });

    run_case("boundary_mesh_membership_internal_face", []() {
        Mesh m = make_unit_cube();
        // Only the boundary-membership contract is under test here: a
        // nonnegative neighbour makes this a non-boundary face.
        m.ownership().set_neighbour(0, 1);
        EXPECT_FALSE(m.boundary().is_consistent_with_mesh(m));
    });

    run_case("boundary_mesh_membership_out_of_range", []() {
        Mesh m = make_unit_cube();
        m.boundary().patch(0).face_ids[0] = static_cast<FaceIndex>(m.n_faces());
        EXPECT_FALSE(m.boundary().is_consistent_with_mesh(m));
    });

    run_case("validate_mesh_full_ok", []() {
        Mesh m = make_unit_cube();
        auto report = validate_mesh(m);
        if (!report.ok) {
            for (const auto& e : report.errors) {
                std::fprintf(stderr, "  error: %s\n", e.c_str());
            }
        }
        EXPECT_TRUE(report.ok);
    });

    run_case("validate_mesh_detects_bad_geometry", []() {
        Mesh m = make_unit_cube();
        m.points().set(0, std::nan(""), 0.0, 0.0);
        auto report = validate_mesh(m);
        EXPECT_FALSE(report.ok);
    });

    run_case("validate_mesh_reports_quality", []() {
        Mesh m = make_unit_cube();
        auto report = validate_mesh(m);
        EXPECT_TRUE(report.max_skewness >= 0.0);
        EXPECT_TRUE(report.surface_closure_error >= 0.0);
        // Pour un cube parfait, la fermeture des surfaces est ~0.
        EXPECT_NEAR(report.surface_closure_error, 0.0, 1e-12);
        // Le volume du cube unité est 1.0.
        EXPECT_NEAR(report.max_cell_volume, 1.0, 1e-12);
        EXPECT_NEAR(report.min_cell_volume, 1.0, 1e-12);
    });

    return run_all();
}
```

### `tests/unit/test_wall_distance_advection_residual.cpp`

SHA256 `f22dedd6b2f355b4582dbef8e738726205522c39b88d44dc1b84ab754299e069`

```cpp
#include "cfdx/physics/wall_distance.h"
#include <cmath>
#include <iostream>
using namespace cfdx::physics;
int main() {
    int failed=0;
    const auto expect=[&](bool ok,const char* message){if(!ok){++failed;std::cerr<<message<<'\n';}};
    // Manufactured stationary ridge for the actual discrete transport rule.
    // Opposing derivatives cancel in U. With h=1/4, gamma=1/20 and
    // d_E=d_W=1/4, gamma*(2*d_P-d_E-d_W)/h^2=1 gives d_P=7/8.
    auto g=make_wall_distance_grid(5,3,3,{0,-1,-1},{.25,1,1},
        [](const WallDistanceVec3& p){return p.x<=0 || p.x>=1;});
    std::vector<double> d(g.points.size(),0.0);
    std::vector<unsigned char> fixed(g.points.size(),1);
    for(std::size_t id=0;id<d.size();++id)
        if(!g.solid[id]) d[id]=g.points[id].x==.5 ? .875 : .25;
    const auto centre=g.index(2,1,1);fixed[centre]=0;
    const double stationary=wall_distance_pde_residual_inf(
        WallDistanceMethod::ADVECTION_DIFFUSION,d,g,&fixed,.05);
    expect(stationary<1e-13,"stationary transport row was measured with another operator");
    // Zero field must retain the unit source even when U is undefined/zero.
    std::fill(d.begin(),d.end(),0.0);
    expect(std::abs(wall_distance_pde_residual_inf(
        WallDistanceMethod::ADVECTION_DIFFUSION,d,g,&fixed,.05)-1.0)<1e-13,
        "zero gradient erased the unit source");

    // Fractional wall at x=.1: at P=.25, d=(x-.1)^2 gives d_P=.0225,
    // d_E=.16, delta_W=.15. U=1, transport=.15 and the assembled diffusion
    // Laplacian is -.0225/(.25*.15)+(.16-.0225)/.25^2=1.6.
    // Hence |.15-1-.05*1.6|=.93, versus .95 with Shortley-Weller.
    WallSurface wall;wall.points={{.1,-4,-4},{.1,4,-4},{.1,4,4},{.1,-4,4}};
    wall.triangles={{{0,1,2}},{{0,2,3}}};WallDistanceBvh bvh(wall);
    auto cut=make_wall_distance_grid(5,3,3,{0,-1,-1},{.25,1,1},
        [](const WallDistanceVec3& p){return p.x<.1;});
    std::vector<double> quadratic(cut.points.size(),0.0);
    std::vector<unsigned char> mask(cut.points.size(),1);
    for(std::size_t id=0;id<quadratic.size();++id)
        if(!cut.solid[id])quadratic[id]=std::pow(cut.points[id].x-.1,2);
    mask[cut.index(1,1,1)]=0;
    const double fractional=wall_distance_pde_residual_inf(
        WallDistanceMethod::ADVECTION_DIFFUSION,quadratic,cut,&mask,.05,&bvh);
    expect(std::abs(fractional-.93)<1e-13,"fractional-wall residual changed the assembled diffusion stencil");
    // A stationary fractional-cut row with nonzero U also distinguishes the
    // operators without relying on the old zero-gradient skip: d_P=.2,
    // d_E=.95 gives transport=4/3 and the assembled Laplacian=20/3.
    // Thus transport-gamma*lap=1, while the S-W residual is 1/12.
    for(std::size_t id=0;id<quadratic.size();++id)
        if(!cut.solid[id])quadratic[id]=.2+3.0*(cut.points[id].x-.25);
    const double stationary_cut=wall_distance_pde_residual_inf(
        WallDistanceMethod::ADVECTION_DIFFUSION,quadratic,cut,&mask,.05,&bvh);
    expect(stationary_cut<1e-13,"stationary cut row was certified with a different diffusion operator");
    std::cout<<"stationary_residual="<<stationary<<" fractional_residual="<<fractional
             <<" stationary_cut_residual="<<stationary_cut<<" failed="<<failed<<'\n';
    return failed==0?0:1;
}

```

### `tests/unit/test_wall_distance_godunov_ridges.cpp`

SHA256 `71e1d748ed3a4920b9aebb017044d95d72624e1a4be2e08c5ac6ac22e2caebde`

```cpp
#include "cfdx/physics/wall_distance.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace cfdx::physics;
int main() {
    WallSurface wall;
    for(double x:{0.0,2.0}) {
        const std::size_t b=wall.points.size();
        wall.points.insert(wall.points.end(),{{x,-4,-4},{x,4,-4},{x,4,4},{x,-4,4}});
        wall.triangles.push_back({{b,b+1,b+2}});
        wall.triangles.push_back({{b,b+2,b+3}});
    }
    const auto grid=make_wall_distance_grid(9,3,3,{0,-1,-1},{.25,1,1},
        [](const WallDistanceVec3& p){return p.x<=0.0 || p.x>=2.0;});
    const WallDistanceBvh bvh(wall);
    const WallDistanceCutFaceOffsetCache cache(bvh,grid);
    const auto expected=exact_reference(wall,grid);
    const auto result=compute_wall_distance(WallDistanceMethod::EIKONAL,wall,grid,50);
    const std::size_t ridge=grid.index(4,1,1);
    const auto require=[](bool ok,const char* message){if(!ok)throw std::runtime_error(message);};
    require(std::abs(result.distance[ridge]-1.0)<1e-13,"channel ridge distance is not exact");
    for(const std::size_t i:{2u,4u,6u}) {
        const std::size_t id=grid.index(i,1,1);
        require(std::abs(godunov_gradient_at(expected,grid,id)-1.0)<1e-13,
                "Godunov gradient double-counted opposite directions");
        require(std::abs(godunov_gradient_at(bvh,expected,grid,id)-1.0)<1e-13,
                "BVH Godunov gradient double-counted opposite directions");
        require(std::abs(godunov_gradient_at(bvh,expected,grid,id,&cache)-1.0)<1e-13,
                "cached Godunov gradient is inconsistent");
    }
    const auto metrics=compare_wall_distance(grid,expected,result.distance,.5);
    require(metrics.l2_relative<1e-13,"channel distance differs from independent two-wall oracle");
    require(result.converged && result.residual_inf<1e-13,
            "exact channel distance is incorrectly reported as unconverged");
    std::cout << "Godunov ridge and monotone branches: PASS, L2=" << metrics.l2_relative
              << " residual=" << result.residual_inf << '\n';
}

```

### `tests/unit/test_wall_distance_invalid_results.cpp`

SHA256 `3b077dc0dc1ec7a4aa9fc3440db6ecd67f1fe9f707892ca5e71f7bbaf2f65a00`

```cpp
#include "cfdx/physics/wall_distance.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

using namespace cfdx::physics;

int main() {
    int failures=0, checks=0;
    const auto expect=[&](bool ok,const char* label) {
        ++checks;
        if(!ok) { ++failures; std::cerr << "FAIL: " << label << '\n'; }
    };
    const auto grid=make_wall_distance_grid(4,4,4,{10,0,0},{1,1,1},
        [](const WallDistanceVec3&){ return false; });
    std::vector<double> reference(grid.points.size(),1.0), distance=reference;
    auto metrics=compare_wall_distance(grid,reference,distance,2.0);
    expect(metrics.samples==grid.points.size() && metrics.l2_relative==0.0,
           "valid identical fields retain zero error");
    distance[0]=std::numeric_limits<double>::quiet_NaN();
    metrics=compare_wall_distance(grid,reference,distance,2.0);
    expect(std::isinf(metrics.l2_relative) && std::isinf(metrics.linf_relative),
           "one NaN cannot be excluded from the accuracy verdict");
    distance.assign(grid.points.size(),std::numeric_limits<double>::infinity());
    metrics=compare_wall_distance(grid,reference,distance,2.0);
    expect(metrics.samples==0 && std::isinf(metrics.l2_relative) &&
           std::isinf(metrics.near_wall_l2_relative),
           "a wholly invalid field cannot have zero global or near-wall error");
    distance=reference;
    reference[1]=std::numeric_limits<double>::quiet_NaN();
    metrics=compare_wall_distance(grid,reference,distance,2.0);
    expect(std::isinf(metrics.l2_relative),"a corrupt oracle cannot certify accuracy");
    reference=distance;
    distance[2]=-1.0;
    metrics=compare_wall_distance(grid,reference,distance,2.0);
    expect(std::isinf(metrics.l2_relative),"negative distance is an invalid field");
    distance=reference;
    metrics=compare_wall_distance(grid,reference,distance,0.5);
    expect(std::isinf(metrics.near_wall_l2_relative),
           "an empty near-wall sample band cannot report zero error");

    WallSurface wall;
    wall.points={{0,-4,-4},{0,4,-4},{0,4,4},{0,-4,4}};
    wall.triangles={{{0,1,2}},{{0,2,3}}};
    // The wall is outside the grid and its seed band. Geometric distance is
    // well defined, but propagation methods cannot reach any of these nodes.
    for(const auto method:{WallDistanceMethod::MESH_WAVE,
                           WallDistanceMethod::DIRECTIONAL_MESH_WAVE,
                           WallDistanceMethod::EIKONAL,
                           WallDistanceMethod::HAMILTON_JACOBI,
                           WallDistanceMethod::ADVECTION_DIFFUSION}) {
        const auto result=compute_wall_distance(method,wall,grid,20);
        expect(!result.converged && std::isinf(result.residual_inf),
               "a propagation result with no reachable sample cannot converge");
    }
    auto solid=grid;
    solid.solid.assign(grid.points.size(),1);
    const auto empty_result=compute_wall_distance(WallDistanceMethod::EXACT_GEOMETRIC,wall,solid,20);
    expect(!empty_result.converged,"an all-solid domain is not an accuracy certificate");
    metrics=compare_wall_distance(solid,reference,distance,2.0);
    expect(metrics.samples==0 && std::isinf(metrics.l2_relative),
           "an all-solid comparison has no accuracy evidence");

    distance.assign(grid.points.size(),1.0);
    distance[0]=std::numeric_limits<double>::quiet_NaN();
    std::vector<unsigned char> fixed(grid.points.size(),1);
    expect(std::isinf(wall_distance_pde_residual_inf(WallDistanceMethod::EIKONAL,distance,grid,&fixed)),
           "invalid fixed seeds invalidate the PDE residual");
    distance.assign(grid.points.size(),0.0);
    expect(wall_distance_pde_residual_inf(WallDistanceMethod::ADVECTION_DIFFUSION,distance,grid)>0.99,
           "zero gradient does not erase the unit source residual");
    std::cout << "checks=" << checks << " failures=" << failures << '\n';
    return failures==0 ? 0 : 1;
}

```

### `tests/validation/test_n11_scalar_conservation.cpp`

SHA256 `b7d7b4b9b4e8b9e41e56e602a6a60413432f9e597ada1c4f371e0bbb8f0308bc`

```cpp
#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/conservation_boundedness.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh make_unit_cube_with_x_patches()
{
    Mesh m;
    m.points().resize(8);
    const double p[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    for (std::size_t i = 0; i < 8; ++i)
        m.points().set(i,p[i][0],p[i][1],p[i][2]);

    m.faces().push_face({0,3,2,1}); // x=0
    m.faces().push_face({4,5,6,7}); // x=1
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3});
    m.faces().push_face({1,2,6,5});

    m.ownership().resize(6);
    for (std::size_t f = 0; f < 6; ++f) {
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});

    Patch x0;
    x0.name = "x0";
    x0.type = PatchType::WALL;
    x0.face_ids = {0};
    m.boundary().add_patch(x0);

    Patch x1;
    x1.name = "x1";
    x1.type = PatchType::WALL;
    x1.face_ids = {1};
    m.boundary().add_patch(x1);

    Patch walls;
    walls.name = "walls";
    walls.type = PatchType::WALL;
    walls.face_ids = {2,3,4,5};
    m.boundary().add_patch(walls);

    return m;
}

int main()
{
    run_case("n11_scalar_post_solve_balance_and_bounds", [] {
        const Mesh mesh = make_unit_cube_with_x_patches();
        const auto geometry = build_fv_geometry(mesh);

        Field<double,Location::FACE> mass_flux(mesh.n_faces(),"phi","kg/s",1);
        Field<double,Location::CELL> source(mesh.n_cells(),"source","1/s",1);
        Field<double,Location::CELL> implicit_source(mesh.n_cells(),"Sp","1/s",1);
        Field<double,Location::CELL> scalar(mesh.n_cells(),"scalar","1",1);
        mass_flux.fill(0.0);
        source.fill(0.0);
        implicit_source.fill(0.0);
        scalar.fill(0.0);

        ScalarBoundaryConditions bc;
        bc["x0"] = {ScalarBoundaryType::FIXED_VALUE,0.0,0.0};
        bc["x1"] = {ScalarBoundaryType::FIXED_VALUE,1.0,0.0};
        bc["walls"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

        const auto equation = assemble_scalar_equation(
            mesh, geometry, mass_flux, 1.0, source, implicit_source, bc,
            true, nullptr, nullptr, nullptr, nullptr,
            ConvectionScheme::UPWIND, &scalar);

        Vector solution(mesh.n_cells(), 0.0);
        const auto linear = solve_scalar_equation(
            equation, solution, ScalarSolveControls{200,1e-12,1.0});
        EXPECT_TRUE(linear.status == SolverStatus::CONVERGED);
        for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell)
            scalar(cell) = solution(cell);
        EXPECT_TRUE(std::isfinite(scalar(0)));
        EXPECT_TRUE(scalar(0) >= 0.0 && scalar(0) <= 1.0);

        // Rebuild the physical scalar flux from the accepted post-solve field.
        // This is deliberately independent of the assembled matrix residual.
        const auto reconstructed = reconstruct_scalar_transport_flux(
            mesh, geometry, mass_flux, scalar, 1.0, bc, true,
            ConvectionScheme::UPWIND);

        const auto balance = audit_transport_balance(
            mesh, reconstructed, source);

        EXPECT_TRUE(balance.finite());
        EXPECT_NEAR(balance.residual, 0.0, 1e-12);
        EXPECT_TRUE(balance.normalized_residual < 1e-12);
        EXPECT_TRUE(balance.max_cell_residual < 1e-12);

        const auto bounds = audit_scalar_bounds(scalar, 0.0, 1.0, 0.0);
        EXPECT_TRUE(bounds.bounded());
        EXPECT_TRUE(bounds.nonfinite == 0);
        EXPECT_TRUE(bounds.below == 0);
        EXPECT_TRUE(bounds.above == 0);

        std::cout << "N11_SCALAR_QUALIFICATION: PASS"
                  << " scalar=" << scalar(0)
                  << " balance_residual=" << balance.residual
                  << " balance_normalized=" << balance.normalized_residual
                  << " min=" << bounds.minimum
                  << " max=" << bounds.maximum
                  << "\n";
    });

    run_case("n11_scalar_nonfinite_and_bound_violation_are_hard_failures", [] {
        Field<double,Location::CELL> scalar(3,"scalar","1",1);
        scalar(0) = -1e-3;
        scalar(1) = 0.5;
        scalar(2) = std::numeric_limits<double>::quiet_NaN();

        const auto bounds = audit_scalar_bounds(scalar, 0.0, 1.0, 0.0);
        EXPECT_TRUE(!bounds.bounded());
        EXPECT_TRUE(bounds.below == 1);
        EXPECT_TRUE(bounds.above == 0);
        EXPECT_TRUE(bounds.nonfinite == 1);
    });
    return run_all();
}

```

### `tests/validation/test_n9_s10_pressure_reference_policies.cpp`

SHA256 `c02632151b9fd7e6ab9eebf321b9969780fe4076dff5b6b6e4c49705ca8c9950`

```cpp
#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"

#include <cstddef>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh make_two_cell_channel()
{
    Mesh m;
    m.points().resize(12);
    const double xyz[12][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1},
        {2,0,0},{2,1,0},{2,0,1},{2,1,1}
    };
    for (std::size_t i = 0; i < 12; ++i)
        m.points().set(i, xyz[i][0], xyz[i][1], xyz[i][2]);

    m.faces().push_face({0,3,7,4});
    m.faces().push_face({1,2,6,5});
    m.faces().push_face({8,9,11,10});
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({1,8,10,5});
    m.faces().push_face({2,6,11,9});
    m.faces().push_face({1,2,9,8});
    m.faces().push_face({5,10,11,6});

    m.ownership().resize(11);
    for (std::size_t f = 0; f < 11; ++f) {
        const auto owner = (f == 2 || f >= 7) ? 1 : 0;
        m.ownership().set_owner(f, owner);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(1, 0);
    m.ownership().set_neighbour(1, 1);
    m.cells().push_cell({0,1,3,4,5,6});
    m.cells().push_cell({1,2,7,8,9,10});

    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    wall.face_ids = {0,2,3,4,5,6,7,8,9,10};
    m.boundary().add_patch(wall);
    return m;
}

static void configure_case(IncompressibleSolverControls& c, PressureGaugePolicy policy)
{
    c.algorithm = PressureVelocityAlgorithm::SIMPLE;
    c.pressure_gauge_policy = policy;
    c.pressure_reference_cell = 0;
    c.pressure_reference_value = 3.0;
    c.convergence.max_iterations = 2;
    c.convergence.continuity_tolerance = 1e-12;
    c.convergence.relative_tolerance = 1e-12;
    c.linear_max_iterations = 100;
    c.linear_tolerance = 1e-12;
}

int main()
{
    run_case("n9_s10_reference_cell_policy", [] {
        const Mesh mesh = make_two_cell_channel();
        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(0.0);
        p.fill(7.0);
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        IncompressibleSolverControls c;
        configure_case(c, PressureGaugePolicy::REFERENCE_CELL);
        const auto result = solve_steady_incompressible(mesh, U, p, ubc, pbc, c);
        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(p(0), 3.0, 1e-12);
        EXPECT_NEAR(p(1), 3.0, 1e-12);
    });

    run_case("n9_s10_zero_mean_shift_preserves_pressure_differences", [] {
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        p(0) = 2.0;
        p(1) = 6.0;

        // Gauge removal changes only the additive constant. This is an
        // operator contract, independent of the momentum/continuity solve.
        apply_pressure_gauge_policy(
            p, geometry, PressureGaugePolicy::ZERO_MEAN, 0, 3.0, false);
        EXPECT_NEAR(0.5 * (p(0) + p(1)), 0.0, 1e-12);
        EXPECT_NEAR(p(1) - p(0), 4.0, 1e-12);
    });

    run_case("n9_s10_zero_mean_equilibrium_solver", [] {
        const Mesh mesh = make_two_cell_channel();
        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(0.0);
        // A closed, unforced box at rest has constant pressure. A nonzero
        // initial pressure gradient is not an equilibrium that a steady solve
        // should preserve; that property belongs only to the gauge helper.
        p.fill(7.0);
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        IncompressibleSolverControls c;
        configure_case(c, PressureGaugePolicy::ZERO_MEAN);
        const auto result = solve_steady_incompressible(mesh, U, p, ubc, pbc, c);
        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(0.5 * (p(0) + p(1)), 0.0, 1e-12);
        EXPECT_NEAR(p(1) - p(0), 0.0, 1e-12);
        for (std::size_t cell = 0; cell < U.size(); ++cell)
            for (std::size_t component = 0; component < 3; ++component)
                EXPECT_NEAR(U(cell, component), 0.0, 1e-12);
    });

    run_case("n9_s10_none_policy_preserves_pressure_level", [] {
        const Mesh mesh = make_two_cell_channel();
        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(0.0);
        p.fill(7.0);
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        IncompressibleSolverControls c;
        configure_case(c, PressureGaugePolicy::NONE);
        c.pressure_linear_solver.null_space = NullSpaceModel::Constant;
        const auto result = solve_steady_incompressible(mesh, U, p, ubc, pbc, c);
        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(p(0), 7.0, 1e-12);
        EXPECT_NEAR(p(1), 7.0, 1e-12);
    });
    return run_all();
}

```

### `tests/validation/test_n9_s9_rhie_chow_consistency.cpp`

SHA256 `6126d2c7369228b3e50c1525d756fa0e41f303c23ab27142ac1891deaa5c30ff`

```cpp
#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh make_two_cell_channel()
{
    Mesh m;
    m.points().resize(12);
    const double xyz[12][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1},
        {2,0,0},{2,1,0},{2,0,1},{2,1,1}
    };
    for (std::size_t i = 0; i < 12; ++i)
        m.points().set(i, xyz[i][0], xyz[i][1], xyz[i][2]);

    m.faces().push_face({0,3,7,4});       // x = 0
    m.faces().push_face({1,2,6,5});       // internal
    m.faces().push_face({8,9,11,10});     // x = 2
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({1,8,10,5});
    m.faces().push_face({2,6,11,9});
    m.faces().push_face({1,2,9,8});
    m.faces().push_face({5,10,11,6});

    m.ownership().resize(11);
    for (std::size_t f = 0; f < 11; ++f) {
        const auto owner = (f == 2 || f >= 7) ? 1 : 0;
        m.ownership().set_owner(f, owner);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(1, 0);
    m.ownership().set_neighbour(1, 1);

    m.cells().push_cell({0,1,3,4,5,6});
    m.cells().push_cell({1,2,7,8,9,10});

    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    wall.face_ids = {0,2,3,4,5,6,7,8,9,10};
    m.boundary().add_patch(wall);
    return m;
}

static std::array<std::vector<double>, 3> unit_rAU(std::size_t n)
{
    return {
        std::vector<double>(n, 1.0),
        std::vector<double>(n, 1.0),
        std::vector<double>(n, 1.0)};
}

// Stretched two-cell channel: cell 0 spans x in [0, 0.5], cell 1 spans
// x in [0.5, 2] (stretch ratio 3). Same topology and face winding as
// make_two_cell_channel(); only the x-coordinates of the internal and
// right-boundary points change.
static Mesh make_two_cell_channel_stretched()
{
    Mesh m;
    m.points().resize(12);
    const double xyz[12][3] = {
        {0,0,0},{0.5,0,0},{0.5,1,0},{0,1,0},
        {0,0,1},{0.5,0,1},{0.5,1,1},{0,1,1},
        {2,0,0},{2,1,0},{2,0,1},{2,1,1}
    };
    for (std::size_t i = 0; i < 12; ++i)
        m.points().set(i, xyz[i][0], xyz[i][1], xyz[i][2]);

    m.faces().push_face({0,3,7,4});       // x = 0
    m.faces().push_face({1,2,6,5});       // internal, x = 0.5
    m.faces().push_face({8,9,11,10});     // x = 2
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({1,8,10,5});
    m.faces().push_face({2,6,11,9});
    m.faces().push_face({1,2,9,8});
    m.faces().push_face({5,10,11,6});

    m.ownership().resize(11);
    for (std::size_t f = 0; f < 11; ++f) {
        const auto owner = (f == 2 || f >= 7) ? 1 : 0;
        m.ownership().set_owner(f, owner);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(1, 0);
    m.ownership().set_neighbour(1, 1);

    m.cells().push_cell({0,1,3,4,5,6});
    m.cells().push_cell({1,2,7,8,9,10});

    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    wall.face_ids = {0,2,3,4,5,6,7,8,9,10};
    m.boundary().add_patch(wall);
    return m;
}

static std::array<std::vector<double>, 3> scale_rAU(
    const std::array<std::vector<double>, 3>& rAU, double s)
{
    std::array<std::vector<double>, 3> out{};
    for (std::size_t k = 0; k < 3; ++k) {
        out[k].resize(rAU[k].size());
        for (std::size_t i = 0; i < rAU[k].size(); ++i)
            out[k][i] = s * rAU[k][i];
    }
    return out;
}

// The Rhie-Chow pressure correction of an internal face, extracted as
// phi_interpolated - phi_RC. Both fluxes come from the production paths
// (make_mass_flux / make_rhie_chow_mass_flux), so no formula is duplicated
// here: only the difference is taken.
static double rc_correction(std::size_t face,
    const Mesh& mesh,
    const FvGeometry& geometry,
    const Field<double, Location::CELL>& U,
    const Field<double, Location::CELL>& p,
    const std::array<std::vector<double>, 3>& rAU,
    double rho,
    const VelocityBoundaryConditions& ubc,
    const ScalarBoundaryConditions& pbc)
{
    const auto phi_interp = make_mass_flux(mesh, geometry, U, rho, ubc);
    const auto phi_rc = make_rhie_chow_mass_flux(
        mesh, geometry, U, p, rAU, rho, ubc, pbc);
    return phi_interp(face) - phi_rc(face);
}

static Mesh make_two_cell_channel_oblique()
{
    Mesh m;
    m.points().resize(12);
    constexpr double shear = 0.35;
    const double xyz[12][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1},
        {2,0,0},{2,1,0},{2,0,1},{2,1,1}
    };
    for (std::size_t i = 0; i < 12; ++i)
        m.points().set(i, xyz[i][0] + shear * xyz[i][1],
                       xyz[i][1], xyz[i][2]);

    m.faces().push_face({0,3,7,4});
    m.faces().push_face({1,2,6,5}); // Internal face is non-orthogonal to C_O C_N.
    m.faces().push_face({8,9,11,10});
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({1,8,10,5});
    m.faces().push_face({2,6,11,9});
    m.faces().push_face({1,2,9,8});
    m.faces().push_face({5,10,11,6});

    m.ownership().resize(11);
    for (std::size_t f = 0; f < 11; ++f) {
        const auto owner = (f == 2 || f >= 7) ? 1 : 0;
        m.ownership().set_owner(f, owner);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(1, 0);
    m.ownership().set_neighbour(1, 1);
    m.cells().push_cell({0,1,3,4,5,6});
    m.cells().push_cell({1,2,7,8,9,10});

    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    wall.face_ids = {0,2,3,4,5,6,7,8,9,10};
    m.boundary().add_patch(wall);
    return m;
}

int main()
{
    run_case("n9_s9_rhie_chow_linear_pressure_consistency", [] {
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        Field<double, Location::CELL> grad_p(2, "grad_p", "Pa/m", 3);
        U.set(0, 2.0, 0.0, 0.0);
        U.set(1, 2.0, 0.0, 0.0);
        p(0) = 1.0;
        p(1) = 3.0;

        // Exact linear pressure field: dp/dx = 2 Pa/m.
        grad_p.set(0, 2.0, 0.0, 0.0);
        grad_p.set(1, 2.0, 0.0, 0.0);

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        const auto phi = make_rhie_chow_mass_flux(
            mesh, geometry, U, p, rAU, 1.0, ubc, pbc);
        const auto phi_linear = make_mass_flux(mesh, geometry, U, 1.0, ubc);

        // The Rhie-Chow pressure correction must vanish for a pressure field
        // whose cell gradient is exactly the same as its face pressure jump.
        EXPECT_NEAR(phi(1), phi_linear(1), 1e-12);
    });

    run_case("n9_s9_rhie_chow_pressure_offset_invariance", [] {
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p0(2, "p", "Pa", 1);
        Field<double, Location::CELL> p1(2, "p", "Pa", 1);
        U.fill(0.0);
        p0(0) = -2.0;
        p0(1) = 4.0;
        p1(0) = 1002.0;
        p1(1) = 1008.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        const auto phi0 = make_rhie_chow_mass_flux(
            mesh, geometry, U, p0, rAU, 1.0, ubc, pbc);
        const auto phi1 = make_rhie_chow_mass_flux(
            mesh, geometry, U, p1, rAU, 1.0, ubc, pbc);

        for (std::size_t f = 0; f < mesh.n_faces(); ++f)
            EXPECT_NEAR(phi0(f), phi1(f), 1e-12);
    });

    run_case("n9_s9_rhie_chow_alternating_pressure_mode_is_coupled", [] {
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(0.0);
        p(0) = 1.0;
        p(1) = -1.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        const auto phi = make_rhie_chow_mass_flux(
            mesh, geometry, U, p, rAU, 1.0, ubc, pbc);

        // A collocated cell-centred velocity interpolation alone would produce
        // zero mass flux. Rhie-Chow must couple this alternating pressure mode
        // to the face flux instead of allowing a checkerboard pressure mode to
        // become invisible to continuity.
        EXPECT_TRUE(std::abs(phi(1)) > 1e-12);
        EXPECT_NEAR(phi(1), 2.0, 1e-12);
    });


    run_case("n9_s9_rhie_chow_correction_scales_linearly_with_rAU_majumdar", [] {
        // D^RC = V/a_P^r = alpha_u * V/a_P = alpha_u * D^pc: the production
        // RC term weights the face pressure flux with the RELAXED inverse
        // momentum diagonal. Scaling rAU by alpha_u must therefore scale the
        // RC correction by exactly alpha_u. This is the operator-level
        // statement of the Majumdar relation documented in
        // N9_S9_RHIE_CHOW_CONSISTENCY.md.
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(0.0);
        p(0) = 1.0;
        p(1) = -1.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        const double c1 = rc_correction(
            1, mesh, geometry, U, p, rAU, 1.0, ubc, pbc);
        EXPECT_TRUE(std::abs(c1) > 1e-12);
        for (const double alpha : {0.65, 0.8}) {
            const auto rAU_a = scale_rAU(rAU, alpha);
            const double ca = rc_correction(
                1, mesh, geometry, U, p, rAU_a, 1.0, ubc, pbc);
            EXPECT_NEAR(ca, alpha * c1, 1e-12);
        }
    });

    run_case("n9_s9_rhie_chow_correction_is_linear_in_pressure", [] {
        // The internal RC operator is linear in p (the Green-Gauss gradient
        // and the owner-to-neighbour difference are linear, and the
        // zero-gradient wall treatment preserves linearity):
        // corr(p_a + p_b) = corr(p_a) + corr(p_b).
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> pa(2, "pa", "Pa", 1);
        Field<double, Location::CELL> pb(2, "pb", "Pa", 1);
        Field<double, Location::CELL> psum(2, "psum", "Pa", 1);
        U.fill(0.0);
        pa(0) = 1.0;
        pa(1) = -1.0;      // alternating mode
        pb(0) = 1.0;
        pb(1) = 3.0;       // linear mode
        psum(0) = pa(0) + pb(0);
        psum(1) = pa(1) + pb(1);

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        const double ca = rc_correction(
            1, mesh, geometry, U, pa, rAU, 1.0, ubc, pbc);
        const double cb = rc_correction(
            1, mesh, geometry, U, pb, rAU, 1.0, ubc, pbc);
        const double cs = rc_correction(
            1, mesh, geometry, U, psum, rAU, 1.0, ubc, pbc);
        EXPECT_TRUE(std::abs(ca) > 1e-12);
        EXPECT_TRUE(std::abs(cb) > 1e-12);
        EXPECT_NEAR(cs, ca + cb, 1e-12);
    });

    run_case("n9_s9_rhie_chow_directional_rAU_face_normal_projection", [] {
        // rfn = cx*ex^2 + cy*ey^2 + cz*ez^2: only the face-normal projection
        // of the directional inverse diagonals enters the correction. The
        // internal face normal is +x, so anisotropic rAU in y/z must not
        // change the correction, while scaling the x component scales it
        // exactly.
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(0.0);
        p(0) = 1.0;
        p(1) = -1.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const auto rAU_iso = unit_rAU(mesh.n_cells());
        std::array<std::vector<double>, 3> rAU_aniso{};
        rAU_aniso[0] = std::vector<double>(mesh.n_cells(), 1.0);
        rAU_aniso[1] = std::vector<double>(mesh.n_cells(), 2.0);
        rAU_aniso[2] = std::vector<double>(mesh.n_cells(), 2.0);
        std::array<std::vector<double>, 3> rAU_x2{};
        rAU_x2[0] = std::vector<double>(mesh.n_cells(), 2.0);
        rAU_x2[1] = std::vector<double>(mesh.n_cells(), 1.0);
        rAU_x2[2] = std::vector<double>(mesh.n_cells(), 1.0);

        const double c_iso = rc_correction(
            1, mesh, geometry, U, p, rAU_iso, 1.0, ubc, pbc);
        const double c_aniso = rc_correction(
            1, mesh, geometry, U, p, rAU_aniso, 1.0, ubc, pbc);
        const double c_x2 = rc_correction(
            1, mesh, geometry, U, p, rAU_x2, 1.0, ubc, pbc);
        EXPECT_TRUE(std::abs(c_iso) > 1e-12);
        EXPECT_NEAR(c_aniso, c_iso, 1e-12);
        EXPECT_NEAR(c_x2, 2.0 * c_iso, 1e-12);
    });

    run_case("n9_s9_rhie_chow_stretched_mesh_volume_weighted_dAU", [] {
        // On a stretched mesh with unequal cell volumes and per-cell inverse
        // diagonals, the face coefficient must be the volume-weighted dAU
        // interpolation dAU_f = 0.5*(V_P*rAU_P + V_N*rAU_N) applied along the
        // owner-to-neighbour line:
        // corr = rho * dAU_f * (p_N - p_P)/d * (S_f . e).
        // The expected value is rebuilt from the same FvGeometry the
        // production path uses - no second implementation of the operator.
        const Mesh mesh = make_two_cell_channel_stretched();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(0.0);
        p(0) = 1.0;
        p(1) = -1.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const double rAU0 = 0.5;
        const double rAU1 = 1.7;
        std::array<std::vector<double>, 3> rAU{};
        for (std::size_t k = 0; k < 3; ++k) {
            rAU[k] = std::vector<double>(mesh.n_cells(), 1.0);
            rAU[k][0] = rAU0;
            rAU[k][1] = rAU1;
        }

        const double corr = rc_correction(
            1, mesh, geometry, U, p, rAU, 1.0, ubc, pbc);

        const Vec3 c0 = geometry.cell_centres[0];
        const Vec3 c1 = geometry.cell_centres[1];
        const double dx = c1.x - c0.x;
        const double dy = c1.y - c0.y;
        const double dz = c1.z - c0.z;
        const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
        const double ex = dx / d;
        const double ey = dy / d;
        const double ez = dz / d;
        const Vec3 Sf = geometry.face_area_vectors[1];
        const double Sf_dot_e = Sf.x * ex + Sf.y * ey + Sf.z * ez;
        // The mesh is orthogonal: S_non = 0, so the deferred non-orthogonal
        // part vanishes and the expected correction is exact.
        const double dAU_f = 0.5 * (geometry.cell_volumes[0] * rAU0 +
                                   geometry.cell_volumes[1] * rAU1);
        const double expected =
            1.0 * dAU_f * ((p(1) - p(0)) / d) * Sf_dot_e;

        EXPECT_TRUE(std::abs(geometry.cell_volumes[0] - 0.5) < 1e-12);
        EXPECT_TRUE(std::abs(geometry.cell_volumes[1] - 1.5) < 1e-12);
        EXPECT_NEAR(corr, expected, 1e-12);
    });

    run_case("n9_s9_rhie_chow_uniform_mesh_limit_exact_face_operator", [] {
        // On the uniform unit two-cell channel the production correction
        // must collapse to the analytic face operator
        // rho * V * rAU * (p_N - p_P)/d * (S_f . e) with V0 = V1, and the
        // geometry must be the unit construction the analytic constants
        // assume (V = d = S_f.e = 1).
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(0.0);
        p(0) = 1.0;
        p(1) = 3.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        const double corr = rc_correction(
            1, mesh, geometry, U, p, rAU, 1.0, ubc, pbc);

        const Vec3 c0 = geometry.cell_centres[0];
        const Vec3 c1 = geometry.cell_centres[1];
        const double dx = c1.x - c0.x;
        const double dy = c1.y - c0.y;
        const double dz = c1.z - c0.z;
        const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
        const double ex = dx / d;
        const double ey = dy / d;
        const double ez = dz / d;
        const Vec3 Sf = geometry.face_area_vectors[1];
        const double Sf_dot_e = Sf.x * ex + Sf.y * ey + Sf.z * ez;

        EXPECT_NEAR(geometry.cell_volumes[0], 1.0, 1e-12);
        EXPECT_NEAR(geometry.cell_volumes[1], 1.0, 1e-12);
        EXPECT_NEAR(d, 1.0, 1e-12);
        EXPECT_NEAR(Sf_dot_e, 1.0, 1e-12);

        const double dAU_f = 0.5 * (geometry.cell_volumes[0] + geometry.cell_volumes[1]);
        const double expected =
            1.0 * dAU_f * ((p(1) - p(0)) / d) * Sf_dot_e;
        EXPECT_NEAR(corr, expected, 1e-12);
    });


    run_case("n9_s9_pressure_matrix_coefficient_matches_rhie_chow_on_skew_face", [] {
        // The implicit pressure matrix must use the same centre-line
        // projection as rhie_chow_pressure_flux_internal(). On this skew face
        // d points along x while Sf has a tangential y component. The
        // directional inverse diagonal is intentionally anisotropic so that
        // projecting it on Sf (rather than d) is also detected.
        const Vec3 dvec{1.0, 0.0, 0.0};
        const Vec3 Sf{1.0, 0.4, 0.0};
        const Vec3 dAU_diagonal{1.0, 4.0, 9.0};

        const double coefficient =
            pressure_face_orthogonal_coefficient(Sf, dvec, dAU_diagonal);
        EXPECT_NEAR(coefficient, 1.0, 1e-12);

        const double area = Sf.mag();
        const Vec3 face_normal{Sf.x / area, Sf.y / area, Sf.z / area};
        const double old_face_normal_rfn =
            dAU_diagonal.x * face_normal.x * face_normal.x +
            dAU_diagonal.y * face_normal.y * face_normal.y +
            dAU_diagonal.z * face_normal.z * face_normal.z;
        const double old_area_over_distance = old_face_normal_rfn * area / dvec.mag();

        // This is deliberately not a loose numerical tolerance: the old
        // area/d plus face-normal projection differs by more than 0.5 on this
        // analytic geometry, while the corrected coefficient is exactly 1.
        EXPECT_TRUE(std::abs(old_area_over_distance - coefficient) > 0.4);

        // Orthogonal limit: retain the previous coefficient exactly even for
        // anisotropic momentum response.
        const double orthogonal = pressure_face_orthogonal_coefficient(
            Vec3{1.0, 0.0, 0.0}, dvec, dAU_diagonal);
        EXPECT_NEAR(orthogonal, 1.0, 1e-12);
    });


    run_case("n9_s9_rhie_chow_skewed_face_uses_shared_implicit_coefficient", [] {
        const Mesh mesh = make_two_cell_channel();
        auto geometry = build_fv_geometry(mesh);
        geometry.face_area_vectors[1].y = 0.4;
        // Rotate the centre-line as well: anisotropic response then couples
        // the centre-line and tangential components in global coordinates.
        geometry.cell_centres[1].y += 0.3;

        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        Field<double, Location::CELL> grad_p(2, "grad_p", "Pa/m", 3);
        p(0) = 0.0;
        p(1) = 0.6;
        // A linear tangential pressure gradient with anisotropic dAU exposes
        // whether the deferred term preserves the full face flux.
        grad_p.set(0, 0.0, 2.0, 0.0);
        grad_p.set(1, 0.0, 2.0, 0.0);

        std::array<std::vector<double>, 3> rAU{
            std::vector<double>(mesh.n_cells(), 1.0),
            std::vector<double>(mesh.n_cells(), 4.0),
            std::vector<double>(mesh.n_cells(), 9.0)};
        const Vec3 dvec = geometry.cell_centres[1] - geometry.cell_centres[0];
        const Vec3 Sf = geometry.face_area_vectors[1];
        const Vec3 dAU_face{
            0.5 * (geometry.cell_volumes[0] * rAU[0][0] +
                   geometry.cell_volumes[1] * rAU[0][1]),
            0.5 * (geometry.cell_volumes[0] * rAU[1][0] +
                   geometry.cell_volumes[1] * rAU[1][1]),
            0.5 * (geometry.cell_volumes[0] * rAU[2][0] +
                   geometry.cell_volumes[1] * rAU[2][1])};
        const double d = dvec.mag();
        const Vec3 e{dvec.x / d, dvec.y / d, dvec.z / d};
        const double orthogonal_area = Sf.dot(e);
        const Vec3 gp{0.0, 2.0, 0.0};
        const double rfn = dAU_face.x * e.x * e.x +
                           dAU_face.y * e.y * e.y +
                           dAU_face.z * e.z * e.z;
        const double implicit = pressure_face_orthogonal_coefficient(
            Sf, dvec, dAU_face) * (p(1) - p(0));
        const double full_face_flux = dAU_face.x * gp.x * Sf.x +
                                      dAU_face.y * gp.y * Sf.y +
                                      dAU_face.z * gp.z * Sf.z;
        const double deferred = full_face_flux -
                               rfn * gp.dot(e) * orthogonal_area;
        const double expected = implicit + deferred;
        EXPECT_NEAR(expected, full_face_flux, 1e-12);
        const double actual = rhie_chow_pressure_flux_internal(
            mesh, geometry, 1, p, grad_p, rAU, 1.0);
        EXPECT_NEAR(actual, expected, 1e-12);
    });

    run_case("n9_s9_skew_boundary_deferred_flux_preserves_linear_pressure", [] {
        // A prescribed boundary pressure lies along an oblique centre-to-face
        // vector. The implicit jump plus deferred remainder must reproduce the
        // full anisotropic pressure-response flux for a linear pressure field.
        const Vec3 Sf{1.0, 0.4, 0.0};
        const Vec3 dvec{1.0, 0.3, 0.0};
        const Vec3 dAU{1.0, 4.0, 9.0};
        const Vec3 grad_p{0.0, 2.0, 0.0};
        const double pressure_jump = grad_p.dot(dvec);
        const double implicit = pressure_face_orthogonal_coefficient(
            Sf, dvec, dAU) * pressure_jump;
        const double deferred = pressure_face_nonorthogonal_flux(
            Sf, dvec, dAU, grad_p);
        const double full_flux = dAU.x * grad_p.x * Sf.x +
                                 dAU.y * grad_p.y * Sf.y +
                                 dAU.z * grad_p.z * Sf.z;
        EXPECT_NEAR(implicit + deferred, full_flux, 1e-12);
        EXPECT_TRUE(std::abs(deferred) > 1.0);
    });


    run_case("n9_s9_rhie_chow_oblique_3d_face_operator_and_relaxation_scaling", [] {
        // An affine shear creates a valid 3D hexahedral mesh whose internal
        // face normal is not parallel to the owner-neighbour centre vector.
        // For an exact linear pressure field, the production face operator
        // (#923 contract: implicit orthogonal coefficient plus deferred
        // non-orthogonal flux) must return the full anisotropic directional
        // flux rho * sum_d dAU_d * grad_d(p) * Sf_d, with
        // dAU_d = 0.5*(V_P*rAU_d[P] + V_N*rAU_d[N]). The complete pressure
        // flux must then scale linearly with the relaxed inverse momentum
        // diagonal; no extra alpha_u factor belongs in Rhie-Chow.
        const Mesh mesh = make_two_cell_channel_oblique();
        const auto geometry = build_fv_geometry(mesh);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        Field<double, Location::CELL> grad_p(2, "grad_p", "Pa/m", 3);
        const Vec3 exact_gradient{2.0, -1.25, 0.75};
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            const Vec3 centre = geometry.cell_centres[c];
            p(c) = exact_gradient.dot(centre);
            grad_p.set(c, exact_gradient.x, exact_gradient.y, exact_gradient.z);
        }

        const Vec3 Sf = geometry.face_area_vectors[1];
        const Vec3 dvec = geometry.cell_centres[1] - geometry.cell_centres[0];
        const double d = dvec.mag();
        const Vec3 e{dvec.x / d, dvec.y / d, dvec.z / d};
        const double orthogonal_area = Sf.dot(e);
        const Vec3 Snon = Sf - e * orthogonal_area;
        EXPECT_TRUE(Snon.mag() > 1e-3);
        EXPECT_TRUE(std::abs(exact_gradient.dot(Snon)) > 1e-3);

        std::array<std::vector<double>, 3> rAU{
            std::vector<double>{0.7, 1.3},
            std::vector<double>{1.1, 0.6},
            std::vector<double>{0.9, 1.8}};
        constexpr double rho = 1.7;
        const double actual = rhie_chow_pressure_flux_internal(
            mesh, geometry, 1, p, grad_p, rAU, rho);

        const double cx = 0.5 * (geometry.cell_volumes[0] * rAU[0][0] +
                                 geometry.cell_volumes[1] * rAU[0][1]);
        const double cy = 0.5 * (geometry.cell_volumes[0] * rAU[1][0] +
                                 geometry.cell_volumes[1] * rAU[1][1]);
        const double cz = 0.5 * (geometry.cell_volumes[0] * rAU[2][0] +
                                 geometry.cell_volumes[1] * rAU[2][1]);
        // Oracle adapted to the #923 operator contract: the implicit
        // orthogonal coefficient plus the deferred non-orthogonal flux
        // equals the full directional tensor contraction for a linear
        // pressure field.
        const double expected = rho * (cx * exact_gradient.x * Sf.x +
                                       cy * exact_gradient.y * Sf.y +
                                       cz * exact_gradient.z * Sf.z);
        EXPECT_NEAR(actual, expected, 1e-12);

        constexpr double alpha_u = 0.4;
        const auto relaxed_rAU = scale_rAU(rAU, alpha_u);
        const double relaxed_actual = rhie_chow_pressure_flux_internal(
            mesh, geometry, 1, p, grad_p, relaxed_rAU, rho);
        EXPECT_NEAR(relaxed_actual, alpha_u * actual, 1e-12);
    });

    run_case("n9_s9_flux_velocity_mismatch_detects_known_inconsistency", [] {
        // Sensitivity proof for the restored independent audit channel:
        // phi_U (plain interpolation of the corrected velocity) and phi_RC
        // (authoritative Rhie-Chow flux) are two independent paths. For the
        // checkerboard pressure mode the channels MUST disagree by the known
        // Rhie-Chow correction (the channel is not blind); for the consistent
        // linear pressure field the internal-face mismatch MUST read zero
        // (the channel does not report an artifact when the identity holds).
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        const auto rAU = unit_rAU(mesh.n_cells());

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        U.fill(0.0);
        Field<double, Location::CELL> p_check(2, "p", "Pa", 1);
        p_check(0) = 1.0;
        p_check(1) = -1.0;
        const auto phi_rc_check = make_rhie_chow_mass_flux(
            mesh, geometry, U, p_check, rAU, 1.0, ubc, pbc);
        const auto phi_u = make_mass_flux(mesh, geometry, U, 1.0, ubc);
        const double mismatch_check = std::abs(phi_rc_check(1) - phi_u(1));
        EXPECT_TRUE(mismatch_check > 1e-12);

        Field<double, Location::CELL> p_linear(2, "p", "Pa", 1);
        p_linear(0) = 1.0;
        p_linear(1) = 3.0;
        const auto phi_rc_linear = make_rhie_chow_mass_flux(
            mesh, geometry, U, p_linear, rAU, 1.0, ubc, pbc);
        EXPECT_NEAR(phi_rc_linear(1), phi_u(1), 1e-12);
    });

    run_case("n9_s9_flux_velocity_mismatch_boundary_channel_sees_expected_correction", [] {
        // The boundary channel must not be blind either: on a fixed-pressure
        // boundary face the two reconstruction paths differ by exactly the
        // Rhie-Chow boundary correction
        // rho*(implicit_coeff*(p_bc - p_P) + deferred), evaluated with the
        // owner-cell dAU = V*rAU and the owner-centre -> face-centre
        // distance, using the same production coefficient helpers.
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        U.fill(0.0);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        p(0) = 1.0;
        p(1) = 3.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::FIXED_VALUE, 0.5, 0.0};

        const auto rAU = unit_rAU(mesh.n_cells());
        const auto phi_rc = make_rhie_chow_mass_flux(
            mesh, geometry, U, p, rAU, 1.0, ubc, pbc);
        const auto phi_u = make_mass_flux(mesh, geometry, U, 1.0, ubc);

        const std::size_t f = 0;  // x = 0 plane, owner = cell 0
        const std::size_t o = mesh.ownership().owner(f);
        const Vec3 Sf = geometry.face_area_vectors[f];
        const Vec3 dvec = geometry.face_centres[f] - geometry.cell_centres[o];
        const auto grad_p = gauss_gradient_with_boundary(p, mesh, geometry, pbc);
        const Vec3 dAU_boundary{
            geometry.cell_volumes[o] * rAU[0][o],
            geometry.cell_volumes[o] * rAU[1][o],
            geometry.cell_volumes[o] * rAU[2][o]};
        const double implicit_coefficient =
            pressure_face_orthogonal_coefficient(Sf, dvec, dAU_boundary);
        const Vec3 gp{
            grad_p.component_data(0)[o],
            grad_p.component_data(1)[o],
            grad_p.component_data(2)[o]};
        const double deferred =
            pressure_face_nonorthogonal_flux(Sf, dvec, dAU_boundary, gp);
        const double expected_correction =
            1.0 * (implicit_coefficient * (0.5 - p(o)) + deferred);

        EXPECT_TRUE(std::abs(expected_correction) > 1e-12);
        EXPECT_NEAR(phi_u(f), 0.0, 1e-12);
        EXPECT_NEAR(phi_rc(f), -expected_correction, 1e-12);
        EXPECT_NEAR(phi_rc(f) - phi_u(f), -expected_correction, 1e-12);
    });

    run_case("n9_s9_rhie_chow_identity_detects_stale_state", [] {
        // The identity channel re-derives the authoritative flux from the
        // current state: re-deriving from the same state reads exactly zero,
        // while a stored flux built from a different (stale) state MUST be
        // detected as a nonzero identity error.
        const Mesh mesh = make_two_cell_channel();
        const auto geometry = build_fv_geometry(mesh);

        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        U.fill(0.0);
        Field<double, Location::CELL> p_a(2, "p", "Pa", 1);
        Field<double, Location::CELL> p_b(2, "p", "Pa", 1);
        p_a(0) = 1.0; p_a(1) = -1.0;
        p_b(0) = 4.0; p_b(1) = -2.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        const auto rAU = unit_rAU(mesh.n_cells());

        const auto phi_stored = make_rhie_chow_mass_flux(
            mesh, geometry, U, p_a, rAU, 1.0, ubc, pbc);
        const auto phi_rederived_same = make_rhie_chow_mass_flux(
            mesh, geometry, U, p_a, rAU, 1.0, ubc, pbc);
        const auto phi_rederived_stale = make_rhie_chow_mass_flux(
            mesh, geometry, U, p_b, rAU, 1.0, ubc, pbc);

        double identity_same = 0.0;
        double identity_stale = 0.0;
        for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
            identity_same = std::max(
                identity_same, std::abs(phi_stored(f) - phi_rederived_same(f)));
            identity_stale = std::max(
                identity_stale, std::abs(phi_stored(f) - phi_rederived_stale(f)));
        }
        EXPECT_TRUE(identity_same == 0.0);
        EXPECT_TRUE(identity_stale > 1e-12);
    });

    return run_all();
}

```

### `validation/run_su2_naca0012.py`

SHA256 `bb5101afebca97162f5da1cbcba6133ac9b070b58ee267e066dab193d500eec0`

```python
#!/usr/bin/env python3
"""Qualify reading the SU2 NACA0012 Euler fixture, with explicit scope.

This fixture is not a CFDX aerodynamic validation.  Its compressible Euler
physics and JST discretization have no production path in the current
incompressible executable.  The faithful converter must keep reporting
unsupported numerical selections instead of silently replacing them.
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "tests" / "data" / "su2"
sys.path.insert(0, str(ROOT / "src" / "cfdx" / "python"))

from cfdx.io.adapters.su2 import Su2Adapter  # noqa: E402
from cfdx.io.converter import convert  # noqa: E402
from cfdx.io.gap_analysis import GapAnalysisReport  # noqa: E402
from cfdx.io.interfaces import ConversionResult, SourceInfo  # noqa: E402


def qualify_reader(directory: Path) -> dict:
    """Retain separate adapter, canonical-conversion and physics verdicts."""
    directory.mkdir(parents=True, exist_ok=True)
    shutil.copy2(DATA_DIR / "mesh_NACA0012_inv.su2", directory / "mesh.su2")
    shutil.copy2(DATA_DIR / "inv_NACA0012_basic.cfg", directory / "config.cfg")
    result = ConversionResult(source=SourceInfo())
    read_ok = Su2Adapter().convert(directory, result)
    setup = result.setup
    checks = {
        "reader_completed": bool(read_ok),
        "no_blocking_reader_finding": not result.gap_report.has_blocking(),
        "points": result.mesh is not None and result.mesh["points"].shape == (5233, 3),
        "triangles": result.mesh is not None and len(result.mesh["cells"].get("triangle", [])) == 10216,
        "physics": setup is not None and setup.physics_model == "compressible_euler",
        "energy": setup is not None and setup.energy_model == "compressible",
        "boundaries": setup is not None and {bc.patch_name for bc in setup.boundary_conditions} == {"airfoil", "farfield"},
        "pressure": setup is not None and setup.initial_condition.pressure == 101325.0,
        "temperature": setup is not None and setup.initial_condition.temperature == 288.15,
    }
    # Keep the stricter converter result visible.  A reader PASS must not
    # be interpreted as acceptance of unsupported JST or compressible flow.
    canonical = convert(directory)
    boundary_by_name = {bc.patch_name: bc for bc in setup.boundary_conditions} if setup else {}
    airfoil = boundary_by_name.get("airfoil")
    farfield = boundary_by_name.get("farfield")
    velocity = setup.initial_condition if setup else None
    physical_metadata_checks = {
        "euler_wall_slip": airfoil is not None and airfoil.value_type.value == "wall_slip",
        "initial_velocity_vector_matches_speed": velocity is not None and math.isclose(
            math.sqrt(sum(value * value for value in velocity.velocity_vector)),
            velocity.velocity, rel_tol=1e-12, abs_tol=1e-12),
        "farfield_velocity_matches_initial_speed": velocity is not None and farfield is not None and math.isclose(
            math.sqrt(sum(value * value for value in farfield.velocity_vector)) if farfield.velocity_vector
            else farfield.velocity_magnitude,
            velocity.velocity, rel_tol=1e-12, abs_tol=1e-12),
    }
    return {
        "scope": "SU2 NACA0012 input-reader qualification",
        "reader_passed": all(checks.values()),
        "reader_checks": checks,
        "physical_metadata_checks": physical_metadata_checks,
        "boundary_initial_condition_equivalence_passed": all(physical_metadata_checks.values()),
        "reader_findings": json.loads(GapAnalysisReport(result.gap_report).to_json()),
        "canonical_conversion_passed": canonical.success,
        "canonical_conversion_findings": json.loads(GapAnalysisReport(canonical.gap_analysis).to_json()),
        "production_supported": False,
        "production_reason": "The CFDX production executable supports an incompressible carrier; this fixture declares compressible Euler at Mach 0.8.",
        "physical_validation": "NOT_RUN",
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        help="directory for input copies and the reader report")
    parser.add_argument("--reader-only", action="store_true",
                        help="request only format-reader checks; permit success without a flow solve")
    args = parser.parse_args(argv)
    directory = args.output_dir or Path(tempfile.mkdtemp(prefix="cfdx_naca0012_reader_"))
    report = qualify_reader(directory)
    report_path = directory / "reader_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("=== SU2 NACA0012 input-reader qualification ===")
    print("Reader:", "PASS" if report["reader_passed"] else "FAIL")
    print("Boundary/initial-condition equivalence:",
          "PASS" if report["boundary_initial_condition_equivalence_passed"] else "FAIL")
    for check, passed in report["physical_metadata_checks"].items():
        if not passed:
            print(f"  MISMATCH: {check}")
    print("Canonical conversion:", "PASS" if report["canonical_conversion_passed"] else "BLOCKED")
    for finding in report["canonical_conversion_findings"]["findings"]:
        if finding["severity"] == "unsupported_blocking":
            print(f"  BLOCKING: {finding['feature']} - {finding['detail']}")
    print("Production: UNSUPPORTED for this compressible Euler fixture")
    print("Physical validation: NOT_RUN; no Cd, Cl or Cp prediction is asserted")
    print(f"Reader evidence: {report_path}")
    if not report["reader_passed"]:
        return 1
    if args.reader_only:
        return 0
    # Preserve a failing default verdict for legacy callers that expected
    # this command to validate the full aerodynamic run.
    print("Full aerodynamic validation is unavailable. Use --reader-only for the explicitly limited reader scope.")
    return 2


if __name__ == "__main__":
    sys.exit(main())

```

### `validation/run_su2_naca0012.sh`

SHA256 `8afc0d26fd7fad38a73b15b7b407c0fd50441d5de449674c867c8ae360fdcec0`

```bash
#!/bin/bash
# Input-reader qualification only; the Python entry point records its scope.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
exec "${PYTHON:-python3}" "$ROOT/validation/run_su2_naca0012.py" "$@"

```

## Traçabilité

Contenu extrait sans changement de l’archive `CFDX_preuves_et_correctifs_2026-10-25.zip` : `audit_output/cfdx_audit_candidate.patch` et `candidate_sources/`. Ce document fournit **le code existant** et signale les corrections **encore absentes** ; il ne les présente pas comme implémentées.
