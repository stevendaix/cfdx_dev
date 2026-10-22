# CFDX Workbench evolution

## Objective

Make the Workbench a context-driven engineering environment rather than a fixed collection of dock widgets. The viewport should remain the primary workspace; panels and actions should appear according to the task, loaded data, and solver state.

## Design principles

1. **One application model, multiple views.** Qt widgets are adapters over application/session state. They must not become independent owners of case, result, or execution data.
2. **Task-oriented workspaces.** Setup, Mesh, Run, and Analyze are named layout presets. Switching tasks changes visibility and focus, not the lifetime or identity of panels.
3. **Capability-driven UI.** An action is enabled only when its required inputs exist. For example, contour requires a loaded dataset and a scalar field; glyph requires a vector field. Explain missing prerequisites inline instead of presenting dead controls.
4. **Composable visualization.** Keep the loaded source dataset separate from derived filters. Contour, slice, clip, threshold, and glyph operations should form a small pipeline whose nodes can be reordered, enabled, edited, or removed without re-reading the source.
5. **Schema-driven engineering forms.** Physics, materials, boundaries, initialization, and numerics should be rendered from typed descriptors, units, constraints, and dependencies rather than model-specific widget code.
6. **Persistent and testable UX contracts.** Stable object names, workspace IDs, action IDs, and state transitions are public UI contracts covered by headless Qt tests and screenshot artifacts.
7. **No silent numerical changes.** UI refactoring must not alter solver settings, validation thresholds, or numerical behavior without a separate engineering change and evidence.

## Proposed information architecture

| Workspace | Primary task | Visible context |
|---|---|---|
| Setup | Configure physics and the case | Workflow, mesh context, case setup, properties |
| Mesh | Inspect geometry and boundary patches | Workflow, mesh browser, properties |
| Run | Launch and monitor execution | Workflow, setup context, properties, run center |
| Analyze | Explore result frames and derived views | Results, probes, properties, run/console context |

The workspace presets are starting points, not hard constraints: users should still be able to move docks and save a personal layout.

## Innovative capabilities to introduce incrementally

### 1. Command palette

A keyboard-first command palette (for example, Ctrl+K) should search actions and panels, expose shortcuts, and show why unavailable actions are disabled. Commands should invoke existing application actions rather than duplicating business logic.

### 2. Capability registry

Describe each UI capability with a stable ID, category, required inputs, enabled predicate, action, and optional panel. The same registry can drive menus, toolbars, command palette, contextual menus, and automated tests. This prevents separate UI surfaces from drifting out of sync.

### 3. Non-destructive visualization pipeline

Represent visualization as a source dataset plus a list/tree of filter nodes. Each node records its operation and parameters; rendering is a projection of the pipeline. Start with contour and slice, then add clip, threshold, vector glyphs, and derived quantities. Preserve source data and allow users to toggle or remove filters.

### 4. Contextual inspector

Replace the generic property form with an inspector that reflects the selected object: boundary patch, material, physics model, probe, result frame, or visualization filter. Editors should come from shared typed descriptors with units, valid ranges, validation messages, and explicit Apply/Revert behavior.

### 5. Run center as an engineering cockpit

Show solver state, iteration/time, residual and monitor trends, warnings, and output in one panel. Keep execution controls state-aware, and preserve a searchable log rather than replacing the entire output on each update.

### 6. Extensible panel and tool contracts

Panels and post-processing tools should register through stable descriptors, not require adding another permanent dock or hard-coded toolbar branch. New physics or analysis operations can then be added with localized implementation and tests.

## Current architectural observations

- `workbench.py` is the primary Workbench composition root and owns a number of dock definitions and toolbar actions.
- `gui.py` still contains the earlier `CFDXMainWindow` implementation while `launch()` delegates to `CFDXWorkbenchWindow`. Keep the legacy API working during migration, but avoid adding new product features to both shells. Plan a compatibility layer or deprecation only after import and API consumers are inventoried.
- `gui_3d.py` currently exposes imperative operations that clear and replace the plotter contents. A visualization pipeline should be introduced as a separate layer rather than embedding filter orchestration in the main window.
- Several style rules are inline in widget construction. Move toward a shared token/theme module once the workspace structure stabilizes.

## Migration sequence

1. **Contextual workspaces** — add Setup/Mesh/Run/Analyze presets, persist selection, and test dock visibility. (Current increment.)
2. **UI contracts and theme** — centralize stable action IDs, capability requirements, spacing/colors/typography, and reusable widget styles.
3. **Capability-driven actions** — derive toolbar/menu/command-palette availability from selected data and solver state.
4. **Visualization pipeline** — separate source data, filter nodes, and renderer state; add regression tests for add/update/remove operations.
5. **Schema-driven inspector** — reuse typed descriptors across setup forms and selected-object properties.
6. **Legacy shell consolidation** — inventory public imports and tests, route supported entry points to the Workbench, and deprecate duplicate shell behavior without breaking imports.

## Acceptance criteria for each increment

- Focused tests cover state transitions and visible UI contracts.
- GUI CI captures a screenshot artifact for review.
- Production solver E2E remains green when the change touches Workbench execution.
- Release/Debug and fast guards remain green.
- No tolerance inflation, skipped tests, or hidden numerical fallbacks.
