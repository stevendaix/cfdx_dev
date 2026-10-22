# CFDX Scientific Workbench — Visual Design System

Status: proposal for incremental implementation
Parent: #426
Related: #1037 (keyboard command palette)
Scope: visual language and UI contracts; no solver behavior changes

## 1. Design intent

CFDX should feel like a modern scientific instrument: calm, precise, dense where useful, and visually focused on the simulation rather than on decorative chrome. The design should be recognizable without imitating a web dashboard or a generic IDE.

The guiding idea is **Scientific Workbench**: a dark graphite workspace, a restrained cyan/teal interaction color, a high-contrast 3D viewport, and data displays that privilege legibility over ornament.

This document defines a direction, not a claim that the current GUI already implements it. Implementation must begin with an audit of existing Qt styles, object names, dock widgets, icons, and screenshot baselines.

## 2. Visual principles

1. **Viewport first.** The mesh or result field is the primary visual object. Toolbars and panels should frame it, not compete with it.
2. **Information hierarchy over decoration.** Use spacing, typography, alignment, and contrast before adding borders, gradients, or shadows.
3. **State is semantic.** Ready, running, paused, failed, stale, unsupported, and validated states must use consistent text/icon/color combinations. Color alone must never convey meaning.
4. **Stable interaction grammar.** Selection, hover, focus, disabled, warning, error, and modified states look and behave consistently across all panels.
5. **Progressive disclosure.** Show common controls directly; put advanced solver and rendering options in contextual sections without hiding their availability or support status.
6. **Composable, not hard-coded.** Theme tokens and reusable widgets should be independent of case data, solver state, and numerical algorithms.
7. **Light mode is a first-class future option.** Avoid embedding literal dark colors in widgets; semantic tokens should make a second theme possible without rewriting components.

## 3. Proposed palette

The values below are starting tokens for visual review, not final accessibility-certified values.

| Token | Proposed value | Purpose |
| --- | --- | --- |
| `surface.canvas` | `#10151B` | Main application background |
| `surface.panel` | `#171E26` | Dock and sidebar surfaces |
| `surface.raised` | `#202A35` | Menus, popovers, elevated controls |
| `surface.input` | `#0D1217` | Editable fields and terminal surfaces |
| `border.subtle` | `#2B3744` | Low-emphasis separators |
| `text.primary` | `#E7EDF3` | Primary labels and values |
| `text.secondary` | `#A7B4C2` | Supporting labels and metadata |
| `text.muted` | `#758495` | Inactive metadata; never critical values |
| `accent.primary` | `#39C6D8` | Focus, selection, active navigation |
| `accent.secondary` | `#65D6B0` | Secondary positive action / active data |
| `status.success` | `#65C995` | Completed / healthy |
| `status.warning` | `#E8B45B` | Warning / attention needed |
| `status.error` | `#F07878` | Failure / blocking diagnostic |
| `status.info` | `#7DAEFF` | Informational state |

Rules:
- Verify text and control contrast against WCAG 2.2 AA where applicable; adjust tokens based on measured contrast, not taste alone.
- Do not use the accent color for every button. Reserve the strongest emphasis for the primary action and current focus/selection.
- Do not encode CFD scalar-field palettes using the UI status palette. Scientific colormaps are a separate domain with their own perceptual and quantitative requirements.
- Use color-blind-friendly plot and status combinations; pair color with labels, icons, or line styles.

## 4. Typography and density

- Use the platform UI font for controls and a consistent tabular-numeric treatment for residuals, iteration counts, cell counts, and timings.
- Establish a small type scale: 12 px metadata, 13–14 px standard controls, 16 px section titles, and 20–24 px page/workspace headings. Respect platform DPI scaling rather than hard-coding physical sizes.
- Use compact spacing for engineering data, but keep clickable targets comfortable and keyboard focus visible.
- Align numeric columns by decimal point where practical; keep units visible and separate from editable values.
- Avoid all-caps labels except short, low-emphasis panel headers where they improve scanning.

## 5. Workbench composition

The target layout remains the existing dockable Workbench, not a new navigation rewrite:

- **Top command strip:** project identity, save/validation, run controls, and command search.
- **Workflow/navigation rail:** case tree and setup/run/results navigation with explicit status.
- **Central viewport:** mesh/result visualization with minimal overlay controls and clear selection feedback.
- **Context panel:** properties for the selected domain object; group by meaning, not by implementation class.
- **Bottom analysis region:** convergence monitors, console, logs, and diagnostics as switchable docks.
- **Status bar:** authoritative run state, mesh scale/counts, iteration/time, and the most relevant residual summary.

Keep panels dockable and hideable. Preserve user layouts, but provide a deterministic default layout and a safe fallback if saved state is incompatible. Avoid multiple simultaneous navigation systems that repeat the same hierarchy.

## 6. Interaction details

- Use one icon family and consistent stroke/size. Icons must have accessible tooltips and text alternatives for critical operations.
- Keyboard focus must be visible and must not be confused with selected-object highlighting.
- The command palette (Ctrl+K) should be the extensible entry point for discoverable actions; each command should expose its label, category, shortcut, enabled state, and reason when unavailable where the action model supports it.
- Use subtle transitions only for short state changes; respect reduced-motion preferences where available. Never animate solver progress in a way that implies numerical progress not reported by the backend.
- Destructive or irreversible actions require explicit wording and appropriate confirmation.
- Empty states should explain the next useful action; errors should explain what failed and how to recover, without inventing unsupported remediation.

## 7. Scientific visualization

- Keep the viewport background neutral and avoid decorative grids that obscure geometry.
- Make scalar bars, units, limits, and field names readable and consistent.
- Separate source datasets from derived display objects: contour, slice, glyph, and probe settings should not overwrite authoritative result fields.
- Provide clear visual distinctions between selected geometry, boundary-condition assignment, mesh-quality overlays, and result fields.
- The UI theme must not silently change scientific colormaps, interpolation, clipping, range selection, or numerical display precision.
- When data are unavailable, stale, or not qualified, communicate that state explicitly instead of showing a plausible-looking placeholder.

## 8. Theme architecture

Implement semantic tokens in one theme module/resource rather than scattering hex values across widgets. The theme layer should provide:
- palette tokens;
- typography and spacing constants;
- shared styles for buttons, inputs, tree views, tables, dock titles, tabs, menus, and status badges;
- icon sizing and semantic status mapping;
- a stable theme application entry point;
- tests that ensure the optional GUI can still be imported without Qt installed.

Qt stylesheets may be used for controls, but should not be the sole source of layout or interaction logic. Prefer reusable widget classes and shared palette roles where Qt supports them. Avoid styling by fragile widget ancestry or incidental text content.

The theme must not own or duplicate application state. It consumes presentation state and semantic status from the application layer.

## 9. Incremental delivery plan

### D0 — Audit and baseline
- Inventory existing QSS, palette overrides, widget object names, icons, and duplicated style literals.
- Capture a stable screenshot of the current Workbench at a documented viewport size and DPI.
- Identify current theme-related tests and screenshot artifact workflow.
- Record gaps before changing visuals.

### D1 — Tokens and shell
- Add centralized semantic tokens and apply them only to the main shell, docks, navigation, and status bar.
- Preserve existing behavior and dock persistence.
- Add focused Qt tests for theme application, keyboard focus visibility where testable, and optional-dependency import behavior.
- Attach before/after screenshots from the same environment.

### D2 — Common controls
- Migrate buttons, inputs, tables, tree views, menus, and dialogs to reusable styles.
- Verify disabled, selected, hover, focus, warning, and error states.
- Avoid broad find-and-replace edits across all GUI modules.

### D3 — Visualization and data density
- Refine viewport overlays, scalar bars, monitor plots, and diagnostic tables.
- Keep scientific colormaps independent from UI colors.
- Verify representative mesh, scalar, vector, and empty-result states.

### D4 — Preferences and qualification
- Add theme preference only after the centralized token layer is stable.
- Consider light mode as a separate follow-up, not as a second set of ad hoc overrides.
- Run headless Qt tests, GUI CI, and controlled screenshot review at the exact PR head.

## 10. Acceptance criteria

- One authoritative semantic token set; no new scattered hard-coded theme colors.
- Existing Workbench, keyboard actions, dock behavior, and case/run/results semantics remain unchanged.
- GUI remains an optional dependency; importing the core package does not require PySide6/PyVista.
- Theme changes do not modify solver behavior, numerical settings, scientific colormaps, or validation gates.
- Automated GUI tests cover core theme application and at least representative enabled/disabled/selected states.
- Before/after screenshots use the same viewport size and environment.
- Release/Debug and relevant GUI checks pass on the exact PR head; skipped validation-total jobs are reported as skipped, not passed.

## 11. Out of scope

- Rewriting the Workbench or its application architecture.
- Moving CFD semantics, solver state, or numerical algorithms into widgets.
- Introducing decorative animations that obscure run state.
- Changing numerical result ranges or default scientific colormaps to match the UI brand.
- Claiming visual or accessibility qualification without screenshot and contrast evidence.
