# Next-generation Workbench prototype and theme proposal

Status: interactive UX prototype for review; not production integration.
Parent architecture: existing CFDX Scientific Workbench and `Application` facade.
Scope: interaction model, workspace composition, and visual theme. No solver or numerical changes.

## Why this is not a GUI replacement

The prototype is a new presentation concept, not a second application model. The existing PySide6 Workbench remains the supported GUI while the concept is evaluated. Its application facade, CFDXSession, commands/events, capability/precondition registry, persistence, renderer adapters, and numerical workflows remain the authoritative implementation.

The target architecture is:

```
React / TypeScript workspace (candidate production frontend)
    ├── Dockable personal workspaces and contextual panels
    ├── Workflow graph and command palette
    ├── Viewport(s), synchronized comparison and probes
    └── Run cockpit, logs, terminal and resource monitors
                 │ typed commands / events
                 ▼
       Existing Python Application facade
                 │
       CFDXSession / capability registry
                 │
       persistence · solver runner · results
                 │
        C++20 solver / HDF5 / VTU
```

This HTML prototype is intentionally dependency-free so reviewers can open it directly from the repository. It is not React, does not connect to CFDX, and uses illustrative geometry, values and status. The Run button never launches a solver.

## UX ideas demonstrated

- Viewport-first workspace with minimal permanent chrome.
- Contextual floating workflow and selection inspector panels.
- Setup / Mesh / Run / Analyze workspace contexts.
- Command palette (Ctrl+K).
- Run cockpit with residual-history visual hierarchy.
- Two-view synchronized comparison concept.
- Dark and light theme switch using semantic tokens.

The workflow card is a compact stage navigator in this first slice, not yet a full React Flow node editor. Dragging panels demonstrates spatial composition; persistent docking, resizing handles, tabs, detached windows and saved layouts are production follow-up work, not claimed as complete here.

## Theme direction

The theme extends the existing proposal in `docs/gui/SCIENTIFIC_DESIGN_SYSTEM.md`, rather than inventing a competing palette.

### Design principles

1. Scientific viewport remains visually dominant.
2. High information density without dashboard clutter.
3. Semantic color tokens; no scattered literal colors in production widgets.
4. Status never communicated by color alone.
5. Dark and light are first-class themes.
6. Scientific scalar maps are separate from application status colors.
7. Keyboard focus, readable units and tabular numeric values are mandatory.

### Semantic token families

| Family | Tokens | Usage |
| --- | --- | --- |
| Surfaces | `surface.canvas`, `surface.panel`, `surface.raised`, `surface.input` | Workspace, panels, menus, editors |
| Text | `text.primary`, `text.secondary`, `text.muted` | Labels and metadata |
| Borders | `border.subtle`, `border.strong`, `focus.ring` | Separation and keyboard focus |
| Interaction | `accent.primary`, `accent.secondary`, `selection.background` | Active context, selected objects |
| Status | `status.success`, `status.warning`, `status.error`, `status.info` | Workflow and execution state |
| Data | independent scientific colormap registry | Temperature, pressure, velocity, residuals |

The prototype uses a graphite/cyan dark theme and a neutral cool light theme. The final palette must be measured against WCAG 2.2 AA contrast targets for text and controls; these prototype values are not accessibility certification. Theme changes must never modify scientific colormap range, interpolation, units, or precision.

### Theme behavior to qualify

- Persist user theme preference and support system preference.
- Keep semantic tokens identical across light/dark; only token values change.
- Audit all existing Qt widgets and screenshots for hard-coded colors.
- Validate focus visibility, disabled state, selected state, warnings and errors in both themes.
- Run automated contrast checks and manual high-DPI review.
- Keep scientific palettes perceptually appropriate and independent of UI theme.

## Migration and compatibility

Do not delete or disable existing GUI functionality to adopt this direction.

1. Review this prototype and freeze the interaction principles.
2. Inventory current GUI features and map each to an existing facade command/event or identify a missing backend contract.
3. Build a frontend shell behind an optional feature flag; initially keep PySide6 as the default.
4. Connect one vertical slice: open existing case → inspect mesh → select region → inspect properties. Use the existing facade and capability checks.
5. Add run monitoring and results browsing only after parity tests prove identical backend semantics.
6. Add docking persistence, workflow graph, terminal/resource panels and synchronized multi-view as independent increments.
7. Keep old and new GUI paths available until feature parity and regression evidence justify a change in default.

No feature should be marked available merely because a frontend control exists. Enabled state and explanatory reasons must come from backend capabilities and preconditions.

## Review checklist

- [ ] Is viewport-first spatial interaction the right direction?
- [ ] Are floating contextual panels preferable to fixed permanent sidebars?
- [ ] Does the dark/light theme fit CFDX and remain legible for long sessions?
- [ ] Which workflow stages deserve graph nodes rather than a compact stage navigator?
- [ ] Should the first production vertical slice target mesh inspection or run monitoring?
- [ ] Confirm that the existing GUI remains available throughout migration.

## Validation boundary

This PR adds only static HTML/CSS/JavaScript and documentation. It does not alter solver code, numerical algorithms, HDF5 schemas, PySide6 production widgets, or CI qualification gates. The demo data are illustrative and must not be used for engineering decisions.
