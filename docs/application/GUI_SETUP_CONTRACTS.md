# GUI setup contracts

Issue #168 extends the GUI foundation from #112 without creating a second numerical model.

The headless application layer exposes typed parameters, stable mesh selection handles, machine-readable diagnostics and case validation shared by CLI, TUI and GUI. The C++ Mesh, Patch, Field and physics structures remain authoritative; Python provides adapters and orchestration.

Validation separates blocking errors from warnings. The next GUI PR can therefore build a real mesh browser and patch editor without embedding static CFD objects in widgets.
