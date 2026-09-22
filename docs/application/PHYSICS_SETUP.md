# Structured physics setup

The application schema mirrors the existing C++ control families rather than inventing GUI-only physics. The first schema covers incompressible flow, energy and turbulence, including the turbulence enum used by the C++ transport controls: LAMINAR, KEPSILON, SST, SMAGORINSKY and DES.

Fields carry a type, optional unit and default. Dependencies are explicit so the GUI can hide or disable dependent controls without changing solver semantics.
