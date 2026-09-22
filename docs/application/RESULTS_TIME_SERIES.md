# Time-series results model

The GUI no longer needs to treat one VTU file as the whole result set. `discover_result_series()` deterministically discovers supported VTK-family files in a results directory and exposes frames through a headless `ResultSeries` model.

The model does not guess physical time when the filename contains no numeric timestep. Such a frame remains valid but has `time=None`. Rendering and field enumeration remain renderer responsibilities.
