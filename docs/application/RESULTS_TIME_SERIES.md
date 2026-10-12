# Time-series results model

The GUI no longer needs to treat one VTU file as the whole result set. `discover_result_series()` deterministically discovers supported VTK-family files in a results directory and exposes frames through a headless `ResultSeries` model.

The model does not guess physical time when the filename contains no numeric timestep. Such a frame remains valid but has `time=None`. Rendering and field enumeration remain renderer responsibilities.

## Frame completeness

`ResultFrame.complete` reports whether every readability check that discovery actually performed succeeded. It is never an unverified assumption:

- without `inspect_fields`, the frame must be non-empty and carry the format banner of a supported format (`<?xml` for `.vtu`/`.vtp`/`.pvtu`, `# vtk DataFile Version` for legacy `.vtk`). This rejects empty, truncated and non-VTK files, but cannot judge a file whose body was cut after a valid header.
- with `inspect_fields`, a VTK reader is additionally asked. `pv.read` does not raise on a structurally broken file; it logs, warns and returns an empty dataset, so an empty dataset is the readability signal.

A frame that fails a check is retained in the series and flagged, not dropped, so the GUI can label it `[incomplete]` instead of silently hiding it.
