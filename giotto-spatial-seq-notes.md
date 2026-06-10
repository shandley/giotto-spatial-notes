# Giotto Suite for Spatial Transcriptomics

Giotto is a comprehensive R package (and now also has Python bindings via GiottoSuite/GiottoClass) for spatial transcriptomics analysis.

## Supported platforms

### Sequencing-based
- Visium, CytAssist Visium (10x)
- Stereo-seq, Seq-Scope, Spatial Genomics

### Imaging-based (RNA)
- Xenium (10x), Vizgen MERSCOPE
- MERFISH, seqFISH, seqFISH+, osmFISH

### Imaging-based (protein)
- Nanostring CosMx, CODEX

### Other
- Resolve Biosciences
- scRNA-seq (for integration)

## Data format support

### Technology-specific constructors

| Function | Platform | Key inputs |
|---|---|---|
| `createGiottoVisiumObject()` | Visium / CytAssist | SpaceRanger output dir, `.h5`, `.csv` coords, `.json` scalefactors, PNG image |
| `createGiottoXeniumObject()` | Xenium | Xenium output dir, `.csv` or `.parquet`, `.h5` cell-feature matrix |
| `createGiottoCosMxObject()` | CosMx | Output dir; modes: `all`, `subcellular`, or `aggregate` |
| `createGiottoObjectSubcellular()` | Generic subcellular | Transcript coords + polygon boundaries |

### Generic `createGiottoObject()` inputs

Expression matrices:
- File paths: `.csv`, `.tsv`, `.txt`
- R objects: `matrix`, `data.table`, `DelayedMatrix`
- Named list for multi-modal data

Spatial coordinates:
- File paths: `.csv`, `.tsv`, `.txt` with x/y/(z) columns
- `data.table` objects

### Format notes
- Visium: supports both SpaceRanger directory layout and direct `.h5` + coord file paths
- Xenium: `.parquet` format supported in addition to `.csv` (useful for large datasets)
- Vizgen: dedicated tutorial using its polygon/coordinate file structure
- Anything without a platform constructor can be loaded generically with a count matrix + coordinate file

## Core analysis capabilities

- Spatial variable genes (SpatialDE, SPARK, bespoke methods)
- Cell-type deconvolution (RCTD, CARD, Spotlight integration)
- Spatial co-expression networks
- Cell-proximity enrichment (neighbor interaction analysis)
- Niche clustering (combining expression + spatial context)
- Spatially-aware dimensionality reduction

## Architecture

Giotto stores data in a `giotto` S4 object with slots for expression, spatial coordinates, cell metadata, spatial networks, and multi-modal layers. Since v3/GiottoSuite, it's modularized into `GiottoClass`, `GiottoUtils`, `GiottoVisuals`.

## Things worth knowing before you start

- R-native; the Python interface exists but is less complete
- Spatial network construction (Delaunay, kNN by distance) is a prerequisite for most spatial stats
- Large datasets (Xenium, Visium HD) benefit from their image-backed (`terra`/`SpatRaster`) object approach rather than loading everything into memory
- Interop with Seurat/SCE is reasonably good via conversion functions

## Open questions

- What platform is your spatial data from?
- Starting from raw data (FASTQ/images) or already have count matrices + coordinates?
- Specific analysis goal (cell typing, spatial domains, ligand-receptor, etc.)?
- R or Python preference?
