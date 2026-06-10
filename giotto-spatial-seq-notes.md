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
- **BMKMANU S3000** (BMKGene) — see loading notes below

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

## BMKMANU S3000 (BMKGene) — loading into Giotto

No dedicated constructor exists, but loading is straightforward because the S3000 uses the
same DNB array technology and GEM/GEF output format as BGI's Stereo-seq. The proprietary
BSTMatrix pipeline produces files that are structurally identical to SAW (Stereo-seq Analysis
Workflow) outputs.

### BSTMatrix output formats

**GEF (recommended, HDF5-based)**
- Square bin: `sample.gef` — contains `/geneExp/binN/` groups for each bin resolution
- Cell bin: `sample.cellbin.gef` — contains `/cellBin/` group with per-cell expression + polygon boundaries

**GEM (plain text, tab-delimited)**

Square bin GEM columns:
```
geneID  geneName  x  y  MIDCount  ExonCount
```

Cell bin GEM columns (adds CellID):
```
geneID  geneName  x  y  MIDCount  ExonCount  CellID
```

Header lines are commented with `#` (FileFormat, BinType, BinSize, OffsetX/Y, etc.)

### Loading paths

**Path 1: GEF file — use `gefToGiotto()` (Giotto v3.2.0+, recommended)**

```r
library(Giotto)

# Cell-level data (cell segmentation was run)
bmk_go <- gefToGiotto(
  gef_file = "path/to/sample.cellbin.gef",
  bin_size = "cellbin"
)

# Square bin at a specific resolution (e.g. bin50 = 50 × 3.5 µm ≈ 175 µm bins)
bmk_go <- gefToGiotto(
  gef_file = "path/to/sample.gef",
  bin_size = "bin50"   # also: bin20, bin100, bin200 depending on what BSTMatrix produced
)
```

**Path 2: h5ad export — use `anndataToGiotto()`**

BSTViewer can export to AnnData (.h5ad). Load with:

```r
bmk_go <- anndataToGiotto("path/to/bstviewer_export.h5ad")
```

**Path 3: Manual GEM parsing — use `createGiottoObject()`**

Use this when you only have the GEM text file.

```r
library(Giotto)
library(data.table)
library(Matrix)

# Read GEM (skip # header lines)
gem <- fread("path/to/sample.cellbin.gem", skip = "geneID")

# --- Cell bin GEM (has CellID column) ---

# Aggregate UMI counts per cell × gene
gem_agg <- gem[, .(MIDCount = sum(MIDCount)), by = .(CellID, geneID)]

# Build sparse matrix (genes × cells)
cells   <- unique(gem_agg$CellID)
genes   <- unique(gem_agg$geneID)
i_idx   <- match(gem_agg$geneID, genes)
j_idx   <- match(gem_agg$CellID, cells)
count_mat <- sparseMatrix(
  i = i_idx, j = j_idx, x = gem_agg$MIDCount,
  dims = c(length(genes), length(cells)),
  dimnames = list(genes, as.character(cells))
)

# Get cell centroids from transcript positions
cell_coords <- gem[, .(sdimx = mean(x), sdimy = mean(y)), by = CellID]
setnames(cell_coords, "CellID", "cell_ID")
cell_coords[, cell_ID := as.character(cell_ID)]

bmk_go <- createGiottoObject(
  expression   = count_mat,
  spatial_locs = cell_coords
)

# --- Square bin GEM (no CellID) ---
# Each unique (x, y) pair is a spot

gem_agg <- gem[, .(MIDCount = sum(MIDCount)), by = .(x, y, geneID)]
gem_agg[, spot_ID := paste0("x", x, "y", y)]

spots <- unique(gem_agg$spot_ID)
genes <- unique(gem_agg$geneID)
count_mat <- sparseMatrix(
  i = match(gem_agg$geneID, genes),
  j = match(gem_agg$spot_ID, spots),
  x = gem_agg$MIDCount,
  dims = c(length(genes), length(spots)),
  dimnames = list(genes, spots)
)

spot_coords <- gem_agg[, .(sdimx = x[1], sdimy = y[1]), by = spot_ID]
setnames(spot_coords, "spot_ID", "cell_ID")

bmk_go <- createGiottoObject(
  expression   = count_mat,
  spatial_locs = spot_coords
)
```

### After loading — standard Giotto workflow

```r
# Build spatial network (required before spatial stats)
bmk_go <- createSpatialNetwork(bmk_go, method = "kNN", k = 6)

# Filter and normalize
bmk_go <- filterGiotto(bmk_go,
  expression_threshold = 1,
  feat_det_in_min_cells = 5,
  min_det_feats_per_cell = 100
)
bmk_go <- normalizeGiotto(bmk_go, scalefactor = 5000)

# Dimension reduction + clustering
bmk_go <- calculateHVF(bmk_go)
bmk_go <- runPCA(bmk_go)
bmk_go <- runUMAP(bmk_go, dimensions_to_use = 1:30)
bmk_go <- createNearestNetwork(bmk_go, dimensions_to_use = 1:30, k = 12)
bmk_go <- doLeidenCluster(bmk_go, resolution = 0.5)

# Spatially variable genes
bmk_go <- binSpect(bmk_go)
```

### Published precedent

- Poplar root paper (PMC11540759): BSTMatrix v2.2g → Seurat v4.3.0.1 with SCTransform
- Artemisinin paper (Nature Comms 2025): BMKMANU S3000 → BSTMatrix → downstream analysis
- ENS Lyon Spatial-Cell-ID group: GitLab repo with `CreateBmkObject.R` and `BSTMatrix_v2.3.r`
  (https://gitbio.ens-lyon.fr/spatial-cell-id/bmkmanu_s1000)

### Key notes

- `gefToGiotto()` requires Giotto v3.2.0 or newer
- BSTMatrix bin sizes mirror SAW: bin1, bin20, bin50, bin100, bin200 (unit = spot diameter × bin)
- At 3.5 µm spot diameter, bin50 ≈ 175 µm (roughly Visium-scale), bin20 ≈ 70 µm (sub-Visium)
- Cell bin data requires that cell segmentation was run in BSTMatrix (uses Cellpose internally)
- The ENS Lyon repo (S1000 but same pipeline) confirms Seurat v4.3.0 as the intended downstream tool;
  Giotto loads identically via `gefToGiotto()` since it uses the same GEF format as Stereo-seq
