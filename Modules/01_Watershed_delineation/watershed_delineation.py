# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 11:30:38 2026

@author: acer
"""
import sys
from pathlib import Path

# ─────────────────────────────────────────────────────────────
# Path Configuration (Makes script portable across computers)
# ─────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
OUTPUTS_DIR = PROJECT_ROOT / "outputs"
DEM_DIR = OUTPUTS_DIR / "dem"

# Add script directory to path so local modules can be imported
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

# Create output directories if they don't exist
DEM_DIR.mkdir(parents=True, exist_ok=True)

# Import local modules
from hydrology_tool import run_hydrology_workflow
from dem_downloader import build_config, download_dem

# ─────────────────────────────────────────────────────────────
# DEM_Download
# ─────────────────────────────────────────────────────────────

if __name__ == "__main__":

    cfg = build_config(
        west = 85.15, south = 27.5, east = 85.6, north = 27.83,
        dem_source="SRTM",      # SRTM | NASADEM | FABDEM | MERIT | ALOS | COPERNICUS
        opentopo_key="3e9fdba0886e5b6926b25485e8f111a4",
        work_dir=str(DEM_DIR),
    )

    dem = download_dem(cfg)

# import matplotlib.pyplot as plt
# import rasterio

# # Path to the downloaded DEM file
# dem_path = cfg["work_dir"] / "ALOS_DEM.tif"

# # Read the DEM data
# with rasterio.open(dem_path) as src:
#     dem_data = src.read(1)
#     transform = src.transform
#     # Get bounds for correct plotting extent
#     extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]

# # Create a figure and an axes object
# fig, ax = plt.subplots(figsize=(10, 8))

# # Display the DEM data as an image (heatmap)
# im = ax.imshow(dem_data, cmap='terrain', origin='upper', extent=extent)

# # Add a colorbar
# plt.colorbar(im, ax=ax, label='Elevation (m)')

# # Set title and labels
# ax.set_title(f'Digital Elevation Model ({cfg["dem_source"]})')
# ax.set_xlabel('Longitude')
# ax.set_ylabel('Latitude')

# # Show the plot
# plt.show()


results = run_hydrology_workflow(
    # ------------------------------------------------------------------ #
    # REQUIRED                                                             #
    # ------------------------------------------------------------------ #

    # Path to your input DEM (GeoTIFF or any GDAL-supported format)
    dem_path=str(DEM_DIR / "SRTM_DEM.tif"),

    # Outlet coordinate(s) in WGS84 (lon, lat) order
    # Single outlet:
    outlet_coords=(85.29356864039967, 27.658033127349174),
    # Multiple outlets (uncomment to use — enables Mode 3 subwatersheds):
    # outlet_coords = [
    #     (85.29356864039967, 27.658033127349174),
    #     (85.36315667973906, 27.669570166448718),  # Bishnumati Khola outlet (Bagmati confluence zone)
    #     (85.36157906231159, 27.667862066149464),  # Dhobi Khola outlet
    #     (85.35535507697145, 27.675512377431993),  # Manohara Khola outlet
    #     (85.34155984447189, 27.682171704309937),  # Hanumante Khola outlet (Bhaktapur side)
    #     (85.32831139677776, 27.691040844193257),  # Kodku Khola outlet
    #     (85.31891111402587, 27.696261798413875),  # Nakkhu Khola outlet
    #     (85.29959843664501, 27.69388498115798),   # Balkhu Khola outlet
    #     (85.29419325690164, 27.69618114752559),   # Upper Bagmati / Shivapuri outlet zone
    #     (85.29676733076668, 27.68600479969747),
    #     (85.29943396774114, 27.665478778024568),  # Tukucha Khola outlet
    #     ],
    # # OUTPUT                                                               #
    # ------------------------------------------------------------------ #

    # Root directory for all outputs
    # Creates subfolders: terrain/ hydrology/ watershed/ streams/
    #                     subwatersheds_auto/ subwatersheds_points/
    #                     statistics/ figures/ logs/ preprocessed/
    output_dir=str(OUTPUTS_DIR),
    
    # ── MODE ──────────────────────────────────────────────────────────────
    # "main"   → main watershed only
    # "auto"   → main watershed + auto subwatersheds (via stream links)
    # "points" → main watershed + user-defined sub-watersheds from extra outlets
    # "hybrid" → all three combined (default)
    mode="main",
    
    # ── STREAM EXTRACTION ─────────────────────────────────────────────────
    stream_threshold=None,               # Manual flow-acc threshold (cells). None = auto
    auto_percentile=99,                # Percentile used when threshold=None (0–100)
                                         # Lower → more/finer streams; higher → fewer/main channels

    # ── OUTLET SNAPPING ───────────────────────────────────────────────────
    snap_distance_cells=60,              # Search radius (pixels) to snap outlet to nearest stream
                                         # Increase if outlet misses the channel

    # ── DEM PREPROCESSING ─────────────────────────────────────────────────
    preprocessing_method="breach",       # "breach" (default, preferred) or "fill"
                                         # breach = least-cost depression removal (preserves topo)
                                         # fill  = raises all depressions to spill point

    # ── COORDINATE SYSTEM ─────────────────────────────────────────────────
    src_epsg=4326,                       # EPSG of your outlet_coords (default: 4326 = WGS84)
                                         # Change if your coords are in a projected CRS

    # ── OUTPUT QUALITY ────────────────────────────────────────────────────
    dpi=300,                             # Figure resolution (150 = good quality, 300 = publication)

    # ── DEBUGGING ─────────────────────────────────────────────────────────
    verbose=False,                       # True = WhiteboxTools prints progress to stdout
)