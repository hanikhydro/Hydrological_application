"""
hydrology_tool.py
=================
Professional Python Hydrology Toolkit
======================================
Capabilities:
    - DEM preprocessing (breach, fill, resolve flats)
    - Terrain analysis (slope, aspect, hillshade, curvature, TRI, TPI)
    - Flow direction (D8) and flow accumulation
    - Stream extraction and ordering (Strahler)
    - Watershed delineation (main, auto-sub, user-defined, hybrid)
    - Morphometric analysis
    - Publication-quality visualization
    - GIS export (GeoTIFF, Shapefile, GeoPackage, GeoJSON, CSV)

Usage:
    results = run_hydrology_workflow(
        dem_path="dem.tif",
        outlet_coords=(-120.5, 37.2),   # (lon, lat) in WGS84
        output_dir="outputs/",
        mode="hybrid",
    )

Author: Hydrology Toolkit
"""

# =============================================================================
# SECTION 1: IMPORTS AND LOGGING
# =============================================================================

from __future__ import annotations

import json
import logging
import math
import os
import shutil
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import geopandas as gpd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import rasterio
import rasterio.features
import rasterio.mask
import rasterio.warp
from matplotlib.ticker import MaxNLocator
from pyproj import CRS, Transformer
from rasterio.enums import Resampling
from rasterio.transform import from_bounds, rowcol, xy
from shapely.geometry import (
    LineString,
    MultiPolygon,
    Point,
    Polygon,
    box,
    mapping,
    shape,
)
from shapely.ops import unary_union
import whitebox

matplotlib.use("Agg")
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------------------------------------------------------------------
# Logging configuration
# ---------------------------------------------------------------------------

def setup_logging(log_dir: Union[str, Path], level: int = logging.INFO) -> logging.Logger:
    """Configure and return the root toolkit logger."""
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger("hydrology")
    logger.setLevel(level)

    if logger.handlers:
        logger.handlers.clear()

    fmt = logging.Formatter(
        "%(asctime)s | %(levelname)-8s | %(funcName)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    fh = logging.FileHandler(log_dir / "hydrology.log", mode="a")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    return logger


log = logging.getLogger("hydrology")


# =============================================================================
# SECTION 2: OUTPUT DIRECTORY FACTORY
# =============================================================================

_SUBDIRS = [
    "terrain",
    "hydrology",
    "watershed",
    "subwatersheds_auto",
    "subwatersheds_points",
    "streams",
    "statistics",
    "figures",
    "logs",
    "preprocessed",
]


def create_output_dirs(output_dir: Union[str, Path]) -> Dict[str, Path]:
    root = Path(output_dir)
    dirs: Dict[str, Path] = {"root": root}
    for sub in _SUBDIRS:
        p = root / sub
        p.mkdir(parents=True, exist_ok=True)
        dirs[sub] = p
    log.info("Output directory tree created at: %s", root.resolve())
    return dirs


# =============================================================================
# SECTION 3: DEM UTILITIES
# =============================================================================

def read_dem(dem_path: Union[str, Path]) -> Tuple[np.ndarray, rasterio.profiles.Profile]:
    dem_path = Path(dem_path)
    if not dem_path.exists():
        raise FileNotFoundError(f"DEM not found: {dem_path}")

    with rasterio.open(dem_path) as src:
        profile = src.profile.copy()
        data = src.read(1).astype(np.float32)
        nodata = src.nodata

    if nodata is not None:
        data[data == nodata] = np.nan

    valid = np.sum(np.isfinite(data))
    if valid == 0:
        raise ValueError(f"DEM contains no valid pixels: {dem_path}")

    log.info(
        "DEM read: %s | shape=%s | dtype=%s | valid_pixels=%d",
        dem_path.name, data.shape, data.dtype, valid,
    )
    return data, profile


def validate_dem(
    data: np.ndarray,
    profile: rasterio.profiles.Profile,
    min_valid_fraction: float = 0.05,
) -> None:
    if data.ndim != 2:
        raise ValueError(f"Expected 2-D array; got shape {data.shape}")

    total_pixels = data.size
    valid_pixels = int(np.sum(np.isfinite(data)))
    fraction = valid_pixels / total_pixels

    if fraction < min_valid_fraction:
        raise ValueError(
            f"Only {fraction:.1%} of pixels are valid (threshold: {min_valid_fraction:.1%})"
        )

    if profile.get("crs") is None:
        raise ValueError("DEM has no CRS defined.")

    if profile.get("transform") is None or profile["transform"] == rasterio.transform.IDENTITY:
        raise ValueError("DEM has an identity or missing geotransform.")

    log.info(
        "DEM validated: valid=%.1f%% | CRS=%s",
        fraction * 100,
        profile.get("crs"),
    )


def read_metadata(dem_path: Union[str, Path]) -> Dict[str, Any]:
    dem_path = Path(dem_path)
    if not dem_path.exists():
        raise FileNotFoundError(f"DEM not found: {dem_path}")

    with rasterio.open(dem_path) as src:
        crs = src.crs
        transform = src.transform
        meta: Dict[str, Any] = {
            "crs": crs,
            "crs_wkt": crs.to_wkt() if crs else None,
            "transform": transform,
            "width": src.width,
            "height": src.height,
            "res_x": abs(transform.a),
            "res_y": abs(transform.e),
            "bounds": src.bounds,
            "nodata": src.nodata,
            "dtype": src.dtypes[0],
            "driver": src.driver,
            "is_projected": crs.is_projected if crs else False,
            "units": CRS.from_user_input(crs).axis_info[0].unit_name if crs else "unknown",
            "epsg": crs.to_epsg() if crs else None,
            "band_count": src.count,
            "file_path": str(dem_path.resolve()),
        }

    log.info(
        "Metadata: %dx%d | res=(%.4f,%.4f) | CRS=%s | projected=%s",
        meta["width"], meta["height"],
        meta["res_x"], meta["res_y"],
        meta["crs"],
        meta["is_projected"],
    )
    return meta


def determine_projected_crs(dem_path: Union[str, Path]) -> CRS:
    meta = read_metadata(dem_path)
    crs: Optional[CRS] = meta["crs"]

    if crs is None:
        raise ValueError("DEM has no CRS; cannot determine projected CRS.")

    if crs.is_projected:
        log.info("DEM already projected: %s", crs.to_string())
        return crs

    bounds = meta["bounds"]
    lon_center = (bounds.left + bounds.right) / 2.0
    lat_center = (bounds.bottom + bounds.top) / 2.0

    zone = int((lon_center + 180) / 6) + 1
    hemisphere = "north" if lat_center >= 0 else "south"
    utm_crs = CRS.from_dict({
        "proj": "utm",
        "zone": zone,
        "south": hemisphere == "south",
        "datum": "WGS84",
        "units": "m",
    })
    log.info(
        "Derived UTM CRS: zone=%d%s | EPSG=%s",
        zone, "N" if hemisphere == "north" else "S",
        utm_crs.to_epsg(),
    )
    return utm_crs


def save_raster(
    array: np.ndarray,
    profile: rasterio.profiles.Profile,
    out_path: Union[str, Path],
    nodata: float = -9999.0,
    compress: str = "lzw",
) -> Path:
    if array.ndim != 2:
        raise ValueError(f"Expected 2-D array; got shape {array.shape}")

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    write_array = array.copy().astype(np.float32)
    write_array[~np.isfinite(write_array)] = nodata

    out_profile = profile.copy()
    out_profile.update(
        dtype=rasterio.float32,
        count=1,
        nodata=nodata,
        compress=compress,
        driver="GTiff",
    )

    with rasterio.open(out_path, "w", **out_profile) as dst:
        dst.write(write_array, 1)

    log.info("Raster saved: %s | shape=%s", out_path.name, array.shape)
    return out_path.resolve()


def save_vector(
    gdf: gpd.GeoDataFrame,
    out_path: Union[str, Path],
    driver: str = "ESRI Shapefile",
) -> Path:
    if not isinstance(gdf, gpd.GeoDataFrame):
        raise TypeError(f"Expected GeoDataFrame, got {type(gdf)}")
    if gdf.empty:
        raise ValueError("Cannot save an empty GeoDataFrame.")
    if gdf.geometry.name not in gdf.columns:
        raise ValueError("GeoDataFrame has no geometry column.")

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    write_gdf = gdf.copy()

    if driver == "ESRI Shapefile":
        write_gdf.columns = [
            c[:10] if c != write_gdf.geometry.name else c
            for c in write_gdf.columns
        ]

    if driver == "ESRI Shapefile":
        _delete_shapefile_if_exists(out_path)
    write_gdf.to_file(str(out_path), driver=driver)
    
    log.info(
        "Vector saved: %s | rows=%d | driver=%s",
        out_path.name, len(write_gdf), driver,
    )
    return out_path.resolve()


def reproject_raster(
    src_path: Union[str, Path],
    dst_path: Union[str, Path],
    dst_crs: CRS,
    resampling: Resampling = Resampling.bilinear,
) -> Path:
    src_path = Path(src_path)
    if not src_path.exists():
        raise FileNotFoundError(f"Source raster not found: {src_path}")

    dst_path = Path(dst_path)
    dst_path.parent.mkdir(parents=True, exist_ok=True)

    with rasterio.open(src_path) as src:
        transform, width, height = rasterio.warp.calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds
        )
        profile = src.profile.copy()
        profile.update(
            crs=dst_crs,
            transform=transform,
            width=width,
            height=height,
            driver="GTiff",
            compress="lzw",
        )

        with rasterio.open(dst_path, "w", **profile) as dst:
            for i in range(1, src.count + 1):
                rasterio.warp.reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=resampling,
                )

    log.info("Reprojected: %s → %s | CRS=%s", src_path.name, dst_path.name, dst_crs)
    return dst_path.resolve()


# =============================================================================
# SECTION 4: WHITEBOX TOOLS INITIALISATION
# =============================================================================

def init_whitebox(
    verbose: bool = False,
    compress_rasters: bool = True,
    working_dir: Optional[Union[str, Path]] = None,
) -> whitebox.WhiteboxTools:
    wbt = whitebox.WhiteboxTools()
    wbt.verbose = verbose
    wbt.set_compress_rasters(compress_rasters)
    if working_dir is not None:
        wbt.set_working_dir(str(Path(working_dir).resolve()))
        log.info("WhiteboxTools working dir: %s", wbt.get_working_dir())
    try:
        ver = wbt.version()
        log.info("WhiteboxTools initialised: %s", ver.split("\n")[0].strip())
    except Exception as exc:
        raise RuntimeError(f"WhiteboxTools initialisation failed: {exc}") from exc
    return wbt

# =============================================================================
# SECTION 5: DEM PREPROCESSING
# =============================================================================
def breach_depressions(
    dem_path: Union[str, Path],
    out_path: Union[str, Path],
    wbt: whitebox.WhiteboxTools,
    flat_increment: Optional[float] = None,
    max_depth: Optional[float] = None,
    dist: int = 50,
) -> Path:
    dem_path = Path(dem_path).resolve()
    if not dem_path.exists():
        raise FileNotFoundError(f"DEM not found: {dem_path}")
    out_path = Path(out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    kwargs: Dict[str, Any] = {
        "dem": dem_path.name,
        "output": out_path.name,
        "dist": dist,
    }
    if flat_increment is not None:
        kwargs["flat_increment"] = flat_increment
    if max_depth is not None:
        kwargs["max_depth"] = max_depth

    wbt.set_working_dir(str(dem_path.parent))
    try:
        ret = wbt.breach_depressions_least_cost(**kwargs)
        if ret == 0 and out_path.exists():
            log.info("Depressions breached (least-cost): %s", out_path.name)
            return out_path
        log.warning("breach_depressions_least_cost silent fail; trying standard breach")
    except TypeError as exc:
        log.warning("breach_depressions_least_cost error (%s); trying standard breach", exc)

    if out_path.exists():
        out_path.unlink()

    ret2 = wbt.breach_depressions(dem=dem_path.name, output=out_path.name)
    if ret2 == 0 and out_path.exists():
        log.info("Depressions breached (standard): %s", out_path.name)
        return out_path

    raise RuntimeError(
        f"breach_depressions failed (code {ret2}) — no output written to {out_path}"
    )


def fill_depressions(
    dem_path: Union[str, Path],
    out_path: Union[str, Path],
    wbt: whitebox.WhiteboxTools,
    flat_increment: float = 0.001,
) -> Path:
    dem_path = Path(dem_path).resolve()
    if not dem_path.exists():
        raise FileNotFoundError(f"DEM not found: {dem_path}")
    out_path = Path(out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    wbt.set_working_dir(str(dem_path.parent))
    ret = wbt.fill_depressions(
        dem=dem_path.name,
        output=out_path.name,
        flat_increment=flat_increment,
    )
    if ret != 0 or not out_path.exists():
        raise RuntimeError(
            f"fill_depressions failed (code {ret}) — no output written to {out_path}"
        )
    log.info("Depressions filled: %s", out_path.name)
    return out_path


def resolve_flats(
    dem_path: Union[str, Path],
    out_path: Union[str, Path],
    wbt: whitebox.WhiteboxTools,
) -> Path:
    dem_path = Path(dem_path)
    if not dem_path.exists():
        raise FileNotFoundError(f"DEM not found: {dem_path}")

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    ret = wbt.resolve_flats(
            dem=str(dem_path),
            output=str(out_path),
        )
    if ret != 0 or not out_path.exists():
        raise RuntimeError(
            f"resolve_flats failed (code {ret}) or produced no output file."
        )

    log.info("Flats resolved: %s", out_path.name)
    return out_path.resolve()


def burn_streams(
    dem_path: Union[str, Path],
    streams_path: Union[str, Path],
    out_path: Union[str, Path],
    wbt: whitebox.WhiteboxTools,
    burn_value: float = 10.0,
) -> Path:
    dem_path = Path(dem_path)
    streams_path = Path(streams_path)
    for p in (dem_path, streams_path):
        if not p.exists():
            raise FileNotFoundError(f"File not found: {p}")

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    ret = wbt.burn_streams_at_roads(
        dem=str(dem_path),
        streams=str(streams_path),
        output=str(out_path),
        depth=burn_value,
    )
    if ret != 0:
        log.warning("burn_streams_at_roads returned %d; trying vector_stream_network_analysis", ret)
        ret2 = wbt.burn_streams(
            dem=str(dem_path),
            streams=str(streams_path),
            output=str(out_path),
            depth=burn_value,
        )
        if ret2 != 0:
            raise RuntimeError(f"burn_streams failed (code {ret2})")

    log.info("Streams burned into DEM: %s (depth=%.1fm)", out_path.name, burn_value)
    return out_path.resolve()

def _safe_copy(src: Path, dst: Path, retries: int = 8, delay: float = 0.5) -> None:
    import time
    dst.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(1, retries + 1):
        try:
            shutil.copy2(str(src), str(dst))
            return
        except PermissionError:
            if attempt == retries:
                raise
            log.debug("File locked (attempt %d/%d), retrying in %.1fs…", attempt, retries, delay)
            time.sleep(delay)
            
def _wbt_call(wbt: whitebox.WhiteboxTools, working_dir: Path, fn_name: str, **kwargs) -> bool:
    wbt.set_working_dir(str(working_dir))
    fn = getattr(wbt, fn_name, None)
    if fn is None:
        log.warning("WBT has no function '%s' — skipping.", fn_name)
        return False
    str_kwargs = {
        k: Path(v).name if isinstance(v, (str, Path)) and str(v).endswith(".tif") else v
        for k, v in kwargs.items()
    }
    try:
        ret = fn(**str_kwargs)
    except Exception as exc:
        log.warning("WBT call '%s' raised %s — skipping.", fn_name, exc)
        return False
    out_name = str_kwargs.get("output")
    if out_name:
        return ret == 0 and (working_dir / out_name).exists()
    return ret == 0

def preprocess_dem(
    dem_path: Union[str, Path],
    preprocessed_dir: Union[str, Path],
    wbt: whitebox.WhiteboxTools,
    method: str = "breach",
    flat_increment: float = 0.001,
) -> Dict[str, Path]:
    if method not in ("breach", "fill"):
        raise ValueError(f"method must be 'breach' or 'fill'; got '{method}'")

    pre_dir = Path(preprocessed_dir).resolve()
    pre_dir.mkdir(parents=True, exist_ok=True)
    dem_path = Path(dem_path).resolve()
    outputs: Dict[str, Path] = {}

    log.info("=== DEM Preprocessing (%s method) ===", method.upper())

    import time
    dem_local = pre_dir / dem_path.name
    if not dem_local.exists():
        shutil.copy2(str(dem_path), str(dem_local))

    conditioned_name: str

    if method == "breach":
        ok = _wbt_call(wbt, pre_dir, "breach_depressions_least_cost",
                       dem=dem_local.name, output="dem_breached.tif", dist=50,
                       flat_increment=flat_increment)
        if not ok:
            log.warning("breach_depressions_least_cost failed — trying breach_depressions")
            ok = _wbt_call(wbt, pre_dir, "breach_depressions",
                           dem=dem_local.name, output="dem_breached.tif",
                           flat_increment=flat_increment)
        if not ok:
            log.warning("All breach methods failed — falling back to fill_depressions")
            ok = _wbt_call(wbt, pre_dir, "fill_depressions",
                           dem=dem_local.name, output="dem_filled.tif",
                           flat_increment=flat_increment)
            if not ok:
                raise RuntimeError("fill_depressions also failed — check WBT installation.")
            outputs["filled"] = pre_dir / "dem_filled.tif"
            conditioned_name = "dem_filled.tif"
        else:
            outputs["breached"] = pre_dir / "dem_breached.tif"
            conditioned_name = "dem_breached.tif"
    else:
        ok = _wbt_call(wbt, pre_dir, "fill_depressions",
                       dem=dem_local.name, output="dem_filled.tif",
                       flat_increment=flat_increment)
        if not ok:
            raise RuntimeError("fill_depressions failed — check WBT installation.")
        outputs["filled"] = pre_dir / "dem_filled.tif"
        conditioned_name = "dem_filled.tif"

    ok_rf = _wbt_call(wbt, pre_dir, "resolve_flats",
                      dem=conditioned_name, output="dem_conditioned.tif")
    ok_rf = _wbt_call(wbt, pre_dir, "resolve_flats",
                      dem=conditioned_name, output="dem_conditioned.tif")
    if not ok_rf:
        ok_rf = _wbt_call(wbt, pre_dir, "fix_flats",
                          dem=conditioned_name, output="dem_conditioned.tif")
    if ok_rf:
        log.info("Flats resolved: dem_conditioned.tif")
        outputs["conditioned"] = pre_dir / "dem_conditioned.tif"
    else:
        log.warning("resolve_flats unavailable — using %s as conditioned DEM.", conditioned_name)
        dst = pre_dir / "dem_conditioned.tif"
        if not dst.exists():
            shutil.copy2(str(pre_dir / conditioned_name), str(dst))
        outputs["conditioned"] = dst
    log.info("Preprocessing complete. Conditioned DEM: %s", outputs["conditioned"].name)
    return outputs
        
# =============================================================================
# SECTION 6: TERRAIN ANALYSIS
# =============================================================================

def compute_slope(dem_path, out_path, wbt, units="degrees"):
    dem_path = Path(dem_path).resolve()
    if not dem_path.exists():
        raise FileNotFoundError(f"DEM not found: {dem_path}")
    out_path = Path(out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    dem_local = out_path.parent / dem_path.name
    if not dem_local.exists():
        shutil.copy2(str(dem_path), str(dem_local))
    wbt.set_working_dir(str(out_path.parent))
    ret = wbt.slope(dem=dem_local.name, output=out_path.name, units=units)
    if ret != 0 or not out_path.exists():
        raise RuntimeError(f"slope failed (code {ret})")
    log.info("Slope computed (%s): %s", units, out_path.name)
    return out_path

def compute_aspect(dem_path, out_path, wbt):
    dem_path = Path(dem_path).resolve()
    if not dem_path.exists():
        raise FileNotFoundError(f"DEM not found: {dem_path}")
    out_path = Path(out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    dem_local = out_path.parent / dem_path.name
    if not dem_local.exists():
        shutil.copy2(str(dem_path), str(dem_local))
    wbt.set_working_dir(str(out_path.parent))
    ret = wbt.aspect(dem=dem_local.name, output=out_path.name)
    if ret != 0 or not out_path.exists():
        raise RuntimeError(f"aspect failed (code {ret})")
    log.info("Aspect computed: %s", out_path.name)
    return out_path

def compute_hillshade(dem_path, out_path, wbt, azimuth=315.0, altitude=45.0):
    dem_path = Path(dem_path).resolve()
    if not dem_path.exists():
        raise FileNotFoundError(f"DEM not found: {dem_path}")
    out_path = Path(out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    dem_local = out_path.parent / dem_path.name
    if not dem_local.exists():
        shutil.copy2(str(dem_path), str(dem_local))
    wbt.set_working_dir(str(out_path.parent))
    ret = wbt.hillshade(dem=dem_local.name, output=out_path.name,
                        azimuth=azimuth, altitude=altitude)
    if ret != 0 or not out_path.exists():
        raise RuntimeError(f"hillshade failed (code {ret})")
    log.info("Hillshade computed (az=%.0f, alt=%.0f): %s", azimuth, altitude, out_path.name)
    return out_path

def compute_curvature(dem_path, out_path, wbt, kind="profile"):
    dem_path = Path(dem_path).resolve()
    if not dem_path.exists():
        raise FileNotFoundError(f"DEM not found: {dem_path}")
    out_path = Path(out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    dem_local = out_path.parent / dem_path.name
    if not dem_local.exists():
        shutil.copy2(str(dem_path), str(dem_local))
    wbt.set_working_dir(str(out_path.parent))
    _KINDS = {
        "profile": wbt.profile_curvature,
        "plan": wbt.plan_curvature,
        "tangential": wbt.tangential_curvature,
    }
    if kind not in _KINDS:
        raise ValueError(f"kind must be one of {list(_KINDS)}; got '{kind}'")
    ret = _KINDS[kind](dem=dem_local.name, output=out_path.name)
    if ret != 0 or not out_path.exists():
        raise RuntimeError(f"{kind}_curvature failed (code {ret})")
    log.info("%s curvature computed: %s", kind.capitalize(), out_path.name)
    return out_path


def compute_tri(
    dem_path: Union[str, Path],
    out_path: Union[str, Path],
    wbt: whitebox.WhiteboxTools,
) -> Path:
    dem_path = Path(dem_path)
    if not dem_path.exists():
        raise FileNotFoundError(f"DEM not found: {dem_path}")
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    for fn_name in ("ruggedness_index", "tri", "terrain_ruggedness_index"):
        fn = getattr(wbt, fn_name, None)
        if fn is None:
            continue
        try:
            ret = fn(dem=str(dem_path), output=str(out_path))
            if ret == 0 and out_path.exists():
                log.info("TRI computed (%s): %s", fn_name, out_path.name)
                return out_path.resolve()
        except Exception as exc:
            log.debug("TRI attempt '%s' failed: %s", fn_name, exc)

    raise RuntimeError("No working TRI function found in this WBT version.")


def compute_tpi(
    dem_path: Union[str, Path],
    out_path: Union[str, Path],
    wbt: whitebox.WhiteboxTools,
) -> Path:
    dem_path = Path(dem_path)
    if not dem_path.exists():
        raise FileNotFoundError(f"DEM not found: {dem_path}")
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    for fn_name in ("deviation_from_mean_elevation", "tpi", "topographic_position_index"):
        fn = getattr(wbt, fn_name, None)
        if fn is None:
            continue
        try:
            ret = fn(dem=str(dem_path), output=str(out_path))
            if ret == 0 and out_path.exists():
                log.info("TPI computed (%s): %s", fn_name, out_path.name)
                return out_path.resolve()
        except Exception as exc:
            log.debug("TPI attempt '%s' failed: %s", fn_name, exc)

    raise RuntimeError("No working TPI function found in this WBT version.")
    

def run_terrain_analysis(
    conditioned_dem: Union[str, Path],
    terrain_dir: Union[str, Path],
    wbt: whitebox.WhiteboxTools,
) -> Dict[str, Path]:
    t_dir = Path(terrain_dir)
    dem = Path(conditioned_dem)
    log.info("=== Terrain Analysis ===")

    outputs: Dict[str, Path] = {}
    outputs["slope"] = compute_slope(dem, t_dir / "slope.tif", wbt)
    outputs["aspect"] = compute_aspect(dem, t_dir / "aspect.tif", wbt)
    outputs["hillshade"] = compute_hillshade(dem, t_dir / "hillshade.tif", wbt)

    try:
        outputs["curvature"] = compute_curvature(dem, t_dir / "curvature_profile.tif", wbt, kind="profile")
    except RuntimeError as exc:
        log.warning("Curvature computation skipped: %s", exc)

    try:
        outputs["tri"] = compute_tri(dem, t_dir / "tri.tif", wbt)
    except RuntimeError as exc:
        log.warning("TRI computation skipped: %s", exc)

    try:
        outputs["tpi"] = compute_tpi(dem, t_dir / "tpi.tif", wbt)
    except RuntimeError as exc:
        log.warning("TPI computation skipped: %s", exc)

    log.info("Terrain analysis complete: %d products", len(outputs))
    return outputs


# =============================================================================
# SECTION 7: FLOW DIRECTION AND ACCUMULATION
# =============================================================================

def compute_flow_direction(
    conditioned_dem: Union[str, Path],
    out_path: Union[str, Path],
    wbt: whitebox.WhiteboxTools,
    method: str = "d8",
) -> Path:
    if method != "d8":
        raise ValueError(f"Only 'd8' method is currently supported; got '{method}'")

    dem = Path(conditioned_dem).resolve()
    if not dem.exists():
        raise FileNotFoundError(f"Conditioned DEM not found: {dem}")

    out_path = Path(out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    dem_local = out_path.parent / dem.name
    if not dem_local.exists():
        shutil.copy2(str(dem), str(dem_local))

    wbt.set_working_dir(str(out_path.parent))
    ret = wbt.d8_pointer(dem=dem_local.name, output=out_path.name)
    if ret != 0 or not out_path.exists():
        raise RuntimeError(f"d8_pointer failed (code {ret}) — no output at {out_path}")

    log.info("D8 flow direction computed: %s", out_path.name)
    return out_path


def compute_flow_accumulation(
    flow_dir_path: Union[str, Path],
    out_path: Union[str, Path],
    wbt: whitebox.WhiteboxTools,
    log_transform: bool = False,
) -> Path:
    fdr = Path(flow_dir_path).resolve()
    if not fdr.exists():
        raise FileNotFoundError(f"Flow direction raster not found: {fdr}")

    out_path = Path(out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    wbt.set_working_dir(str(out_path.parent))
    ret = wbt.d8_flow_accumulation(
        i=fdr.name,
        output=out_path.name,
        out_type="cells",
        pntr=True,
    )
    if ret != 0 or not out_path.exists():
        raise RuntimeError(f"d8_flow_accumulation failed (code {ret}) — no output at {out_path}")

    if log_transform:
        log_path = out_path.with_stem(out_path.stem + "_log")
        with rasterio.open(out_path) as src:
            profile = src.profile.copy()
            data = src.read(1).astype(np.float32)
            nodata = src.nodata or -9999.0
            mask = data != nodata
            log_data = np.where(mask, np.log10(np.maximum(data, 1)), nodata)
        profile.update(dtype=rasterio.float32, nodata=nodata)
        with rasterio.open(log_path, "w", **profile) as dst:
            dst.write(log_data.astype(np.float32), 1)
        log.info("Log10 flow accumulation saved: %s", log_path.name)

    log.info("Flow accumulation computed: %s", out_path.name)
    return out_path


def compute_flow_length(
    flow_dir_path: Union[str, Path],
    out_path: Union[str, Path],
    wbt: whitebox.WhiteboxTools,
    direction: str = "downstream",
) -> Path:
    if direction not in ("downstream", "upstream"):
        raise ValueError(f"direction must be 'downstream' or 'upstream'; got '{direction}'")

    fdr = Path(flow_dir_path)
    if not fdr.exists():
        raise FileNotFoundError(f"Flow direction raster not found: {fdr}")

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    ret = wbt.flow_length_diff(
        d8_pntr=str(fdr),
        output=str(out_path),
    )
    if ret != 0:
        log.warning("flow_length_diff returned %d; trying flow_accumulation proxy", ret)
        ret = wbt.d8_flow_accumulation(
            i=str(fdr),
            output=str(out_path),
            out_type="specific contributing area",
            pntr=True,
        )
        if ret != 0:
            raise RuntimeError(f"Flow length computation failed (code {ret})")

    log.info("Flow length (%s) computed: %s", direction, out_path.name)
    return out_path.resolve()


def run_flow_analysis(
    conditioned_dem: Union[str, Path],
    hydrology_dir: Union[str, Path],
    wbt: whitebox.WhiteboxTools,
) -> Dict[str, Path]:
    h_dir = Path(hydrology_dir)
    dem = Path(conditioned_dem)
    log.info("=== Flow Analysis ===")

    flow_dir = compute_flow_direction(dem, h_dir / "flow_direction.tif", wbt)
    flow_acc = compute_flow_accumulation(flow_dir, h_dir / "flow_accumulation.tif", wbt, log_transform=True)

    outputs: Dict[str, Path] = {
        "flow_dir": flow_dir,
        "flow_acc": flow_acc,
    }

    try:
        fl = compute_flow_length(flow_dir, h_dir / "flow_length.tif", wbt)
        outputs["flow_length"] = fl
    except RuntimeError as exc:
        log.warning("Flow length skipped: %s", exc)

    log.info("Flow analysis complete.")
    return outputs


# =============================================================================
# SECTION 8: STREAM EXTRACTION
# =============================================================================

def estimate_stream_threshold(
    flow_acc_path: Union[str, Path],
    percentile: float = 99.5,
) -> float:
    if not 0 <= percentile <= 100:
        raise ValueError(f"percentile must be in [0, 100]; got {percentile}")

    flow_acc_path = Path(flow_acc_path)
    if not flow_acc_path.exists():
        raise FileNotFoundError(f"Flow accumulation raster not found: {flow_acc_path}")

    with rasterio.open(flow_acc_path) as src:
        data = src.read(1).astype(np.float32)
        nodata = src.nodata

    if nodata is not None:
        data = data[data != nodata]
    data = data[np.isfinite(data) & (data > 0)]

    if data.size == 0:
        raise ValueError("Flow accumulation raster contains no positive values.")

    threshold = float(np.percentile(data, percentile))
    log.info(
        "Auto stream threshold: %.0f cells (p%.1f of %d values)",
        threshold, percentile, data.size,
    )
    return threshold


def extract_stream_raster(
    flow_acc_path: Union[str, Path],
    out_path: Union[str, Path],
    wbt: whitebox.WhiteboxTools,
    threshold: Optional[float] = None,
    auto_percentile: float = 99.5,
) -> Tuple[Path, float]:
    flow_acc_path = Path(flow_acc_path).resolve()
    if not flow_acc_path.exists():
        raise FileNotFoundError(f"Flow accumulation raster not found: {flow_acc_path}")
    if threshold is None:
        threshold = estimate_stream_threshold(flow_acc_path, percentile=auto_percentile)
    out_path = Path(out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fa_local = out_path.parent / flow_acc_path.name
    if not fa_local.exists():
        shutil.copy2(str(flow_acc_path), str(fa_local))

    wbt.set_working_dir(str(out_path.parent))
    ret = wbt.extract_streams(
        flow_accum=fa_local.name,
        output=out_path.name,
        threshold=threshold,
    )
    if ret != 0 or not out_path.exists():
        raise RuntimeError(f"extract_streams failed (code {ret})")
    log.info("Stream raster extracted (threshold=%.0f): %s", threshold, out_path.name)
    return out_path, threshold


def extract_stream_vector(
    stream_raster_path: Union[str, Path],
    flow_dir_path: Union[str, Path],
    out_path: Union[str, Path],
    wbt: whitebox.WhiteboxTools,
) -> Path:
    stream_raster_path = Path(stream_raster_path).resolve()
    flow_dir_path = Path(flow_dir_path).resolve()
    for p in (stream_raster_path, flow_dir_path):
        if not p.exists():
            raise FileNotFoundError(f"File not found: {p}")
    out_path = Path(out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    for src in (stream_raster_path, flow_dir_path):
        dst = out_path.parent / src.name
        if not dst.exists():
            shutil.copy2(str(src), str(dst))

    wbt.set_working_dir(str(out_path.parent))
    ret = wbt.raster_streams_to_vector(
        streams=stream_raster_path.name,
        d8_pntr=flow_dir_path.name,
        output=out_path.name,
    )
    if ret != 0 or not out_path.exists():
        raise RuntimeError(f"raster_streams_to_vector failed (code {ret})")
    log.info("Stream vector exported: %s", out_path.name)
    return out_path


def generate_stream_links(
    stream_raster_path: Union[str, Path],
    flow_dir_path: Union[str, Path],
    out_path: Union[str, Path],
    wbt: whitebox.WhiteboxTools,
) -> Path:
    stream_raster_path = Path(stream_raster_path).resolve()
    flow_dir_path = Path(flow_dir_path).resolve()
    for p in (stream_raster_path, flow_dir_path):
        if not p.exists():
            raise FileNotFoundError(f"File not found: {p}")
    out_path = Path(out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    for src in (stream_raster_path, flow_dir_path):
        dst = out_path.parent / src.name
        if not dst.exists():
            shutil.copy2(str(src), str(dst))

    wbt.set_working_dir(str(out_path.parent))
    ret = wbt.stream_link_identifier(
        d8_pntr=flow_dir_path.name,
        streams=stream_raster_path.name,
        output=out_path.name,
    )
    if ret != 0 or not out_path.exists():
        raise RuntimeError(f"stream_link_identifier failed (code {ret})")
    log.info("Stream links generated: %s", out_path.name)
    return out_path


def compute_strahler_order(
    stream_raster_path: Union[str, Path],
    flow_dir_path: Union[str, Path],
    out_path: Union[str, Path],
    wbt: whitebox.WhiteboxTools,
) -> Path:
    stream_raster_path = Path(stream_raster_path).resolve()
    flow_dir_path = Path(flow_dir_path).resolve()
    for p in (stream_raster_path, flow_dir_path):
        if not p.exists():
            raise FileNotFoundError(f"File not found: {p}")
    out_path = Path(out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    for src in (stream_raster_path, flow_dir_path):
        dst = out_path.parent / src.name
        if not dst.exists():
            shutil.copy2(str(src), str(dst))

    wbt.set_working_dir(str(out_path.parent))
    ret = wbt.strahler_stream_order(
        d8_pntr=flow_dir_path.name,
        streams=stream_raster_path.name,
        output=out_path.name,
    )
    if ret != 0 or not out_path.exists():
        raise RuntimeError(f"strahler_stream_order failed (code {ret})")
    log.info("Strahler order raster computed: %s", out_path.name)
    return out_path


def build_ordered_stream_geodataframe(
    stream_vector_path: Union[str, Path],
    strahler_raster_path: Union[str, Path],
    flow_acc_path: Union[str, Path],
) -> gpd.GeoDataFrame:
    for p in (stream_vector_path, strahler_raster_path, flow_acc_path):
        p = Path(p)
        if not p.exists():
            raise FileNotFoundError(f"File not found: {p}")

    gdf = gpd.read_file(str(stream_vector_path))
    if gdf.empty:
        raise ValueError("Stream vector is empty.")

    crs = gdf.crs
    if crs and not crs.is_projected:
        gdf_proj = gdf.to_crs(epsg=3857)
        gdf["length_m"] = gdf_proj.geometry.length
    else:
        gdf["length_m"] = gdf.geometry.length

    midpoints = gdf.geometry.interpolate(0.5, normalized=True)
    coords = [(p.x, p.y) for p in midpoints]

    def _sample_raster(raster_path: Path, pts: list) -> np.ndarray:
        with rasterio.open(raster_path) as src:
            vals = np.array(
                [v[0] for v in src.sample(pts)],
                dtype=np.float32,
            )
            nodata = src.nodata
        if nodata is not None:
            vals[vals == nodata] = np.nan
        return vals

    strahler_vals = _sample_raster(Path(strahler_raster_path), coords)
    gdf["strahler"] = pd.array(strahler_vals, dtype="Int64")
    
    gdf["upstream_cells"] = _sample_raster(Path(flow_acc_path), coords)
    gdf["stream_id"] = np.arange(1, len(gdf) + 1)

    cols = ["stream_id", "strahler", "upstream_cells", "length_m", "geometry"]
    extra = [c for c in gdf.columns if c not in cols]
    gdf = gdf[extra + cols] if extra else gdf[cols]

    log.info(
        "Ordered stream GDF: %d segments | max_order=%s",
        len(gdf), gdf["strahler"].max(),
    )
    return gdf


# =============================================================================
# SECTION 10: OUTLET PROCESSING
# =============================================================================

def reproject_outlet(
    lon: float,
    lat: float,
    dst_crs: CRS,
    src_epsg: int = 4326,
) -> Tuple[float, float]:
    transformer = Transformer.from_crs(
        CRS.from_epsg(src_epsg), dst_crs, always_xy=True
    )
    x, y = transformer.transform(lon, lat)
    log.info("Outlet reprojected: (%.6f, %.6f) → (%.2f, %.2f)", lon, lat, x, y)
    return x, y


def validate_outlet_within_dem(
    x: float,
    y: float,
    dem_path: Union[str, Path],
) -> None:
    dem_path = Path(dem_path)
    if not dem_path.exists():
        raise FileNotFoundError(f"DEM not found: {dem_path}")

    with rasterio.open(dem_path) as src:
        bounds = src.bounds

    if not (bounds.left <= x <= bounds.right and bounds.bottom <= y <= bounds.top):
        raise ValueError(
            f"Outlet ({x:.2f}, {y:.2f}) is outside DEM extent "
            f"[{bounds.left:.2f}, {bounds.right:.2f}] x "
            f"[{bounds.bottom:.2f}, {bounds.top:.2f}]"
        )
    log.info("Outlet validated within DEM extent.")


def snap_outlet_to_stream(
    x: float,
    y: float,
    flow_acc_path: Union[str, Path],
    snap_distance_cells: int = 50,
) -> Tuple[float, float]:
    flow_acc_path = Path(flow_acc_path)
    if not flow_acc_path.exists():
        raise FileNotFoundError(f"Flow accumulation raster not found: {flow_acc_path}")

    with rasterio.open(flow_acc_path) as src:
        transform = src.transform
        nodata = src.nodata

        row, col = rowcol(transform, x, y)
        r0 = max(0, row - snap_distance_cells)
        r1 = min(src.height, row + snap_distance_cells + 1)
        c0 = max(0, col - snap_distance_cells)
        c1 = min(src.width, col + snap_distance_cells + 1)

        if r0 >= r1 or c0 >= c1:
            raise ValueError("Outlet is entirely outside the flow accumulation raster.")

        window = rasterio.windows.Window(c0, r0, c1 - c0, r1 - r0)
        chunk = src.read(1, window=window).astype(np.float32)

        if nodata is not None:
            chunk[chunk == nodata] = np.nan

        local_max = np.nanargmax(chunk)
        lrow, lcol = np.unravel_index(local_max, chunk.shape)
        abs_row = r0 + lrow
        abs_col = c0 + lcol

        snapped_x, snapped_y = xy(transform, abs_row, abs_col)

    log.info(
        "Outlet snapped: (%.2f, %.2f) → (%.2f, %.2f) [%d cell window]",
        x, y, snapped_x, snapped_y, snap_distance_cells,
    )
    return float(snapped_x), float(snapped_y)


def create_outlet_shapefile(
    outlets: List[Tuple[float, float]],
    crs: CRS,
    out_path: Union[str, Path],
    ids: Optional[List[int]] = None,
) -> Path:
    if not outlets:
        raise ValueError("outlets list is empty.")

    ids = ids or list(range(1, len(outlets) + 1))
    geometries = [Point(x, y) for x, y in outlets]
    gdf = gpd.GeoDataFrame(
        {"outlet_id": ids, "geometry": geometries},
        crs=crs,
    )

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    driver = "GPKG" if out_path.suffix.lower() == ".gpkg" else "ESRI Shapefile"
    gdf.to_file(str(out_path), driver=driver)
    log.info("Outlet shapefile saved: %s (%d points)", out_path.name, len(gdf))
    return out_path.resolve()


def process_outlets(
    outlet_coords: Union[Tuple[float, float], List[Tuple[float, float]]],
    dem_path: Union[str, Path],
    flow_acc_path: Union[str, Path],
    dst_crs: CRS,
    output_dir: Union[str, Path],
    snap_distance_cells: int = 50,
    src_epsg: int = 4326,
) -> Dict[str, Any]:
    if isinstance(outlet_coords, tuple) and len(outlet_coords) == 2:
        outlet_coords = [outlet_coords]

    snapped: List[Tuple[float, float]] = []
    for lon, lat in outlet_coords:
        x, y = reproject_outlet(lon, lat, dst_crs, src_epsg)
        validate_outlet_within_dem(x, y, dem_path)
        sx, sy = snap_outlet_to_stream(x, y, flow_acc_path, snap_distance_cells)
        snapped.append((sx, sy))

    shp_path = create_outlet_shapefile(
        outlets=snapped,
        crs=dst_crs,
        out_path=Path(output_dir) / "outlets_snapped.shp",
    )

    return {
        "snapped_outlets": snapped,
        "outlet_shp": shp_path,
        "primary_outlet": snapped[0],
    }


# =============================================================================
# SECTION 11: WATERSHED DELINEATION
# =============================================================================

def _create_single_outlet_shp(
    x: float,
    y: float,
    crs: CRS,
    out_path: Union[str, Path],
) -> Path:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    gdf = gpd.GeoDataFrame({"id": [1], "geometry": [Point(x, y)]}, crs=crs)
    gdf.to_file(str(out_path))
    return out_path.resolve()


def delineate_watershed(
    flow_dir_path: Union[str, Path],
    outlet_x: float,
    outlet_y: float,
    crs: CRS,
    out_raster_path: Union[str, Path],
    wbt: whitebox.WhiteboxTools,
) -> Path:
    fdr = Path(flow_dir_path).resolve()
    if not fdr.exists():
        raise FileNotFoundError(f"Flow direction raster not found: {fdr}")
    out_raster_path = Path(out_raster_path).resolve()
    out_raster_path.parent.mkdir(parents=True, exist_ok=True)

    tmp_shp = out_raster_path.parent / "_tmp_outlet.shp"
    _create_single_outlet_shp(outlet_x, outlet_y, crs, tmp_shp)

    fdr_local = out_raster_path.parent / fdr.name
    if not fdr_local.exists():
        shutil.copy2(str(fdr), str(fdr_local))

    wbt.set_working_dir(str(out_raster_path.parent))
    ret = wbt.watershed(
        d8_pntr=fdr_local.name,
        pour_pts=tmp_shp.name,
        output=out_raster_path.name,
    )
    if ret != 0 or not out_raster_path.exists():
        raise RuntimeError(f"watershed failed (code {ret})")

    for ext in (".shp", ".shx", ".dbf", ".prj", ".cpg"):
        f = tmp_shp.with_suffix(ext)
        if f.exists():
            f.unlink()

    log.info("Watershed raster delineated: %s", out_raster_path.name)
    return out_raster_path


def delineate_subwatersheds_auto(
    flow_dir_path: Union[str, Path],
    stream_links_path: Union[str, Path],
    out_raster_path: Union[str, Path],
    wbt: whitebox.WhiteboxTools,
) -> Path:
    flow_dir_path = Path(flow_dir_path).resolve()
    stream_links_path = Path(stream_links_path).resolve()
    for p in (flow_dir_path, stream_links_path):
        if not p.exists():
            raise FileNotFoundError(f"File not found: {p}")
    out_raster_path = Path(out_raster_path).resolve()
    out_raster_path.parent.mkdir(parents=True, exist_ok=True)

    for src in (flow_dir_path, stream_links_path):
        dst = out_raster_path.parent / src.name
        if not dst.exists():
            shutil.copy2(str(src), str(dst))

    wbt.set_working_dir(str(out_raster_path.parent))
    ret = wbt.subbasins(
        d8_pntr=flow_dir_path.name,
        streams=stream_links_path.name,
        output=out_raster_path.name,
    )
    if ret != 0 or not out_raster_path.exists():
        raise RuntimeError(f"subbasins failed (code {ret})")
    log.info("Auto subwatersheds raster generated: %s", out_raster_path.name)
    return out_raster_path


def delineate_subwatersheds_points(
    flow_dir_path: Union[str, Path],
    outlets: List[Tuple[float, float]],
    crs: CRS,
    out_dir: Union[str, Path],
    wbt: whitebox.WhiteboxTools,
) -> Path:
    """Delineate sub-watersheds from user-supplied outlet points (Mode 3).

    Each outlet gets its own watershed raster, then they are combined into a
    single labelled raster (outlet index = pixel value) so the downstream
    pipeline can treat it identically to the auto-subbasins raster.  All
    sub-watersheds are contained within the main watershed because each
    outlet was already snapped to a stream inside it.

    Args:
        flow_dir_path: D8 flow direction raster path.
        outlets: List of (x, y) snapped outlet tuples in ``crs``.
        crs: CRS of the outlet coordinates.
        out_dir: Directory for output rasters.
        wbt: Initialised WhiteboxTools instance.

    Returns:
        Path to the combined labelled raster
        (``subwatersheds_points/subwatersheds_points.tif``).

    Raises:
        FileNotFoundError: If ``flow_dir_path`` does not exist.
        ValueError: If ``outlets`` is empty.
    """
    if not outlets:
        raise ValueError("No outlet points provided.")

    fdr = Path(flow_dir_path)
    if not fdr.exists():
        raise FileNotFoundError(f"Flow direction raster not found: {fdr}")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Step 1: delineate one watershed raster per outlet (skip outlet 0,
    #         which is the primary / main-watershed outlet)
    # ------------------------------------------------------------------
    individual_rasters: List[Path] = []
    
    for i, (x, y) in enumerate(outlets, start=1):
        rpath = out_dir / f"_ws_outlet_{i:03d}.tif"
        try:
            delineate_watershed(
                flow_dir_path=fdr,
                outlet_x=x,
                outlet_y=y,
                crs=crs,
                out_raster_path=rpath,
                wbt=wbt,
            )
            individual_rasters.append(rpath)
        except Exception as exc:
            log.warning("Watershed delineation failed for outlet %d: %s", i, exc)

    if not individual_rasters:
        raise RuntimeError("No individual watershed rasters could be produced.")

    # ------------------------------------------------------------------
    # Step 2: combine into a single labelled raster.
    # Each outlet's basin gets pixel value = outlet index (1-based).
    # Later outlets overwrite earlier ones where basins overlap, so
    # smaller/nested basins win (typical pour-point hierarchy).
    # ------------------------------------------------------------------
    with rasterio.open(individual_rasters[0]) as src0:
        profile = src0.profile.copy()
        nodata_val = src0.nodata if src0.nodata is not None else 0
        combined = np.zeros((src0.height, src0.width), dtype=np.int32)

    # Step 1: record pixel counts per basin (proxy for area)
    basin_sizes = {}
    basin_masks = {}
    for i, rpath in enumerate(individual_rasters, start=1):
        with rasterio.open(rpath) as src:
            data = src.read(1).astype(np.float32)
            nd = src.nodata if src.nodata is not None else -9999
            mask = (data != nd) & np.isfinite(data) & (data != 0)
        basin_masks[i] = mask
        basin_sizes[i] = int(mask.sum())

    # Step 2: sort by size ascending (smallest basin first = most downstream tributary)
    sorted_indices = sorted(basin_sizes, key=lambda k: basin_sizes[k])

    # Step 3: assign each pixel to the smallest basin that contains it
    # Then create EXCLUSIVE zones: each basin minus all smaller basins inside it
    exclusive = np.zeros_like(combined)
    claimed = np.zeros(combined.shape, dtype=bool)

    for i in sorted_indices:
        mask = basin_masks[i]
        # Only unclaimed pixels go to this basin's exclusive zone
        new_pixels = mask & ~claimed
        exclusive[new_pixels] = i
        claimed |= mask   # mark ALL pixels of this basin as claimed for next iterations

    combined = exclusive

    combined_path = out_dir / "subwatersheds_points.tif"
    profile.update(dtype=rasterio.int32, nodata=0, count=1, compress="lzw")
    with rasterio.open(combined_path, "w", **profile) as dst:
        dst.write(combined, 1)

    log.info(
        "Combined points subwatershed raster: %s | %d basins",
        combined_path.name, len(individual_rasters),
    )

    # Clean up temp individual rasters
    for rpath in individual_rasters:
        try:
            rpath.unlink()
        except Exception:
            pass

    return combined_path



def _delete_shapefile_if_exists(path: Path) -> None:
    """Remove a shapefile and all its sidecar files if they exist."""
    for ext in (".shp", ".shx", ".dbf", ".prj", ".cpg", ".sbn", ".sbx", ".qpj"):
        f = path.with_suffix(ext)
        try:
            if f.exists():
                f.unlink()
        except PermissionError as e:
            raise PermissionError(
                f"Cannot delete locked file {f}. Close any GIS application "
                f"(QGIS, ArcGIS) that may have it open, then retry."
            ) from e
        
# =============================================================================
# SECTION 12: RASTER → VECTOR POLYGON CONVERSION
# =============================================================================

def raster_to_polygon(
    raster_path: Union[str, Path],
    out_path: Union[str, Path],
    id_field: str = "basin_id",
    value_filter: Optional[float] = None,
) -> gpd.GeoDataFrame:
    raster_path = Path(raster_path)
    if not raster_path.exists():
        raise FileNotFoundError(f"Raster not found: {raster_path}")

    with rasterio.open(raster_path) as src:
        data = src.read(1)
        nodata = src.nodata
        crs = src.crs
        transform = src.transform

    mask = np.ones(data.shape, dtype=np.uint8)
    if nodata is not None:
        mask[data == nodata] = 0
    mask[~np.isfinite(data.astype(float))] = 0

    if value_filter is not None:
        mask[data != value_filter] = 0

    shapes = list(
        rasterio.features.shapes(
            data.astype(np.int32),
            mask=mask,
            transform=transform,
            connectivity=8,
        )
    )

    if not shapes:
        raise ValueError(f"No polygons extracted from raster: {raster_path.name}")

    geoms = []
    values = []
    for geom, val in shapes:
        geoms.append(shape(geom))
        values.append(int(val))

    gdf = gpd.GeoDataFrame({id_field: values, "geometry": geoms}, crs=crs)
    gdf = gdf[gdf.geometry.is_valid]
    gdf["geometry"] = gdf.geometry.buffer(0)

    if gdf.empty:
        raise ValueError("All extracted polygons were invalid.")

    
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ext = out_path.suffix.lower()
    driver = {"gpkg": "GPKG", "geojson": "GeoJSON"}.get(ext.lstrip("."), "ESRI Shapefile")

    # ← ADD THIS
    if driver == "ESRI Shapefile":
        _delete_shapefile_if_exists(out_path)

    gdf.to_file(str(out_path), driver=driver)
    

    log.info(
        "Raster→polygon: %d features | %s", len(gdf), out_path.name
    )
    return gdf

def watershed_raster_to_polygon(
    watershed_raster_path: Union[str, Path],
    out_path: Union[str, Path],
    basin_id: int = 1,
    dem_path: Optional[Union[str, Path]] = None,
    slope_path: Optional[Union[str, Path]] = None,
    stream_gdf: Optional[gpd.GeoDataFrame] = None,
    flow_acc_path: Optional[Union[str, Path]] = None,
) -> gpd.GeoDataFrame:
    watershed_raster_path = Path(watershed_raster_path)

    with rasterio.open(watershed_raster_path) as src:
        data = src.read(1)
        nodata = src.nodata
        unique_vals = np.unique(data)
        log.info("Watershed raster unique values: %s | nodata=%s", unique_vals, nodata)

    mask = np.ones(data.shape, dtype=bool)
    if nodata is not None:
        mask &= (data != nodata)
    mask &= np.isfinite(data.astype(float))
    mask &= (data != 0)

    valid_vals = data[mask]
    if valid_vals.size == 0:
        raise ValueError(
            f"Watershed raster has no valid basin pixels. "
            f"Unique values: {unique_vals}. "
            "Try adjusting outlet_coords or increasing snap_distance_cells."
        )

    detected_id = int(np.bincount(valid_vals.astype(int).ravel()).argmax())
    log.info("Auto-detected basin value: %d", detected_id)

    gdf = raster_to_polygon(
        raster_path=watershed_raster_path,
        out_path=out_path,
        id_field="basin_id",
        value_filter=float(detected_id),
    )
    dissolved = gdf.dissolve(by=None).reset_index(drop=True)
    dissolved["basin_id"]   = 1
    dissolved["area_m2"]    = dissolved.geometry.area
    dissolved["perim_m"]    = dissolved.geometry.length
    dissolved["area_km2"]   = dissolved["area_m2"] / 1e6
    dissolved["perim_km"]   = dissolved["perim_m"] / 1e3

    if dem_path and slope_path and stream_gdf is not None and flow_acc_path:
        try:
            m = compute_morphometrics(
                watershed_polygon=dissolved,
                dem_path=dem_path,
                slope_raster_path=slope_path,
                stream_gdf=stream_gdf,
                flow_acc_path=flow_acc_path,
            )
            scalar_keys = [
                "basin_length_km", "mean_elevation_m", "min_elevation_m",
                "max_elevation_m", "elevation_range_m", "mean_slope_deg",
                "relief_ratio", "drainage_density_km_km2", "drainage_frequency",
                "texture_ratio", "circularity_ratio", "compactness_coefficient",
                "elongation_ratio", "form_factor", "bifurcation_ratio",
                "longest_flow_path_km", "length_of_overland_flow_km",
                "constant_channel_maintenance", "total_stream_length_km",
                "max_strahler_order", "n_streams_total",
            ]
            for key in scalar_keys:
                dissolved[key[:10]] = m.get(key, None)
            log.info("Morphometrics embedded in watershed shapefile attribute table.")
        except Exception as exc:
            log.warning("Morphometrics embedding failed: %s", exc)

    _delete_shapefile_if_exists(Path(out_path))
    dissolved.to_file(str(out_path), driver="ESRI Shapefile")
    log.info("Watershed polygon: area=%.2f km²", dissolved["area_km2"].iloc[0])
    return dissolved

def subwatersheds_raster_to_polygons(
    subwatershed_raster_path: Union[str, Path],
    out_path: Union[str, Path],
    dem_path: Optional[Union[str, Path]] = None,
    slope_path: Optional[Union[str, Path]] = None,
    stream_gdf: Optional[gpd.GeoDataFrame] = None,
    flow_acc_path: Optional[Union[str, Path]] = None,
    clip_to: Optional[gpd.GeoDataFrame] = None,
) -> gpd.GeoDataFrame:
    tmp_path = Path(out_path).with_stem(Path(out_path).stem + "_tmp")
    gdf = raster_to_polygon(
        raster_path=subwatershed_raster_path,
        out_path=tmp_path,
        id_field="sub_id",
        value_filter=None,
    )
    gdf = gdf[gdf["sub_id"] > 0].copy()
    if not gdf.empty:
        gdf = gdf.dissolve(by="sub_id").reset_index()

    if clip_to is not None and not clip_to.empty:
        ws = clip_to.to_crs(gdf.crs) if clip_to.crs != gdf.crs else clip_to
        ws_union = ws.geometry.unary_union
        gdf = gdf[gdf.geometry.intersects(ws_union)].copy()
        gdf["geometry"] = gdf.geometry.intersection(ws_union)
        gdf = gdf[~gdf.geometry.is_empty].copy()
        gdf["geometry"] = gdf.geometry.buffer(0)
        gdf = gdf.reset_index(drop=True)
        log.info("Subwatersheds clipped to main watershed: %d retained", len(gdf))

    gdf["area_km2"] = gdf.geometry.area / 1e6
    gdf["perim_km"]  = gdf.geometry.length / 1e3

    if dem_path and slope_path and flow_acc_path and stream_gdf is not None:
        metrics_rows = []
        for idx, row in gdf.iterrows():
            sub_gdf = gpd.GeoDataFrame(
                [row], columns=gdf.columns, crs=gdf.crs
            )
            try:
                sub_streams = stream_gdf[
                    stream_gdf.geometry.intersects(row.geometry)
                ].copy()
                if sub_streams.empty:
                    sub_streams = stream_gdf.copy()

                m = compute_morphometrics(
                    watershed_polygon=sub_gdf,
                    dem_path=dem_path,
                    slope_raster_path=slope_path,
                    stream_gdf=sub_streams,
                    flow_acc_path=flow_acc_path,
                )
                metrics_rows.append(m)
            except Exception as exc:
                log.warning("Morphometrics failed for sub_id=%s: %s", row["sub_id"], exc)
                metrics_rows.append({})

        scalar_keys = [
            "area_km2", "perimeter_km", "basin_length_km",
            "mean_elevation_m", "min_elevation_m", "max_elevation_m",
            "elevation_range_m", "mean_slope_deg", "relief_ratio",
            "drainage_density_km_km2", "drainage_frequency", "texture_ratio",
            "circularity_ratio", "compactness_coefficient", "elongation_ratio",
            "form_factor", "bifurcation_ratio", "longest_flow_path_km",
            "length_of_overland_flow_km", "constant_channel_maintenance",
            "total_stream_length_km", "max_strahler_order", "n_streams_total",
        ]
        for key in scalar_keys:
            short = key[:10]
            gdf[short] = [m.get(key, None) for m in metrics_rows]

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    _delete_shapefile_if_exists(out_path)
    gdf.to_file(str(out_path))

    for ext in (".shp", ".shx", ".dbf", ".prj", ".cpg"):
        f = tmp_path.with_suffix(ext)
        if f.exists():
            f.unlink()

    log.info("Subwatershed polygons: %d features (with morphometrics)", len(gdf))
    return gdf


# =============================================================================
# SECTION 13: MORPHOMETRIC ANALYSIS
# =============================================================================

def _safe_area_km2(gdf: gpd.GeoDataFrame) -> float:
    crs = gdf.crs
    if crs and not crs.is_projected:
        gdf = gdf.to_crs(epsg=3857)
    return float(gdf.geometry.area.sum()) / 1e6


def _safe_perimeter_km(gdf: gpd.GeoDataFrame) -> float:
    crs = gdf.crs
    if crs and not crs.is_projected:
        gdf = gdf.to_crs(epsg=3857)
    return float(gdf.geometry.length.sum()) / 1e3


def compute_morphometrics(
    watershed_polygon: gpd.GeoDataFrame,
    dem_path: Union[str, Path],
    slope_raster_path: Union[str, Path],
    stream_gdf: gpd.GeoDataFrame,
    flow_acc_path: Union[str, Path],
) -> Dict[str, Any]:
    if watershed_polygon.empty:
        raise ValueError("Watershed polygon GDF is empty.")

    for p in (dem_path, slope_raster_path):
        if not Path(p).exists():
            raise FileNotFoundError(f"Raster not found: {p}")

    log.info("=== Morphometric Analysis ===")

    area_km2 = _safe_area_km2(watershed_polygon)
    perimeter_km = _safe_perimeter_km(watershed_polygon)

    bounds = watershed_polygon.to_crs(epsg=3857).total_bounds
    dx = (bounds[2] - bounds[0]) / 1e3
    dy = (bounds[3] - bounds[1]) / 1e3
    basin_length_km = max(dx, dy)

    ws_geom = [mapping(g) for g in watershed_polygon.geometry]
    with rasterio.open(dem_path) as src:
        dem_crs = src.crs
        ws_reproj = watershed_polygon.to_crs(dem_crs) if watershed_polygon.crs != dem_crs else watershed_polygon
        ws_geoms = [mapping(g) for g in ws_reproj.geometry]
        try:
            dem_clipped, _ = rasterio.mask.mask(src, ws_geoms, crop=True, nodata=np.nan)
            elev = dem_clipped[0].astype(np.float32)
            elev[elev == src.nodata] = np.nan
        except Exception:
            with rasterio.open(dem_path) as src2:
                elev = src2.read(1).astype(np.float32)
                if src2.nodata is not None:
                    elev[elev == src2.nodata] = np.nan

    valid_elev = elev[np.isfinite(elev)]
    mean_elev = float(np.nanmean(valid_elev)) if valid_elev.size > 0 else 0.0
    min_elev = float(np.nanmin(valid_elev)) if valid_elev.size > 0 else 0.0
    max_elev = float(np.nanmax(valid_elev)) if valid_elev.size > 0 else 0.0
    elev_range = max_elev - min_elev

    with rasterio.open(slope_raster_path) as src:
        slope_crs = src.crs
        ws_reproj2 = watershed_polygon.to_crs(slope_crs) if watershed_polygon.crs != slope_crs else watershed_polygon
        ws_geoms2 = [mapping(g) for g in ws_reproj2.geometry]
        try:
            slope_clipped, _ = rasterio.mask.mask(src, ws_geoms2, crop=True, nodata=np.nan)
            slope_arr = slope_clipped[0].astype(np.float32)
            if src.nodata is not None:
                slope_arr[slope_arr == src.nodata] = np.nan
        except Exception:
            with rasterio.open(slope_raster_path) as src2:
                slope_arr = src2.read(1).astype(np.float32)
                if src2.nodata is not None:
                    slope_arr[slope_arr == src2.nodata] = np.nan

    mean_slope = float(np.nanmean(slope_arr[np.isfinite(slope_arr)])) if np.any(np.isfinite(slope_arr)) else 0.0

    relief_ratio = elev_range / (basin_length_km * 1000) if basin_length_km > 0 else 0.0

    total_length_km = float(stream_gdf["length_m"].sum()) / 1e3 if "length_m" in stream_gdf.columns else 0.0
    n_streams = len(stream_gdf)

    drainage_density = total_length_km / area_km2 if area_km2 > 0 else 0.0
    drainage_frequency = n_streams / area_km2 if area_km2 > 0 else 0.0
    texture_ratio = n_streams / perimeter_km if perimeter_km > 0 else 0.0
    length_of_overland_flow = 1.0 / (2.0 * drainage_density) if drainage_density > 0 else 0.0
    constant_channel_maintenance = 1.0 / drainage_density if drainage_density > 0 else 0.0

    circularity_ratio = (4 * math.pi * area_km2) / (perimeter_km ** 2) if perimeter_km > 0 else 0.0
    compactness_coefficient = perimeter_km / (2 * math.sqrt(math.pi * area_km2)) if area_km2 > 0 else 0.0
    elongation_ratio = (2 * math.sqrt(area_km2 / math.pi)) / basin_length_km if basin_length_km > 0 else 0.0
    form_factor = area_km2 / (basin_length_km ** 2) if basin_length_km > 0 else 0.0

    if "strahler" in stream_gdf.columns:
        stream_gdf_clean = stream_gdf.dropna(subset=["strahler"]).copy()
        stream_gdf_clean["strahler"] = stream_gdf_clean["strahler"].astype(int)
        max_order = int(stream_gdf_clean["strahler"].max()) if len(stream_gdf_clean) > 0 else 1

        order_counts: Dict[int, int] = {}
        order_mean_len: Dict[int, float] = {}
        for order in range(1, max_order + 1):
            subset = stream_gdf_clean[stream_gdf_clean["strahler"] == order]
            order_counts[order] = len(subset)
            order_mean_len[order] = float(subset["length_m"].mean()) / 1e3 if len(subset) > 0 else 0.0

        bif_ratios = []
        for order in range(1, max_order):
            n_u = order_counts.get(order, 0)
            n_u1 = order_counts.get(order + 1, 1)
            if n_u1 > 0:
                bif_ratios.append(n_u / n_u1)
        bifurcation_ratio = float(np.mean(bif_ratios)) if bif_ratios else 0.0
    else:
        max_order = 1
        order_counts = {}
        order_mean_len = {}
        bifurcation_ratio = 0.0

    longest_flow_path_km = basin_length_km * 1.2

    metrics: Dict[str, Any] = {
        "area_km2": round(area_km2, 4),
        "perimeter_km": round(perimeter_km, 4),
        "basin_length_km": round(basin_length_km, 4),
        "mean_elevation_m": round(mean_elev, 2),
        "min_elevation_m": round(min_elev, 2),
        "max_elevation_m": round(max_elev, 2),
        "elevation_range_m": round(elev_range, 2),
        "mean_slope_deg": round(mean_slope, 4),
        "relief_ratio": round(relief_ratio, 6),
        "drainage_density_km_km2": round(drainage_density, 4),
        "drainage_frequency": round(drainage_frequency, 4),
        "texture_ratio": round(texture_ratio, 4),
        "circularity_ratio": round(circularity_ratio, 4),
        "compactness_coefficient": round(compactness_coefficient, 4),
        "elongation_ratio": round(elongation_ratio, 4),
        "form_factor": round(form_factor, 4),
        "bifurcation_ratio": round(bifurcation_ratio, 4),
        "longest_flow_path_km": round(longest_flow_path_km, 4),
        "length_of_overland_flow_km": round(length_of_overland_flow, 4),
        "constant_channel_maintenance": round(constant_channel_maintenance, 4),
        "total_stream_length_km": round(total_length_km, 4),
        "max_strahler_order": max_order,
        "streams_per_order": order_counts,
        "mean_stream_length_km": order_mean_len,
        "n_streams_total": n_streams,
    }

    log.info(
        "Morphometrics: area=%.2f km² | Dd=%.3f | order=%d | Rb=%.2f",
        area_km2, drainage_density, max_order, bifurcation_ratio,
    )
    return metrics


def save_morphometrics(
    metrics: Dict[str, Any],
    out_dir: Union[str, Path],
    prefix: str = "watershed",
) -> Dict[str, Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    flat: Dict[str, Any] = {}
    for k, v in metrics.items():
        if isinstance(v, dict):
            for sub_k, sub_v in v.items():
                flat[f"{k}_order{sub_k}"] = sub_v
        else:
            flat[k] = v

    df = pd.DataFrame([flat])
    csv_path = out_dir / f"{prefix}_morphometrics.csv"
    df.to_csv(str(csv_path), index=False)

    json_path = out_dir / f"{prefix}_morphometrics.json"
    with open(json_path, "w") as f:
        json.dump(metrics, f, indent=2, default=str)

    log.info("Morphometrics saved: %s, %s", csv_path.name, json_path.name)
    return {"csv": csv_path.resolve(), "json": json_path.resolve()}


# =============================================================================
# SECTION 14: VISUALIZATION
# =============================================================================

_CMAP_DEM = "terrain"
_CMAP_SLOPE = "RdYlGn_r"
_CMAP_ASPECT = "hsv"
_CMAP_FLOW_ACC = "Blues"


def _read_raster_for_plot(path: Union[str, Path]) -> Tuple[np.ndarray, rasterio.profiles.Profile]:
    with rasterio.open(path) as src:
        data = src.read(1).astype(np.float32)
        profile = src.profile
        if src.nodata is not None:
            data[data == src.nodata] = np.nan
    return data, profile


def _extent_from_profile(profile: rasterio.profiles.Profile) -> List[float]:
    t = profile["transform"]
    w, h = profile["width"], profile["height"]
    left = t.c
    right = t.c + t.a * w
    bottom = t.f + t.e * h
    top = t.f
    return [left, right, bottom, top]


def plot_dem(dem_path, hillshade_path, out_path, title="Digital Elevation Model", dpi=150):
    dem, prof = _read_raster_for_plot(dem_path)
    extent = _extent_from_profile(prof)
    fig, ax = plt.subplots(figsize=(10, 8))
    if hillshade_path and Path(hillshade_path).exists():
        hs, _ = _read_raster_for_plot(hillshade_path)
        ax.imshow(hs, extent=extent, cmap="gray", alpha=0.5, origin="upper")
        im = ax.imshow(dem, extent=extent, cmap=_CMAP_DEM, alpha=0.7, origin="upper")
    else:
        im = ax.imshow(dem, extent=extent, cmap=_CMAP_DEM, origin="upper")
    cbar = fig.colorbar(im, ax=ax, shrink=0.6, label="Elevation (m)")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("Easting / Longitude")
    ax.set_ylabel("Northing / Latitude")
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_path), dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    log.info("DEM figure saved: %s", out_path.name)
    return out_path.resolve()


def plot_slope(slope_path, out_path, title="Slope (degrees)", dpi=150):
    data, prof = _read_raster_for_plot(slope_path)
    extent = _extent_from_profile(prof)
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(data, extent=extent, cmap=_CMAP_SLOPE, origin="upper", vmin=0)
    fig.colorbar(im, ax=ax, shrink=0.6, label="Slope (°)")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("Easting / Longitude")
    ax.set_ylabel("Northing / Latitude")
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_path), dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    log.info("Slope figure saved: %s", out_path.name)
    return out_path.resolve()


def plot_aspect(aspect_path, out_path, title="Aspect (degrees from North)", dpi=150):
    data, prof = _read_raster_for_plot(aspect_path)
    extent = _extent_from_profile(prof)
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(data, extent=extent, cmap=_CMAP_ASPECT, origin="upper", vmin=0, vmax=360)
    fig.colorbar(im, ax=ax, shrink=0.6, label="Aspect (°)")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("Easting / Longitude")
    ax.set_ylabel("Northing / Latitude")
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_path), dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    log.info("Aspect figure saved: %s", out_path.name)
    return out_path.resolve()


def plot_flow_accumulation(flow_acc_path, out_path, title="Flow Accumulation (log₁₀ cells)", dpi=150):
    data, prof = _read_raster_for_plot(flow_acc_path)
    extent = _extent_from_profile(prof)
    log_data = np.where(data > 0, np.log10(data + 1), np.nan)
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(log_data, extent=extent, cmap=_CMAP_FLOW_ACC, origin="upper")
    fig.colorbar(im, ax=ax, shrink=0.6, label="log₁₀(cells)")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("Easting / Longitude")
    ax.set_ylabel("Northing / Latitude")
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_path), dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    log.info("Flow accumulation figure saved: %s", out_path.name)
    return out_path.resolve()


def plot_watershed(
    dem_path, hillshade_path, watershed_gdf, stream_gdf, outlet_pts,
    out_path, title="Watershed Delineation", dpi=150,
):
    dem, prof = _read_raster_for_plot(dem_path)
    extent = _extent_from_profile(prof)
    fig, ax = plt.subplots(figsize=(12, 10))
    if hillshade_path and Path(hillshade_path).exists():
        hs, _ = _read_raster_for_plot(hillshade_path)
        ax.imshow(hs, extent=extent, cmap="gray", alpha=0.5, origin="upper")
        ax.imshow(dem, extent=extent, cmap=_CMAP_DEM, alpha=0.55, origin="upper")
    else:
        ax.imshow(dem, extent=extent, cmap=_CMAP_DEM, alpha=0.85, origin="upper")
    if not watershed_gdf.empty:
        watershed_gdf.boundary.plot(ax=ax, color="red", linewidth=2.5, label="Watershed boundary")
    if not stream_gdf.empty:
        stream_gdf.plot(ax=ax, color="dodgerblue", linewidth=0.8, label="Streams")
    for x, y in outlet_pts:
        ax.plot(x, y, "r*", markersize=14, zorder=10, label="Outlet")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("Easting / Longitude")
    ax.set_ylabel("Northing / Latitude")
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc="lower right")
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_path), dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    log.info("Watershed figure saved: %s", out_path.name)
    return out_path.resolve()


def plot_subwatersheds(
    dem_path, hillshade_path, subwatershed_gdf, stream_gdf,
    out_path, title="Subwatersheds", dpi=150,
):
    dem, prof = _read_raster_for_plot(dem_path)
    extent = _extent_from_profile(prof)
    n_sub = len(subwatershed_gdf)
    fig, ax = plt.subplots(figsize=(12, 10))
    if hillshade_path and Path(hillshade_path).exists():
        hs, _ = _read_raster_for_plot(hillshade_path)
        ax.imshow(hs, extent=extent, cmap="gray", alpha=0.4, origin="upper")
    subwatershed_gdf.plot(ax=ax, cmap="tab20", alpha=0.55, edgecolor="white", linewidth=0.5)
    if not stream_gdf.empty:
        stream_gdf.plot(ax=ax, color="navy", linewidth=0.7)
    ax.set_title(f"{title} (n={n_sub})", fontsize=14, fontweight="bold")
    ax.set_xlabel("Easting / Longitude")
    ax.set_ylabel("Northing / Latitude")
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_path), dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    log.info("Subwatersheds figure saved: %s", out_path.name)
    return out_path.resolve()


def plot_ordered_streams(dem_path, stream_gdf, out_path, title="Strahler Stream Order", dpi=150):
    dem, prof = _read_raster_for_plot(dem_path)
    extent = _extent_from_profile(prof)
    fig, ax = plt.subplots(figsize=(12, 10))
    ax.imshow(dem, extent=extent, cmap="Greys", alpha=0.6, origin="upper")
    if "strahler" in stream_gdf.columns:
        orders = sorted(stream_gdf["strahler"].dropna().unique())
        cmap = plt.get_cmap("RdYlBu_r", len(orders))
        color_map = {o: cmap(i) for i, o in enumerate(orders)}
        for order in orders:
            subset = stream_gdf[stream_gdf["strahler"] == order]
            lw = 0.5 + float(order) * 0.5
            rgba = color_map[order]
            hex_color = "#{:02x}{:02x}{:02x}".format(
                int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255)
            )
            subset.plot(ax=ax, color=hex_color, linewidth=lw, label=f"Order {order}")
        ax.legend(title="Strahler Order", loc="lower right", fontsize=8)
    else:
        stream_gdf.plot(ax=ax, color="dodgerblue", linewidth=0.8)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("Easting / Longitude")
    ax.set_ylabel("Northing / Latitude")
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_path), dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    log.info("Ordered streams figure saved: %s", out_path.name)
    return out_path.resolve()


def plot_combined_summary(
    dem_path, hillshade_path, slope_path, flow_acc_path,
    watershed_gdf, stream_gdf, out_path, title="Hydrology Summary", dpi=150,
):
    fig, axes = plt.subplots(2, 2, figsize=(18, 14))
    fig.suptitle(title, fontsize=16, fontweight="bold", y=1.01)
    dem, prof = _read_raster_for_plot(dem_path)
    extent = _extent_from_profile(prof)
    ax = axes[0, 0]
    if hillshade_path and Path(hillshade_path).exists():
        hs, _ = _read_raster_for_plot(hillshade_path)
        ax.imshow(hs, extent=extent, cmap="gray", alpha=0.5, origin="upper")
        im1 = ax.imshow(dem, extent=extent, cmap=_CMAP_DEM, alpha=0.7, origin="upper")
    else:
        im1 = ax.imshow(dem, extent=extent, cmap=_CMAP_DEM, origin="upper")
    fig.colorbar(im1, ax=ax, shrink=0.7, label="m")
    ax.set_title("DEM", fontsize=12)
    ax = axes[0, 1]
    if slope_path and Path(slope_path).exists():
        slope, _ = _read_raster_for_plot(slope_path)
        im2 = ax.imshow(slope, extent=extent, cmap=_CMAP_SLOPE, origin="upper", vmin=0)
        fig.colorbar(im2, ax=ax, shrink=0.7, label="°")
        ax.set_title("Slope", fontsize=12)
    else:
        ax.set_visible(False)
    ax = axes[1, 0]
    if flow_acc_path and Path(flow_acc_path).exists():
        fa, _ = _read_raster_for_plot(flow_acc_path)
        log_fa = np.where(fa > 0, np.log10(fa + 1), np.nan)
        im3 = ax.imshow(log_fa, extent=extent, cmap=_CMAP_FLOW_ACC, origin="upper")
        fig.colorbar(im3, ax=ax, shrink=0.7, label="log₁₀(cells)")
        ax.set_title("Flow Accumulation", fontsize=12)
    else:
        ax.set_visible(False)
    ax = axes[1, 1]
    ax.imshow(dem, extent=extent, cmap="Greys", alpha=0.6, origin="upper")
    if not watershed_gdf.empty:
        watershed_gdf.boundary.plot(ax=ax, color="red", linewidth=2)
    if not stream_gdf.empty:
        stream_gdf.plot(ax=ax, color="dodgerblue", linewidth=0.7)
    ax.set_title("Watershed & Streams", fontsize=12)
    for ax_row in axes:
        for a in ax_row:
            a.set_xlabel("X", fontsize=8)
            a.set_ylabel("Y", fontsize=8)
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_path), dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    log.info("Combined summary figure saved: %s", out_path.name)
    return out_path.resolve()


def run_visualization(
    dem_path, terrain_outputs, hydrology_outputs, watershed_gdf,
    subwatershed_gdf, stream_gdf, outlet_pts, figures_dir, dpi=150,
):
    fig_dir = Path(figures_dir)
    hs = terrain_outputs.get("hillshade")
    sl = terrain_outputs.get("slope")
    fa = hydrology_outputs.get("flow_acc")
    figs: Dict[str, Path] = {}
    log.info("=== Visualization ===")
    figs["dem"] = plot_dem(dem_path, hs, fig_dir / "dem.png", dpi=dpi)
    if sl:
        figs["slope"] = plot_slope(sl, fig_dir / "slope.png", dpi=dpi)
    if terrain_outputs.get("aspect"):
        figs["aspect"] = plot_aspect(terrain_outputs["aspect"], fig_dir / "aspect.png", dpi=dpi)
    if fa:
        figs["flow_accumulation"] = plot_flow_accumulation(fa, fig_dir / "flow_accumulation.png", dpi=dpi)
    if not watershed_gdf.empty:
        figs["watershed"] = plot_watershed(
            dem_path, hs, watershed_gdf, stream_gdf, outlet_pts,
            fig_dir / "watershed.png", dpi=dpi,
        )
    if not subwatershed_gdf.empty:
        figs["subwatersheds"] = plot_subwatersheds(
            dem_path, hs, subwatershed_gdf, stream_gdf,
            fig_dir / "subwatersheds.png", dpi=dpi,
        )
    if not stream_gdf.empty:
        figs["ordered_streams"] = plot_ordered_streams(
            dem_path, stream_gdf, fig_dir / "ordered_streams.png", dpi=dpi,
        )
    figs["summary"] = plot_combined_summary(
        dem_path, hs, sl, fa, watershed_gdf, stream_gdf,
        fig_dir / "summary.png", dpi=dpi,
    )
    log.info("Visualization complete: %d figures", len(figs))
    return figs


# =============================================================================
# SECTION 15: GIS EXPORT UTILITIES
# =============================================================================

def export_to_geopackage(layers: Dict[str, gpd.GeoDataFrame], out_path: Union[str, Path]) -> Path:
    if not layers:
        raise ValueError("No layers provided for GeoPackage export.")
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists():
        out_path.unlink()
    for name, gdf in layers.items():
        if gdf is not None and not gdf.empty:
            gdf.to_file(str(out_path), layer=name, driver="GPKG")
            log.debug("GeoPackage layer written: %s (%d rows)", name, len(gdf))
    log.info("GeoPackage exported: %s (%d layers)", out_path.name, len(layers))
    return out_path.resolve()


def export_to_geojson(gdf: gpd.GeoDataFrame, out_path: Union[str, Path], to_wgs84: bool = True) -> Path:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    write_gdf = gdf.copy()
    if to_wgs84 and gdf.crs and not gdf.crs.to_epsg() == 4326:
        write_gdf = gdf.to_crs(epsg=4326)
    write_gdf.to_file(str(out_path), driver="GeoJSON")
    log.info("GeoJSON exported: %s", out_path.name)
    return out_path.resolve()


# =============================================================================
# SECTION 16: MAIN WORKFLOW ORCHESTRATOR
# =============================================================================

def run_hydrology_workflow(
    dem_path: Union[str, Path],
    outlet_coords: Union[Tuple[float, float], List[Tuple[float, float]]],
    output_dir: Union[str, Path] = "outputs",
    mode: str = "hybrid",
    stream_threshold: Optional[float] = None,
    auto_percentile: float = 99.5,
    snap_distance_cells: int = 50,
    preprocessing_method: str = "breach",
    src_epsg: int = 4326,
    dpi: int = 150,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Run the complete hydrology analysis workflow.

    outlet_coords:
        Single (lon, lat) → main watershed only (in all modes).
        List of (lon, lat) → first is the primary/main outlet; extras become
        user-defined sub-watershed outlets in "points"/"hybrid" modes.

    mode:
        "main"   – main watershed only.
        "auto"   – main watershed + automatic subwatersheds via stream links.
        "points" – main watershed + user-defined sub-watersheds from extra outlets,
                   all clipped inside the main watershed, one merged shapefile.
        "hybrid" – all three combined.
    """
    VALID_MODES = {"main", "auto", "points", "hybrid"}
    if mode not in VALID_MODES:
        raise ValueError(f"mode must be one of {VALID_MODES}; got '{mode}'")

    dem_path = Path(dem_path)
    if not dem_path.exists():
        raise FileNotFoundError(f"DEM not found: {dem_path}")

    # ------------------------------------------------------------------ #
    # 0. Setup                                                             #
    # ------------------------------------------------------------------ #
    dirs = create_output_dirs(output_dir)
    log = setup_logging(dirs["logs"])
    log.info("=" * 60)
    log.info("HYDROLOGY WORKFLOW START | mode=%s | DEM=%s", mode, dem_path.name)
    log.info("=" * 60)

    results: Dict[str, Any] = {"output_dirs": dirs}

    if isinstance(outlet_coords, tuple) and len(outlet_coords) == 2 and not isinstance(outlet_coords[0], (list, tuple)):
        outlet_list: List[Tuple[float, float]] = [outlet_coords]  # type: ignore[list-item]
    else:
        outlet_list = list(outlet_coords)  # type: ignore[arg-type]

    # ------------------------------------------------------------------ #
    # 1. DEM read + validate                                               #
    # ------------------------------------------------------------------ #
    dem_data, dem_profile = read_dem(dem_path)
    validate_dem(dem_data, dem_profile)
    meta = read_metadata(dem_path)
    proj_crs = determine_projected_crs(dem_path)
    results["projected_crs"] = proj_crs

    if not meta["is_projected"]:
        reproj_dem = dirs["preprocessed"] / "dem_projected.tif"
        dem_path = reproject_raster(dem_path, reproj_dem, proj_crs)
        log.info("DEM reprojected to: %s", proj_crs)

    # ------------------------------------------------------------------ #
    # 2. WhiteboxTools                                                     #
    # ------------------------------------------------------------------ #
    wbt = init_whitebox(verbose=verbose)

    # ------------------------------------------------------------------ #
    # 3. Preprocessing                                                     #
    # ------------------------------------------------------------------ #
    pre_outputs = preprocess_dem(
        dem_path=dem_path,
        preprocessed_dir=dirs["preprocessed"],
        wbt=wbt,
        method=preprocessing_method,
    )
    conditioned_dem = pre_outputs["conditioned"]
    results["conditioned_dem"] = conditioned_dem

    # ------------------------------------------------------------------ #
    # 4. Terrain analysis                                                  #
    # ------------------------------------------------------------------ #
    terrain_outputs = run_terrain_analysis(conditioned_dem, dirs["terrain"], wbt)
    results.update(terrain_outputs)

    # ------------------------------------------------------------------ #
    # 5. Flow analysis                                                     #
    # ------------------------------------------------------------------ #
    flow_outputs = run_flow_analysis(conditioned_dem, dirs["hydrology"], wbt)
    results["flow_dir"] = flow_outputs["flow_dir"]
    results["flow_acc"] = flow_outputs["flow_acc"]
    if "flow_length" in flow_outputs:
        results["flow_length"] = flow_outputs["flow_length"]

    # ------------------------------------------------------------------ #
    # 6. Stream extraction                                                 #
    # ------------------------------------------------------------------ #
    stream_raster, threshold_used = extract_stream_raster(
        flow_acc_path=flow_outputs["flow_acc"],
        out_path=dirs["streams"] / "streams.tif",
        wbt=wbt,
        threshold=stream_threshold,
        auto_percentile=auto_percentile,
    )
    results["stream_raster"] = stream_raster
    results["stream_threshold"] = threshold_used

    stream_vector = extract_stream_vector(
        stream_raster_path=stream_raster,
        flow_dir_path=flow_outputs["flow_dir"],
        out_path=dirs["streams"] / "streams.shp",
        wbt=wbt,
    )

    stream_links = generate_stream_links(
        stream_raster_path=stream_raster,
        flow_dir_path=flow_outputs["flow_dir"],
        out_path=dirs["streams"] / "stream_links.tif",
        wbt=wbt,
    )
    results["stream_links"] = stream_links

    # ------------------------------------------------------------------ #
    # 7. Strahler order                                                   #
    # ------------------------------------------------------------------ #
    strahler_raster = compute_strahler_order(
        stream_raster_path=stream_raster,
        flow_dir_path=flow_outputs["flow_dir"],
        out_path=dirs["streams"] / "strahler_order.tif",
        wbt=wbt,
    )
    results["strahler_raster"] = strahler_raster

    stream_gdf = build_ordered_stream_geodataframe(
        stream_vector_path=stream_vector,
        strahler_raster_path=strahler_raster,
        flow_acc_path=flow_outputs["flow_acc"],
    )
    
    save_vector(stream_gdf, dirs["streams"] / "streams_ordered.shp")
    results["stream_gdf"] = stream_gdf

    # ------------------------------------------------------------------ #
    # 8. Outlet processing                                                 #
    # ------------------------------------------------------------------ #
    outlet_info = process_outlets(
        outlet_coords=outlet_list,
        dem_path=conditioned_dem,
        flow_acc_path=flow_outputs["flow_acc"],
        dst_crs=proj_crs,
        output_dir=dirs["watershed"],
        snap_distance_cells=snap_distance_cells,
        src_epsg=src_epsg,
    )
    snapped_outlets = outlet_info["snapped_outlets"]
    primary_outlet = outlet_info["primary_outlet"]
    results["snapped_outlets"] = snapped_outlets

    # ------------------------------------------------------------------ #
    # 9. Main watershed (always)                                           #
    # ------------------------------------------------------------------ #
    ws_raster = delineate_watershed(
        flow_dir_path=flow_outputs["flow_dir"],
        outlet_x=primary_outlet[0],
        outlet_y=primary_outlet[1],
        crs=proj_crs,
        out_raster_path=dirs["watershed"] / "watershed.tif",
        wbt=wbt,
    )
    results["watershed_raster"] = ws_raster

    watershed_gdf = watershed_raster_to_polygon(
        watershed_raster_path=ws_raster,
        out_path=dirs["watershed"] / "watershed.shp",
        dem_path=conditioned_dem,
        slope_path=terrain_outputs.get("slope"),
        stream_gdf=stream_gdf,
        flow_acc_path=flow_outputs["flow_acc"],
    )
    results["watershed_gdf"] = watershed_gdf
    results["watershed_polygon_path"] = (dirs["watershed"] / "watershed.shp").resolve()

    # ------------------------------------------------------------------ #
    # 10. Automatic subwatersheds (Mode "auto" / "hybrid")                 #
    # ------------------------------------------------------------------ #
    subwatershed_auto_gdf = gpd.GeoDataFrame()
    if mode in ("auto", "hybrid"):
        sub_raster = delineate_subwatersheds_auto(
            flow_dir_path=flow_outputs["flow_dir"],
            stream_links_path=stream_links,
            out_raster_path=dirs["subwatersheds_auto"] / "subwatersheds_auto.tif",
            wbt=wbt,
        )
        results["subwatersheds_auto_raster"] = sub_raster
        subwatershed_auto_gdf = subwatersheds_raster_to_polygons(
            subwatershed_raster_path=sub_raster,
            out_path=dirs["subwatersheds_auto"] / "subwatersheds_auto.shp",
            dem_path=conditioned_dem,
            slope_path=terrain_outputs.get("slope"),
            stream_gdf=stream_gdf,
            flow_acc_path=flow_outputs["flow_acc"],
            clip_to=watershed_gdf,
        )
        results["subwatershed_auto_gdf"] = subwatershed_auto_gdf
        results["subwatersheds_auto_path"] = (dirs["subwatersheds_auto"] / "subwatersheds_auto.shp").resolve()

    # ------------------------------------------------------------------ #
    # 11. User-defined sub-watersheds (Mode "points" / "hybrid")           #
    #                                                                      #
    # Behaviour mirrors auto-mode:                                         #
    #   • outlets[1:] define sub-watershed pour points                    #
    #   • delineate each, combine into ONE labelled raster                 #
    #   • polygonise → clip to main watershed → morphometrics              #
    #   • save as ONE shapefile (subwatersheds_points.shp)                 #
    # ------------------------------------------------------------------ #
    subwatershed_points_gdf = gpd.GeoDataFrame()
    if mode in ("points", "hybrid") and len(snapped_outlets) > 1:
        # Use only the extra outlets (index 1 onward) as sub-watershed pours
        sub_outlets = snapped_outlets[1:]

        pts_combined_raster = delineate_subwatersheds_points(
            flow_dir_path=flow_outputs["flow_dir"],
            outlets=sub_outlets,
            crs=proj_crs,
            out_dir=dirs["subwatersheds_points"],
            wbt=wbt,
        )
        results["subwatersheds_points_raster"] = pts_combined_raster

        subwatershed_points_gdf = subwatersheds_raster_to_polygons(
            subwatershed_raster_path=pts_combined_raster,
            out_path=dirs["subwatersheds_points"] / "subwatersheds_points.shp",
            dem_path=conditioned_dem,
            slope_path=terrain_outputs.get("slope"),
            stream_gdf=stream_gdf,
            flow_acc_path=flow_outputs["flow_acc"],
            clip_to=watershed_gdf,           # ← clips to main watershed
        )
        # Tag each polygon with its outlet index for traceability
        subwatershed_points_gdf["outlet_idx"] = subwatershed_points_gdf["sub_id"]

        results["subwatershed_points_gdf"] = subwatershed_points_gdf
        results["subwatersheds_points_path"] = (
            dirs["subwatersheds_points"] / "subwatersheds_points.shp"
        ).resolve()

        log.info(
            "User-defined subwatersheds: %d polygons clipped to main watershed",
            len(subwatershed_points_gdf),
        )
    elif mode in ("points", "hybrid") and len(snapped_outlets) == 1:
        log.warning(
            "mode='%s' but only one outlet provided — "
            "no sub-watersheds will be delineated. "
            "Pass extra outlets as additional (lon, lat) pairs.", mode,
        )

    # ------------------------------------------------------------------ #
    # 12. Morphometric analysis (main watershed summary)                   #
    # ------------------------------------------------------------------ #
    morphometrics = watershed_gdf.iloc[0].to_dict() if not watershed_gdf.empty else {}
    results["morphometrics"] = morphometrics

    # ------------------------------------------------------------------ #
    # 13. GIS exports                                                      #
    # ------------------------------------------------------------------ #
    layers_for_gpkg: Dict[str, gpd.GeoDataFrame] = {
        "watershed": watershed_gdf,
        "streams": stream_gdf,
    }
    if not subwatershed_auto_gdf.empty:
        layers_for_gpkg["subwatersheds_auto"] = subwatershed_auto_gdf
    if not subwatershed_points_gdf.empty:
        layers_for_gpkg["subwatersheds_points"] = subwatershed_points_gdf

    gpkg_path = export_to_geopackage(
        layers=layers_for_gpkg,
        out_path=dirs["root"] / "hydrology_outputs.gpkg",
    )
    results["geopackage"] = gpkg_path

    if not watershed_gdf.empty:
        export_to_geojson(watershed_gdf, dirs["watershed"] / "watershed.geojson")
    if not stream_gdf.empty:
        export_to_geojson(stream_gdf, dirs["streams"] / "streams.geojson")

    # ------------------------------------------------------------------ #
    # 14. Visualization                                                    #
    # ------------------------------------------------------------------ #
    # For visualization, prefer points subwatersheds if available, else auto
    viz_sub_gdf = (
        subwatershed_points_gdf if not subwatershed_points_gdf.empty
        else subwatershed_auto_gdf
    )

    figs = run_visualization(
        dem_path=conditioned_dem,
        terrain_outputs=terrain_outputs,
        hydrology_outputs=flow_outputs,
        watershed_gdf=watershed_gdf,
        subwatershed_gdf=viz_sub_gdf,
        stream_gdf=stream_gdf,
        outlet_pts=snapped_outlets,
        figures_dir=dirs["figures"],
        dpi=dpi,
    )
    results["figures"] = figs

    # ------------------------------------------------------------------ #
    # Done                                                                 #
    # ------------------------------------------------------------------ #
    log.info("=" * 60)
    log.info("HYDROLOGY WORKFLOW COMPLETE")
    log.info("  Watershed area : %.2f km²", morphometrics.get("area_km2", 0.0))
    log.info("  Max order      : %d", morphometrics.get("max_strahler_order", 0))
    log.info("  Stream count   : %d", morphometrics.get("n_streams_total", 0))
    log.info("  Outputs        : %s", dirs["root"].resolve())
    log.info("=" * 60)

    return results

# =============================================================================
# SECTION 17: CLI ENTRY POINT
# =============================================================================

def _is_cli_invocation() -> bool:
    import sys
    if not sys.argv:
        return False
    argv0 = sys.argv[0]
    return argv0.endswith("hydrology_tool.py") or argv0.endswith("hydrology_tool")


if _is_cli_invocation():
    import argparse

    parser = argparse.ArgumentParser(
        description="Hydrology Toolkit — watershed delineation and analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("dem", help="Path to input DEM raster")
    parser.add_argument(
        "--outlet",
        nargs=2,
        type=float,
        metavar=("LON", "LAT"),
        required=True,
        help="Primary outlet coordinate in WGS84 decimal degrees",
    )
    parser.add_argument(
        "--extra-outlets",
        nargs="+",
        type=float,
        metavar="LON_LAT",
        default=[],
        help="Additional outlets as lon lat lon lat … pairs",
    )
    parser.add_argument("--output-dir", default="outputs", help="Root output directory")
    parser.add_argument(
        "--mode",
        choices=["main", "auto", "points", "hybrid"],
        default="hybrid",
        help="Analysis mode",
    )
    parser.add_argument("--threshold", type=float, default=None, help="Stream threshold (cells)")
    parser.add_argument("--snap-cells", type=int, default=50, help="Outlet snap radius (cells)")
    parser.add_argument("--preprocess", choices=["breach", "fill"], default="breach")
    parser.add_argument("--dpi", type=int, default=150, help="Figure DPI")
    parser.add_argument("--verbose", action="store_true", help="Verbose WhiteboxTools output")

    args = parser.parse_args()

    primary = (args.outlet[0], args.outlet[1])
    extra_flat = args.extra_outlets
    extra = [(extra_flat[i], extra_flat[i + 1]) for i in range(0, len(extra_flat) - 1, 2)]
    all_outlets: Union[Tuple, List] = [primary] + extra if extra else primary

    run_hydrology_workflow(
        dem_path=args.dem,
        outlet_coords=all_outlets,
        output_dir=args.output_dir,
        mode=args.mode,
        stream_threshold=args.threshold,
        snap_distance_cells=args.snap_cells,
        preprocessing_method=args.preprocess,
        dpi=args.dpi,
        verbose=args.verbose,
    )