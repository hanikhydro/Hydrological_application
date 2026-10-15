# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 12:48:23 2026

@author: acer
"""

# ============================================================
# DEM Downloader
# Supports:
#   - SRTM
#   - NASADEM
#   - FABDEM
#   - MERIT
#   - ALOS
#   - COPERNICUS
# ============================================================

# pip install rasterio requests planetary-computer pystac-client elevation whitebox


import os
import shutil
import time
import requests
import rasterio
from pathlib import Path
from rasterio.merge import merge as rio_merge

import planetary_computer as pc
import pystac_client
import elevation

# ------------------------------------------------------------
# DEM Sources
# ------------------------------------------------------------

OPENTOPO_PRODUCT_CODES = {
    "SRTM": "SRTMGL1",
    "NASADEM": "NASADEM",
    "FABDEM": "FABDEM",
    "MERIT": "MERIT",
    "ALOS": "AW3D30",
    "COPERNICUS": "COP30",
}

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

def build_config(
    west,
    south,
    east,
    north,
    dem_source="COPERNICUS",
    opentopo_key="demoapikeyot2022",
    work_dir="./dem_download",
):
    """
    Create configuration for DEM download using a bounding box.

    Parameters
    ----------
    west, south, east, north : float
        Bounding box coordinates (EPSG:4326)

    dem_source : str
        SRTM | NASADEM | FABDEM | MERIT | ALOS | COPERNICUS

    opentopo_key : str
        OpenTopography API key

    work_dir : str
        Directory to save downloaded DEM

    Returns
    -------
    dict
    """

    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    return {
        "bbox": (west, south, east, north),
        "dem_source": dem_source.upper(),
        "opentopo_key": opentopo_key,
        "work_dir": work_dir,
    }

# ------------------------------------------------------------
# OpenTopography Download
# ------------------------------------------------------------

def download_from_opentopo(source, bbox, outfile, api_key):

    west, south, east, north = bbox

    product = OPENTOPO_PRODUCT_CODES[source]

    url = (
        "https://portal.opentopography.org/API/globaldem"
        f"?demtype={product}"
        f"&south={south}"
        f"&north={north}"
        f"&west={west}"
        f"&east={east}"
        f"&outputFormat=GTiff"
        f"&API_Key={api_key}"
    )

    print("Downloading from OpenTopography...")

    r = requests.get(url, stream=True, timeout=300)
    r.raise_for_status()

    with open(outfile, "wb") as f:
        for chunk in r.iter_content(1024 * 1024):
            if chunk:
                f.write(chunk)

# ------------------------------------------------------------
# Planetary Computer (Copernicus)
# ------------------------------------------------------------

def download_copernicus_pc(bbox, outfile):

    west, south, east, north = bbox

    catalog = pystac_client.Client.open(
        "https://planetarycomputer.microsoft.com/api/stac/v1",
        modifier=pc.sign_inplace,
    )

    search = catalog.search(
        collections=["cop-dem-glo-30"],
        bbox=[west, south, east, north],
    )

    items = list(search.items())

    if len(items) == 0:
        raise RuntimeError("No Copernicus tiles found.")

    temp_files = []

    for i, item in enumerate(items):

        item = pc.sign(item)

        asset = next(iter(item.assets.values()))

        tif = outfile.parent / f"tile_{i}.tif"

        r = requests.get(asset.href, stream=True)

        with open(tif, "wb") as f:
            for chunk in r.iter_content(1024 * 1024):
                if chunk:
                    f.write(chunk)

        temp_files.append(str(tif))

    if len(temp_files) == 1:

        shutil.copy(temp_files[0], outfile)

    else:

        srcs = [rasterio.open(i) for i in temp_files]

        mosaic, transform = rio_merge(srcs)

        meta = srcs[0].meta.copy()

        meta.update(
            height=mosaic.shape[1],
            width=mosaic.shape[2],
            transform=transform,
        )

        with rasterio.open(outfile, "w", **meta) as dst:
            dst.write(mosaic)

        for s in srcs:
            s.close()

# ------------------------------------------------------------
# elevation package fallback
# ------------------------------------------------------------

def download_srtm_elevation(bbox, outfile):

    west, south, east, north = bbox

    tmp = outfile.parent / "tmp.tif"

    elevation.clip(
        bounds=(west, south, east, north),
        output=str(tmp),
        product="SRTM3",
    )

    shutil.copy(tmp, outfile)

# ------------------------------------------------------------
# Validation
# ------------------------------------------------------------

def validate_dem(path):

    with rasterio.open(path) as src:

        arr = src.read(1)

        nodata = src.nodata

        if nodata is not None:
            arr = arr[arr != nodata]

        print("Elevation range:")
        print(arr.min(), "to", arr.max(), "m")

        print("Raster size:")
        print(src.width, "x", src.height)

# ------------------------------------------------------------
# Main Download Function
# ------------------------------------------------------------

def download_dem(cfg):

    outfile = cfg["work_dir"] / f"{cfg["dem_source"]}_DEM.tif"

    source = cfg["dem_source"]

    print(f"\nDownloading {source} DEM")

    start = time.time()

    if source == "COPERNICUS":

        try:
            download_copernicus_pc(cfg["bbox"], outfile)
        except:
            download_from_opentopo(
                "COPERNICUS",
                cfg["bbox"],
                outfile,
                cfg["opentopo_key"],
            )

    else:

        try:
            download_from_opentopo(
                source,
                cfg["bbox"],
                outfile,
                cfg["opentopo_key"],
            )

        except:

            print("Primary download failed.")
            print("Trying SRTM fallback...")

            try:
                download_from_opentopo(
                    "SRTM",
                    cfg["bbox"],
                    outfile,
                    cfg["opentopo_key"],
                )

            except:

                print("Trying elevation package...")
                download_srtm_elevation(cfg["bbox"], outfile)

    validate_dem(outfile)

    print(f"\nFinished in {time.time()-start:.1f} sec")

    print("Saved to:", outfile)

    return outfile