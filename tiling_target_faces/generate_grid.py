from __future__ import annotations

from pathlib import Path
from typing import Optional
import numpy as np
import xarray as xr

try:
    import pygeohash as pgh
except Exception as e:
    raise ImportError(
        "pygeohash is required. Install with: pip install pygeohash"
    ) from e


def _format_stretch_factor(sf: float) -> str:
    """
    Format stretch factor like the Bash version:
    1.5 -> '1d50' (i.e., two decimals, then '.' -> 'd').
    """
    s = f"{sf:.2f}"
    return s.replace(".", "d")


def get_gridspec_prefix(
    cs_res,
    stretch_grid,
    stretch_factor: float = 1.0,
    target_lat: float = -90.0,
    target_lon: float = 170.0,
    geohash_precision=12
) -> str:
    """
    Build the GridSpec prefix.

    Examples
    --------
    >>> get_gridspec_prefix(48, False)
    'c48'
    >>> get_gridspec_prefix(36, True, 1.5, 32.0, -103.0)
    'c36_s1d50_t<geohash>'
    """
    if not isinstance(cs_res, int) or cs_res <= 0:
        raise ValueError("cs_res must be a positive integer.")

    if not stretch_grid:
        return f"c{cs_res}"

    sf_formatted = _format_stretch_factor(float(stretch_factor))
    target_lon = ((float(target_lon) + 180) % 360) - 180

    geoh = pgh.encode(float(target_lat), target_lon, precision=geohash_precision)
    return f"c{cs_res}_s{sf_formatted}_t{geoh}"


def generate_grid_from_gridspec(
    cs_res: int,
    output_file: str | Path,
    stretch_grid: bool,
    stretch_factor: float = 1.0,
    target_lat: float = -90.0,
    target_lon: float = 170.0,
    workdir: str | Path = ".",
    prefix: Optional[str] = None,
    geohash_precision: int = 12,
    validate_inputs: bool = True,
    save: bool = True,
) -> xr.Dataset:
    """
    Create a cubed-sphere (CS or SGCS) grid NetCDF by combining 6 tile files.

    Parameters
    ----------
    cs_res : int
        Cubed-sphere resolution (e.g., 48, 90, 180).
    output_file : str | Path
        Path to write the combined NetCDF (e.g., 'c48_grid.nc').
    stretch_grid : bool
        Whether this is a stretched grid (SGCS).
    stretch_factor : float, default 1.0
        Stretch factor when `stretch_grid=True`.
    target_lat, target_lon : float
        Target point for stretch (lat in [-90,90], lon usually in [-180,180] or [0,360)).
    workdir : str | Path, default '.'
        Directory where tile files are located (and where to check the gridspec file).
    prefix : str | None
        If provided, use this exact prefix. Otherwise it will be generated from args.
    geohash_precision : int, default 7
        Precision for the geohash in the prefix if `prefix` is not provided.
    validate_inputs : bool, default True
        If True, raise with helpful messages when inputs/files are missing.
    save : bool, default True
        If True, writes to `output_file`. Returns the Dataset either way.

    Returns
    -------
    xr.Dataset
        Dataset with variables: area, corner_lons, corner_lats
        and coordinate-like variables: lats, lons.

    Notes
    -----
    - Expects per-tile files named: '{prefix}.tile{1..6}.nc'
      containing variables 'lats', 'lons', and 'area'.
    - Also expects a '{prefix}_gridspec.nc' file to exist, matching your prior check.
    - Longitudes are mapped to [0, 360) to match your original behavior.
    """
    workdir = Path(workdir)
    output_file = Path(output_file)

    if prefix is None:
        prefix = get_gridspec_prefix(
            cs_res,
            stretch_grid,
            stretch_factor,
            target_lat,
            target_lon,
            geohash_precision=geohash_precision,
        )

    if validate_inputs:
        if not isinstance(cs_res, int) or cs_res <= 0:
            raise ValueError("cs_res must be a positive integer.")

    # Preallocate arrays
    lon = np.zeros((6, cs_res, cs_res), dtype=np.float64)
    lat = np.zeros((6, cs_res, cs_res), dtype=np.float64)
    lon_b = np.zeros((6, cs_res + 1, cs_res + 1), dtype=np.float64)
    lat_b = np.zeros((6, cs_res + 1, cs_res + 1), dtype=np.float64)
    area = np.zeros((6, cs_res, cs_res), dtype=np.float64)

    # Check gridspec file presence (as in your original script)
    gridspec_path = workdir / f"{prefix}_gridspec.nc"
    if validate_inputs and not gridspec_path.is_file():
        raise FileNotFoundError(
            f"ERROR: {gridspec_path} does not exist. "
            f"Use GridSpec to generate it first!"
        )

    # Read each tile
    for i in range(6):
        tile_idx = i + 1
        tile_path = workdir / f"{prefix}.tile{tile_idx}.nc"
        if validate_inputs and not tile_path.is_file():
            raise FileNotFoundError(
                f"ERROR: Tile file {tile_path} does not exist!"
            )

        ds = xr.open_dataset(tile_path)
        try:
            # Expecting shapes: lats/lons -> (2*cs_res+1, 2*cs_res+1), area -> (cs_res, cs_res)
            super_lats = ds["lats"].values
            super_lons = ds["lons"].values
            tile_area = ds["area"].values

            # Sample the supergrid corners/centers like your Bash/Python snippet
            lat_b[i] = super_lats[::2, ::2]
            lon_b[i] = super_lons[::2, ::2]
            lat[i] = super_lats[1::2, 1::2]
            lon[i] = super_lons[1::2, 1::2]
            area[i] = tile_area
        finally:
            ds.close()

    # Adjust longitudes to [0, 360)
    lon = np.where(lon < 0.0, lon + 360.0, lon)
    lon_b = np.where(lon_b < 0.0, lon_b + 360.0, lon_b)

    # Build xarray structures
    area_da = xr.DataArray(
        area,
        dims=["nf", "Ydim", "Xdim"],
        coords={
            "lats": (["nf", "Ydim", "Xdim"], lat),
            "lons": (["nf", "Ydim", "Xdim"], lon),
        },
        attrs=dict(units="m2", long_name="Surface area of each grid box"),
    )

    corner_lons_da = xr.DataArray(
        lon_b,
        dims=["nf", "YCdim", "XCdim"],
        attrs=dict(units="degrees_east", long_name="Longitude"),
    )

    corner_lats_da = xr.DataArray(
        lat_b,
        dims=["nf", "YCdim", "XCdim"],
        attrs=dict(units="degrees_north", long_name="Latitude"),
    )

    data = xr.Dataset(
        data_vars={
            "area": area_da,
            "corner_lons": corner_lons_da,
            "corner_lats": corner_lats_da,
        }
    )

    # Attach coordinate metadata (to match your original)
    data["lats"].attrs = dict(units="degrees_north", long_name="Latitude")
    data["lons"].attrs = dict(units="degrees_east", long_name="Longitude")

    # Global attrs for stretched grids
    if stretch_grid:
        data.attrs["STRETCH_FACTOR"] = np.float32(stretch_factor)
        data.attrs["TARGET_LAT"] = np.float32(target_lat)
        # Keep lon in [0,360) if provided that way
        t_lon = float(target_lon)
        data.attrs["TARGET_LON"] = np.float32(t_lon if t_lon >= 0 else (t_lon % 360))

    if save:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        data.to_netcdf(output_file)
        print(f"Combined c{cs_res} cubed-sphere grid data saved to {output_file}")

    return data


# -------- Optional CLI shim (mirrors your Bash usage) --------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Create CS/SGCS grid NetCDF by combining 6 tile files."
    )
    parser.add_argument("CS_RES", type=int, help="Cubed-sphere resolution, e.g., 48")
    parser.add_argument("OUTPUT_FILE", type=Path, help="Path to output NetCDF")
    parser.add_argument(
        "STRETCH_GRID",
        type=str,
        help="True/False (case-insensitive) for stretched grid",
    )
    parser.add_argument("--stretch-factor", type=float, default=1.0)
    parser.add_argument("--target-lat", type=float, default=-90.0)
    parser.add_argument("--target-lon", type=float, default=170.0)
    parser.add_argument("--workdir", type=Path, default=Path("."))
    parser.add_argument("--prefix", type=str, default=None)
    parser.add_argument("--geohash-precision", type=int, default=7)
    parser.add_argument("--no-save", action="store_true", help="Do not write NetCDF")

    args = parser.parse_args()
    stretch_grid_bool = args.STRETCH_GRID.lower() == "true"

    if args.prefix is None:
        pref = get_gridspec_prefix(
            args.CS_RES,
            stretch_grid_bool,
            args.stretch_factor,
            args.target_lat,
            args.target_lon,
            geohash_precision=args.geohash_precision,
        )
    else:
        pref = args.prefix

    generate_grid_from_gridspec(
        args.CS_RES,
        args.OUTPUT_FILE,
        stretch_grid_bool,
        args.stretch_factor,
        args.target_lat,
        args.target_lon,
        workdir=args.workdir,
        prefix=pref,
        geohash_precision=args.geohash_precision,
        save=not args.no_save,
    )