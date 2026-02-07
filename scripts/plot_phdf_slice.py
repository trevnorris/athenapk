#!/usr/bin/env python3
"""Quicklook 2D slice plot for Parthenon phdf cell-centered datasets."""

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


def assemble_uniform_level_field(dataset, logical_locations, levels, component):
    """Assemble block-local data into a global array for single-level runs."""
    data = np.asarray(dataset)
    if data.ndim != 5:
        raise RuntimeError(f"Expected 5D dataset (blocks, comps, z, y, x), got {data.shape}")
    if component < 0 or component >= data.shape[1]:
        raise RuntimeError(f"Component {component} out of bounds for dataset with {data.shape[1]} comps")

    levels = np.asarray(levels)
    if not np.all(levels == levels[0]):
        raise RuntimeError("This quicklook currently supports single-level outputs only")

    logical = np.asarray(logical_locations)
    if logical.shape[1] != 3:
        raise RuntimeError(f"Unexpected LogicalLocations shape: {logical.shape}")

    nb, _, nz, ny, nx = data.shape
    if logical.shape[0] != nb:
        raise RuntimeError("LogicalLocations block count mismatch")

    max_lx = int(np.max(logical[:, 0]))
    max_ly = int(np.max(logical[:, 1]))
    max_lz = int(np.max(logical[:, 2]))

    global_field = np.zeros(((max_lz + 1) * nz, (max_ly + 1) * ny, (max_lx + 1) * nx), dtype=data.dtype)
    for b in range(nb):
        lx, ly, lz = [int(v) for v in logical[b]]
        i0 = lx * nx
        j0 = ly * ny
        k0 = lz * nz
        global_field[k0 : k0 + nz, j0 : j0 + ny, i0 : i0 + nx] = data[b, component]
    return global_field


def assemble_axis_coords(volume_locations, logical_locations, global_shape):
    logical = np.asarray(logical_locations)
    x_blk = np.asarray(volume_locations["x"])
    y_blk = np.asarray(volume_locations["y"])
    z_blk = np.asarray(volume_locations["z"])

    x = np.zeros(global_shape[2], dtype=x_blk.dtype)
    y = np.zeros(global_shape[1], dtype=y_blk.dtype)
    z = np.zeros(global_shape[0], dtype=z_blk.dtype)

    nb = logical.shape[0]
    nx = x_blk.shape[1]
    ny = y_blk.shape[1]
    nz = z_blk.shape[1]
    for b in range(nb):
        lx, ly, lz = [int(v) for v in logical[b]]
        x[lx * nx : (lx + 1) * nx] = x_blk[b]
        y[ly * ny : (ly + 1) * ny] = y_blk[b]
        z[lz * nz : (lz + 1) * nz] = z_blk[b]
    return x, y, z


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--phdf", required=True, help="Path to .phdf file")
    p.add_argument("--dataset", default="em4d_a", help="Cell-centered dataset name")
    p.add_argument("--component", type=int, default=0, help="Component index in dataset")
    p.add_argument("--axis", choices=["x", "y", "z"], default="z", help="Normal axis of slice")
    p.add_argument("--index", type=int, default=None, help="Slice index (default: center)")
    p.add_argument("--output", default="phdf_slice.png", help="Output PNG file")
    p.add_argument("--cmap", default="viridis", help="Matplotlib colormap")
    p.add_argument("--title", default=None, help="Custom plot title")
    p.add_argument("--absolute", action="store_true", help="Plot absolute value")
    args = p.parse_args()

    phdf = Path(args.phdf)
    if not phdf.exists():
        raise SystemExit(f"phdf file not found: {phdf}")

    with h5py.File(phdf, "r") as f:
        if args.dataset not in f:
            raise SystemExit(f"dataset '{args.dataset}' not found in {phdf}")
        if "LogicalLocations" not in f or "Levels" not in f or "VolumeLocations" not in f:
            raise SystemExit("phdf missing required mesh metadata groups/datasets")

        field = assemble_uniform_level_field(
            f[args.dataset], f["LogicalLocations"], f["Levels"], args.component
        )
        x, y, z = assemble_axis_coords(f["VolumeLocations"], f["LogicalLocations"], field.shape)

    if args.axis == "z":
        idx = field.shape[0] // 2 if args.index is None else args.index
        img = field[idx, :, :]
        x_axis, y_axis = x, y
        xlabel, ylabel = "x", "y"
    elif args.axis == "y":
        idx = field.shape[1] // 2 if args.index is None else args.index
        img = field[:, idx, :]
        x_axis, y_axis = x, z
        xlabel, ylabel = "x", "z"
    else:
        idx = field.shape[2] // 2 if args.index is None else args.index
        img = field[:, :, idx]
        x_axis, y_axis = y, z
        xlabel, ylabel = "y", "z"

    if idx < 0:
        raise SystemExit("slice index must be non-negative")
    if args.axis == "z" and idx >= field.shape[0]:
        raise SystemExit(f"z slice index out of range: {idx}")
    if args.axis == "y" and idx >= field.shape[1]:
        raise SystemExit(f"y slice index out of range: {idx}")
    if args.axis == "x" and idx >= field.shape[2]:
        raise SystemExit(f"x slice index out of range: {idx}")

    if args.absolute:
        img = np.abs(img)

    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    extent = [float(np.min(x_axis)), float(np.max(x_axis)), float(np.min(y_axis)), float(np.max(y_axis))]
    im = ax.imshow(img, origin="lower", cmap=args.cmap, extent=extent, aspect="auto")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if args.title is None:
        title = f"{args.dataset}[comp={args.component}] slice {args.axis}={idx}"
    else:
        title = args.title
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label=args.dataset)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=150)
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
