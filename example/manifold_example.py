#!/usr/bin/env python3
# coding: utf-8
# @file      manifold_example.py
# @author    Ignacio Vizzo     [ivizzo@uni-bonn.de]
#
# Copyright (c) 2020 Ignacio Vizzo, all rights reserved
import os

import click
import numpy as np
import open3d as o3d

import manifold


@click.command()
@click.argument("filename", type=click.Path(exists=True))
@click.option("--depth", type=int, default=8)
def main(filename, depth):
    """Convert a triangular mesh to a watertight model."""
    filename = os.path.abspath(filename)
    model_name = filename.split(".")[0]
    extension = filename.split(".")[1]
    out_filename = model_name + "_manifold." + extension

    print("Reading", filename)
    in_mesh = o3d.io.read_triangle_mesh(filename)

    input_vertices = np.asarray(in_mesh.vertices)
    input_triangles = np.asarray(in_mesh.triangles)
    input_is_manifold = manifold.is_manifold(input_vertices, input_triangles)

    # Create the ManifoldProcessor object
    processor = manifold.Processor(input_vertices, input_triangles)

    print("Converting {} to a manifold".format(filename))
    output_vertices, output_triangles = processor.get_manifold_mesh(depth)
    output_is_manifold = manifold.is_manifold(output_vertices, output_triangles)

    print("Saving results to", out_filename)
    out_mesh = o3d.geometry.TriangleMesh(
        o3d.utility.Vector3dVector(output_vertices),
        o3d.utility.Vector3iVector(output_triangles),
    )
    o3d.io.write_triangle_mesh(out_filename, out_mesh)

    print("Results:")
    print("Input  mesh: {} is manifold? {}".format(in_mesh, input_is_manifold))
    print("Output mesh: {} is manifold? {}".format(out_mesh, output_is_manifold))


if __name__ == "__main__":
    main()
