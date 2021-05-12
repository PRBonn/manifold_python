import os
import unittest

import numpy as np
import open3d as o3d

import manifold


def file_relative_path(path):
    """Returns full path relative to this given file."""
    file_path = os.path.realpath(__file__)
    file_dir = os.path.dirname(file_path)
    return os.path.abspath(os.path.join(file_dir, path))


class ManifoldTest(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Use the Delta Airplane used as exaple in the original implementation
        self.model_mesh_path = file_relative_path("../example/input.obj")
        self.input_mesh_model = o3d.io.read_triangle_mesh(self.model_mesh_path)
        self.input_vertices = np.asarray(self.input_mesh_model.vertices)
        self.input_triangles = np.asarray(self.input_mesh_model.triangles)

    def test_is_manifold(self):
        self.assertFalse(
            manifold.is_manifold(self.input_vertices, self.input_triangles)
        )

    def test_manifold_processing(self):
        processor = manifold.Processor(self.input_vertices, self.input_triangles)
        vertices, triangles = processor.get_manifold_mesh()
        self.assertTrue(len(vertices) > 0)
        self.assertTrue(len(triangles) > 0)
        self.assertTrue(manifold.is_manifold(vertices, triangles))

    def test_wrong_input_types(self):
        with self.assertRaises(AssertionError):
            vertices = self.model_mesh_path
            triangles = self.input_mesh_model.vertices
            processor = manifold.Processor(vertices, triangles)
        with self.assertRaises(AssertionError):
            vertices = self.input_vertices.astype(np.float32)
            processor = manifold.Processor(vertices, triangles)
        with self.assertRaises(AssertionError):
            triangles = self.input_triangles.astype(np.int64)
            processor = manifold.Processor(vertices, triangles)
