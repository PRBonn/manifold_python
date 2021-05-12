from typing import Tuple

import numpy as np

from .pybind import manifold_pybind


class Processor:
    """Wrapper class around the low level C++ python bindings.

    This wrapper us to better interact with numpy arrays to make the
    client code more simple and readable. I tried to fix this from the
    C++ API and, while it's doable(check Open3D) the resulting pybind11
    code looks more dense and harder to mantain. I tought it was better
    to just wrap the numpy arrays to custom vector types and use this
    class as a proxy to acess those custom types, namely VectorEigen3d
    and VectorEigen3i. The user will only need to interact using
    numpy.ndarrays.
    """

    def __init__(self, vertices: np.ndarray, triangles: np.ndarray):
        """Creates a Manifold Mesh Processor."""
        assert isinstance(vertices, np.ndarray), "vertices must by np.ndarray(n, 3)"
        assert isinstance(triangles, np.ndarray), "triangles must by np.ndarray(n, 3)"
        assert vertices.dtype == np.float64, "vertices dtype must be np.float64"
        assert triangles.dtype == np.int32, "triangles dtype must be np.int32"
        self.processor = manifold_pybind._Processor(
            manifold_pybind._VectorEigen3d(vertices),
            manifold_pybind._VectorEigen3i(triangles),
        )

    def get_manifold_mesh(self, depth=8) -> Tuple[np.ndarray, np.ndarray]:
        """Obtain the manifold mesh."""
        vertices, triangles = self.processor._get_manifold_mesh(depth)
        return np.asarray(vertices), np.asarray(triangles)


def is_manifold(vertices: np.ndarray, triangles: np.ndarray) -> bool:
    """Check wheter the specified mesh by its vertices and triangles is
    watertight or not."""
    assert isinstance(vertices, np.ndarray), "vertices must by np.ndarray(n, 3)"
    assert isinstance(triangles, np.ndarray), "triangles must by np.ndarray(n, 3)"
    assert vertices.dtype == np.float64, "vertices dtype must be np.float64"
    assert triangles.dtype == np.int32, "triangles dtype must be np.int32"
    return manifold_pybind._Processor._is_manifold(
        manifold_pybind._VectorEigen3d(vertices),
        manifold_pybind._VectorEigen3i(triangles),
    )
