"""
Unit tests for the geometry program
"""

import geometry
import math


class TestGeometry:

    def test_rectangle_slopes(self):
        assert (2.747477419454621, 2.747477419454621, -0.3639702342662024,
                -0.3639702342662024) == geometry.rectangle_slopes(math.radians(70))

    def test_intersect(self):
        assert 10.5 == geometry.intersect(.5, 3, 12)

    def test_find_point_on_line(self):
        assert (5.47213595499958, 3.23606797749979) == geometry.find_point_on_line(1, 1, 5, .5)
