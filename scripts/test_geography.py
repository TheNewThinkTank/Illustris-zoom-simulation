"""
Unit tests for the geography program
"""

import geography
import math


class TestGeography:

    def test_rectangle_slopes(self):
        assert (2.747477419454621, 2.747477419454621, -0.3639702342662024,
                -0.3639702342662024) == geography.rectangle_slopes(math.radians(70))

    def test_intersect(self):
        assert 10.5 == geography.intersect(.5, 3, 12)

    def test_find_point_on_line(self):
        assert (5.47213595499958, 3.23606797749979) == geography.find_point_on_line(1, 1, 5, .5)
