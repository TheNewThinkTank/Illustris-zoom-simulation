"""
Unit tests for the geometry program
"""

# Standard library
import math
from typing import Final, Tuple

# Local application
import geometry


def test_rectangle_slopes():
    angle: Final = math.radians(70)
    slopes: Tuple = (2.747477419454621, 2.747477419454621, -0.3639702342662024, -0.3639702342662024)
    assert geometry.rectangle_slopes(angle) == slopes, f"Should be {slopes}"


def test_intersect():
    intersection: int = 10.5
    assert geometry.intersect(.5, 3, 12) == intersection, f"Should be {intersection}"


def test_find_point_on_line():
    point_on_line: Tuple = (5.47213595499958, 3.23606797749979)
    assert geometry.find_point_on_line(1, 1, 5, .5) == point_on_line, f"Should be {point_on_line}"

