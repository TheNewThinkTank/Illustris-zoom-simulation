"""
Unit tests for the geometry program
"""

# Standard library
import math

# Local application
import geometry


def test_rectangle_slopes():
    angle = math.radians(70)
    slopes = (2.747477419454621, 2.747477419454621, -0.3639702342662024, -0.3639702342662024)
    assert slopes == geometry.rectangle_slopes(angle), f"Should be {slopes}"


def test_intersect():
    assert 10.5 == geometry.intersect(.5, 3, 12)


def test_find_point_on_line():
    point_on_line = (5.47213595499958, 3.23606797749979)
    assert point_on_line == geometry.find_point_on_line(1, 1, 5, .5), f"Should be {point_on_line}"


if __name__ == "__main__":
    test_rectangle_slopes()
    test_intersect()
    test_find_point_on_line()
    print("Everything passed")
