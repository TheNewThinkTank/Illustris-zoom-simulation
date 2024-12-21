
import math
from typing import Final
import src.rectangle as rectangle


def test_rectangle_slopes():
    angle: Final = math.radians(70)
    slopes: tuple = (2.747477419454621,
                     2.747477419454621,
                     -0.3639702342662024,
                     -0.3639702342662024
                     )
    assert rectangle.rectangle_slopes(angle) == slopes, f"Should be {slopes}"


def test_intersect():
    intersection: int = 10.5
    assert rectangle.intersect(.5, 3, 12) == intersection, f"Should be {intersection}"


def test_find_point_on_line():
    point_on_line: tuple = (5.47213595499958, 3.23606797749979)
    assert rectangle.find_point_on_line(1,
                                       1,
                                       5,
                                       .5) == point_on_line, \
                                        f"Should be {point_on_line}"
