# python
import pytest
from main import isInParallelogram

@pytest.mark.parametrize("coord,x,a,y,expected", [
    # inside examples
    ([1, 1], 5, 2, 3, True),
    ([0, 0], 5, 2, 3, True),        # vertex
    ([5 + 2, 3], 5, 2, 3, True),    # opposite vertex (x+a, y)

    # on-edge examples
    ([2, 0], 5, 2, 3, True),        # on bottom edge
    ([1, int((3 / 2) * 1)], 5, 2, 3, True),  # on slanted edge (approx)
    ([3, 3], 5, 2, 3, True),

    # outside examples
    ([-1, 1], 5, 2, 3, False),
    ([10, 10], 5, 2, 3, False),
    ([3, 4], 5, 2, 3, False)
])
def test_is_in_parallelogram(coord, x, a, y, expected):
    assert isInParallelogram(coord, x, a, y) == expected


import pytest
import numpy as np
from main import distToDisc

@pytest.mark.parametrize("point,disc,expected", [
    # Point at disc center (inside)
    ((0, 0), ((0, 0), 1), 0),
    # Point on circumference (distance 0)
    ((1, 0), ((0, 0), 1), 0),
    # Point inside disc
    ((0.5, 0), ((0, 0), 1), 0),
    # Point outside disc
    ((2, 0), ((0, 0), 1), 1),
    # Point far outside
    ((5, 0), ((0, 0), 1), 4),
    # Non-origin center
    ((3, 4), ((1, 2), 1), np.sqrt((3-1)**2 + (4-2)**2) - 1),
    # Negative coordinates
    ((-1, -1), ((0, 0), 1), np.sqrt(2) - 1),
])
def test_distToDisc(point, disc, expected):
    result = distToDisc(point, disc)
    assert abs(result - expected) < 1e-6  # Floating point tolerance

# Add to main_test.py
import pytest
import numpy as np
from main import min_distance_to_any_disc

@pytest.mark.parametrize("point,discs,expected", [
    # Point at center of one disc (distance 0)
    ((0, 0), [((0, 0), 1)], 0),
    # Point on circumference of one disc (distance 0)
    ((1, 0), [((0, 0), 1)], 0),
    # Point inside one disc
    ((0.5, 0), [((0, 0), 1)], 0),
    # Point outside single disc
    ((3, 0), [((0, 0), 1)], 2),
    # Point between two discs (closer to first)
    ((1.5, 0), [((0, 0), 1), ((4, 0), 1)], 0.5),
    # Point inside one disc, outside another
    ((0.5, 0), [((0, 0), 1), ((5, 0), 1)], 0),
    # Point outside all discs (multiple discs)
    ((10, 10), [((0, 0), 1), ((3, 0), 1), ((0, 3), 1)], np.sqrt(10**2 + 7**2) - 1),
    # Point equidistant to two discs
    ((2, 0), [((0, 0), 1), ((4, 0), 1)], 1),
    # Empty discs list (edge case - will cause error, so omit or handle separately)
    # Point with negative coordinates
    ((-5, -5), [((0, 0), 1), ((2, 2), 0.5)], np.sqrt(5**2 + 5**2) - 1),
    # Multiple discs, point inside smallest
    ((5, 5), [((0, 0), 2), ((5, 5), 0.5), ((10, 10), 3)], 0),
])
def test_min_distance_to_any_disc(point, discs, expected):
    result = min_distance_to_any_disc(point, discs)
    assert abs(result - expected) < 1e-6  # Floating point tolerance

# Add to main_test.py
import pytest
from main import is_disc_fully_in_cell

@pytest.mark.parametrize("disc,qx,qa,qy,expected", [
    # Disc fully inside cell (center at middle, small radius)
    (((2.5, 1.5), 0.5), 5, 2, 3, True),
    # Disc at vertex (0,0) - partially outside
    (((0, 0), 1.0), 5, 2, 3, False),
    # Disc at vertex (qx,0) - partially outside
    (((5, 0), 1.0), 5, 2, 3, False),
    # Disc at vertex (qa,qy) - partially outside
    (((2, 3), 1.0), 5, 2, 3, False),
    # Disc at vertex (qx+qa,qy) - partially outside
    (((7, 3), 1.0), 5, 2, 3, False),
    # Very small disc at center (always inside)
    (((3, 1.5), 0.1), 5, 2, 3, True),
    # Large disc covering entire cell (extends outside)
    (((3, 1.5), 5.0), 5, 2, 3, False),
    # Disc touching left edge but inside
    (((1.0, 1.5), 0.5), 5, 2, 3, False),
    # Disc touching right edge but inside
    (((6.0, 1.5), 0.5), 5, 2, 3, True),
    # Disc touching bottom edge but inside
    (((3, 0.5), 0.4), 5, 2, 3, True),
    # Disc touching top edge but inside
    (((3, 2.5), 0.4), 5, 2, 3, True),
    # Disc extending beyond left slanted edge
    (((0.5, 1.0), 0.8), 5, 2, 3, False),
    # Disc extending beyond right slanted edge
    (((6.5, 2.0), 0.8), 5, 2, 3, False),
    # Rectangle case (qa=0): disc fully inside
    (((2, 1.5), 0.8), 5, 0, 3, True),
    # Rectangle case (qa=0): disc extending beyond right edge
    (((4.5, 1.5), 0.8), 5, 0, 3, False),
])
def test_is_disc_fully_in_cell(disc, qx, qa, qy, expected):
    result = is_disc_fully_in_cell(disc, qx, qa, qy)
    assert result == expected





# Add to main_test.py
import pytest
from main import genIntPointsInCell, isInParallelogram


@pytest.mark.parametrize("qx,qa,qy,description", [
    # Rectangle cases (qa=0)
    (3, 0, 2, "Simple rectangle 3x2"),
    (1, 0, 1, "Unit square"),
    (5, 0, 3, "Rectangle 5x3"),

    # Parallelogram cases with positive slant
    (3, 2, 2, "Parallelogram with qa=2"),
    (5, 2, 3, "Parallelogram 5+2 by 3"),
    (4, 3, 4, "Parallelogram 4+3 by 4"),

    # Edge cases
    (1, 1, 1, "Minimal parallelogram"),
    (10, 5, 8, "Larger parallelogram"),
])
def test_genIntPointsInCell(qx, qa, qy, description):
    """Test that genIntPointsInCell returns correct integer points inside parallelogram"""
    result = genIntPointsInCell(qx, qa, qy)

    # All returned points must be integers
    assert all(isinstance(x, int) and isinstance(y, int) for x, y in result), \
        f"{description}: Non-integer coordinates found"

    # All returned points must be inside parallelogram
    assert all(isInParallelogram((x, y), qx, qa, qy) for x, y in result), \
        f"{description}: Points outside parallelogram found"

    # Verify corner points are included
    assert (0, 0) in result, f"{description}: Missing vertex (0, 0)"
    assert (qx, 0) in result, f"{description}: Missing vertex (qx, 0)"
    assert (qa, qy) in result, f"{description}: Missing vertex (qa, qy)"
    assert (qx + qa, qy) in result, f"{description}: Missing vertex (qx+qa, qy)"

    # Check no duplicate points
    assert len(result) == len(set(result)), \
        f"{description}: Duplicate points found"

    # Verify minimum count (at least 4 vertices)
    assert len(result) >= 4, \
        f"{description}: Less than 4 points (should include at least vertices)"

    # For rectangles (qa=0), verify exact count
    if qa == 0:
        expected_count = (qx + 1) * (qy + 1)
        assert len(result) == expected_count, \
            f"{description}: Rectangle should have {expected_count} points, got {len(result)}"


@pytest.mark.parametrize("qx,qa,qy,point,should_include", [
    # Test specific points that should/shouldn't be included
    (5, 2, 3, (3, 1), True),  # Middle of parallelogram
    (5, 2, 3, (0, 0), True),  # Bottom-left vertex
    (5, 2, 3, (5, 0), True),  # Bottom-right vertex
    (5, 2, 3, (2, 3), True),  # Top-left vertex
    (5, 2, 3, (7, 3), True),  # Top-right vertex
    (5, 2, 3, (-1, 1), False),  # Outside left
    (5, 2, 3, (8, 3), False),  # Outside right
    (5, 2, 3, (3, 4), False),  # Above top edge
    (5, 2, 3, (3, -1), False),  # Below bottom edge
])
def test_genIntPointsInCell_specific_points(qx, qa, qy, point, should_include):
    """Test that specific points are correctly included or excluded"""
    result = genIntPointsInCell(qx, qa, qy)

    if should_include:
        assert point in result, \
            f"Point {point} should be inside parallelogram ({qx}, {qa}, {qy})"
    else:
        assert point not in result, \
            f"Point {point} should be outside parallelogram ({qx}, {qa}, {qy})"


def test_genIntPointsInCell_exhaustive_check():
    """Exhaustively verify no points are missed or incorrectly included"""
    qx, qa, qy = 5, 2, 3
    result = genIntPointsInCell(qx, qa, qy)

    # Manually check all candidate points
    expected = []
    for x in range(-1, qx + qa + 2):  # Check beyond bounds
        for y in range(-1, qy + 2):
            if isInParallelogram((x, y), qx, qa, qy):
                expected.append((x, y))

    assert sorted(result) == sorted(expected), \
        f"Mismatch between expected and actual points.\n" \
        f"Missing: {set(expected) - set(result)}\n" \
        f"Extra: {set(result) - set(expected)}"
