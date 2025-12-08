import numpy as np
from scipy.optimize import minimize

# import matplotlib.pyplot as plt

"""
Parallelogram Unit Cell:
                 (a,y)             (x+a,y)
                  ----------------------
                /                     /
               /                     /
              ----------------------
              (0,0)               (x,0)
"""


def isInParallelogram(px, py, x, a, y):
    """Check if point (px, py) is inside the parallelogram defined by (0,0), (x,0), (a,y), (x+a,y)"""
    if px < 0 or px > x + a or py < 0 or py > y:
        return False
    # Check if point is below the line from (0,0) to (a,y)
    if py > (y / a) * px:
        return False
    # Check if point is below the line from (x,0) to (x+a)
    if py > (y / a) * (px - x):
        return False
    return True


def distToDisc(point, center, radius):
    """Calculate minimum distance from point to disc edge"""
    px, py = point
    cx, cy = center
    dist_to_center = np.sqrt((px - cx) ** 2 + (py - cy) ** 2)
    return max(0, dist_to_center - radius)


def genIntPointsInCell(x, a, y):
    """Generate all integer grid points within the unit cell"""
    points = []
    # Find bounds for integer points in the parallelogram
    min_i = 0
    max_i = x + a
    min_j = 0
    max_j = y

    for i in range(min_i, max_i + 1):
        for j in range(min_j, max_j + 1):
            if isInParallelogram(i, j, x, a, y):
                points.append((i, j))
    return points

def is_disc_fully_in_cell(center, radius, x, a, y):
    cx, cy = center
    # Check if all points on the circumference are within the cell
    # Test points at angles 0, 90, 180, 270 degrees
    test_points = [(cx + radius, cy), (cx - radius, cy),
                   (cx, cy + radius), (cx, cy - radius)]
    return all(isInParallelogram(px, py, x, a, y) for px, py in test_points)

def are_discs_overlapping(center1, radius1, center2, radius2):
    dist = np.sqrt((center1[0] - center2[0])**2 + (center1[1] - center2[1])**2)
    return dist < radius1 + radius2

#TODO
def calc_mean_points_per_disc(integer_points, num_discs):
    return len(integer_points) / num_discs if num_discs > 0 else 0

def handle_boundary_disc(center, radius, x, a, y):
    # executed after determining that the disc is not fully in the cell
    cx, cy = center
    translated_copies = []
    # Check if disc extends beyond right or left boundaries
    if cx + radius > (a*cy/y)+x:
        translated_copies.append((cx - (x + a), cy))
    if cx - radius < (a*cy/y):
        translated_copies.append((cx + (x + a), cy))
    # Check if disc extends beyond top or bottom boundaries
    if cy + radius > y:
        translated_copies.append((cx, cy - y))
    if cy - radius < 0:
        translated_copies.append((cx, cy + y))
    return translated_copies

def min_distance_to_any_disc(point, discs):
    return min(distToDisc(point, center, radius)
               for center, radius in discs)






"""
TODO:


- create an objective function which can be used to maximize the mean number of points per disc
- create a method which calculates the mean number of points per disc

- create a method which can return the minimum distance between a point and any disc in the cell


DONE:
- create a method which can determine whether the disc is fully in the cell or not
- create a method which can determine whether two discs are overlapping or not
- create a method which is executed after the placing of a disc to check if the disc is fully in the cell and if not (given that it is not the discs centred 
at the vertices) and to duplicate and translate the disc so that the extruding part of the disc is would be on the other side of the cell
"""



