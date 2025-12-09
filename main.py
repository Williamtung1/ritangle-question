import numpy as np
from scipy.optimize import minimize

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib.patches as mpatches

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
    dist = np.sqrt((center1[0] - center2[0]) ** 2 + (center1[1] - center2[1]) ** 2)
    return dist < radius1 + radius2

def P(integer_points, num_discs):
    return len(integer_points) / num_discs if num_discs > 0 else 0

def Q(discs, x, a, y):
    """
    Fractional coverage of the parallelogram by the circles.

    Parameters:
    - discs: iterable of (center, radius) where center is (cx, cy)
    - x, a, y: parallelogram parameters (area = x * y)

    Deduplication: copies that differ by integer multiples of (x + a, 0) or (0, y)
    are considered the same disc and counted only once.
    """
    width = x + a
    height = y
    cell_area = abs(x * y)

    unique = set()
    for center, radius in discs:
        cx, cy = center
        # map center into the fundamental cell [0, width) x [0, height)
        nx = cx % width
        ny = cy % height
        # round to avoid floating point equality issues
        key = (round(nx, 9), round(ny, 9), round(float(radius), 9))
        unique.add(key)

    total_disc_area = sum(np.pi * (r ** 2) for (_, _, r) in unique)
    return total_disc_area / cell_area



def handle_boundary_disc(center, radius, x, a, y):
    # executed after determining that the disc is not fully in the cell
    cx, cy = center
    translated_copies = []
    # Check if disc extends beyond right or left boundaries
    if cx + radius > (a * cy / y) + x:
        translated_copies.append((cx - (x + a), cy))
    if cx - radius < (a * cy / y):
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

def all_points_covered(integer_points, discs):
    """Check if all integer points are covered by at least one disc"""
    for point in integer_points:
        if min_distance_to_any_disc(point, discs) > 0:
            return False
    return True


def optimize_three_radii_for_P(x, a, y, r1):
    """
    Find optimal r2 and r3 that maximize P (points per disc).

    Parameters:
    - x, a, y: parallelogram parameters
    - r1: fixed radius at vertices

    Returns:
    - (r2, r3, optimal_discs): optimal radii and resulting disc placement
    """
    integer_points = genIntPointsInCell(x, a, y)

    def objective(radii):
        r2, r3 = radii

        # Enforce ordering constraint: r1 > r2 > r3 > 0
        if r2 >= r1 or r3 >= r2 or r3 <= 0.1:
            return 1e10  # Large penalty for invalid ordering

        # Place discs with fixed radii
        discs = place_discs_with_fixed_radii(x, a, y, r1, r2, r3)

        # Check if all points are covered
        if not all_points_covered(integer_points, discs):
            return 1e10  # Penalty for incomplete coverage

        # Calculate P and return negative (since we minimize)
        num_discs = len(discs)
        p_value = P(integer_points, num_discs)
        return -p_value  # Maximize P = minimize -P

    # Initial guess: r2 = 0.7*r1, r3 = 0.4*r1
    initial_guess = [r1 * 0.7, r1 * 0.4]

    # Bounds: r3 < r2 < r1
    bounds = [(0.1, r1 - 0.01), (0.1, r1 * 0.7)]

    # Optimize
    result = minimize(objective, initial_guess, method='Nelder-Mead', bounds=bounds)

    optimal_r2, optimal_r3 = result.x
    optimal_discs = place_discs_with_fixed_radii(x, a, y, r1, optimal_r2, optimal_r3)

    return optimal_r2, optimal_r3, optimal_discs


def place_discs_with_fixed_radii(x, a, y, r1, r2, r3):
    """
    Place discs with three fixed radii using greedy strategy.

    Parameters:
    - x, a, y: parallelogram parameters
    - r1, r2, r3: three radii (r1 > r2 > r3)

    Returns:
    - discs: list of (center, radius) tuples
    """
    integer_points = genIntPointsInCell(x, a, y)
    discs = []

    # Step 1: Place r1 discs at vertices
    vertices = [(0, 0), (x, 0), (a, y), (x + a, y)]
    for vertex in vertices:
        discs.append((vertex, r1))

    uncovered = [p for p in integer_points if min_distance_to_any_disc(p, discs) > 0]

    # Step 2: Greedily place r2 discs
    while uncovered:
        candidate = max(uncovered, key=lambda p: min_distance_to_any_disc(p, discs))

        # Try to place r2 disc
        if can_place_disc(candidate, r2, discs):
            discs.append((candidate, r2))
            if not is_disc_fully_in_cell(candidate, r2, x, a, y):
                copies = handle_boundary_disc(candidate, r2, x, a, y)
                discs.extend((copy, r2) for copy in copies)
            uncovered = [p for p in uncovered if min_distance_to_any_disc(p, discs) > 0]
        else:
            break  # Cannot place r2, move to r3

    # Step 3: Fill remaining gaps with r3 discs
    while uncovered:
        candidate = max(uncovered, key=lambda p: min_distance_to_any_disc(p, discs))

        if can_place_disc(candidate, r3, discs):
            discs.append((candidate, r3))
            if not is_disc_fully_in_cell(candidate, r3, x, a, y):
                copies = handle_boundary_disc(candidate, r3, x, a, y)
                discs.extend((copy, r3) for copy in copies)
            uncovered = [p for p in uncovered if min_distance_to_any_disc(p, discs) > 0]
        else:
            # Cannot place disc without overlap
            uncovered.remove(candidate)

    return discs


def can_place_disc(center, radius, existing_discs):
    """Check if disc at center with radius can be placed without overlap"""
    for existing_center, existing_radius in existing_discs:
        dist = np.sqrt((center[0] - existing_center[0])**2 +
                      (center[1] - existing_center[1])**2)
        if dist < radius + existing_radius:
            return False
    return True


def render_discs(discs, x, a, y, integer_points):
    """
    Render the parallelogram cell with discs and integer points.
    """
    # use a square-ish figure so aspect locking behaves predictably
    fig, ax = plt.subplots(figsize=(8, 8))

    # Draw parallelogram boundary
    parallelogram = mpatches.Polygon(
        [(0, 0), (x, 0), (x + a, y), (a, y)],
        fill=False,
        edgecolor='black',
        linewidth=2,
        label='Cell boundary'
    )
    ax.add_patch(parallelogram)

    # Get unique radii for color coding
    unique_radii = sorted(set(r for _, r in discs), reverse=True)
    colors = ['red', 'blue', 'green']  # r1, r2, r3
    radius_to_color = {r: colors[i] for i, r in enumerate(unique_radii)}

    # Draw discs
    for center, radius in discs:
        cx, cy = center
        # Only draw discs whose centers are in the fundamental cell
        if 0 <= cx <= x + a and 0 <= cy <= y:
            circle = Circle(
                (cx, cy),
                radius,
                fill=False,
                edgecolor=radius_to_color.get(radius, 'gray'),
                linewidth=1.5,
                alpha=0.7
            )
            ax.add_patch(circle)

    # Draw integer points
    int_x = [p[0] for p in integer_points]
    int_y = [p[1] for p in integer_points]
    ax.scatter(int_x, int_y, c='black', s=20, marker='o',
               label='Integer points', zorder=5)

    # Draw vertices
    vertices = [(0, 0), (x, 0), (a, y), (x + a, y)]
    vert_x = [v[0] for v in vertices]
    vert_y = [v[1] for v in vertices]
    ax.scatter(vert_x, vert_y, c='purple', s=100, marker='s',
               label='Vertices', zorder=6)

    # Create legend for radii
    legend_elements = [
        mpatches.Patch(color='black', label='Cell boundary'),
        mpatches.Patch(color='purple', label='Vertices')
    ]
    for i, r in enumerate(unique_radii):
        legend_elements.append(
            mpatches.Patch(color=colors[i],
                           label=f'r{i + 1} = {r:.3f}')
        )
    legend_elements.append(
        mpatches.Patch(color='black', label='Integer points')
    )

    ax.legend(handles=legend_elements, loc='upper right')

    # Set axis properties
    x_min, x_max = -0.5, x + a + 0.5
    y_min, y_max = -0.5, y + 0.5
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    # prefer set_box_aspect to enforce equal data-unit scaling (width : height = xspan : yspan)
    xspan = x_max - x_min
    yspan = y_max - y_min
    try:
        ax.set_box_aspect((xspan, yspan))
    except Exception:
        # fallback for older matplotlib versions
        ax.set_aspect('equal', adjustable='box')

    ax.grid(True, alpha=0.3)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title(f'Optimal Disc Placement (P = {P(integer_points, len(discs)):.4f})',
                 fontsize=14)

    plt.tight_layout()
    plt.show()


x, a, y = 10, 4, 6  # Example parallelogram
r1 = 3   # Fixed vertex radius

optimal_r2, optimal_r3, optimal_discs = optimize_three_radii_for_P(x, a, y, r1)

print(f"Optimal r2: {optimal_r2:.4f}")
print(f"Optimal r3: {optimal_r3:.4f}")
print(f"Number of discs: {len(optimal_discs)}")
print(f"The discs are: {optimal_discs}")
render_discs(optimal_discs, x, a, y, genIntPointsInCell(x, a, y))
print(f"P (points per disc): {P(genIntPointsInCell(x, a, y), len(optimal_discs)):.4f}")
"""
TODO:
- create a greedy algorithm with the constraint that no discs are overlapping and that all 
integer points are covered and to set the two other radii given the one centred at the vertices,
that the two other radii will be smaller than the one at the vertices, using the min_distance_to_any_disc method
calculate the largest possible radius for one of the disc in the insides of the cell and fill the rest of the cell with discs of that radius and with the disc with the smallest radius as well

DONE:

- create a method which can return the minimum distance between a point and any disc in the cell
- create a method to check if all the integers points are covered by at least one disc
- create a method which can determine whether the disc is fully in the cell or not
- create a method which can determine whether two discs are overlapping or not
- create a method which is executed after the placing of a disc to check if the disc is fully in the cell and if not (given that it is not the discs centred 
at the vertices) and to duplicate and translate the disc so that the extruding part of the disc is would be on the other side of the cell
- create a method which calculates the mean number of points per disc
- create an objective function which can be used to maximize the mean number of points per disc
"""
