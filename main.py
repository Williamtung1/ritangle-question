import csv

import numpy as np

"""
Parallelogram Unit Cell:
                 (qa,qy)             (qx+qa,qy)
                  ----------------------
                /                     /
               /                     /
              ----------------------
              (0,0)               (qx,0)
"""


def isInParallelogram(coord, qx, qa, qy):
    px, py = coord
    flag = False

    # Points must be within vertical bounds
    if 0 <= py <= qy and 0 <= px:
        # Points in the middle vertical strip (no slant)
        if qa <= px <= qx:
            flag = True
        # Left region: points must be above the line from (0,0) to (qa, qy)
        elif px < qa:
            flag = (py <= (qy / qa) * px)
        # Right region: points must be above the line from (qx, 0) to (qx+qa, qy)
        elif px > qx:
            flag = (py >= (qy / qa) * (px - qx))
    return flag


def distToDisc(point, disc):  # shortest distance from a point to the circumference of a disc
    px, py = point
    (cx, cy), radius = disc
    dist = np.sqrt((px - cx) ** 2 + (py - cy) ** 2)
    return max(0, dist - radius)


def min_distance_to_any_disc(point, discs):
    distances = []
    for center, radius in discs:
        d = distToDisc(point, (center, radius))
        distances.append(d)

    # return the smallest distance (0 if point is inside/on any disc)
    return min(distances)


def genIntPointsInCell(qx, qa, qy):
    points = []
    for i in range(qx + qa + 1):
        for j in range(qy + 1):
            if isInParallelogram((i, j), qx, qa, qy):
                points.append((i, j))
    return points


def is_disc_fully_in_cell(disc, qx, qa, qy):
    (cx, cy), radius = disc
    angles = [i * np.pi / 4 for i in range(8)]  # 0, 45, 90, ..., 315 degrees
    test_points = [(cx + radius * np.cos(angle), cy + radius * np.sin(angle))
                   for angle in angles]
    return all(isInParallelogram((px, py), qx, qa, qy) for px, py in test_points)


def are_discs_overlapping(center1, radius1, center2, radius2):
    dist = np.sqrt((center1[0] - center2[0]) ** 2 + (center1[1] - center2[1]) ** 2)
    return dist < radius1 + radius2


def P(integer_points, num_discs):
    return len(integer_points) / num_discs if num_discs > 0 else 0


def Q(discs, qx, qa, qy):
    """
    Fractional coverage of the parallelogram by the circles.

    Parameters:
    - discs: iterable of (center, radius) where center is (cx, cy)
    - x, a, y: parallelogram parameters (area = x * y)

    Deduplication: copies that differ by integer multiples of (x + a, 0) or (0, y)
    are considered the same disc and counted only once.
    """
    width = qx + qa
    height = qy
    cell_area = abs(qx * qy)

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


def handle_boundary_disc(center, radius, qx, qa, qy):
    # executed after determining that the disc is not fully in the cell
    cx, cy = center
    translated_copies = []
    # Check if disc extends beyond right or left boundaries
    if cx + radius > (qa * cy / qy) + qx:
        translated_copies.append((cx - (qx + qa), cy))
    if cx - radius < (qa * cy / qy):
        translated_copies.append((cx + (qx + qa), cy))
    # Check if disc extends beyond top or bottom boundaries
    if cy + radius > qy:
        translated_copies.append((cx, cy - qy))
    if cy - radius < 0:
        translated_copies.append((cx, cy + qy))
    return translated_copies


def all_points_covered(integer_points, discs):
    """Check if all integer points are covered by at least one disc"""
    arr = []
    for point in integer_points:
        if min_distance_to_any_disc(point, discs) > 0:
            arr.append(point)
    return arr

def can_place_disc(center, radius, existing_discs):
    """Check if disc at center with radius can be placed without overlap"""
    for existing_center, existing_radius in existing_discs:
        dist = np.sqrt((center[0] - existing_center[0]) ** 2 +
                       (center[1] - existing_center[1]) ** 2)
        if dist < radius + existing_radius:
            return False
    return True

import numpy as np
from itertools import combinations
from scipy.optimize import minimize, differential_evolution

def check_all_discs_non_overlapping(discs):
    """Check if any pair of discs overlaps"""
    for (c1, r1), (c2, r2) in combinations(discs, 2):
        dist = np.sqrt((c1[0] - c2[0])**2 + (c1[1] - c2[1])**2)
        if dist < r1 + r2 - 1e-9:  # small tolerance
            return False
    return True

def get_all_discs_with_copies(discs, qx, qa, qy):
    """Get all discs including periodic boundary copies"""
    all_discs = list(discs)
    for center, radius in discs:
        copies = handle_boundary_disc(center, radius, qx, qa, qy)
        for copy_center in copies:
            all_discs.append((copy_center, radius))
    return all_discs

def point_covered_by_disc(point, center, radius):
    """Check if point is within or on disc boundary"""
    px, py = point
    cx, cy = center
    dist = np.sqrt((px - cx)**2 + (py - cy)**2)
    return dist <= radius + 1e-9

def all_integer_points_covered(integer_points, discs, qx, qa, qy):
    """Check if all integer points are covered by discs (including periodic copies)"""
    all_discs = get_all_discs_with_copies(discs, qx, qa, qy)
    for point in integer_points:
        covered = False
        for center, radius in all_discs:
            if point_covered_by_disc(point, center, radius):
                covered = True
                break
        if not covered:
            return False
    return True

def greedy_disc_placement(qx, qa, qy, radii, integer_points):
    """
    Greedy algorithm to place discs covering all integer points.

    Strategy:
    1. Place vertex discs first (largest radius)
    2. Find uncovered points
    3. Place discs to cover maximum uncovered points per disc
    4. Repeat until all points covered
    """
    r1, r2, r3 = sorted(radii, reverse=True)  # r1 >= r2 >= r3
    discs = []

    # Place vertex discs (using largest radius)
    vertices = [(0, 0), (qx, 0), (qa, qy), (qx + qa, qy)]
    for v in vertices:
        if isInParallelogram(v, qx, qa, qy) or v in [(0, 0), (qx, 0), (qa, qy), (qx + qa, qy)]:
            discs.append((v, r1))

    def get_uncovered_points():
        all_discs = get_all_discs_with_copies(discs, qx, qa, qy)
        uncovered = []
        for p in integer_points:
            if not any(point_covered_by_disc(p, c, r) for c, r in all_discs):
                uncovered.append(p)
        return uncovered

    def can_place_at(center, radius):
        """Check if disc can be placed without overlapping existing discs"""
        all_existing = get_all_discs_with_copies(discs, qx, qa, qy)
        new_disc = [(center, radius)]
        new_with_copies = get_all_discs_with_copies(new_disc, qx, qa, qy)

        for nc, nr in new_with_copies:
            for ec, er in all_existing:
                dist = np.sqrt((nc[0] - ec[0])**2 + (nc[1] - ec[1])**2)
                if dist < nr + er - 1e-9:
                    return False
        return True

    def count_points_covered(center, radius, uncovered):
        """Count how many uncovered points would be covered by this disc"""
        count = 0
        disc = [(center, radius)]
        disc_with_copies = get_all_discs_with_copies(disc, qx, qa, qy)
        for p in uncovered:
            if any(point_covered_by_disc(p, c, r) for c, r in disc_with_copies):
                count += 1
        return count

    # Iteratively place discs to cover remaining points
    max_iterations = 1000
    iteration = 0

    while iteration < max_iterations:
        uncovered = get_uncovered_points()
        if not uncovered:
            break

        best_placement = None
        best_score = 0

        # Try placing discs centered at or near uncovered points
        for radius in [r2, r3, r1]:  # Try smaller radii first for efficiency
            for p in uncovered:
                # Try exact point
                candidates = [p]
                # Try nearby positions for better coverage
                for dx in np.linspace(-radius/2, radius/2, 5):
                    for dy in np.linspace(-radius/2, radius/2, 5):
                        candidates.append((p[0] + dx, p[1] + dy))

                for candidate in candidates:
                    if not isInParallelogram(candidate, qx, qa, qy):
                        continue
                    if can_place_at(candidate, radius):
                        score = count_points_covered(candidate, radius, uncovered)
                        if score > best_score:
                            best_score = score
                            best_placement = (candidate, radius)

        if best_placement is None:
            # No valid placement found, try smaller radius at uncovered point
            for p in uncovered:
                for radius in [r3, r2, r1]:
                    if can_place_at(p, radius):
                        best_placement = (p, radius)
                        break
                if best_placement:
                    break

        if best_placement is None:
            break  # Cannot cover all points with current configuration

        discs.append(best_placement)
        iteration += 1

    return discs

def evaluate_configuration(params, integer_points_count):
    """
    Evaluate a configuration and return negative P (for minimization).

    params: [qx, qa, qy, r1, r2, r3]
    """
    qx, qa, qy, r1, r2, r3 = params
    qx, qa, qy = int(round(qx)), int(round(qa)), int(round(qy))

    if qx < 1 or qa < 0 or qy < 1:
        return float('inf')

    integer_points = genIntPointsInCell(qx, qa, qy)
    if len(integer_points) == 0:
        return float('inf')

    radii = [r1, r2, r3]
    discs = greedy_disc_placement(qx, qa, qy, radii, integer_points)

    if not all_integer_points_covered(integer_points, discs, qx, qa, qy):
        return float('inf')

    all_discs = get_all_discs_with_copies(discs, qx, qa, qy)
    if not check_all_discs_non_overlapping(discs):
        return float('inf')

    num_discs = len(discs)
    if num_discs == 0:
        return float('inf')

    p_value = len(integer_points) / num_discs
    return -p_value  # Negative for minimization

def optimize_packing(qx_range=(1, 20), qa_range=(0, 10), qy_range=(1, 20),
                     r_range=(0.5, 10), max_iter=100):
    """
    Optimize the disc packing to maximize P.

    Returns: dict with optimal parameters and results
    """
    bounds = [
        qx_range,  # qx
        qa_range,  # qa
        qy_range,  # qy
        r_range,   # r1
        r_range,   # r2
        r_range,   # r3
    ]

    best_result = None
    best_p = 0

    # Use differential evolution for global optimization
    result = differential_evolution(
        lambda p: evaluate_configuration(p, None),
        bounds,
        maxiter=max_iter,
        seed=42,
        polish=False,
        workers=1
    )

    if result.fun != float('inf'):
        qx, qa, qy = int(round(result.x[0])), int(round(result.x[1])), int(round(result.x[2]))
        r1, r2, r3 = result.x[3], result.x[4], result.x[5]

        integer_points = genIntPointsInCell(qx, qa, qy)
        discs = greedy_disc_placement(qx, qa, qy, [r1, r2, r3], integer_points)

        best_result = {
            'qx': qx,
            'qa': qa,
            'qy': qy,
            'r1': r1,
            'r2': r2,
            'r3': r3,
            'discs': discs,
            'num_discs': len(discs),
            'num_integer_points': len(integer_points),
            'P': len(integer_points) / len(discs) if discs else 0,
            'Q': Q(discs, qx, qa, qy)
        }

    return best_result


if __name__ == "__main__":
    output_file = 'optimization_results.csv'

    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write header
        writer.writerow(['x', 'a', 'y', 'r1', 'optimal_r2', 'optimal_r3', 'num_discs', 'discs', 'P_value', 'Q_value'])

        for x in range(1, 11):
            for a in range(1, 11):
                for y in range(1, 11):
                    for r1 in range(1, 11):
                        if r1 < x / 2 and r1 < y / 2:"""
                            try:
                                optimal_r2, optimal_r3, optimal_discs = optimize_three_radii_for_P(x, a, y, r1)

                                num_discs = len(optimal_discs)
                                discs = optimal_discs
                                integer_points = genIntPointsInCell(x, a, y)
                                p_value = P(integer_points, num_discs)
                                q_value = Q(optimal_discs, x, a, y)

                                # Write row to CSV
                                writer.writerow(
                                    [x, a, y, r1, optimal_r2, optimal_r3, num_discs, discs, p_value, q_value])

                            except Exception as e:
                                print(f"Error for x={x}, a={a}, y={y}, r1={r1}: {e}")"""

    print(f"\nResults saved to {output_file}")

"""
TODO:
- create a greedy algorithm with the constraint that no discs are overlapping and that all 
integer points are covered and to set the two other radii given the one centred at the vertices,
that the two other radii will be smaller than the one at the vertices, using the min_distance_to_any_disc method
calculate the largest possible radius for one of the disc in the insides of the cell and fill the rest of the cell with discs of that radius and with the disc with the smallest radius as well

- create an algorithm that can check if is allowed to place a disc at at a certain point with a certain radius without overlapping with other discs
- 




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
