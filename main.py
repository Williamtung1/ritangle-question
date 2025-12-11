import math
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
            if qa < qx:
                flag = True
            else:
                flag = ((qy / qa) * px >= py >= (qy / qa) * (px - qx))
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

def genIntPointsInCell(qx, qa, qy):
    points = []
    for i in range(math.floor(qx + qa + 1)):
        for j in range(math.floor(qy + 1)):
            if isInParallelogram((i, j), qx, qa, qy):
                points.append((i, j))
    return points

def isDiscFullyInCell(disc, qx, qa, qy):
    (cx, cy), radius = disc
    angles = [i * np.pi / 4 for i in range(8)]  # 0, 45, 90, ..., 315 degrees
    test_points = [(cx + radius * np.cos(angle), cy + radius * np.sin(angle))
                   for angle in angles]
    return all(isInParallelogram((px, py), qx, qa, qy) for px, py in test_points)

def areTwoDiscsOverlapping(center1, radius1, center2, radius2, tolerance=1e-6):
    dist = np.sqrt((center1[0] - center2[0]) ** 2 + (center1[1] - center2[1]) ** 2)
    return dist < (radius1 + radius2 - tolerance)

def checkIfAllDiscsNotOverlapping(discs):
    for i in range(len(discs)):
        for j in range(i + 1, len(discs)):
            center1, radius1 = discs[i]
            center2, radius2 = discs[j]
            if areTwoDiscsOverlapping(center1, radius1, center2, radius2):
                return False
    return True

def P(integer_points, num_discs, qx, qa, qy):
    counter = 0
    for px, py in integer_points:
        if py == (qy / qa) * px and py >= (qy / qa) * (px - qx):
            counter += 0.25
        elif py == 0 or py == qy or py == (qy / qa) * px or py >= (qy / qa) * (px - qx):
            counter += 0.5
        else:
            counter += 1

    return counter / num_discs

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

def dupeBoundaryDisc(disc, qx, qa, qy):
    # executed after determining that the disc is not fully in the cell
    (cx, cy), radius = disc
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

def dupeAllBoundaryDiscs(discs, qx, qa, qy):
    """Get all discs including periodic boundary copies"""
    all_discs = list(discs)
    for disc in discs:
        if not isDiscFullyInCell(disc, qx, qa, qy):
            copy = dupeBoundaryDisc(disc, qx, qa, qy)
            all_discs.append((copy, disc[1]))
    return all_discs

def minDistanceToAnyDisc(point, discs):
    distances = []
    for center, radius in discs:
        d = distToDisc(point, (center, radius))
        distances.append(d)

    # return the smallest distance (0 if point is inside/on any disc)
    return min(distances)

def areAllPointsCovered(integer_points, discs):
    """Check if all integer points are covered by at least one disc"""
    arr = []
    for point in integer_points:
        if minDistanceToAnyDisc(point, discs) > 0:
            arr.append(point)
    return arr

def canPlaceDisc(disc, existing_discs, tolerance=1e-6):
    center, radius = disc
    """Check if disc at center with radius can be placed without overlap"""
    for existing_center, existing_radius in existing_discs:
        dist = np.sqrt((center[0] - existing_center[0]) ** 2 +
                       (center[1] - existing_center[1]) ** 2)
        if dist < (radius + existing_radius - tolerance):
            return False
    return True

def isPointCoveredByDisc(point, disc):
    """Check if point is within or on disc boundary"""
    center, radius = disc
    px, py = point
    cx, cy = center
    dist = np.sqrt((px - cx) ** 2 + (py - cy) ** 2)
    return dist <= radius + 1e-9

def areAllIntPointsCovered(integer_points, discs, qx, qa, qy):
    """Check if all integer points are covered by discs (including periodic copies)"""
    all_discs = dupeAllBoundaryDiscs(discs, qx, qa, qy)
    for point in integer_points:
        covered = False
        for center, radius in all_discs:
            if isPointCoveredByDisc(point, (center, radius)):
                covered = True
                break
        if not covered:
            return False
    return True

def countDiscsInCell(discs, qx, qa, qy):
    """
    Count the effective number of discs in the unit cell.
    Vertex discs count as 0.25 each (shared by 4 cells)
    Edge discs count as 0.5 each (shared by 2 cells)
    Interior discs count as 1.0 each
    """
    vertices = [(0, 0), (qx, 0), (qa, qy), (qx + qa, qy)]
    tolerance = 1e-6
    total_count = 0.0

    for center, radius in discs:
        cx, cy = center

        # Check if disc is at a vertex
        is_vertex = any(abs(cx - vx) < tolerance and abs(cy - vy) < tolerance
                        for vx, vy in vertices)
        if is_vertex:
            total_count += 0.25
            continue

        # Check if disc is on an edge (not at vertex)
        on_edge = False
        if abs(cy) < tolerance or abs(cy - qy) < tolerance:  # top or bottom edge
            on_edge = True
        elif abs(cx) < tolerance and abs(cy - (qy / qa) * cx) < tolerance:  # left edge
            on_edge = True
        elif abs(cx - qx) < tolerance and abs(cy - (qy / qa) * (cx - qx)) < tolerance:  # right edge
            on_edge = True

        if on_edge:
            total_count += 0.5
        else:
            total_count += 1.0

    return total_count

def placeInteriorDiscs_v5(qx, qa, qy, r1, r2, r3, max_iterations=100):
    """
    Greedy algorithm v5: Recursive search with backtracking.
    This attempts to solve the placement problem using a depth-first search,
    allowing it to backtrack from decisions that lead to dead ends.
    """
    # Initial state
    vertex_discs = [
        ((0, 0), r1), ((qx, 0), r1), ((qa, qy), r1), ((qx + qa, qy), r1)
    ]
    integer_points = genIntPointsInCell(qx, qa, qy)
    available_radii = sorted([r2, r3], reverse=True)

    # Memoization cache to avoid re-computing for the same state
    memo = {}

    def solve(current_discs):
        # Create a hashable key for the current state
        state_key = tuple(sorted(current_discs))
        if state_key in memo:
            return memo[state_key]

        uncovered_points = [p for p in integer_points if minDistanceToAnyDisc(p, current_discs) > 1e-9]

        if not uncovered_points:
            return True, current_discs  # Success

        # If we have too many discs, it's likely a dead end path
        if len(current_discs) > len(integer_points):  # Heuristic limit
            return False, []

        # Focus on covering the first uncovered point
        point_to_cover = uncovered_points[0]
        p_np = np.array(point_to_cover)

        candidate_discs = set()

        # --- Generate all possible candidates to cover this one point ---
        # Strategy 1: Center on point
        for r_new in available_radii:
            candidate = (point_to_cover, r_new)
            if canPlaceDisc(candidate, current_discs):
                candidate_discs.add(candidate)

        # Strategy 2: Tangent placements
        for existing_disc in current_discs:
            c_exist, r_exist = existing_disc
            c_exist_np = np.array(c_exist)
            for r_new in available_radii:
                d = np.linalg.norm(p_np - c_exist_np)
                if d > r_new + (r_exist + r_new) or d < abs(r_new - (r_exist + r_new)):
                    continue

                a = (d ** 2 - (r_exist + r_new) ** 2 + r_new ** 2) / (2 * d)
                h_sq = r_new ** 2 - a ** 2
                if h_sq < 0: continue
                h = np.sqrt(h_sq)

                p2 = p_np + a * (c_exist_np - p_np) / d
                perp_vec = np.array([-(c_exist_np[1] - p_np[1]), c_exist_np[0] - p_np[0]]) / d

                c1 = tuple(p2 + h * perp_vec)
                c2 = tuple(p2 - h * perp_vec)

                for c_new in [c1, c2]:
                    candidate = (c_new, r_new)
                    if canPlaceDisc(candidate, current_discs):
                        candidate_discs.add(candidate)

        # --- Try candidates recursively ---
        # Sort candidates to try the most promising ones first
        sorted_candidates = sorted(
            candidate_discs,
            key=lambda disc: sum(1 for p in uncovered_points if isPointCoveredByDisc(p, disc)),
            reverse=True
        )

        for candidate in sorted_candidates:
            new_discs = current_discs + [candidate]
            # Recursive call
            success, final_discs = solve(new_discs)
            if success:
                memo[state_key] = (True, final_discs)
                return True, final_discs

        # If no candidate leads to a solution, this path is a failure
        memo[state_key] = (False, [])
        return False, []

    # Initial call to the recursive solver
    success, final_discs = solve(list(vertex_discs))
    return final_discs, success

def evaluateConfiguration(qx, qa, qy, r1, r2, r3, k=0.01, verbose=False):
    """
    Evaluate a configuration and return the objective value F = P - k*Q

    Returns:
        - F: objective value (or large negative penalty if invalid)
        - discs: list of placed discs
        - P_val: points per disc
        - Q_val: area coverage fraction
        - success: whether all constraints are satisfied
    """
    PENALTY = -1e6

    # Basic validity checks
    if qx <= 0 or qa <= 0 or qy <= 0 or r1 <= 0 or r2 <= 0 or r3 <= 0:
        return PENALTY, [], 0, 0, False

    if qx < 1 or qy < 1:  # Cell too small
        return PENALTY, [], 0, 0, False

    # Check that r1, r2, r3 are different (use tolerance)
    radii = sorted([r1, r2, r3])
    if radii[1] - radii[0] < 0.01 or radii[2] - radii[1] < 0.01:
        if verbose:
            print("Radii not sufficiently different")
        return PENALTY, [], 0, 0, False

    # Check that vertex discs don't overlap
    vertices = [(0, 0), (qx, 0), (qa, qy), (qx + qa, qy)]
    for i in range(len(vertices)):
        for j in range(i + 1, len(vertices)):
            v1, v2 = vertices[i], vertices[j]
            dist = np.sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2)
            if dist < 2 * r1 - 1e-6:
                if verbose:
                    print(f"Vertex discs overlap: r1={r1:.3f}, dist={dist:.3f} between {v1} and {v2}")
                return PENALTY, [], 0, 0, False

    try:
        # Place discs using the new v5 algorithm
        discs, success = placeInteriorDiscs_v5(qx, qa, qy, r1, r2, r3)

        if not success:
            if verbose:
                print("Failed to cover all points")
            return PENALTY, discs, 0, 0, False

        # Check for overlaps
        if not checkIfAllDiscsNotOverlapping(discs):
            if verbose:
                print("Discs overlap")
            return PENALTY, discs, 0, 0, False

        # Calculate objective
        integer_points = genIntPointsInCell(qx, qa, qy)
        num_discs = countDiscsInCell(discs, qx, qa, qy)

        if num_discs == 0:
            return PENALTY, discs, 0, 0, False

        P_val = len(integer_points) / num_discs
        Q_val = Q(discs, qx, qa, qy)
        F = P_val - k * Q_val

        if verbose:
            print(f"Valid config: P={P_val:.4f}, Q={Q_val:.4f}, F={F:.4f}, num_discs={num_discs:.2f}")

        return F, discs, P_val, Q_val, True

    except Exception as e:
        if verbose:
            print(f"Error: {e}")
        return PENALTY, [], 0, 0, False

def optimizationObjective(params, k=0.01):
    """
    Objective function for scipy.optimize
    Takes params = [qx, qa, qy, r1, r2, r3]
    Returns negative F (for minimization)
    """
    qx, qa, qy, r1, r2, r3 = params
    F, discs, P_val, Q_val, success = evaluateConfiguration(qx, qa, qy, r1, r2, r3, k)
    return -F  # Negate because we minimize

def runOptimization(k=0.01, method='differential_evolution', max_iter=50):
    """
    Run the optimization to find best disc packing configuration.

    Parameters:
        k: small constant for objective F = P - k*Q
        method: 'differential_evolution', 'basin_hopping', or 'multistart'
        max_iter: maximum iterations

    Returns:
        best_params, best_F, best_discs, best_P, best_Q
    """
    from scipy.optimize import differential_evolution, minimize
    import time

    print(f"Starting optimization with method: {method}")
    print(f"Objective: F = P - {k}*Q")
    print("-" * 60)

    # Parameter bounds: [qx, qa, qy, r1, r2, r3]
    bounds = [
        (3, 12),  # qx
        (0.5, 8),  # qa
        (3, 12),  # qy
        (1.5, 6),  # r1 (vertex discs - larger)
        (0.8, 4),  # r2 (medium discs)
        (0.5, 3),  # r3 (smallest discs)
    ]

    start_time = time.time()

    if method == 'differential_evolution':
        result = differential_evolution(
            lambda x: optimizationObjective(x, k),
            bounds,
            maxiter=max_iter,
            popsize=15,
            mutation=(0.5, 1.5),
            recombination=0.7,
            seed=42,
            workers=1,
            disp=True,
            polish=True
        )
        best_params = result.x
        best_F = -result.fun

    elif method == 'multistart':
        best_result = None
        best_F = -np.inf
        best_params = None

        n_starts = max_iter
        for i in range(n_starts):
            # Random initial guess
            x0 = [np.random.uniform(b[0], b[1]) for b in bounds]

            try:
                result = minimize(
                    lambda x: optimizationObjective(x, k),
                    x0,
                    method='L-BFGS-B',
                    bounds=bounds,
                    options={'maxiter': 100}
                )

                if -result.fun > best_F:
                    best_F = -result.fun
                    best_params = result.x
                    best_result = result
                    print(f"Start {i + 1}/{n_starts}: F = {best_F:.4f}")
            except:
                pass

        if best_params is None:
            print("All starts failed!")
            return None, None, None, None, None

    else:
        raise ValueError(f"Unknown method: {method}")

    elapsed = time.time() - start_time

    # Evaluate best configuration
    qx, qa, qy, r1, r2, r3 = best_params
    F, discs, P_val, Q_val, success = evaluateConfiguration(qx, qa, qy, r1, r2, r3, k, verbose=True)

    print("-" * 60)
    print(f"Optimization completed in {elapsed:.2f}s")
    print(f"Best configuration found:")
    print(f"  Cell: qx={qx:.4f}, qa={qa:.4f}, qy={qy:.4f}")
    print(f"  Radii: r1={r1:.4f}, r2={r2:.4f}, r3={r3:.4f}")
    print(f"  Objective F = {F:.4f}")
    print(f"  Points per disc (P) = {P_val:.4f}")
    print(f"  Area coverage (Q) = {Q_val:.4f}")
    print(f"  Number of discs = {len(discs)}")
    print(f"  Valid solution: {success}")

    return best_params, F, discs, P_val, Q_val

def saveResults(params, F, discs, P_val, Q_val, filename='optimization_results.csv'):
    """Save optimization results to CSV file"""
    import csv
    import os

    qx, qa, qy, r1, r2, r3 = params

    # Read existing data
    file_exists = os.path.exists(filename)

    with open(filename, 'a', newline='') as f:
        writer = csv.writer(f)
        if not file_exists or os.path.getsize(filename) == 0:
            writer.writerow(['qx', 'qa', 'qy', 'r1', 'r2', 'r3', 'num_discs', 'P_value', 'Q_value', 'F_value', 'discs'])

        discs_str = str(discs)
        writer.writerow([qx, qa, qy, r1, r2, r3, len(discs), P_val, Q_val, F, discs_str])

    print(f"Results saved to {filename}")

if __name__ == "__main__":
    # Run optimization
    print("=" * 60)
    print("DISC PACKING OPTIMIZATION")
    print("=" * 60)

    # You can adjust k (smaller k prioritizes P more, larger k considers Q more)
    k_value = 0.01

    # Run optimization with differential evolution (best for global optimization)
    best_params, best_F, best_discs, best_P, best_Q = runOptimization(
        k=k_value,
        method='differential_evolution',
        max_iter=50  # Increase for better results (but slower)
    )

    if best_params is not None:
        # Save results
        saveResults(best_params, best_F, best_discs, best_P, best_Q)

        # Visualize the result
        print("\nGenerating visualization...")
        from display import display_discs

        qx, qa, qy, r1, r2, r3 = best_params
        display_discs(best_discs, qx, qa, qy)

    print("\nOptimization complete!")
