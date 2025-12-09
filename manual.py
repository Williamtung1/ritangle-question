import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

"""
Parallelogram Unit Cell:
                 (qa,qy)             (qx+qa,qy)
                 ______________________
                /                     /
               /                     /
              -----------------------
            (0,0)               (qx,0)
"""
##### Parameters #####
inputArr = [((0, 0), 3), ((6, 0), 3), ((9, 5), 3), ((3, 5), 3)]
qx = 6
qa = 3
qy = 5
overlap_check = True
all_points_covered_check = True
#######################


def display_discs(discs, qx, qa, qy):
    fig, ax = plt.subplots(figsize=(10, 8))

    # Draw parallelogram boundary
    vertices = [
        (0, 0),
        (qx, 0),
        (qx + qa, qy),
        (qa, qy),
        (0, 0)  # close the shape
    ]
    xs, ys = zip(*vertices)
    ax.plot(xs, ys, 'k-', linewidth=2, label='Cell boundary')

    # Draw discs
    colors = plt.cm.Set3(np.linspace(0, 1, len(discs)))
    for idx, (center, radius) in enumerate(discs):
        circle = patches.Circle(center, radius, fill=True, alpha=0.3,
                                edgecolor='black', linewidth=1.5,
                                facecolor=colors[idx])
        ax.add_patch(circle)
        # Mark center
        ax.plot(center[0], center[1], 'ko', markersize=4)

    from main import genIntPointsInCell
    int_points = genIntPointsInCell(qx, qa, qy)
    px, py = zip(*int_points)
    ax.plot(px, py, 'r.', markersize=8, label='Integer points')

    # Configure axes
    ax.set_aspect('equal')
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)

    ax.grid(True, alpha=0.3)
    ax.set_title(f'Parallelogram Cell (qx={qx}, qa={qa}, qy={qy})', fontsize=14)
    ax.legend(loc='upper right')

    # Set plot limits with some padding
    margin = max(qx, qa, qy) * 0.1
    ax.set_xlim(-margin, qx + qa + margin)
    ax.set_ylim(-margin, qy + margin)

    plt.tight_layout()
    plt.show()


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


# from main import are_discs_overlapping, P, Q, handle_boundary_disc, all_points_covered, genIntPointsInCell, is_disc_fully_in_cell
# from display import display_discs



all_int_points = genIntPointsInCell(qx, qa, qy)
#arrWithPartialBoundaryDisc = []

if overlap_check:
    for i in range(len(inputArr)):
        for j in range(i + 1, len(inputArr)):
            center1, radius1 = inputArr[i]
            center2, radius2 = inputArr[j]
            if are_discs_overlapping(center1, radius1, center2, radius2):
                print(f"Overlap detected between disc {i} and disc {j}")
                break

if all_points_covered_check:
    arr = all_points_covered(all_int_points, inputArr)
    if not arr:
        print("All the points in the unit cell is covered")
    else:
        print(f"List of uncovered points: {arr}")
"""
for disc in inputArr:
    if not is_disc_fully_in_cell(disc, qx, qa, qy):
        arrWithPartialBoundaryDisc.append(handle_boundary_disc(disc[0], disc[1], qx, qa, qy))
"""
P = P(all_int_points, len(inputArr))
Q = Q(inputArr, qx, qa, qy)

print(f"Parameters: qx: {qx},qy: {qy},qa: {qa}, Number of discs: {len(inputArr)}, Discs: {inputArr}")
print(f"P: {P}, Q: {Q}")
display_discs(inputArr, qx, qa, qy)
