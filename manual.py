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

qx = 10
qa = 2
qy = 10
r1 = 5
r2 = 3
r3 = 2
inputArr = [((0,0),r1),((qx,0),r1),((qa,qy),r1),((qx+qa,qy),r1)]
overlap_check = False
all_points_covered_check = True
#######################

from display import display_discs
from main import areTwoDiscsOverlapping, P, Q, areAllPointsCovered, genIntPointsInCell

all_int_points = genIntPointsInCell(qx, qa, qy)
# arrWithPartialBoundaryDisc = []

if overlap_check:
    for i in range(len(inputArr)):
        for j in range(i + 1, len(inputArr)):
            center1, radius1 = inputArr[i]
            center2, radius2 = inputArr[j]
            if areTwoDiscsOverlapping(center1, radius1, center2, radius2):
                print(f"Overlap detected between disc {i} and disc {j}")
                break

if all_points_covered_check:
    arr = areAllPointsCovered(all_int_points, inputArr)
    if not arr:
        print("All the points in the unit cell is covered")
    else:
        print(f"List of uncovered points: {arr}")
"""
for disc in inputArr:
    if not is_disc_fully_in_cell(disc, qx, qa, qy):
        arrWithPartialBoundaryDisc.append(handle_boundary_disc(disc[0], disc[1], qx, qa, qy))
"""
P = P(all_int_points, len(inputArr), qx, qa, qy)
Q = Q(inputArr, qx, qa, qy)

print(f"Parameters: qx: {qx},qy: {qy},qa: {qa}, Number of discs: {len(inputArr)}, Discs: {inputArr}")
print(f"P: {P}, Q: {Q}")
display_discs(inputArr, qx, qa, qy)
