- Methods
  - an iterator through the different lengths of the radii to give the optimum arrangements and radii of circles
  - a optimizer of configuration of the given radii of the circles to create the optimum unit cell
  - a comparator between two configurations of disks as to which is better
  - a calculator of mean number of points per disc
  - a calculator of the fraction of the area of the plane that is covered by the discs
  - a renderer of a given configuration and a given max size to render it at
  
- Objective
- MAX number of points/disc
- tiebreaker min lowest area coverage
  
----------------------------------

Consider the infinite x-y plane. Your task is to cover all the integer grid points on that plane (i.e. all the points (i,j) where i and j are integers with a repeating pattern of discs. A point on the edge of a disc is regarded as being covered. The discs may touch but not overlap. You may use discs of up to three different sizes. The objective is to use as few discs as possible, i.e. to maximise the mean number of points per disc, P . As a tie-break between solutions with the same P, the solution with the lowest area coverage Q(the fraction of the area of the plane that is covered by the discs) will be preferred. You need to construct a unit cell for your solution. The unit cell must cover the plane using translations only; repeating identical copies of it must generate a set of whole discs that cover the grid points in the whole plane. If rotations of an original shape are required, attach the rotated shape to the original to form a unit cell that needs translations only.

  

Calculate P for the unit cell by adding up the fractions **of** discs included in the unit cell to get a total number of discs. Count points on an edge as 1/2, points on a vertex that is shared by three unit cells as 1/3, points on a 4-way vertex as 1/4, etc.

-------------------------------

**Paraphrasing the problem:**
- cover all the integer grid points on an infinite x-y plane
- a point on the edge of a disc is considered covered
- discs may touch but may not overlap                                     v
- use 3 different sizes of discs in the whole plane 
- Objective: maximise the mean number of points/disc                      v
- Tiebreaker: Minimises the fraction of the plane that is covered by the discs         v
- the unit cell must cover the plane using translation only

**Continuous Variables to be optimised:**
- length of the three radii
- value of a,x,y, ie the cell size and dimensions
- the arrangement of the discs inside the unit cell/in between the vertices of the unit cell(which is the one i am struggling the most with)

**Conditions:**
Covering: Every integer point (x,y) must be covered by at least one disc, center (xc,yc).
∣p−c∣≤ R for every integer point p

Packing (Non-Overlap): The distance between any two disc centers ci and cj must be greater than or equal to the sum of their radii.
∣ci−cj∣ ≥ Ri + Rj (considering periodic boundaries)

The shapes can only be parallelograms

the vertices of the unit cell would be (0,0), (x,0), (a,y) and (x+a, y), where all the vertices would be the centres of the circles

the radius is the same for the four circles centred the vertices of the unit cell

**The objective function F :**
Maximising
F = P - k * Q
where p is the mean number of points/disc
q is the fraction of the plane covered
k is tiny constant

**Thoughts:**
- i think i should implement the calculation of p first and then iterate through a list of dummy data with all Variables constant except from one and get the curve to see approx how they look
- the q func as well
- idea is to iterate through the arr of circles(with the vertices of the unit cell(perhaps as attributes of UnitCell)), check for each if they are fully in the parallelograms or not, if not, (and if they are not on the edges either), find some universal way to integrate them and find the proportion of the circle being in the 
cell itself
- with the given vertices of the unit cell and the equations of the lines and circles, it should be fairly straightforward to determine how many points are in the enclosed area
- would the circles on the inside be smaller than the ones centred on the vertices of the unit cell?




**TODO:**
- to mathematically define the rest of the 3 constraints
- learn about Optimization via Simulated Annealing and the Meta-Heuristic Optimization approach