import numpy as np
import warnings
from collections import Counter

from matplotlib.colors import ColorConverter
from matplotlib.pyplot import gca

from collections import Counter
from matplotlib.patches import Circle, PathPatch
from matplotlib.path import Path
from matplotlib.colors import ColorConverter
from matplotlib.pyplot import gca


from scipy.optimize import brentq

tol = 1e-10


def compute_venn2_areas(diagram_areas, normalize_to=1.0):
    '''
    The list of venn areas is given as 3 values, corresponding to venn diagram areas in the following order:
     (Ab, aB, AB)  (i.e. last element corresponds to the size of intersection A&B&C).
    The return value is a list of areas (A, B, AB), such that the total area is normalized
    to normalize_to. If total area was 0, returns (1e-06, 1e-06, 0.0)

    Assumes all input values are nonnegative (to be more precise, all areas are passed through and abs() function)
    >>> compute_venn2_areas((1, 1, 0))
    (0.5, 0.5, 0.0)
    >>> compute_venn2_areas((0, 0, 0))
    (1e-06, 1e-06, 0.0)
    >>> compute_venn2_areas((1, 1, 1), normalize_to=3)
    (2.0, 2.0, 1.0)
    >>> compute_venn2_areas((1, 2, 3), normalize_to=6)
    (4.0, 5.0, 3.0)
    '''
    # Normalize input values to sum to 1
    areas = np.array(np.abs(diagram_areas), float)
    total_area = np.sum(areas)
    if np.abs(total_area) < tol:
        warnings.warn("Both circles have zero area")
        return (1e-06, 1e-06, 0.0)
    else:
        areas = areas / total_area * normalize_to
        return (areas[0] + areas[2], areas[1] + areas[2], areas[2])


def solve_venn2_circles(venn_areas):
    '''
    Given the list of "venn areas" (as output from compute_venn2_areas, i.e. [A, B, AB]),
    finds the positions and radii of the two circles.
    The return value is a tuple (coords, radii), where coords is a 2x2 array of coordinates and
    radii is a 2x1 array of circle radii.

    Assumes the input values to be nonnegative and not all zero.
    In particular, the first two values must be positive.

    >>> c, r = solve_venn2_circles((1, 1, 0))
    >>> np.round(r, 3)
    array([ 0.564,  0.564])
    >>> c, r = solve_venn2_circles(compute_venn2_areas((1, 2, 3)))
    >>> np.round(r, 3)
    array([ 0.461,  0.515])
    '''
    (A_a, A_b, A_ab) = list(map(float, venn_areas))
    r_a, r_b = np.sqrt(A_a / np.pi), np.sqrt(A_b / np.pi)
    radii = np.array([r_a, r_b])
    if A_ab > tol:
        # Nonzero intersection
        coords = np.zeros((2, 2))
        coords[1][0] = find_distance_by_area(radii[0], radii[1], A_ab)
    else:
        # Zero intersection
        coords = np.zeros((2, 2))
        coords[1][0] = radii[0] + radii[1] + max(np.mean(radii) * 1.1, 0.2)   # The max here is needed for the case r_a = r_b = 0
    coords = normalize_by_center_of_mass(coords, radii)
    return (coords, radii)


def compute_venn2_regions(centers, radii):
    '''
    Returns a triple of VennRegion objects, describing the three regions of the diagram, corresponding to sets
    (Ab, aB, AB)

    >>> centers, radii = solve_venn2_circles((1, 1, 0.5))
    >>> regions = compute_venn2_regions(centers, radii)
    '''
    A = VennCircleRegion(centers[0], radii[0])
    B = VennCircleRegion(centers[1], radii[1])
    Ab, AB = A.subtract_and_intersect_circle(B.center, B.radius)
    aB, _ = B.subtract_and_intersect_circle(A.center, A.radius)
    return (Ab, aB, AB)


def compute_venn2_colors(set_colors):
    '''
    Given two base colors, computes combinations of colors corresponding to all regions of the venn diagram.
    returns a list of 3 elements, providing colors for regions (10, 01, 11).

    >>> compute_venn2_colors(('r', 'g'))
    (array([ 1.,  0.,  0.]), array([ 0. ,  0.5,  0. ]), array([ 0.7 ,  0.35,  0.  ]))
    '''
    ccv = ColorConverter()
    base_colors = [np.array(ccv.to_rgb(c)) for c in set_colors]
    return (base_colors[0], base_colors[1], mix_colors(base_colors[0], base_colors[1]))


def compute_venn2_subsets(a, b):
    '''
    Given two set or Counter objects, computes the sizes of (a & ~b, b & ~a, a & b).
    Returns the result as a tuple.

    >>> compute_venn2_subsets(set([1,2,3,4]), set([2,3,4,5,6]))
    (1, 2, 3)
    >>> compute_venn2_subsets(Counter([1,2,3,4]), Counter([2,3,4,5,6]))
    (1, 2, 3)
    >>> compute_venn2_subsets(Counter([]), Counter([]))
    (0, 0, 0)
    >>> compute_venn2_subsets(set([]), set([]))
    (0, 0, 0)
    >>> compute_venn2_subsets(set([1]), set([]))
    (1, 0, 0)
    >>> compute_venn2_subsets(set([1]), set([1]))
    (0, 0, 1)
    >>> compute_venn2_subsets(Counter([1]), Counter([1]))
    (0, 0, 1)
    >>> compute_venn2_subsets(set([1,2]), set([1]))
    (1, 0, 1)
    >>> compute_venn2_subsets(Counter([1,1,2,2,2]), Counter([1,2,3,3]))
    (3, 2, 2)
    >>> compute_venn2_subsets(Counter([1,1,2]), Counter([1,2,2]))
    (1, 1, 2)
    >>> compute_venn2_subsets(Counter([1,1]), set([]))
    Traceback (most recent call last):
    ...
    ValueError: Both arguments must be of the same type
    '''
    if not (type(a) == type(b)):
        raise ValueError("Both arguments must be of the same type")
    set_size = len if type(a) != Counter else lambda x: sum(x.values())   # We cannot use len to compute the cardinality of a Counter
    return (set_size(a - b), set_size(b - a), set_size(a & b))


def venn2_circles(subsets, normalize_to=1.0, alpha=1.0, color='black', linestyle='solid', linewidth=2.0, ax=None, **kwargs):
    '''
    Plots only the two circles for the corresponding Venn diagram.
    Useful for debugging or enhancing the basic venn diagram.
    parameters ``subsets``, ``normalize_to`` and ``ax`` are the same as in venn2()
    ``kwargs`` are passed as-is to matplotlib.patches.Circle.
    returns a list of three Circle patches.

    >>> c = venn2_circles((1, 2, 3))
    >>> c = venn2_circles({'10': 1, '01': 2, '11': 3}) # Same effect
    >>> c = venn2_circles([set([1,2,3,4]), set([2,3,4,5,6])]) # Also same effect
    '''
    if isinstance(subsets, dict):
        subsets = [subsets.get(t, 0) for t in ['10', '01', '11']]
    elif len(subsets) == 2:
        subsets = compute_venn2_subsets(*subsets)
    areas = compute_venn2_areas(subsets, normalize_to)
    centers, radii = solve_venn2_circles(areas)

    if ax is None:
        ax = gca()
    prepare_venn_axes(ax, centers, radii)
    result = []
    for (c, r) in zip(centers, radii):
        circle = Circle(c, r, alpha=alpha, edgecolor=color, facecolor='none', linestyle=linestyle, linewidth=linewidth, **kwargs)
        ax.add_patch(circle)
        result.append(circle)
    return result


class VennBase(object):
    def __init__(self):
        pass


class Venn2(object):
    def __init__(self, data):
        self.venn2(data)

    def compute_venn2_areas(self, diagram_areas, normalize_to=1.0):
        '''
        The list of venn areas is given as 3 values, corresponding to venn diagram areas in the following order:
         (Ab, aB, AB)  (i.e. last element corresponds to the size of intersection A&B&C).
        The return value is a list of areas (A, B, AB), such that the total area is normalized
        to normalize_to. If total area was 0, returns (1e-06, 1e-06, 0.0)

        Assumes all input values are nonnegative (to be more precise, all areas are passed through and abs() function)
        >>> compute_venn2_areas((1, 1, 0))
        (0.5, 0.5, 0.0)
        >>> compute_venn2_areas((0, 0, 0))
        (1e-06, 1e-06, 0.0)
        >>> compute_venn2_areas((1, 1, 1), normalize_to=3)
        (2.0, 2.0, 1.0)
        >>> compute_venn2_areas((1, 2, 3), normalize_to=6)
        (4.0, 5.0, 3.0)
        '''
        # Normalize input values to sum to 1
        areas = np.array(np.abs(diagram_areas), float)
        total_area = np.sum(areas)
        if np.abs(total_area) < tol:
            warnings.warn("Both circles have zero area")
            return (1e-06, 1e-06, 0.0)
        else:
            areas = areas / total_area * normalize_to
            return (areas[0] + areas[2], areas[1] + areas[2], areas[2])

    def venn2(self, subsets, set_labels=('A', 'B'), set_colors=('r', 'g'), 
              alpha=0.4, normalize_to=1.0, ax=None, subset_label_formatter=None):
        '''Plots a 2-set area-weighted Venn diagram.
        The subsets parameter can be one of the following:
         - A list (or a tuple) containing two set objects.
         - A dict, providing sizes of three diagram regions.
           The regions are identified via two-letter binary codes ('10', '01', and '11'), hence a valid set could look like:
           {'10': 10, '01': 20, '11': 40}. Unmentioned codes are considered to map to 0.
         - A list (or a tuple) with three numbers, denoting the sizes of the regions in the following order:
           (10, 01, 11)
    
        ``set_labels`` parameter is a list of two strings - set labels. Set it to None to disable set labels.
        The ``set_colors`` parameter should be a list of two elements, specifying the "base colors" of the two circles.
        The color of circle intersection will be computed based on those.
    
        The ``normalize_to`` parameter specifies the total (on-axes) area of the circles to be drawn. Sometimes tuning it (together
        with the overall fiture size) may be useful to fit the text labels better.
        The return value is a ``VennDiagram`` object, that keeps references to the ``Text`` and ``Patch`` objects used on the plot
        and lets you know the centers and radii of the circles, if you need it.
    
        The ``ax`` parameter specifies the axes on which the plot will be drawn (None means current axes).
    
        The ``subset_label_formatter`` parameter is a function that can be passed to format the labels
        that describe the size of each subset.
    
        >>> from matplotlib_venn import *
        >>> v = venn2(subsets={'10': 1, '01': 1, '11': 1}, set_labels = ('A', 'B'))
        >>> c = venn2_circles(subsets=(1, 1, 1), linestyle='dashed')
        >>> v.get_patch_by_id('10').set_alpha(1.0)
        >>> v.get_patch_by_id('10').set_color('white')
        >>> v.get_label_by_id('10').set_text('Unknown')
        >>> v.get_label_by_id('A').set_text('Set A')

        You can provide sets themselves rather than subset sizes:
        >>> v = venn2(subsets=[set([1,2]), set([2,3,4,5])], set_labels = ('A', 'B'))
        >>> c = venn2_circles(subsets=[set([1,2]), set([2,3,4,5])], linestyle='dashed')
        >>> print("%0.2f" % (v.get_circle_radius(1)/v.get_circle_radius(0)))
        1.41
        '''
        if isinstance(subsets, dict):
            subsets = [subsets.get(t, 0) for t in ['10', '01', '11']]
        elif len(subsets) == 2:
            subsets = compute_venn2_subsets(*subsets)

        if subset_label_formatter is None:
            subset_label_formatter = str

        areas = self.compute_venn2_areas(subsets, normalize_to)
        centers, radii = solve_venn2_circles(areas)
        regions = compute_venn2_regions(centers, radii)
        colors = compute_venn2_colors(set_colors)

        if ax is None:
            ax = gca()
        prepare_venn_axes(ax, centers, radii)

        # Create and add patches and subset labels
        patches = [r.make_patch() for r in regions]
        for (p, c) in zip(patches, colors):
            if p is not None:
                p.set_facecolor(c)
                p.set_edgecolor('none')
                p.set_alpha(alpha)
                ax.add_patch(p)
        label_positions = [r.label_position() for r in regions]
        subset_labels = [ax.text(lbl[0], lbl[1], subset_label_formatter(s), va='center', ha='center') if lbl is not None else None for (lbl, s) in zip(label_positions, subsets)]

        # Position set labels
        if set_labels is not None:
            padding = np.mean([r * 0.1 for r in radii])
            label_positions = [centers[0] + np.array([0.0, - radii[0] - padding]),
                           centers[1] + np.array([0.0, - radii[1] - padding])]
            labels = [ax.text(pos[0], pos[1], txt, size='large', ha='right', va='top') for (pos, txt) in zip(label_positions, set_labels)]
            labels[1].set_ha('left')
        else:
            labels = None
        return VennDiagram(patches, subset_labels, labels, centers, radii)


class Venn3(object):
    def __init__(self, data):
        venn3(data)


class Venn4(object):
    def __init__(self):
        pass

class Venn5(object):
    def __init__(self):
        pass

class Venn6(object):
    def __init__(self):
        pass


class Venn(object):
    def __init__(self):
        pass



def venn2_unweighted(subsets, set_labels=('A', 'B'), set_colors=('r', 'g'), alpha=0.4, normalize_to=1.0, subset_areas=(1, 1, 1), ax=None, subset_label_formatter=None):
    '''
    The version of venn2 without area-weighting.
    It is implemented as a wrapper around venn2. Namely, venn2 is invoked as usual, but with all subset areas
    set to 1. The subset labels are then replaced in the resulting diagram with the provided subset sizes.

    The parameters are all the same as that of venn2.
    In addition there is a subset_areas parameter, which specifies the actual subset areas.
    (it is (1, 1, 1) by default. You are free to change it, within reason).
    '''
    v = venn2(subset_areas, set_labels, set_colors, alpha, normalize_to, ax)
    # Now rename the labels
    if subset_label_formatter is None:
        subset_label_formatter = str
    subset_ids = ['10', '01', '11']
    if isinstance(subsets, dict):
        subsets = [subsets.get(t, 0) for t in subset_ids]
    elif len(subsets) == 2:
        subsets = compute_venn2_subsets(*subsets)
    for n, id in enumerate(subset_ids):
        lbl = v.get_label_by_id(id)
        if lbl is not None:
            lbl.set_text(subset_label_formatter(subsets[n]))
    return v


def venn3_unweighted(subsets, set_labels=('A', 'B', 'C'), set_colors=('r', 'g', 'b'), alpha=0.4, normalize_to=1.0, subset_areas=(1, 1, 1, 1, 1, 1, 1), ax=None, subset_label_formatter=None):
    '''
    The version of venn3 without area-weighting.
    It is implemented as a wrapper around venn3. Namely, venn3 is invoked as usual, but with all subset areas
    set to 1. The subset labels are then replaced in the resulting diagram with the provided subset sizes.

    The parameters are all the same as that of venn2.
    In addition there is a subset_areas parameter, which specifies the actual subset areas.
    (it is (1, 1, 1, 1, 1, 1, 1) by default. You are free to change it, within reason).
    '''
    v = venn3(subset_areas, set_labels, set_colors, alpha, normalize_to, ax)
    # Now rename the labels
    if subset_label_formatter is None:
        subset_label_formatter = str
    subset_ids = ['100', '010', '110', '001', '101', '011', '111']
    if isinstance(subsets, dict):
        subsets = [subsets.get(t, 0) for t in subset_ids]
    elif len(subsets) == 3:
        subsets = compute_venn3_subsets(*subsets)
    for n, id in enumerate(subset_ids):
        lbl = v.get_label_by_id(id)
        if lbl is not None:
            lbl.set_text(subset_label_formatter(subsets[n]))
    return v



def point_in_circle(pt, center, radius):
    '''
    Returns true if a given point is located inside (or on the border) of a circle.

    >>> point_in_circle((0, 0), (0, 0), 1)
    True
    >>> point_in_circle((1, 0), (0, 0), 1)
    True
    >>> point_in_circle((1, 1), (0, 0), 1)
    False
    '''
    d = np.linalg.norm(np.asarray(pt) - np.asarray(center))
    return d <= radius


def box_product(v1, v2):
    '''Returns a determinant |v1 v2|. The value is equal to the signed area of a parallelogram built on v1 and v2.
    The value is positive is v2 is to the left of v1.

    >>> box_product((0.0, 1.0), (0.0, 1.0))
    0.0
    >>> box_product((1.0, 0.0), (0.0, 1.0))
    1.0
    >>> box_product((0.0, 1.0), (1.0, 0.0))
    -1.0
    '''
    return v1[0]*v2[1] - v1[1]*v2[0]


def circle_intersection_area(r, R, d):
    '''
    Formula from: http://mathworld.wolfram.com/Circle-CircleIntersection.html
    Does not make sense for negative r, R or d

    >>> circle_intersection_area(0.0, 0.0, 0.0)
    0.0
    >>> circle_intersection_area(1.0, 1.0, 0.0)
    3.1415...
    >>> circle_intersection_area(1.0, 1.0, 1.0)
    1.2283...
    '''
    if np.abs(d) < tol:
        minR = np.min([r, R])
        return np.pi * minR**2
    if np.abs(r - 0) < tol or np.abs(R - 0) < tol:
        return 0.0
    d2, r2, R2 = float(d**2), float(r**2), float(R**2)
    arg = (d2 + r2 - R2) / 2 / d / r
    arg = np.max([np.min([arg, 1.0]), -1.0])  # Even with valid arguments, the above computation may result in things like -1.001
    A = r2 * np.arccos(arg)
    arg = (d2 + R2 - r2) / 2 / d / R
    arg = np.max([np.min([arg, 1.0]), -1.0])
    B = R2 * np.arccos(arg)
    arg = (-d + r + R) * (d + r - R) * (d - r + R) * (d + r + R)
    arg = np.max([arg, 0])
    C = -0.5 * np.sqrt(arg)
    return A + B + C


def circle_line_intersection(center, r, a, b):
    '''
    Computes two intersection points between the circle centered at <center> and radius <r> and a line given by two points a and b.
    If no intersection exists, or if a==b, None is returned. If one intersection exists, it is repeated in the answer.

    >>> circle_line_intersection(np.array([0.0, 0.0]), 1, np.array([-1.0, 0.0]), np.array([1.0, 0.0]))
    array([[ 1.,  0.],
           [-1.,  0.]])
    >>> abs(np.round(circle_line_intersection(np.array([1.0, 1.0]), np.sqrt(2), np.array([-1.0, 1.0]), np.array([1.0, -1.0])), 6))
    array([[ 0.,  0.],
           [ 0.,  0.]])
    '''
    s = b - a
    # Quadratic eqn coefs
    A = np.linalg.norm(s)**2
    if abs(A) < tol:
        return None
    B = 2 * np.dot(a - center, s)
    C = np.linalg.norm(a - center)**2 - r**2
    disc = B**2 - 4 * A * C
    if disc < 0.0:
        return None
    t1 = (-B + np.sqrt(disc)) / 2.0 / A
    t2 = (-B - np.sqrt(disc)) / 2.0 / A
    return np.array([a + t1 * s, a + t2 * s])


def find_distance_by_area(r, R, a, numeric_correction=0.0001):
    '''
    Solves circle_intersection_area(r, R, d) == a for d numerically (analytical solution seems to be too ugly to pursue).
    Assumes that a < pi * min(r, R)**2, will fail otherwise.

    The numeric correction parameter is used whenever the computed distance is exactly (R - r) (i.e. one circle must be inside another).
    In this case the result returned is (R-r+correction). This helps later when we position the circles and need to ensure they intersect.

    >>> find_distance_by_area(1, 1, 0, 0.0)
    2.0
    >>> round(find_distance_by_area(1, 1, 3.1415, 0.0), 4)
    0.0
    >>> d = find_distance_by_area(2, 3, 4, 0.0)
    >>> d
    3.37...
    >>> round(circle_intersection_area(2, 3, d), 10)
    4.0
    >>> find_distance_by_area(1, 2, np.pi)
    1.0001
    '''
    if r > R:
        r, R = R, r
    if np.abs(a) < tol:
        return float(r + R)
    if np.abs(min([r, R])**2 * np.pi - a) < tol:
        return np.abs(R - r + numeric_correction)
    return brentq(lambda x: circle_intersection_area(r, R, x) - a, R - r, R + r)


def circle_circle_intersection(C_a, r_a, C_b, r_b):
    '''
    Finds the coordinates of the intersection points of two circles A and B.
    Circle center coordinates C_a and C_b, should be given as tuples (or 1x2 arrays).
    Returns a 2x2 array result with result[0] being the first intersection point (to the right of the vector C_a -> C_b)
    and result[1] being the second intersection point.

    If there is a single intersection point, it is repeated in output.
    If there are no intersection points or an infinite number of those, None is returned.

    >>> circle_circle_intersection([0, 0], 1, [1, 0], 1) # Two intersection points
    array([[ 0.5      , -0.866...],
           [ 0.5      ,  0.866...]])
    >>> circle_circle_intersection([0, 0], 1, [2, 0], 1) # Single intersection point (circles touch from outside)
    array([[ 1.,  0.],
           [ 1.,  0.]])
    >>> circle_circle_intersection([0, 0], 1, [0.5, 0], 0.5) # Single intersection point (circles touch from inside)
    array([[ 1.,  0.],
           [ 1.,  0.]])
    >>> circle_circle_intersection([0, 0], 1, [0, 0], 1) is None # Infinite number of intersections (circles coincide)
    True
    >>> circle_circle_intersection([0, 0], 1, [0, 0.1], 0.8) is None # No intersections (one circle inside another)
    True
    >>> circle_circle_intersection([0, 0], 1, [2.1, 0], 1) is None # No intersections (one circle outside another)
    True
    '''
    C_a, C_b = np.asarray(C_a, float), np.asarray(C_b, float)
    v_ab = C_b - C_a
    d_ab = np.linalg.norm(v_ab)
    if np.abs(d_ab) < tol:  # No intersection points or infinitely many of them (circle centers coincide)
        return None
    cos_gamma = (d_ab**2 + r_a**2 - r_b**2) / 2.0 / d_ab / r_a

    if abs(cos_gamma) > 1.0 + tol/10: # Allow for a tiny numeric tolerance here too (always better to be return something instead of None, if possible)
        return None         # No intersection point (circles do not touch)
    if (cos_gamma > 1.0):
        cos_gamma = 1.0
    if (cos_gamma < -1.0):
        cos_gamma = -1.0

    sin_gamma = np.sqrt(1 - cos_gamma**2)
    u = v_ab / d_ab
    v = np.array([-u[1], u[0]])
    pt1 = C_a + r_a * cos_gamma * u - r_a * sin_gamma * v
    pt2 = C_a + r_a * cos_gamma * u + r_a * sin_gamma * v
    return np.array([pt1, pt2])


def vector_angle_in_degrees(v):
    '''
    Given a vector, returns its elevation angle in degrees (-180..180).

    >>> vector_angle_in_degrees([1, 0])
    0.0
    >>> vector_angle_in_degrees([1, 1])
    45.0
    >>> vector_angle_in_degrees([0, 1])
    90.0
    >>> vector_angle_in_degrees([-1, 1])
    135.0
    >>> vector_angle_in_degrees([-1, 0])
    180.0
    >>> vector_angle_in_degrees([-1, -1])
    -135.0
    >>> vector_angle_in_degrees([0, -1])
    -90.0
    >>> vector_angle_in_degrees([1, -1])
    -45.0
    '''
    return np.arctan2(v[1], v[0]) * 180 / np.pi


def normalize_by_center_of_mass(coords, radii):
    '''
    Given coordinates of circle centers and radii, as two arrays,
    returns new coordinates array, computed such that the center of mass of the
    three circles is (0, 0).

    >>> normalize_by_center_of_mass(np.array([[0.0, 0.0], [2.0, 0.0], [1.0, 3.0]]), np.array([1.0, 1.0, 1.0]))
    array([[-1., -1.],
           [ 1., -1.],
           [ 0.,  2.]])
    >>> normalize_by_center_of_mass(np.array([[0.0, 0.0], [2.0, 0.0], [1.0, 2.0]]), np.array([1.0, 1.0, np.sqrt(2.0)]))
    array([[-1., -1.],
           [ 1., -1.],
           [ 0.,  1.]])
    '''
    # Now find the center of mass.
    radii = radii**2
    sum_r = np.sum(radii)
    if sum_r < tol:
        return coords
    else:
        return coords - np.dot(radii, coords) / np.sum(radii)

'''

The current logic of drawing the venn diagram is the following:
 - Position the circles.
 - Compute the regions of the diagram based on circles
 - Compute the position of the label within each region.
 - Create matplotlib PathPatch or Circle objects for each of the regions.

This module contains functionality necessary for the second, third and fourth steps of this process.

Note that the regions of an up to 3-circle Venn diagram may be of the following kinds:
 - No region
 - A circle
 - A 2, 3 or 4-arc "poly-arc-gon".  (I.e. a polygon with up to 4 vertices, that are connected by circle arcs)
 - A set of two 3-arc-gons.

We create each of the regions by starting with a circle, and then either intersecting or subtracting the second and the third circles.
The classes below implement the region representation, the intersection/subtraction procedures and the conversion to matplotlib patches.
In addition, each region type has a "label positioning" procedure assigned.
'''


class VennRegionException(Exception):
    pass


class VennRegion(object):
    '''
    This is a superclass of a Venn diagram region, defining the interface that has to be supported by the different region types.
    '''
    def subtract_and_intersect_circle(self, center, radius):
        '''
        Given a circular region, compute two new regions:
        one obtained by subtracting the circle from this region, and another obtained by intersecting the circle with the region.

        In all implementations it is assumed that the circle to be subtracted is not completely within
        the current region without touching its borders, i.e. it will not form a "hole" when subtracted.

        Arguments:
           center (tuple):  A two-element tuple-like, representing the coordinates of the center of the circle.
           radius (float):  A nonnegative number, the radius of the circle.

        Returns:
           a list with two elements - the result of subtracting the circle, and the result of intersecting with the circle.
        '''
        raise NotImplementedError("Method not implemented")


    def label_position(self):
        '''Compute the position of a label for this region and return it as a 1x2 numpy array (x, y).
        May return None if label is not applicable.'''
        raise NotImplementedError("Method not implemented")

    def size(self):
        '''Return a number, representing the size of the region. It is not important that the number would be a precise
        measurement, as long as sizes of various regions can be compared to choose the largest one.'''
        raise NotImplementedError("Method not implemented")

    def make_patch(self):
        '''Create a matplotlib patch object, corresponding to this region. May return None if no patch has to be created.'''
        raise NotImplementedError("Method not implemented")

    def verify(self):
        '''Self-verification routine for purposes of testing. Raises a VennRegionException if some inconsistencies of internal representation
        are discovered.'''
        raise NotImplementedError("Method not implemented")


class VennEmptyRegion(VennRegion):
    '''
    An empty region. To save some memory, returns [self, self] on the subtract_and_intersect_circle operation.
    It is possible to create an empty region with a non-None label position, by providing it in the constructor.

    >>> v = VennEmptyRegion()
    >>> [a, b] = v.subtract_and_intersect_circle((1,2), 3)
    >>> assert a == v and b == v
    >>> assert v.label_position() is None
    >>> assert v.size() == 0
    >>> assert v.make_patch() is None
    >>> assert v.is_empty()
    >>> v = VennEmptyRegion((0, 0))
    >>> v.label_position()
    array([ 0.,  0.])
    '''
    def __init__(self, label_pos = None):
        self.label_pos = None if label_pos is None else np.asarray(label_pos, float)
    def subtract_and_intersect_circle(self, center, radius):
        return [self, self]
    def size(self):
        return 0
    def label_position(self):
        return self.label_pos
    def make_patch(self):
        return None
    def is_empty(self):  # We use this in tests as an equivalent of isinstance(VennEmptyRegion)
        return True
    def verify(self):
        pass


class VennCircleRegion(VennRegion):
    '''
    A circle-shaped region.

    >>> vcr = VennCircleRegion((0, 0), 1)
    >>> vcr.size()
    3.1415...
    >>> vcr.label_position()
    array([ 0.,  0.])
    >>> vcr.make_patch()
    <matplotlib.patches.Circle object at ...>
    >>> sr, ir = vcr.subtract_and_intersect_circle((0.5, 0), 1)
    >>> assert abs(sr.size() + ir.size() - vcr.size()) < tol
    '''

    def __init__(self, center, radius):
        self.center = np.asarray(center, float)
        self.radius = abs(radius)
        if (radius < -tol):
            raise VennRegionException("Circle with a negative radius is invalid")

    def subtract_and_intersect_circle(self, center, radius):
        '''Will throw a VennRegionException if the circle to be subtracted is completely inside and not touching the given region.'''

        # Check whether the target circle intersects us
        center = np.asarray(center, float)
        d = np.linalg.norm(center - self.center)
        if d > (radius + self.radius - tol):
            return [self, VennEmptyRegion()] # The circle does not intersect us
        elif d < tol:
            if radius > self.radius - tol:
                # We are completely covered by that circle or we are the same circle
                return [VennEmptyRegion(), self]
            else:
                # That other circle is inside us and smaller than us - we can't deal with it
                raise VennRegionException("Invalid configuration of circular regions (holes are not supported).")
        else:
            # We *must* intersect the other circle. If it is not the case, then it is inside us completely,
            # and we'll complain.
            intersections = circle_circle_intersection(self.center, self.radius, center, radius)

            if intersections is None:
                raise VennRegionException("Invalid configuration of circular regions (holes are not supported).")
            elif np.all(abs(intersections[0] - intersections[1]) < tol) and self.radius < radius:
                # There is a single intersection point (i.e. we are touching the circle),
                # the circle to be subtracted is not outside of us (this was checked before), and is larger than us.
                # This is a particular corner case that is not dealt with correctly by the general-purpose code below and must
                # be handled separately
                return [VennEmptyRegion(), self]
            else:
                # Otherwise the subtracted region is a 2-arc-gon
                # Before we need to convert the intersection points as angles wrt each circle.
                a_1 = vector_angle_in_degrees(intersections[0] - self.center)
                a_2 = vector_angle_in_degrees(intersections[1] - self.center)
                b_1 = vector_angle_in_degrees(intersections[0] - center)
                b_2 = vector_angle_in_degrees(intersections[1] - center)

                # We must take care of the situation where the intersection points happen to be the same
                if (abs(b_1 - b_2) < tol):
                    b_1 = b_2 - tol/2
                if (abs(a_1 - a_2) < tol):
                    a_2 = a_1 + tol/2

                # The subtraction is a 2-arc-gon [(AB, B-), (BA, A+)]
                s_arc1 = Arc(center, radius, b_1, b_2, False)
                s_arc2 = Arc(self.center, self.radius, a_2, a_1, True)
                subtraction = VennArcgonRegion([s_arc1, s_arc2])

                # .. and the intersection is a 2-arc-gon [(AB, A+), (BA, B+)]
                i_arc1 = Arc(self.center, self.radius, a_1, a_2, True)
                i_arc2 = Arc(center, radius, b_2, b_1, True)
                intersection = VennArcgonRegion([i_arc1, i_arc2])
                return [subtraction, intersection]

    def size(self):
        '''
        Return the area of the circle

        >>> VennCircleRegion((0, 0), 1).size()
        3.1415...
        >>> VennCircleRegion((0, 0), 2).size()
        12.56637...
        '''
        return np.pi * self.radius**2;

    def label_position(self):
        '''
        The label should be positioned in the center of the circle

        >>> VennCircleRegion((0, 0), 1).label_position()
        array([ 0.,  0.])
        >>> VennCircleRegion((-1.2, 3.4), 1).label_position()
        array([-1.2,  3.4])
        '''
        return self.center

    def make_patch(self):
        '''
        Returns the corresponding circular patch.

        >>> patch = VennCircleRegion((1, 2), 3).make_patch()
        >>> patch
        <matplotlib.patches.Circle object at ...>
        >>> patch.center, patch.radius
        (array([ 1.,  2.]), 3.0)
        '''
        return Circle(self.center, self.radius)

    def verify(self):
        pass


class VennArcgonRegion(VennRegion):
    '''
    A poly-arc region.
    Note that we essentially only support 2, 3 and 4 arced regions,
    whereas intersections and subtractions only work for 2-arc regions.
    '''

    def __init__(self, arcs):
        '''
        Create a poly-arc region given a list of Arc objects.
        The arcs list must be of length 2, 3 or 4.
        The arcs must form a closed polygon, i.e. the last point of each arc must be the first point of the next arc.
        The vertices of a 3 or 4-arcgon must be listed in a CCW order. Arcs must not intersect.

        This is not verified in the constructor, but a special verify() method can be used to check
        for validity.
        '''
        self.arcs = arcs

    def verify(self):
        '''
        Verify the correctness of the region arcs. Throws an VennRegionException if verification fails
        (or any other exception if it happens during verification).
        '''
        # Verify size of arcs list
        if (len(self.arcs) < 2):
            raise VennRegionException("At least two arcs needed in a poly-arc region")
        if (len(self.arcs) > 4):
            raise VennRegionException("At most 4 arcs are supported currently for poly-arc regions")

        TRIG_TOL = 100*tol  # We need to use looser tolerance level here because conversion to angles and back is prone to large errors.
        # Verify connectedness of arcs
        for i in range(len(self.arcs)):
            if not np.all(self.arcs[i-1].end_point() - self.arcs[i].start_point() < TRIG_TOL):
                raise VennRegionException("Arcs of an poly-arc-gon must be connected via endpoints")

        # Verify that arcs do not cross-intersect except at endpoints
        for i in range(len(self.arcs)-1):
            for j in range(i+1, len(self.arcs)):
                ips = self.arcs[i].intersect_arc(self.arcs[j])
                for ip in ips:
                    if not (np.all(abs(ip - self.arcs[i].start_point()) < TRIG_TOL) or np.all(abs(ip - self.arcs[i].end_point()) < TRIG_TOL)):
                        raise VennRegionException("Arcs of a poly-arc-gon may only intersect at endpoints")

                if len(ips) != 0 and (i - j) % len(self.arcs) > 1 and (j - i) % len(self.arcs) > 1:
                    # Two non-consecutive arcs intersect. This is in general not good, but
                    # may occasionally happen when all arcs inbetween have length 0.
                    pass # raise VennRegionException("Non-consecutive arcs of a poly-arc-gon may not intersect")

        # Verify that vertices are ordered so that at each point the direction along the polyarc changes towards the left.
        # Note that this test only makes sense for polyarcs obtained using circle intersections & subtractions.
        # A "flower-like" polyarc may have its vertices ordered counter-clockwise yet the direction would turn to the right at each of them.
        for i in range(len(self.arcs)):
            prev_arc = self.arcs[i-1]
            cur_arc = self.arcs[i]
            if box_product(prev_arc.direction_vector(prev_arc.to_angle), cur_arc.direction_vector(cur_arc.from_angle)) < -tol:
                raise VennRegionException("Arcs must be ordered so that the direction at each vertex changes counter-clockwise")

    def subtract_and_intersect_circle(self, center, radius):
        '''
        Circle subtraction / intersection only supported by 2-gon regions, otherwise a VennRegionException is thrown.
        In addition, such an exception will be thrown if the circle to be subtracted is completely within the region and forms a "hole".

        The result may be either a VennArcgonRegion or a VennMultipieceRegion (the latter happens when the circle "splits" a crescent in two).
        '''
        if len(self.arcs) != 2:
            raise VennRegionException("Circle subtraction and intersection with poly-arc regions is currently only supported for 2-arc-gons.")

        # In the following we consider the 2-arc-gon case.
        # Before we do anything, we check for a special case, where the circle of interest is one of the two circles forming the arcs.
        # In this case we can determine the answer quite easily.
        matching_arcs = [a for a in self.arcs if a.lies_on_circle(center, radius)]
        if len(matching_arcs) != 0:
            # If the circle matches a positive arc, the result is [empty, self], otherwise [self, empty]
            return [VennEmptyRegion(), self] if matching_arcs[0].direction else [self, VennEmptyRegion()]

        # Consider the intersection points of the circle with the arcs.
        # If any of the intersection points corresponds exactly to any of the arc's endpoints, we will end up with
        # a lot of messy special cases (as if the usual situation is not messy enough, eh).
        # To avoid that, we cheat by slightly increasing the circle's radius until this is not the case any more.
        center = np.asarray(center)
        illegal_intersections = [a.start_point() for a in self.arcs]
        while True:
            valid = True
            intersections = [a.intersect_circle(center, radius) for a in self.arcs]
            for ints in intersections:
                for pt in ints:
                    for illegal_pt in illegal_intersections:
                        if np.all(abs(pt - illegal_pt) < tol):
                            valid = False
            if valid:
                break
            else:
                radius += tol

        center = np.asarray(center)
        intersections = [a.intersect_circle(center, radius) for a in self.arcs]

        if len(intersections[0]) == 0 and len(intersections[1]) == 0:
            # Case I
            if point_in_circle(self.arcs[0].start_point(), center, radius):
                # Case I.a)
                return [VennEmptyRegion(), self]
            else:
                # Case I.b)
                return [self, VennEmptyRegion()]
        elif len(intersections[0]) == 2 and len(intersections[1]) == 2:
            # Case II. a) or b)
            case_II_a = not point_in_circle(self.arcs[0].start_point(), center, radius)

            a1 = self.arcs[0].subarc_between_points(None, intersections[0][0])
            a2 = Arc(center, radius,
                     vector_angle_in_degrees(intersections[0][0] - center),
                     vector_angle_in_degrees(intersections[1][1] - center),
                     not case_II_a)
            a2.fix_360_to_0()
            a3 = self.arcs[1].subarc_between_points(intersections[1][1], None)
            piece1 = VennArcgonRegion([a1, a2, a3])

            b1 = self.arcs[1].subarc_between_points(None, intersections[1][0])
            b2 = Arc(center, radius,
                     vector_angle_in_degrees(intersections[1][0] - center),
                     vector_angle_in_degrees(intersections[0][1] - center),
                     not case_II_a)
            b2.fix_360_to_0()
            b3 = self.arcs[0].subarc_between_points(intersections[0][1], None)
            piece2 = VennArcgonRegion([b1, b2, b3])

            subtraction = VennMultipieceRegion([piece1, piece2])

            c1 = self.arcs[0].subarc(a1.to_angle, b3.from_angle)
            c2 = b2.reversed()
            c3 = self.arcs[1].subarc(b1.to_angle, a3.from_angle)
            c4 = a2.reversed()
            intersection = VennArcgonRegion([c1, c2, c3, c4])

            return [subtraction, intersection] if case_II_a else [intersection, subtraction]
        else:
            # Case III. Yuck.
            if len(intersections[0]) == 0 or len(intersections[1]) == 0:
                # Case III.a)
                x = 0 if len(intersections[0]) != 0 else 1
                y = 1 - x
                if len(intersections[x]) != 2:
                    warnings.warn("Numeric precision error during polyarc intersection, case IIIa. Expect wrong results.")
                    intersections[x] = [intersections[x][0], intersections[x][0]]  # This way we'll at least produce some result, although it will probably be wrong
                if not point_in_circle(self.arcs[0].start_point(), center, radius):
                    # Case III.a.1)
                    #   result_subtraction = {X from start to i, circle i to j (direction = negative), X j to end, Y}
                    a1 = self.arcs[x].subarc_between_points(None, intersections[x][0])
                    a2 = Arc(center, radius,
                             vector_angle_in_degrees(intersections[x][0] - center),
                             vector_angle_in_degrees(intersections[x][1] - center),
                             False)
                    a3 = self.arcs[x].subarc_between_points(intersections[x][1], None)
                    a4 = self.arcs[y]
                    subtraction = VennArcgonRegion([a1, a2, a3, a4])

                    #   result_intersection = {X i to j, circle j to i (direction = positive)}
                    b1 = self.arcs[x].subarc(a1.to_angle, a3.from_angle)
                    b2 = a2.reversed()
                    intersection = VennArcgonRegion([b1, b2])

                    return [subtraction, intersection]
                else:
                    # Case III.a.2)
                    #   result_subtraction = {X i to j, circle j to i negative}
                    a1 = self.arcs[x].subarc_between_points(intersections[x][0], intersections[x][1])
                    a2 = Arc(center, radius,
                             vector_angle_in_degrees(intersections[x][1] - center),
                             vector_angle_in_degrees(intersections[x][0] - center),
                             False)
                    subtraction = VennArcgonRegion([a1, a2])

                    #   result_intersection = {X 0 to i, circle i to j positive, X j to end, Y}
                    b1 = self.arcs[x].subarc(None, a1.from_angle)
                    b2 = a2.reversed()
                    b3 = self.arcs[x].subarc(a1.to_angle, None)
                    b4 = self.arcs[y]
                    intersection = VennArcgonRegion([b1, b2, b3, b4])

                    return [subtraction, intersection]
            else:
                # Case III.b)
                if len(intersections[0]) == 2 or len(intersections[1]) == 2:
                    warnings.warn("Numeric precision error during polyarc intersection, case IIIb. Expect wrong results.")

                # One of the arcs must start outside the circle, call it x
                x = 0 if not point_in_circle(self.arcs[0].start_point(), center, radius) else 1
                y = 1 - x

                a1 = self.arcs[x].subarc_between_points(None, intersections[x][0])
                a2 = Arc(center, radius,
                         vector_angle_in_degrees(intersections[x][0] - center),
                         vector_angle_in_degrees(intersections[y][0] - center), False)
                a3 = self.arcs[y].subarc_between_points(intersections[y][0], None)
                subtraction = VennArcgonRegion([a1, a2, a3])

                b1 = self.arcs[x].subarc(a1.to_angle, None)
                b2 = self.arcs[y].subarc(None, a3.from_angle)
                b3 = a2.reversed()
                intersection = VennArcgonRegion([b1, b2, b3])
                return [subtraction, intersection]

    def label_position(self):
        # Position the label right inbetween the midpoints of the arcs
        midpoints = [a.mid_point() for a in self.arcs]
        # For two-arc regions take the usual average
        # For more than two arcs, use arc lengths as the weights.
        if len(self.arcs) == 2:
            return np.mean(midpoints, 0)
        else:
            lengths = [a.length_degrees() for a in self.arcs]
            avg = np.sum([mp * l for (mp, l) in zip(midpoints, lengths)], 0)
            return avg / np.sum(lengths)

    def size(self):
        '''Return the area of the patch.

        The area can be computed using the standard polygon area formula + signed segment areas of each arc.
        '''
        polygon_area = 0
        for a in self.arcs:
            polygon_area += box_product(a.start_point(), a.end_point())
        polygon_area /= 2.0
        return polygon_area + sum([a.sign * a.segment_area() for a in self.arcs])

    def make_patch(self):
        '''
        Retuns a matplotlib PathPatch representing the current region.
        '''
        path = [self.arcs[0].start_point()]
        for a in self.arcs:
            if a.direction:
                vertices = Path.arc(a.from_angle, a.to_angle).vertices
            else:
                vertices = Path.arc(a.to_angle, a.from_angle).vertices
                vertices = vertices[np.arange(len(vertices) - 1, -1, -1)]
            vertices = vertices * a.radius + a.center
            path = path + list(vertices[1:])
        codes = [1] + [4] * (len(path) - 1)  # NB: We could also add a CLOSEPOLY code (and a random vertex) to the end
        return PathPatch(Path(path, codes))


class VennMultipieceRegion(VennRegion):
    '''
    A region containing several pieces.
    In principle, any number of pieces is supported,
    although no more than 2 should ever be needed in a 3-circle Venn diagram.
    Although subtraction/intersection are straightforward to implement we do
    not need those for matplotlib-venn, we raise exceptions in those methods.
    '''

    def __init__(self, pieces):
        '''
        Create a multi-piece region from a list of VennRegion objects.
        The list may be empty or contain a single item (although those regions can be converted to a
        VennEmptyRegion or a single region of the necessary type.
        '''
        self.pieces = pieces

    def label_position(self):
        '''
        Find the largest region and position the label in that.
        '''
        reg_sizes = [(r.size(), r) for r in self.pieces]
        reg_sizes.sort()
        return reg_sizes[-1][1].label_position()

    def size(self):
        return sum([p.size() for p in self.pieces])

    def make_patch(self):
        '''Currently only works if all the pieces are Arcgons.
           In this case returns a multiple-piece path. Otherwise throws an exception.'''
        paths = [p.make_patch().get_path() for p in self.pieces]
        vertices = np.concatenate([p.vertices for p in paths])
        codes = np.concatenate([p.codes for p in paths])
        return PathPatch(Path(vertices, codes))

    def verify(self):
        for p in self.pieces:
            p.verify()



class Arc(object):
    '''
    A representation of a directed circle arc.
    Essentially it is a namedtuple(center, radius, from_angle, to_angle, direction) with a bunch of helper methods
    for measuring arc lengths and intersections.

    The from_angle and to_angle of an arc must be represented in degrees.
    The direction is a boolean, with True corresponding to counterclockwise (positive) direction, and False - clockwise (negative).
    For convenience, the class defines a "sign" property, which is +1 if direction = True and -1 otherwise.
    '''

    def __init__(self, center, radius, from_angle, to_angle, direction):
        '''Raises a ValueError if radius is negative.

        >>> a = Arc((0, 0), -1, 0, 0, True)
        Traceback (most recent call last):
        ...
        ValueError: Arc's radius may not be negative
        >>> a = Arc((0, 0), 0, 0, 0, True)
        >>> a = Arc((0, 0), 1, 0, 0, True)
        '''
        self.center = np.asarray(center)
        self.radius = float(radius)
        if radius < 0.0:
            raise ValueError("Arc's radius may not be negative")
        self.from_angle = float(from_angle)
        self.to_angle = float(to_angle)
        self.direction = direction
        self.sign = 1 if direction else -1

    def length_degrees(self):
        '''Computes the length of the arc in degrees.
        The length computation corresponds to what you would expect if you would draw the arc using matplotlib taking direction into account.

        >>> Arc((0,0), 1, 0, 0, True).length_degrees()
        0.0
        >>> Arc((0,0), 2, 0, 0, False).length_degrees()
        0.0

        >>> Arc((0,0), 3, 0, 1, True).length_degrees()
        1.0
        >>> Arc((0,0), 4, 0, 1, False).length_degrees()
        359.0

        >>> Arc((0,0), 5, 0, 360, True).length_degrees()
        360.0
        >>> Arc((0,0), 6, 0, 360, False).length_degrees()
        0.0

        >>> Arc((0,0), 7, 0, 361, True).length_degrees()
        360.0
        >>> Arc((0,0), 8, 0, 361, False).length_degrees()
        359.0

        >>> Arc((0,0), 9, 10, -10, True).length_degrees()
        340.0
        >>> Arc((0,0), 10, 10, -10, False).length_degrees()
        20.0

        >>> Arc((0,0), 1, 10, 5, True).length_degrees()
        355.0
        >>> Arc((0,0), 1, -10, -5, False).length_degrees()
        355.0
        >>> Arc((0,0), 1, 180, -180, True).length_degrees()
        0.0
        >>> Arc((0,0), 1, 180, -180, False).length_degrees()
        360.0
        >>> Arc((0,0), 1, -180, 180, True).length_degrees()
        360.0
        >>> Arc((0,0), 1, -180, 180, False).length_degrees()
        0.0
        >>> Arc((0,0), 1, 175, -175, True).length_degrees()
        10.0
        >>> Arc((0,0), 1, 175, -175, False).length_degrees()
        350.0
        '''
        d_angle = self.sign * (self.to_angle - self.from_angle)
        if (d_angle > 360):
            return 360.0
        elif (d_angle < 0):
            return d_angle % 360.0
        else:
            return abs(d_angle)   # Yes, abs() is needed, otherwise we get the weird "-0.0" output in the doctests

    def length_radians(self):
        '''Returns the length of the arc in radians.

        >>> Arc((0,0), 1, 0, 0, True).length_radians()
        0.0
        >>> Arc((0,0), 2, 0, 360, True).length_radians()
        6.283...
        >>> Arc((0,0), 6, -18, 18, True).length_radians()
        0.6283...
        '''
        return self.length_degrees() * np.pi / 180.0

    def length(self):
        '''Returns the actual length of the arc.

        >>> Arc((0,0), 2, 0, 360, True).length()
        12.566...
        >>> Arc((0,0), 2, 90, 360, False).length()
        3.1415...
        >>> Arc((0,0), 0, 90, 360, True).length()
        0.0
        '''
        return self.radius * self.length_radians()

    def sector_area(self):
        '''Returns the area of the corresponding arc sector.

        >>> Arc((0,0), 2, 0, 360, True).sector_area()
        12.566...
        >>> Arc((0,0), 2, 0, 36, True).sector_area()
        1.2566...
        >>> Arc((0,0), 2, 0, 36, False).sector_area()
        11.3097...
        '''
        return self.radius**2 / 2 * self.length_radians()

    def segment_area(self):
        '''Returns the area of the corresponding arc segment.

        >>> Arc((0,0), 2, 0, 360, True).segment_area()
        12.566...
        >>> Arc((0,0), 2, 0, 180, True).segment_area()
        6.283...
        >>> Arc((0,0), 2, 0, 90, True).segment_area()
        1.14159...
        >>> Arc((0,0), 2, 0, 90, False).segment_area()
        11.42477796...
        >>> Arc((0,0), 2, 0, 0, False).segment_area()
        0.0
        >>> Arc((0, 9), 1, 89.99, 90, False).segment_area()
        3.1415...
        '''
        theta = self.length_radians()
        return self.radius**2 / 2 * (theta - np.sin(theta))

    def angle_as_point(self, angle):
        '''
        Converts a given angle in degrees to the point coordinates on the arc's circle.
        Inverse of point_to_angle.

        >>> Arc((1, 1), 1, 0, 0, True).angle_as_point(0)
        array([ 2.,  1.])
        >>> Arc((1, 1), 1, 0, 0, True).angle_as_point(90)
        array([ 1.,  2.])
        >>> Arc((1, 1), 1, 0, 0, True).angle_as_point(-270)
        array([ 1.,  2.])
        '''
        angle_rad = angle * np.pi / 180.0
        return self.center + self.radius * np.array([np.cos(angle_rad), np.sin(angle_rad)])

    def start_point(self):
        '''
        Returns a 2x1 numpy array with the coordinates of the arc's start point.

        >>> Arc((0, 0), 1, 0, 0, True).start_point()
        array([ 1.,  0.])
        >>> Arc((0, 0), 1, 45, 0, True).start_point()
        array([ 0.707...,  0.707...])
        '''
        return self.angle_as_point(self.from_angle)

    def end_point(self):
        '''
        Returns a 2x1 numpy array with the coordinates of the arc's end point.

        >>> np.all(Arc((0, 0), 1, 0, 90, True).end_point() - np.array([0, 1]) < tol)
        True
        '''
        return self.angle_as_point(self.to_angle)

    def mid_point(self):
        '''
        Returns the midpoint of the arc as a 1x2 numpy array.
        '''
        midpoint_angle = self.from_angle + self.sign*self.length_degrees() / 2
        return self.angle_as_point(midpoint_angle)

    def approximately_equal(self, arc, tolerance=tol):
        '''
        Returns true if the parameters of this arc are within <tolerance> of the parameters of the other arc, and the direction is the same.
        Note that no angle simplification is performed (i.e. some arcs that might be equal in principle are not declared as such
        by this method)

        >>> Arc((0, 0), 10, 20, 30, True).approximately_equal(Arc((tol/2, tol/2), 10+tol/2, 20-tol/2, 30-tol/2, True))
        True
        >>> Arc((0, 0), 10, 20, 30, True).approximately_equal(Arc((0, 0), 10, 20, 30, False))
        False
        >>> Arc((0, 0), 10, 20, 30, True).approximately_equal(Arc((0, 0+tol), 10, 20, 30, True))
        False
        '''
        return self.direction == arc.direction \
                and np.all(abs(self.center - arc.center) < tolerance) and abs(self.radius - arc.radius) < tolerance \
                and abs(self.from_angle - arc.from_angle) < tolerance and abs(self.to_angle - arc.to_angle) < tolerance

    def point_as_angle(self, pt):
        '''
        Given a point located on the arc's circle, return the corresponding angle in degrees.
        No check is done that the point lies on the circle
        (this is essentially a convenience wrapper around _math.vector_angle_in_degrees)

        >>> a = Arc((0, 0), 1, 0, 0, True)
        >>> a.point_as_angle((1, 0))
        0.0
        >>> a.point_as_angle((1, 1))
        45.0
        >>> a.point_as_angle((0, 1))
        90.0
        >>> a.point_as_angle((-1, 1))
        135.0
        >>> a.point_as_angle((-1, 0))
        180.0
        >>> a.point_as_angle((-1, -1))
        -135.0
        >>> a.point_as_angle((0, -1))
        -90.0
        >>> a.point_as_angle((1, -1))
        -45.0
        '''
        return vector_angle_in_degrees(np.asarray(pt) - self.center)

    def contains_angle_degrees(self, angle):
        '''
        Returns true, if a point with the corresponding angle (given in degrees) is within the arc.
        Does no tolerance checks (i.e. if the arc is of length 0, you must provide angle == from_angle == to_angle to get a positive answer here)

        >>> a = Arc((0, 0), 1, 0, 0, True)
        >>> assert a.contains_angle_degrees(0)
        >>> assert a.contains_angle_degrees(360)
        >>> assert not a.contains_angle_degrees(1)

        >>> a = Arc((0, 0), 1, 170, -170, True)
        >>> assert not a.contains_angle_degrees(165)
        >>> assert a.contains_angle_degrees(170)
        >>> assert a.contains_angle_degrees(175)
        >>> assert a.contains_angle_degrees(180)
        >>> assert a.contains_angle_degrees(185)
        >>> assert a.contains_angle_degrees(190)
        >>> assert not a.contains_angle_degrees(195)

        >>> assert not a.contains_angle_degrees(-195)
        >>> assert a.contains_angle_degrees(-190)
        >>> assert a.contains_angle_degrees(-185)
        >>> assert a.contains_angle_degrees(-180)
        >>> assert a.contains_angle_degrees(-175)
        >>> assert a.contains_angle_degrees(-170)
        >>> assert not a.contains_angle_degrees(-165)
        >>> assert a.contains_angle_degrees(-170 - 360)
        >>> assert a.contains_angle_degrees(-190 - 360)
        >>> assert a.contains_angle_degrees(170 + 360)
        >>> assert not a.contains_angle_degrees(0)
        >>> assert not a.contains_angle_degrees(100)
        >>> assert not a.contains_angle_degrees(-100)
        '''
        _d = self.sign * (angle - self.from_angle) % 360.0
        return (_d <= self.length_degrees())

    def intersect_circle(self, center, radius):
        '''
        Given a circle, finds the intersection point(s) of the arc with the circle.
        Returns a list of 2x1 numpy arrays. The list has length 0, 1 or 2, depending on how many intesection points there are.
        If the circle touches the arc, it is reported as two intersection points (which are equal).
        Points are ordered along the arc.
        Intersection with the same circle as the arc's own (which means infinitely many points usually) is reported as no intersection at all.

        >>> a = Arc((0, 0), 1, -60, 60, True)
        >>> a.intersect_circle((1, 0), 1)
        [array([ 0.5..., -0.866...]), array([ 0.5...,  0.866...])]
        >>> a.intersect_circle((0.9, 0), 1)
        []
        >>> a.intersect_circle((1,-0.1), 1)
        [array([ 0.586...,  0.810...])]
        >>> a.intersect_circle((1, 0.1), 1)
        [array([ 0.586..., -0.810...])]
        >>> a.intersect_circle((0, 0), 1)  # Infinitely many intersection points
        []
        >>> a.intersect_circle((2, 0), 1)  # Touching point, hence repeated twice
        [array([ 1.,  0.]), array([ 1.,  0.])]

        >>> a = Arc((0, 0), 1, 60, -60, False) # Same arc, different direction
        >>> a.intersect_circle((1, 0), 1)
        [array([ 0.5...,  0.866...]), array([ 0.5..., -0.866...])]

        >>> a = Arc((0, 0), 1, 120, -120, True)
        >>> a.intersect_circle((-1, 0), 1)
        [array([-0.5...,  0.866...]), array([-0.5..., -0.866...])]
        >>> a.intersect_circle((-0.9, 0), 1)
        []
        >>> a.intersect_circle((-1,-0.1), 1)
        [array([-0.586...,  0.810...])]
        >>> a.intersect_circle((-1, 0.1), 1)
        [array([-0.586..., -0.810...])]
        >>> a.intersect_circle((-2, 0), 1)
        [array([-1.,  0.]), array([-1.,  0.])]
        >>> a = Arc((0, 0), 1, -120, 120, False)
        >>> a.intersect_circle((-1, 0), 1)
        [array([-0.5..., -0.866...]), array([-0.5...,  0.866...])]
        '''
        intersections = circle_circle_intersection(self.center, self.radius, center, radius)
        if intersections is None:
            return []

        # Check whether the points lie on the arc and order them accordingly
        _len = self.length_degrees()
        isections = [[self.sign * (self.point_as_angle(pt) - self.from_angle) % 360.0, pt] for pt in intersections]

        # Try to find as many candidate intersections as possible (i.e. +- tol within arc limits)
        # Unless arc's length is 360, interpret intersections just before the arc's starting point as belonging to the starting point.
        if _len < 360.0 - tol:
            for isec in isections:
                if isec[0] > 360.0 - tol:
                    isec[0] = 0.0

        isections = [(a, pt[0], pt[1]) for (a, pt) in isections if a < _len + tol or a > 360 - tol]
        isections.sort()
        return [np.array([b, c]) for (a, b, c) in isections]

    def intersect_arc(self, arc):
        '''
        Given an arc, finds the intersection point(s) of this arc with that.
        Returns a list of 2x1 numpy arrays. The list has length 0, 1 or 2, depending on how many intesection points there are.
        Points are ordered along the arc.
        Intersection with the arc along the same circle (which means infinitely many points usually) is reported as no intersection at all.

        >>> a = Arc((0, 0), 1, -90, 90, True)
        >>> a.intersect_arc(Arc((1, 0), 1, 90, 270, True))
        [array([ 0.5      , -0.866...]), array([ 0.5      ,  0.866...])]
        >>> a.intersect_arc(Arc((1, 0), 1, 90, 180, True))
        [array([ 0.5      ,  0.866...])]
        >>> a.intersect_arc(Arc((1, 0), 1, 121, 239, True))
        []
        >>> a.intersect_arc(Arc((1, 0), 1, 120-tol, 240+tol, True))     # Without -tol and +tol the results differ on different architectures due to rounding (see Debian #813782).
        [array([ 0.5      , -0.866...]), array([ 0.5      ,  0.866...])]
        '''
        intersections = self.intersect_circle(arc.center, arc.radius)
        isections = [pt for pt in intersections if arc.contains_angle_degrees(arc.point_as_angle(pt))]
        return isections

    def subarc(self, from_angle=None, to_angle=None):
        '''
        Creates a sub-arc from a given angle (or beginning of this arc) to a given angle (or end of this arc).
        Verifies that from_angle and to_angle are within the arc and properly ordered.
        If from_angle is None, start of this arc is used instead.
        If to_angle is None, end of this arc is used instead.
        Angles are given in degrees.

        >>> a = Arc((0, 0), 1, 0, 360, True)
        >>> a.subarc(None, None)
        Arc([0.000, 0.000], 1.000,   0.000, 360.000,   True,   degrees=360.000)
        >>> a.subarc(360, None)
        Arc([0.000, 0.000], 1.000,   360.000, 360.000,   True,   degrees=0.000)
        >>> a.subarc(0, None)
        Arc([0.000, 0.000], 1.000,   0.000, 360.000,   True,   degrees=360.000)
        >>> a.subarc(-10, None)
        Arc([0.000, 0.000], 1.000,   350.000, 360.000,   True,   degrees=10.000)
        >>> a.subarc(None, -10)
        Arc([0.000, 0.000], 1.000,   0.000, 350.000,   True,   degrees=350.000)
        >>> a.subarc(1, 359).subarc(2, 358).subarc()
        Arc([0.000, 0.000], 1.000,   2.000, 358.000,   True,   degrees=356.000)
        '''

        if from_angle is None:
            from_angle = self.from_angle
        if to_angle is None:
            to_angle = self.to_angle
        cur_length = self.length_degrees()
        d_new_from = self.sign * (from_angle - self.from_angle)
        if (d_new_from != 360.0):
            d_new_from = d_new_from % 360.0
        d_new_to = self.sign * (to_angle - self.from_angle)
        if (d_new_to != 360.0):
            d_new_to = d_new_to % 360.0
        # Gracefully handle numeric precision issues for zero-length arcs
        if abs(d_new_from - d_new_to) < tol:
            d_new_from = d_new_to
        if d_new_to < d_new_from:
            raise ValueError("Subarc to-angle must be smaller than from-angle.")
        if d_new_to > cur_length + tol:
            raise ValueError("Subarc to-angle must lie within the current arc.")
        return Arc(self.center, self.radius, self.from_angle + self.sign*d_new_from, self.from_angle + self.sign*d_new_to, self.direction)

    def subarc_between_points(self, p_from=None, p_to=None):
        '''
        Given two points on the arc, extract a sub-arc between those points.
        No check is made to verify the points are actually on the arc.
        It is basically a wrapper around subarc(point_as_angle(p_from), point_as_angle(p_to)).
        Either p_from or p_to may be None to denote first or last arc endpoints.

        >>> a = Arc((0, 0), 1, 0, 90, True)
        >>> a.subarc_between_points((1, 0), (np.cos(np.pi/4), np.sin(np.pi/4)))
        Arc([0.000, 0.000], 1.000,   0.000, 45.000,   True,   degrees=45.000)
        >>> a.subarc_between_points(None, None)
        Arc([0.000, 0.000], 1.000,   0.000, 90.000,   True,   degrees=90.000)
        >>> a.subarc_between_points((np.cos(np.pi/4), np.sin(np.pi/4)))
        Arc([0.000, 0.000], 1.000,   45.000, 90.000,   True,   degrees=45.000)
        '''
        a_from = self.point_as_angle(p_from) if p_from is not None else None
        a_to = self.point_as_angle(p_to) if p_to is not None else None
        return self.subarc(a_from, a_to)

    def reversed(self):
        '''
        Returns a copy of this arc, with the direction flipped.

        >>> Arc((0, 0), 1, 0, 360, True).reversed()
        Arc([0.000, 0.000], 1.000,   360.000, 0.000,   False,   degrees=360.000)
        >>> Arc((0, 0), 1, 175, -175, True).reversed()
        Arc([0.000, 0.000], 1.000,   -175.000, 175.000,   False,   degrees=10.000)
        >>> Arc((0, 0), 1, 0, 370, True).reversed()
        Arc([0.000, 0.000], 1.000,   370.000, 0.000,   False,   degrees=360.000)
        '''
        return Arc(self.center, self.radius, self.to_angle, self.from_angle, not self.direction)

    def direction_vector(self, angle):
        '''
        Returns a unit vector, pointing in the arc's movement direction at a given (absolute) angle (in degrees).
        No check is made whether angle lies within the arc's span (the results for angles outside of the arc's span )
        Returns a 2x1 numpy array.

        >>> a = Arc((0, 0), 1, 0, 90, True)
        >>> assert all(abs(a.direction_vector(0) - np.array([0.0, 1.0])) < tol)
        >>> assert all(abs(a.direction_vector(45) - np.array([ -0.70710678, 0.70710678])) < 1e-6)
        >>> assert all(abs(a.direction_vector(90) - np.array([-1.0, 0.0])) < tol)
        >>> assert all(abs(a.direction_vector(135) - np.array([-0.70710678, -0.70710678])) < 1e-6)
        >>> assert all(abs(a.direction_vector(-180) - np.array([0.0, -1.0])) < tol)
        >>> assert all(abs(a.direction_vector(-90) - np.array([1.0, 0.0])) < tol)
        >>> a = a.reversed()
        >>> assert all(abs(a.direction_vector(0) - np.array([0.0, -1.0])) < tol)
        >>> assert all(abs(a.direction_vector(45) - np.array([ 0.70710678, -0.70710678])) < 1e-6)
        >>> assert all(abs(a.direction_vector(90) - np.array([1.0, 0.0])) < tol)
        >>> assert all(abs(a.direction_vector(135) - np.array([0.70710678, 0.70710678])) < 1e-6)
        >>> assert all(abs(a.direction_vector(-180) - np.array([0.0, 1.0])) < tol)
        >>> assert all(abs(a.direction_vector(-90) - np.array([-1.0, 0.0])) < tol)
        '''
        a = angle + self.sign * 90
        a = a * np.pi / 180.0
        return np.array([np.cos(a), np.sin(a)])

    def fix_360_to_0(self):
        '''
        Sometimes we have to create an arc using from_angle and to_angle computed numerically.
        If from_angle == to_angle, it may sometimes happen that a tiny discrepancy will make from_angle > to_angle, and instead of
        getting a 0-length arc we end up with a 360-degree arc.
        Sometimes we know for sure that a 360-degree arc is not what we want, and in those cases
        the problem is easy to fix. This helper method does that. It checks whether from_angle and to_angle are numerically similar,
        and if so makes them equal.

        >>> a = Arc((0, 0), 1, 0, -tol/2, True)
        >>> a
        Arc([0.000, 0.000], 1.000,   0.000, -0.000,   True,   degrees=360.000)
        >>> a.fix_360_to_0()
        >>> a
        Arc([0.000, 0.000], 1.000,   -0.000, -0.000,   True,   degrees=0.000)
        '''
        if abs(self.from_angle - self.to_angle) < tol:
            self.from_angle = self.to_angle

    def lies_on_circle(self, center, radius):
        '''Tests whether the arc circle's center and radius match the given ones within <tol> tolerance.

        >>> a = Arc((0, 0), 1, 0, 0, False)
        >>> a.lies_on_circle((tol/2, tol/2), 1+tol/2)
        True
        >>> a.lies_on_circle((tol/2, tol/2), 1-tol)
        False
        '''
        return np.all(abs(np.asarray(center) - self.center) < tol) and abs(radius - self.radius) < tol

    def __repr__(self):
        return "Arc([%0.3f, %0.3f], %0.3f,   %0.3f, %0.3f,   %s,   degrees=%0.3f)" \
                % (self.center[0], self.center[1], self.radius, self.from_angle, self.to_angle, self.direction, self.length_degrees())


class VennDiagram:
    '''
    A container for a set of patches and patch labels and set labels, which make up the rendered venn diagram.
    This object is returned by a venn2 or venn3 function call.
    '''
    id2idx = {'10': 0, '01': 1, '11': 2,
              '100': 0, '010': 1, '110': 2, '001': 3, '101': 4, '011': 5, '111': 6, 'A': 0, 'B': 1, 'C': 2}

    def __init__(self, patches, subset_labels, set_labels, centers, radii):
        self.patches = patches
        self.subset_labels = subset_labels
        self.set_labels = set_labels
        self.centers = centers
        self.radii = radii

    def get_patch_by_id(self, id):
        '''Returns a patch by a "region id".
           A region id is a string '10', '01' or '11' for 2-circle diagram or a
           string like '001', '010', etc, for 3-circle diagram.'''
        return self.patches[self.id2idx[id]]

    def get_label_by_id(self, id):
        '''
        Returns a subset label by a "region id".
        A region id is a string '10', '01' or '11' for 2-circle diagram or a
        string like '001', '010', etc, for 3-circle diagram.
        Alternatively, if the string 'A', 'B'  (or 'C' for 3-circle diagram) is given, the label of the
        corresponding set is returned (or None).'''
        if len(id) == 1:
            return self.set_labels[self.id2idx[id]] if self.set_labels is not None else None
        else:
            return self.subset_labels[self.id2idx[id]]

    def get_circle_center(self, id):
        '''
        Returns the coordinates of the center of a circle as a numpy array (x,y)
        id must be 0, 1 or 2 (corresponding to the first, second, or third circle).
        This is a getter-only (i.e. changing this value does not affect the diagram)
        '''
        return self.centers[id]

    def get_circle_radius(self, id):
        '''
        Returns the radius of circle id (where id is 0, 1 or 2).
        This is a getter-only (i.e. changing this value does not affect the diagram)
        '''
        return self.radii[id]

    def hide_zeroes(self):
        '''
        Sometimes it makes sense to hide the labels for subsets whose size is zero.
        This utility method does this.
        '''
        for v in self.subset_labels:
            if v is not None and v.get_text() == '0':
                v.set_visible(False)


def mix_colors(col1, col2, col3=None):
    '''
    Mixes two colors to compute a "mixed" color (for purposes of computing
    colors of the intersection regions based on the colors of the sets.
    Note that we do not simply compute averages of given colors as those seem
    too dark for some default configurations. Thus, we lighten the combination up a bit.

    Inputs are (up to) three RGB triples of floats 0.0-1.0 given as numpy arrays.

    >>> mix_colors(np.array([1.0, 0., 0.]), np.array([1.0, 0., 0.])) # doctest: +NORMALIZE_WHITESPACE
    array([ 1.,  0.,  0.])
    >>> mix_colors(np.array([1.0, 1., 0.]), np.array([1.0, 0.9, 0.]), np.array([1.0, 0.8, 0.1])) # doctest: +NORMALIZE_WHITESPACE
    array([ 1. ,  1. , 0.04])
    '''
    if col3 is None:
        mix_color = 0.7 * (col1 + col2)
    else:
        mix_color = 0.4 * (col1 + col2 + col3)
    mix_color = np.min([mix_color, [1.0, 1.0, 1.0]], 0)
    return mix_color


def prepare_venn_axes(ax, centers, radii):
    '''
    Sets properties of the axis object to suit venn plotting. I.e. hides ticks, makes proper xlim/ylim.
    '''
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    min_x = min([centers[i][0] - radii[i] for i in range(len(radii))])
    max_x = max([centers[i][0] + radii[i] for i in range(len(radii))])
    min_y = min([centers[i][1] - radii[i] for i in range(len(radii))])
    max_y = max([centers[i][1] + radii[i] for i in range(len(radii))])
    ax.set_xlim([min_x - 0.1, max_x + 0.1])
    ax.set_ylim([min_y - 0.1, max_y + 0.1])
    ax.set_axis_off()


def compute_venn3_areas(diagram_areas, normalize_to=1.0, _minimal_area=1e-6):
    '''
    The list of venn areas is given as 7 values, corresponding to venn diagram areas in the following order:
     (Abc, aBc, ABc, abC, AbC, aBC, ABC)
    (i.e. last element corresponds to the size of intersection A&B&C).
    The return value is a list of areas (A_a, A_b, A_c, A_ab, A_bc, A_ac, A_abc),
    such that the total area of all circles is normalized to normalize_to.
    If the area of any circle is smaller than _minimal_area, makes it equal to _minimal_area.

    Assumes all input values are nonnegative (to be more precise, all areas are passed through and abs() function)
    >>> compute_venn3_areas((1, 1, 0, 1, 0, 0, 0))
    (0.33..., 0.33..., 0.33..., 0.0, 0.0, 0.0, 0.0)
    >>> compute_venn3_areas((0, 0, 0, 0, 0, 0, 0))
    (1e-06, 1e-06, 1e-06, 0.0, 0.0, 0.0, 0.0)
    >>> compute_venn3_areas((1, 1, 1, 1, 1, 1, 1), normalize_to=7)
    (4.0, 4.0, 4.0, 2.0, 2.0, 2.0, 1.0)
    >>> compute_venn3_areas((1, 2, 3, 4, 5, 6, 7), normalize_to=56/2)
    (16.0, 18.0, 22.0, 10.0, 13.0, 12.0, 7.0)
    '''
    # Normalize input values to sum to 1
    areas = np.array(np.abs(diagram_areas), float)
    total_area = np.sum(areas)
    if np.abs(total_area) < _minimal_area:
        warnings.warn("All circles have zero area")
        return (1e-06, 1e-06, 1e-06, 0.0, 0.0, 0.0, 0.0)
    else:
        areas = areas / total_area * normalize_to
        A_a = areas[0] + areas[2] + areas[4] + areas[6]
        if A_a < _minimal_area:
            warnings.warn("Circle A has zero area")
            A_a = _minimal_area
        A_b = areas[1] + areas[2] + areas[5] + areas[6]
        if A_b < _minimal_area:
            warnings.warn("Circle B has zero area")
            A_b = _minimal_area
        A_c = areas[3] + areas[4] + areas[5] + areas[6]
        if A_c < _minimal_area:
            warnings.warn("Circle C has zero area")
            A_c = _minimal_area

        # Areas of the three intersections (ab, ac, bc)
        A_ab, A_ac, A_bc = areas[2] + areas[6], areas[4] + areas[6], areas[5] + areas[6]

        return (A_a, A_b, A_c, A_ab, A_bc, A_ac, areas[6])


def solve_venn3_circles(venn_areas):
    '''
    Given the list of "venn areas" (as output from compute_venn3_areas, i.e. [A, B, C, AB, BC, AC, ABC]),
    finds the positions and radii of the three circles.
    The return value is a tuple (coords, radii), where coords is a 3x2 array of coordinates and
    radii is a 3x1 array of circle radii.

    Assumes the input values to be nonnegative and not all zero.
    In particular, the first three values must all be positive.

    The overall match is only approximate (to be precise, what is matched are the areas of the circles and the
    three pairwise intersections).

    >>> c, r = solve_venn3_circles((1, 1, 1, 0, 0, 0, 0))
    >>> np.round(r, 3)
    array([ 0.564,  0.564,  0.564])
    >>> c, r = solve_venn3_circles(compute_venn3_areas((1, 2, 40, 30, 4, 40, 4)))
    >>> np.round(r, 3)
    array([ 0.359,  0.476,  0.453])
    '''
    (A_a, A_b, A_c, A_ab, A_bc, A_ac, A_abc) = list(map(float, venn_areas))
    r_a, r_b, r_c = np.sqrt(A_a / np.pi), np.sqrt(A_b / np.pi), np.sqrt(A_c / np.pi)
    intersection_areas = [A_ab, A_bc, A_ac]
    radii = np.array([r_a, r_b, r_c])

    # Hypothetical distances between circle centers that assure
    # that their pairwise intersection areas match the requirements.
    dists = [find_distance_by_area(radii[i], radii[j], intersection_areas[i]) for (i, j) in [(0, 1), (1, 2), (2, 0)]]

    # How many intersections have nonzero area?
    num_nonzero = sum(np.array([A_ab, A_bc, A_ac]) > tol)

    # Handle four separate cases:
    #    1. All pairwise areas nonzero
    #    2. Two pairwise areas nonzero
    #    3. One pairwise area nonzero
    #    4. All pairwise areas zero.

    if num_nonzero == 3:
        # The "generic" case, simply use dists to position circles at the vertices of a triangle.
        # Before we need to ensure that resulting circles can be at all positioned on a triangle,
        # use an ad-hoc fix.
        for i in range(3):
            i, j, k = (i, (i + 1) % 3, (i + 2) % 3)
            if dists[i] > dists[j] + dists[k]:
                a, b = (j, k) if dists[j] < dists[k] else (k, j)
                dists[i] = dists[b] + dists[a]*0.8
                warnings.warn("Bad circle positioning")
        coords = position_venn3_circles_generic(radii, dists)
    elif num_nonzero == 2:
        # One pair of circles is not intersecting.
        # In this case we can position all three circles in a line
        # The two circles that have no intersection will be on either sides.
        for i in range(3):
            if intersection_areas[i] < tol:
                (left, right, middle) = (i, (i + 1) % 3, (i + 2) % 3)
                coords = np.zeros((3, 2))
                coords[middle][0] = dists[middle]
                coords[right][0] = dists[middle] + dists[right]
                # We want to avoid the situation where left & right still intersect
                if coords[left][0] + radii[left] > coords[right][0] - radii[right]:
                    mid = (coords[left][0] + radii[left] + coords[right][0] - radii[right]) / 2.0
                    coords[left][0] = mid - radii[left] - 1e-5
                    coords[right][0] = mid + radii[right] + 1e-5
                break
    elif num_nonzero == 1:
        # Only one pair of circles is intersecting, and one circle is independent.
        # Position all on a line first two intersecting, then the free one.
        for i in range(3):
            if intersection_areas[i] > tol:
                (left, right, side) = (i, (i + 1) % 3, (i + 2) % 3)
                coords = np.zeros((3, 2))
                coords[right][0] = dists[left]
                coords[side][0] = dists[left] + radii[right] + radii[side] * 1.1  # Pad by 10%
                break
    else:
        # All circles are non-touching. Put them all in a sequence
        coords = np.zeros((3, 2))
        coords[1][0] = radii[0] + radii[1] * 1.1
        coords[2][0] = radii[0] + radii[1] * 1.1 + radii[1] + radii[2] * 1.1

    coords = normalize_by_center_of_mass(coords, radii)
    return (coords, radii)


def position_venn3_circles_generic(radii, dists):
    '''
    Given radii = (r_a, r_b, r_c) and distances between the circles = (d_ab, d_bc, d_ac),
    finds the coordinates of the centers for the three circles so that they form a proper triangle.
    The current positioning method puts the center of A and B on a horizontal line y==0,
    and C just below.

    Returns a 3x2 array with circle center coordinates in rows.

    >>> position_venn3_circles_generic((1, 1, 1), (0, 0, 0))
    array([[ 0.,  0.],
           [ 0.,  0.],
           [ 0., -0.]])
    >>> position_venn3_circles_generic((1, 1, 1), (2, 2, 2))
    array([[ 0.        ,  0.        ],
           [ 2.        ,  0.        ],
           [ 1.        , -1.73205081]])
    '''
    (d_ab, d_bc, d_ac) = dists
    (r_a, r_b, r_c) = radii
    coords = np.array([[0, 0], [d_ab, 0], [0, 0]], float)
    C_x = (d_ac**2 - d_bc**2 + d_ab**2) / 2.0 / d_ab if np.abs(d_ab) > tol else 0.0
    C_y = -np.sqrt(d_ac**2 - C_x**2)
    coords[2, :] = C_x, C_y
    return coords


def compute_venn3_regions(centers, radii):
    '''
    Given the 3x2 matrix with circle center coordinates, and a 3-element list (or array) with circle radii [as returned from solve_venn3_circles],
    returns the 7 regions, comprising the venn diagram, as VennRegion objects.

    Regions are returned in order (Abc, aBc, ABc, abC, AbC, aBC, ABC)

    >>> centers, radii = solve_venn3_circles((1, 1, 1, 1, 1, 1, 1))
    >>> regions = compute_venn3_regions(centers, radii)
    '''
    A = VennCircleRegion(centers[0], radii[0])
    B = VennCircleRegion(centers[1], radii[1])
    C = VennCircleRegion(centers[2], radii[2])
    Ab, AB = A.subtract_and_intersect_circle(B.center, B.radius)
    ABc, ABC = AB.subtract_and_intersect_circle(C.center, C.radius)
    Abc, AbC = Ab.subtract_and_intersect_circle(C.center, C.radius)
    aB, _ = B.subtract_and_intersect_circle(A.center, A.radius)
    aBc, aBC = aB.subtract_and_intersect_circle(C.center, C.radius)
    aC, _ = C.subtract_and_intersect_circle(A.center, A.radius)
    abC, _ = aC.subtract_and_intersect_circle(B.center, B.radius)
    return [Abc, aBc, ABc, abC, AbC, aBC, ABC]


def compute_venn3_colors(set_colors):
    '''
    Given three base colors, computes combinations of colors corresponding to all regions of the venn diagram.
    returns a list of 7 elements, providing colors for regions (100, 010, 110, 001, 101, 011, 111).

    >>> compute_venn3_colors(['r', 'g', 'b'])
    (array([ 1.,  0.,  0.]),..., array([ 0.4,  0.2,  0.4]))
    '''
    ccv = ColorConverter()
    base_colors = [np.array(ccv.to_rgb(c)) for c in set_colors]
    return (base_colors[0], base_colors[1], mix_colors(base_colors[0], base_colors[1]), base_colors[2],
            mix_colors(base_colors[0], base_colors[2]), mix_colors(base_colors[1], base_colors[2]), mix_colors(base_colors[0], base_colors[1], base_colors[2]))


def compute_venn3_subsets(a, b, c):
    '''
    Given three set or Counter objects, computes the sizes of (a & ~b & ~c, ~a & b & ~c, a & b & ~c, ....),
    as needed by the subsets parameter of venn3 and venn3_circles.
    Returns the result as a tuple.

    >>> compute_venn3_subsets(set([1,2,3]), set([2,3,4]), set([3,4,5,6]))
    (1, 0, 1, 2, 0, 1, 1)
    >>> compute_venn3_subsets(Counter([1,2,3]), Counter([2,3,4]), Counter([3,4,5,6]))
    (1, 0, 1, 2, 0, 1, 1)
    >>> compute_venn3_subsets(Counter([1,1,1]), Counter([1,1,1]), Counter([1,1,1,1]))
    (0, 0, 0, 1, 0, 0, 3)
    >>> compute_venn3_subsets(Counter([1,1,2,2,3,3]), Counter([2,2,3,3,4,4]), Counter([3,3,4,4,5,5,6,6]))
    (2, 0, 2, 4, 0, 2, 2)
    >>> compute_venn3_subsets(Counter([1,2,3]), Counter([2,2,3,3,4,4]), Counter([3,3,4,4,4,5,5,6]))
    (1, 1, 1, 4, 0, 3, 1)
    >>> compute_venn3_subsets(set([]), set([]), set([]))
    (0, 0, 0, 0, 0, 0, 0)
    >>> compute_venn3_subsets(set([1]), set([]), set([]))
    (1, 0, 0, 0, 0, 0, 0)
    >>> compute_venn3_subsets(set([]), set([1]), set([]))
    (0, 1, 0, 0, 0, 0, 0)
    >>> compute_venn3_subsets(set([]), set([]), set([1]))
    (0, 0, 0, 1, 0, 0, 0)
    >>> compute_venn3_subsets(Counter([]), Counter([]), Counter([1]))
    (0, 0, 0, 1, 0, 0, 0)
    >>> compute_venn3_subsets(set([1]), set([1]), set([1]))
    (0, 0, 0, 0, 0, 0, 1)
    >>> compute_venn3_subsets(set([1,3,5,7]), set([2,3,6,7]), set([4,5,6,7]))
    (1, 1, 1, 1, 1, 1, 1)
    >>> compute_venn3_subsets(Counter([1,3,5,7]), Counter([2,3,6,7]), Counter([4,5,6,7]))
    (1, 1, 1, 1, 1, 1, 1)
    >>> compute_venn3_subsets(Counter([1,3,5,7]), set([2,3,6,7]), set([4,5,6,7]))
    Traceback (most recent call last):
    ...
    ValueError: All arguments must be of the same type
    '''
    if not (type(a) == type(b) == type(c)):
        raise ValueError("All arguments must be of the same type")
    set_size = len if type(a) != Counter else lambda x: sum(x.values())   # We cannot use len to compute the cardinality of a Counter
    return (set_size(a - (b | c)),  # TODO: This is certainly not the most efficient way to compute.
        set_size(b - (a | c)),
        set_size((a & b) - c),
        set_size(c - (a | b)),
        set_size((a & c) - b),
        set_size((b & c) - a),
        set_size(a & b & c))


def venn3_circles(subsets, normalize_to=1.0, alpha=1.0, color='black', linestyle='solid', linewidth=2.0, ax=None, **kwargs):
    '''
    Plots only the three circles for the corresponding Venn diagram.
    Useful for debugging or enhancing the basic venn diagram.
    parameters ``subsets``, ``normalize_to`` and ``ax`` are the same as in venn3()
    kwargs are passed as-is to matplotlib.patches.Circle.
    returns a list of three Circle patches.

        >>> plot = venn3_circles({'001': 10, '100': 20, '010': 21, '110': 13, '011': 14})
        >>> plot = venn3_circles([set(['A','B','C']), set(['A','D','E','F']), set(['D','G','H'])])
    '''
    # Prepare parameters
    if isinstance(subsets, dict):
        subsets = [subsets.get(t, 0) for t in ['100', '010', '110', '001', '101', '011', '111']]
    elif len(subsets) == 3:
        subsets = compute_venn3_subsets(*subsets)

    areas = compute_venn3_areas(subsets, normalize_to)
    centers, radii = solve_venn3_circles(areas)

    if ax is None:
        ax = gca()
    prepare_venn_axes(ax, centers, radii)
    result = []
    for (c, r) in zip(centers, radii):
        circle = Circle(c, r, alpha=alpha, edgecolor=color, facecolor='none', linestyle=linestyle, linewidth=linewidth, **kwargs)
        ax.add_patch(circle)
        result.append(circle)
    return result


def venn3(subsets, set_labels=('A', 'B', 'C'), set_colors=('r', 'g', 'b'), alpha=0.4, normalize_to=1.0, ax=None, subset_label_formatter=None):
    '''Plots a 3-set area-weighted Venn diagram.
    The subsets parameter can be one of the following:
     - A list (or a tuple), containing three set objects.
     - A dict, providing sizes of seven diagram regions.
       The regions are identified via three-letter binary codes ('100', '010', etc), hence a valid set could look like:
       {'001': 10, '010': 20, '110':30, ...}. Unmentioned codes are considered to map to 0.
     - A list (or a tuple) with 7 numbers, denoting the sizes of the regions in the following order:
       (100, 010, 110, 001, 101, 011, 111).

    ``set_labels`` parameter is a list of three strings - set labels. Set it to None to disable set labels.
    The ``set_colors`` parameter should be a list of three elements, specifying the "base colors" of the three circles.
    The colors of circle intersections will be computed based on those.

    The ``normalize_to`` parameter specifies the total (on-axes) area of the circles to be drawn. Sometimes tuning it (together
    with the overall fiture size) may be useful to fit the text labels better.
    The return value is a ``VennDiagram`` object, that keeps references to the ``Text`` and ``Patch`` objects used on the plot
    and lets you know the centers and radii of the circles, if you need it.

    The ``ax`` parameter specifies the axes on which the plot will be drawn (None means current axes).

    The ``subset_label_formatter`` parameter is a function that can be passed to format the labels
    that describe the size of each subset.

    Note: if some of the circles happen to have zero area, you will probably not get a nice picture.

    >>> import matplotlib # (The first two lines prevent the doctest from falling when TCL not installed. Not really necessary in most cases)
    >>> matplotlib.use('Agg')
    >>> from matplotlib_venn import *
    >>> v = venn3(subsets=(1, 1, 1, 1, 1, 1, 1), set_labels = ('A', 'B', 'C'))
    >>> c = venn3_circles(subsets=(1, 1, 1, 1, 1, 1, 1), linestyle='dashed')
    >>> v.get_patch_by_id('100').set_alpha(1.0)
    >>> v.get_patch_by_id('100').set_color('white')
    >>> v.get_label_by_id('100').set_text('Unknown')
    >>> v.get_label_by_id('C').set_text('Set C')

    You can provide sets themselves rather than subset sizes:
    >>> v = venn3(subsets=[set([1,2]), set([2,3,4,5]), set([4,5,6,7,8,9,10,11])])
    >>> print("%0.2f %0.2f %0.2f" % (v.get_circle_radius(0), v.get_circle_radius(1)/v.get_circle_radius(0), v.get_circle_radius(2)/v.get_circle_radius(0)))
    0.24 1.41 2.00
    >>> c = venn3_circles(subsets=[set([1,2]), set([2,3,4,5]), set([4,5,6,7,8,9,10,11])])
    '''
    # Prepare parameters
    if isinstance(subsets, dict):
        subsets = [subsets.get(t, 0) for t in ['100', '010', '110', '001', '101', '011', '111']]
    elif len(subsets) == 3:
        subsets = compute_venn3_subsets(*subsets)

    if subset_label_formatter is None:
        subset_label_formatter = str

    areas = compute_venn3_areas(subsets, normalize_to)
    centers, radii = solve_venn3_circles(areas)
    regions = compute_venn3_regions(centers, radii)
    colors = compute_venn3_colors(set_colors)

    # Remove regions that are too small from the diagram
    MIN_REGION_SIZE = 1e-4
    for i in range(len(regions)):
        if regions[i].size() < MIN_REGION_SIZE and subsets[i] == 0:
            regions[i] = VennEmptyRegion()

    # There is a rare case (Issue #12) when the middle region is visually empty
    # (the positioning of the circles does not let them intersect), yet the corresponding value is not 0.
    # we address it separately here by positioning the label of that empty region in a custom way
    if isinstance(regions[6], VennEmptyRegion) and subsets[6] > 0:
        intersections = [circle_circle_intersection(centers[i], radii[i], centers[j], radii[j]) for (i, j) in [(0, 1), (1, 2), (2, 0)]]
        middle_pos = np.mean([i[0] for i in intersections], 0)
        regions[6] = VennEmptyRegion(middle_pos)

    if ax is None:
        ax = gca()
    prepare_venn_axes(ax, centers, radii)

    # Create and add patches and text
    patches = [r.make_patch() for r in regions]
    for (p, c) in zip(patches, colors):
        if p is not None:
            p.set_facecolor(c)
            p.set_edgecolor('none')
            p.set_alpha(alpha)
            ax.add_patch(p)
    label_positions = [r.label_position() for r in regions]
    subset_labels = [ax.text(lbl[0], lbl[1], subset_label_formatter(s), va='center', ha='center') if lbl is not None else None for (lbl, s) in zip(label_positions, subsets)]

    # Position labels
    if set_labels is not None:
        # There are two situations, when set C is not on the same line with sets A and B, and when the three are on the same line.
        if abs(centers[2][1] - centers[0][1]) > tol:
            # Three circles NOT on the same line
            label_positions = [centers[0] + np.array([-radii[0] / 2, radii[0]]),
                               centers[1] + np.array([radii[1] / 2, radii[1]]),
                               centers[2] + np.array([0.0, -radii[2] * 1.1])]
            labels = [ax.text(pos[0], pos[1], txt, size='large') for (pos, txt) in zip(label_positions, set_labels)]
            labels[0].set_horizontalalignment('right')
            labels[1].set_horizontalalignment('left')
            labels[2].set_verticalalignment('top')
            labels[2].set_horizontalalignment('center')
        else:
            padding = np.mean([r * 0.1 for r in radii])
            # Three circles on the same line
            label_positions = [centers[0] + np.array([0.0, - radii[0] - padding]),
                               centers[1] + np.array([0.0, - radii[1] - padding]),
                               centers[2] + np.array([0.0, - radii[2] - padding])]
            labels = [ax.text(pos[0], pos[1], txt, size='large', ha='center', va='top') for (pos, txt) in zip(label_positions, set_labels)]
    else:
        labels = None
    return VennDiagram(patches, subset_labels, labels, centers, radii)
