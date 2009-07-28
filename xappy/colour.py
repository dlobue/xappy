# Copyright (C) 2009 Lemur Consulting Ltd
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


import collections
import itertools
import marshall

import colormath
import colormath.color_objects
import numpy
import scipy.cluster
import xapian

# in order to decide on the size of the buckets we want to know the
# possible range of value that each coordinate can take.

lab_ranges_from_srgb = ((0, 100), (-128.0, 128.0), (-128.0, 128.0))


def check_in_range(l, a, b):
    r = lab_ranges_from_srgb
    return ( (r[0][0] <= l <= r[0][1]) and
             (r[1][0] <= a <= r[1][1]) and
             (r[2][0] <= b <= r[2][1]) )


def compute_steps(count, source_space='rgb'):    
    if source_space != 'rgb':
        raise NotImplementedError('source space: %s is not supported')
    return ( (float(r[1]) - float(r[0])) / float(count) for
             r in lab_ranges_from_srgb)


def bucket_for_point(l, a, b, l_step, a_step, b_step):
    """ return the coordinates of the least corner of the bucket
    within which the point l, a, b falls, assuming that the space is
    divided up into cuboids of dimensions l_step, a_step, b_step.

    """

    return ((l // l_step) * l_step,
            (a // a_step) * a_step,
            (b // b_step) * b_step)

def bucket_for_rgb_point(r, g, b, steps):
    rbg = colormath.color_objects.RGBColor(r, g, b)
    lab = rbg.convert_to('lab')
    return bucket_for_point(lab.lab_l, lab.lab_a, lab.lab_b, *steps)


def term_for_lab_point(l, a, b, l_step, a_step, b_step):
    bucket_coords = bucket_for_point(l, a, b,
                                     l_step, a_step, b_step)
    return repr(bucket_coords)

def term_for_rgb_point(r, g, b, l_step, a_step, b_step):
    rgb = colormath.color_objects.RGBColor(r, g, b)
    lab = rgb.convert_to('lab')
    assert (lab_ranges_from_srgb[0][0] <= lab.lab_l <= lab_ranges_from_srgb[0][1])
    assert (lab_ranges_from_srgb[1][0] <= lab.lab_a <= lab_ranges_from_srgb[1][1])
    assert (lab_ranges_from_srgb[2][0] <= lab.lab_b <= lab_ranges_from_srgb[2][1])
    
    return term_for_lab_point(lab.lab_l, lab.lab_a, lab.lab_b,
                              l_step, a_step, b_step)

class ColourWeight(xapian.Weight):
    def __init__(self, distance_factor, term):
        self.distance_factor = distance_factor
        self.term = term

    def name(self):
        return "Colour"

    def serialize(self):
        return ""

    def get_sumpart(*args):
        return 1

    def get_maxpart(*args):
        return 1

    def get_sumextra(*args):
        return 0

    def get_maxextra(*args):
        return 0


# this could all be pushed down into numpy, but colormath needs a
# colour object. This is a candidate for speeding up if performance is
# problematic.

def terms_and_weights(colour_freqs, l_step, a_step, b_step):
    weight_dict = collections.defaultdict(float)
    for col, freq, prec in colour_freqs:
        origin = colormath.color_objects.RGBColor(*col).convert_to('lab')
        
        def coord_start(coord, lo, hi):
            return max(coord -  (hi - lo) * prec / 2.0, lo)
        
        lmin = coord_start(origin.lab_l, *lab_ranges_from_srgb[0])
        amin = coord_start(origin.lab_a, *lab_ranges_from_srgb[1])
        bmin = coord_start(origin.lab_b, *lab_ranges_from_srgb[2])

        def coord_end(coord, lo, hi):
            return min(coord + (hi - lo) * prec / 2.0, hi)
        
        lmax = coord_end(origin.lab_l, *lab_ranges_from_srgb[0])
        amax = coord_end(origin.lab_a, *lab_ranges_from_srgb[1])
        bmax = coord_end(origin.lab_b, *lab_ranges_from_srgb[2])

        # we're actually making terms for a cube centred on the origin
        # here - quite possibly a sphere makes more sense
        # if prec is 0 we should just get one term - check!
        l  = lmin
        while l <= lmax:
            a = amin
            while a <= amax:
                b = bmin
                while b <= bmax:
                    term = term_for_lab_point(l, a, b, l_step, a_step, b_step)
                    distance = origin.delta_e(colormath.color_objects.LabColor(l, a, b))
                    weight_dict[term] += freq / (1.0 + distance)
                    b += b_step
                a += a_step
            l += l_step
        return weight_dict

def do_clustering(colour_freqs):
    """ Compute the clusters for the colours in colour_freqs. Generate
    a sequence of colour frequency clusters, the frequency for each
    colour within a cluster being the arithmentic mean of the frequency
    for the cluster.

    """

    colour_list = [x[0] for x in colour_freqs]

    if len(colour_list) <= 1:
        yield colour_freqs
    else:
        # turn the colour data into a numpy array
        colours = numpy.array(colour_list)

        #second parameter may need fiddling with - it's a measure of the
        #inconsistency of the clusters.
        clusters =  scipy.cluster.hierarchy.fclusterdata(colours, 0.5)

        cluster_sums = collections.defaultdict(float)
        cluster_members = collections.defaultdict(list)
    
        def keyf(c):
            return clusters[c[0]]
    
        sfreqs = sorted(enumerate(colour_freqs), key=keyf)
        groups =  itertools.groupby(sfreqs, keyf)
        
        for i, group in groups:
            group = [ x[1] for x in group]
            count = len(group)
            sum_freq = 0.0
            for cf in group:
                sum_freq += cf[1]
                mean_freq = sum_freq / count
            for cf in group:
                cf[1] = mean_freq
            yield group
