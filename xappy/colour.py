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
"""
colour.py utilities for bucketing and clustering colours.

In the database colours are stored as terms representing a 'bucket'.
Each bucket represent a cube in the lab colour space. The limits of
the space that allocate buckets to is given my the conversions from
rgb to lab by the colormath module.

"""
import collections
import operator
import itertools
import marshall

import colormath
import colormath.color_objects
import numpy
import scipy.cluster
import xapian
import xappy

def compute_range_limits(dim=256):
    """ This finds the extremes of the Lab coordinates by iterating
    over all possible rgb colours (assuming `dim` steps in each of the
    rgb axes).

    Warning: this is slow, but it's only needed for checking and
    testing, not for indexing or query generation.
    """
    min_l = min_a = min_b = 10000000
    max_l = max_a = max_b = -10000000
    try:
        gen = itertools.product(xrange(dim), repeat = 3)
    except AttributeError:
        # itertools.product new in 2.6
        gen = ((x,y,z)
               for x in xrange(dim)
               for y in xrange(dim)
               for z in xrange(dim))
    for rgb_coords in gen:
        rgb = colormath.color_objects.RGBColor(*rgb_coords)
        lab = rgb.convert_to('lab')
        min_l = min(min_l, lab.lab_l)
        min_a = min(min_a, lab.lab_a)
        min_b = min(min_b, lab.lab_b)
        max_l = max(max_l, lab.lab_l)
        max_a = max(max_a, lab.lab_a)
        max_b = max(max_b, lab.lab_b)
    return (min_l, max_l), (min_a, max_a), (min_b, max_b)    

# in order to decide on the size of the buckets we want to know the
# possible range of value that each coordinate can take.

#lab_ranges = compute_lab_ranges()

lab_ranges = (
    (0.0, 99.999984533331272),
    (-86.182949405160798, 98.235320176646439),
    (-107.86546414496824, 94.477318179693782))

max_distance = pow(sum(pow(x[1]-x[0], 2) for x in lab_ranges), 0.5)

step_size_cache = {}

def step_sizes(step_count):
    try:
        return step_size_cache[step_count]
    except KeyError:
        sizes = tuple((r[1] - r[0]) / float(step_count) for r in lab_ranges)
        step_size_cache[step_count] = sizes
        return sizes

def rgb2lab(rgb):
    rgb = colormath.color_objects.RGBColor(*rgb)
    return rgb.convert_to('lab').get_value_tuple()

def check_in_range(lab):
    l, a, b = lab
    r = lab_ranges
    return ( (r[0][0] <= l <= r[0][1]) and
             (r[1][0] <= a <= r[1][1]) and
             (r[2][0] <= b <= r[2][1]) )

def compute_steps(count):    
    return ( (float(r[1]) - float(r[0])) / float(count) for
             r in lab_ranges)


def lab2bucket(lab, step_count):
    """ return the indices of the bucket within which the point l, a,
    b falls, assuming that the space is divided up into `step_count`
    steps in each coordinate.

    """
    l, a, b = lab
    l_step, a_step, b_step = step_sizes(step_count)
    return tuple(int(x) for x in ((l - lab_ranges[0][0]) // l_step,
                                  (a - lab_ranges[1][0]) // a_step,
                                  (b - lab_ranges[2][0]) // b_step))

def rgb2bucket(rgb, step_count):
    return lab2bucket(rgb2lab(rgb), step_count)

def encode_bucket(bucket_indices, step_count):
    """ Return a hex string identifying the supplied bucket. Buckets
    are numbered according to their position in the lexicographic
    ordering of their coordinates.
    """
    
    l, a, b = bucket_indices
    position = (l +
                step_count * a +
                pow(step_count, 2) * b)
    return hex(position)

def decode_bucket(bucket_id, step_count):
    """ return the bucket indices of a string encoded with
    encode_bucket"""
    val = int(bucket_id, 16)
    b, rem = divmod(val, pow(step_count,2))
    a, l = divmod(rem, step_count)
    return l, a, b

def lab2term(lab, step_count):
    return encode_bucket(lab2bucket(lab, step_count), step_count)

def rgb2term(rgb, step_count):
    return lab2term(rgb2lab(rgb), step_count)

# synonym
term2lab = decode_bucket

def cluster_coords(coords, coord_fun=None, distance_factor=0.05):
    distance = distance_factor * max_distance
    coord_list = list(coords)
    source = (coord_list if coord_fun is None
              else map(coord_fun, coord_list))
    if len(source) < 2:
        yield coord_list
    else:
        coord_array = numpy.array(source)
    
        clusters = scipy.cluster.hierarchy.fclusterdata(
            coord_array, distance, criterion='distance')

        def keyf(c):
            return clusters[c[0]]

        sfreqs = sorted(enumerate(coord_list), key=keyf)
        groups =  itertools.groupby(sfreqs, keyf)
        for k, group in groups:
            yield map(operator.itemgetter(1), group)

def cluster_terms(terms, step_count, distance_factor=0.05):
    coord_fun = lambda t: term2lab(t, step_count)
    return cluster_coords(
        terms, coord_fun=coord_fun, distance_factor = distance_factor)

def average_weights(terms_and_weights):
    """ terms and weights consists of a iterable, each element of
    which is an itereable of (item, weight) pairs. The weights should
    be numbers. The output is groups of (item, average_weight), pairs
    where the average_weight is the arithmetic mean of the original
    weights in the group.

    """

    for group in terms_and_weights:
        g = list(group)
        average =  sum(iterools.imap(operator.itemgetter(0), g)) / len(g)
        yield [(t, average) for t, _ in g]

# this could all be pushed down into numpy, but colormath needs a
# colour object. This is a candidate for speeding up if performance is
# problematic.

#FIXME: with small precisions we don't get the terms for the target
#colour, but just the one(s) below

def terms_and_weights(colour_freqs, step_count):
    weight_dict = collections.defaultdict(float)
    l_step, a_step, b_step = step_sizes(step_count)
    print l_step, a_step, b_step
    for col, freq, spread in colour_freqs:
        origin = colormath.color_objects.RGBColor(*col).convert_to('lab')
        l, a, b = origin.get_value_tuple()
        def coord_start(coord, lo, hi):
            return max(coord -  (hi - lo) * spread / 2.0, lo)
        
        lmin = coord_start(l, *lab_ranges[0])
        amin = coord_start(a, *lab_ranges[1])
        bmin = coord_start(b, *lab_ranges[2])

        def coord_end(coord, step, lo, hi):
            return min(coord + (hi - lo) * spread / 2.0, hi+step)
        lmax = coord_end(l, l_step, *lab_ranges[0])
        amax = coord_end(a, a_step, *lab_ranges[1])
        bmax = coord_end(b, b_step, *lab_ranges[2])

        print lmin, amin, bmin, lmax, amax, bmax
        # we're actually making terms for a cube centred on the origin
        # here - quite possibly a sphere makes more sense
        # if prec is 0 we should just get one term - check!
        l  = lmin
        while l <= lmax:
            a = amin
            while a <= amax:
                b = bmin
                while b <= bmax:
                    term = lab2term((l, a, b), step_count)
                    distance = origin.delta_e(
                        colormath.color_objects.LabColor(l, a, b))
                    weight_dict[term] += freq / (1.0 + distance)
                    b += b_step
                    print "b", b

                a += a_step
                print "a", a
                            
            l += l_step
            print "l", l
        return weight_dict

def query_colour(sconn, field, colour_freqs, step_count, clustering=False):
    """ Generate a query to find document with similar colours in
    `field` to those specified in `colour_freqs`. `colour_freqs`
    should be at iterable whose members are lists or tuples
    consisting of 3 data. These being (in order) a sequence
    consisting rgb colour coordinates, each in the range 0-255; a
    frequency measure and a precision measure.
    
    If `clustering` is True then individual colours will be grouped
    together into clusters, and the total frequency for the
    cluster used to weight terms for its consituent colours.

    If `clustering` is False then no clustering will be done and each
    frequency is simply used to weight the terms generated from
    that colour.

    In either case each colour will be used to generate terms for
    colours close to that colour, with decreasing weights as the
    distance increases. The number of terms generated is
    controlled by the precision, which indicates the percentage of
    the total range of colour values represented by the
    colour. Note that the higher this value the more terms that
    will be generated, which may affect performance. A value of 0
    means that only one term will be generated. (It is not
    possible to exclude a colour completely with this mechanism -
    simply omit it from `colour_freqs` to achieve this.)
    
    """

    prefix = sconn._field_mappings.get_prefix(field)

    def term_subqs(ts):
        return [xappy.Query(xapian.Query(prefix + term)) * weight for
                term, weight in ts.iteritems()]

    clusters = (cluster_coords(
        colour_freqs, coord_fun=operator.itemgetter(0)) if clustering
                else [colour_freqs])
    
    for cluster in clusters:
        subqs = []
        weighted_terms = terms_and_weights(cluster, step_count)
        subqs.append(xappy.Query.compose(xappy.Query.OP_OR,
                                         term_subqs(weighted_terms)))
    return xappy.Query.compose(xappy.Query.OP_AND, subqs)
