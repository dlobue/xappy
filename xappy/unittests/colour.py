# Copyright (C) 2008 Lemur Consulting Ltd
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
import itertools
from xappytest import *
from xappy.fieldactions import FieldActions
import xappy.colour
import colormath

colourdata = [[(0, 0, 0), 10, 0.01], #black
              [(255, 255, 255), 5,  0.01 ], #white
              [(0, 0, 255), 3, 0.01], #blue
              [(0, 255, 0), 2, 0.01], #green
              [(255, 0, 0), 1, 0.01]  #red
              ]

'''
class ColourTermsTestCase(TestCase):

    def test_within_range(self):
        """ Ensure that all the colour from
        rgb data in the 0-256 range result in LAB colours in the
        expected range.

        """
        min_a = min_l = min_b = 100000
        max_l = max_a = max_b = 0
        for r in xrange(256):
            for g in xrange(256):
                for b in xrange(256):
                    lab = xappy.colour.rgb2bucket((r, g, b), 25)
                    min_l = min(lab[0], min_l)
                    min_a = min(lab[1], min_a)
                    min_b = min(lab[2], min_b)
                    max_l = max(lab[0], max_l)
                    max_a = max(lab[1], max_a)
                    max_b = max(lab[2], max_b)
                    self.assert_(xappy.colour.check_in_range(lab))
        print "%f %f %f %f %f %f" % (min_l, max_l, min_a, max_a, min_b, max_b)
'''
        
class ColourIndexingTestCase(TestCase):

    def pre_test(self):
        self.dbpath = os.path.join(self.tempdir, 'db')
        self.iconn = xappy.IndexerConnection(self.dbpath)

    def post_test(self):
        self.iconn.close()

    def test_basic_indexing(self):
        self.iconn.add_field_action('colour', xappy.FieldActions.COLOUR)
        doc = xappy.UnprocessedDocument()
        for col, freq, spread in colourdata:
            colourterm = xappy.colour.rgb2term(col, 50)
            doc.fields.append(xappy.Field('colour', [(colourterm, freq)]))
        self.iconn.add(doc)

    def test_normalisation(self):
        self.iconn.add_field_action('colour', xappy.FieldActions.COLOUR)
        doc = xappy.UnprocessedDocument()
        for col, freq, spread in colourdata:
            colourterm =  xappy.colour.rgb2term(col, 50)
            doc.fields.append(xappy.Field('colour', [(colourterm, freq)]))
        did = self.iconn.add(doc)
        doc = self.iconn.get_document(did)
        prefix = doc._fieldmappings.get_prefix('colour')
        cumfreq = 0
        for t in doc._doc.termlist():
            if t.term.startswith(prefix):
                cumfreq += t.wdf
        # we may have some rounding errors, so allow a bit of slack
        self.assert_(995 <= cumfreq <= 1005)

class ClusterTestCase(TestCase):

    def test_lab_clusters(self):
        # hopefully two clusters here
        clustercols = [(128, 128, 128),
                       (127, 127, 127),
                       (126, 126, 126),
                       (50, 200, 20)]

        lab_cols = itertools.imap(xappy.colour.rgb2lab, clustercols)

        clusters = xappy.colour.cluster_coords(lab_cols)
        #print list(clusters)
        self.assertEqual(2, len(list(clusters)))

class ColourSearchTestCase(TestCase):

    def pre_test(self):
        self.dbpath = os.path.join(self.tempdir, 'db')
        self.iconn = xappy.IndexerConnection(self.dbpath)
        self.iconn.add_field_action('colour', xappy.FieldActions.COLOUR)
        for col, freq, spread in colourdata:
            doc = xappy.UnprocessedDocument()
            colourterm = xappy.colour.rgb2term(col, 50)
            doc.fields.append(xappy.Field('colour', [(colourterm, freq)]))
            self.iconn.add(doc)
        self.iconn.close()
        self.sconn = xappy.SearchConnection(self.dbpath)

    def post_test(self):
        self.sconn.close()

    def test_basic_search(self):
        """ Ensure that query_colour can be called with
        generating an error.

        """
        
        xappy.colour.query_colour(self.sconn, 'colour', colourdata, 50)
        

    def test_correct_colour(self):
        """ Check that correct document is found when querying for the
        colours that were supplied at indexing time.

        """

        prefix = self.sconn._field_mappings.get_prefix('colour')

        action_params = self.sconn._field_actions['colour']._actions[FieldActions.COLOUR][0]


        print [t.term for t in self.sconn._index.allterms()]
        for colour, frequency, spread in colourdata:
            print colour, frequency, spread
            query = xappy.colour.query_colour(
                self.sconn, 'colour', [[colour, frequency, spread]], 50,
                clustering=False)
            print query
            results = self.sconn.search(query, 0, 10)
            r = results[0]
            d = self.sconn.get_document(r.id)
            terms = d._doc.termlist()
            terms = [x.term for x in terms if x.term.startswith(prefix)]
            terms = set(terms)
            print terms
            colourterm = xappy.colour.rgb2term(colour, 50)
            rgb_term = prefix + colourterm
            print rgb_term
            self.assert_(rgb_term in terms)

if __name__ == '__main__':
    main()
