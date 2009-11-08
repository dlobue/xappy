#!/usr/bin/env python
#
# Copyright (C) 2009 Richard Boulton
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
r"""generic.py: Base cachemanager classes.

"""
__docformat__ = "restructuredtext en"

import cPickle
import UserDict

class CacheManager(object):
    """Base class for caches of precalculated results.

    """
    def prepare_iter_by_docid(self):
        """Prepare to iterate by document ID.

        This does any preparations necessary for iter_by_docid().

        Caches should keep track of whether any changes have been made since
        this was last called, and use the values calculated last time if not.

        """
        # Naive implementation: build up all the data in memory.
        if getattr(self, 'inverted_items', None) is None:
            items = {}
            for queryid in self.iter_queryids():
                for rank, docid in enumerate(self.get_hits(queryid)):
                    items.setdefault(docid, []).append((queryid, rank))
            self.inverted_items = items

    def invalidate_iter_by_docid(self):
        """Invalidate any cached items for the iter_by_docid.

        """
        self.inverted_items = None

    def iter_by_docid(self):
        """Return an iterator which returns all the documents with cached hits,
        in order of document ID, together with the queryids and ranks for those
        queries.

        Returns (docid, <list of (queryid, rank)>) pairs.

        The default implementation is fairly naive, and builds up all the data
        in memory before iterating through it.  Subclasses can override this
        implementation if they wish.

        """
        self.prepare_iter_by_docid()
        items = self.inverted_items
        for docid in sorted(items.keys()):
            yield docid, items[docid]

    def iter_queryids(self):
        """Return an iterator returning all the queryids for which there are
        cached items.

        Doesn't guarantee any ordering on the queryids.

        """
        raise NotImplementedError

    def get_queryid(self, query_str):
        """Get a (numeric) query ID given a query string.

        If the query string is not in the cache, returns None.

        """
        raise NotImplementedError

    def get_or_make_queryid(self, query_str):
        """Get or allocate a (numeric) query ID given a query string.

        If the query string is not in the cache, makes a new query ID and
        return it.  Otherwise, returns the existing query ID.

        """
        raise NotImplementedError

    def get_hits(self, queryid, startrank=0, endrank=None):
        """Get the Xapian document IDs of documents matching a query.

        This returns a list of document IDs, in ascending rank order.

        `queryid` is the numeric ID of the query to look up.

        `startrank` is the rank (starting at 0) of the first result to return.

        `endrank` is the rank at which to stop returning results: ie, 1 more
        than the rank of the last result to return.  If `endrank` is None, return
        all the results.

        """
        raise NotImplementedError

    def set_hits(self, queryid, docids):
        """Set the Xapian document IDs of documents matching a query.

        `queryid` is the numeric ID of the query to look up.

        """
        raise NotImplementedError

    def remove_hits(self, queryid, ranks_and_docids):
        """Remove the hits at given ranks from the cached entry for a query.

        `queryid` is the numeric ID of the query to look up.

        ranks_and_docids is a iterable of (rank, docid) pairs.  The docids in
        these pairs are the docids of the items to remove. The ranks in the
        pairs are _hints_ of the ranks to remove: the hints are allowed to be
        overestimates (but not underestimates).  This will happen when the
        IndexerConnection tries to remove hits from queries which have already
        had some items removed.

        """
        raise NotImplementedError

    def flush(self):
        """Ensure that all changes made to the cache are written to disk.

        """
        raise NotImplementedError

    def close(self):
        """Close the storage for the cache (if this means anything).

        """
        raise NotImplementedError

class KeyValueStoreCacheManager(CacheManager, UserDict.DictMixin):
    """A manager that stores the cached items in chunks in a key-value store.

    Subclasses must implement:

     - __getitem__()
     - __setitem__()
     - __delitem__()
     - flush()
     - close()

    Subclasses may change the serialised representation by implementing:

     - encode()
     - decode()

    or, individually: 

     - encode_int()
     - decode_int()
     - encode_docids()
     - decode_docids()

    """
    encode = staticmethod(lambda x: cPickle.dumps(x, 2))
    decode = staticmethod(cPickle.loads)

    encode_int = encode
    decode_int = decode
    encode_docids = encode
    decode_docids = decode

    def __init__(self, chunksize=None):
        if chunksize is None:
            chunksize = 1000
        self.chunksize=chunksize
        CacheManager.__init__(self)

    def __getitem__(self, key):
        """Get the value for a given key.

        If there is no value for the key, return an empty string.

        This must be implemented by subclasses.

        """
        raise NotImplementedError

    def __setitem__(self, key, value):
        """Set the value for a given key.

        Replaces any existing value.

        This must be implemented by subclasses.

        """
        raise NotImplementedError

    def __delitem__(self, key):
        """Delete the value for a given key.

        If there is no value for the key, do nothing.

        This must be implemented by subclasses.

        """
        raise NotImplementedError

    def keys(self):
        """Iterate through all the keys stored.

        The keys can be returned in any order, as long as each is returned
        exactly once.

        This must be implemented by subclasses.

        """
        raise NotImplementedError

    def iter_queryids(self):
        # Currently, we don't allow sparse queryids, so we can just return an
        # iterator over the range of values.
        # Key 'I' holds the (encoded) value of the next ID to allocate.
        v = self['I']
        if len(v) == 0:
            return iter(())
        return xrange(self.decode_int(v))

    def get_queryid(self, query_str):
        v = self['Q' + query_str]
        if len(v) == 0:
            return None
        return self.decode_int(v)

    def get_or_make_queryid(self, query_str):
        v = self['Q' + query_str]
        if len(v) == 0:
            # Get the next ID to use.
            v = self['I']
            if len(v) == 0:
                thisid = 0
            else:
                thisid = self.decode_int(v)

            self['Q' + query_str] = self.encode_int(thisid)
            self['I'] = self.encode_int(thisid + 1)
            return thisid
        return self.decode_int(v)

    def make_hit_chunk_key(self, queryid, chunk):
        """Make the key for looking up a particular chunk for a given queryid.

        """
        return 'H%d:%d' % (queryid, chunk)

    def get_hits(self, queryid, startrank=0, endrank=None):
        startchunk = int(startrank // self.chunksize)
        if endrank is None:
            endchunk = None
        else:
            endchunk = int(endrank // self.chunksize)

        hits = []
        chunk = startchunk
        startrank_in_chunk = startrank - chunk * self.chunksize
        while endchunk is None or chunk <= endchunk:
            data = self[self.make_hit_chunk_key(queryid, chunk)]
            if len(data) == 0:
                # Chunk doesn't exist: implies that we're at the end
                break
            chunkhits = self.decode_docids(data)

            if endrank is None:
                # Use the whole of the chunk
                if startrank_in_chunk == 0:
                    # Avoid doing a slice if we don't have to.
                    hits.extend(chunkhits)
                else:
                    hits.extend(chunkhits[startrank_in_chunk:])
            else:
                endrank_in_chunk = min(len(chunkhits),
                                       endrank - chunk * self.chunksize)
                hits.extend(chunkhits[startrank_in_chunk:endrank_in_chunk])
            
            startrank_in_chunk = 0
            chunk += 1

        return hits

    def set_hits(self, queryid, docids):
        self.set_hits_internal(queryid, docids, 0)

    def set_hits_internal(self, queryid, docids, chunk):
        """Internal implementation of set_hits.

        Set the hits starting at the chunk number specified by `chunk`.

        """
        # Convert docids into a list - allows iterators to be supplied.
        # Could be done chunk-by-chunk, but it's easier to do this way, and
        # probably fast enough.
        docids = list(docids)

        chunkcount = int(len(docids) // self.chunksize) + 1 + chunk

        # Add the data in chunks
        offset = 0
        while chunk < chunkcount:
            data = self.encode_docids(docids[offset : offset + self.chunksize])
            offset += self.chunksize
            self[self.make_hit_chunk_key(queryid, chunk)] = data
            chunk += 1

        # Ensure that there aren't any chunks for further bits of data
        # remaining
        while True:
            key = self.make_hit_chunk_key(queryid, chunk)
            data = self[key]
            if len(data) == 0:
                break
            del self[key]
            chunk += 1

    def remove_hits(self, queryid, ranks_and_docids):
        if len(ranks_and_docids) == 0:
            return

        # To remove a hit, we need to get all chunks after the one containing
        # the given rank, and update them.

        # First, ensure the ranks are a sorted list - sort in descending order.
        ranks_and_docids = list(ranks_and_docids)
        ranks_and_docids.sort(reverse=True)

        startchunk = int(ranks_and_docids[-1][0] // self.chunksize)
        startrank = startchunk * self.chunksize
        hits = self.get_hits(queryid, startrank)

        unmatched = []
        for rank, docid in ranks_and_docids:
            if ((len(hits) > rank - startrank) and
                (hits[rank - startrank] == docid)):
                del hits[rank - startrank]
            else:
                unmatched.append((rank, docid))

        if len(unmatched) != 0:
            # Get all the hits, so that we can find the mismatched ones.
            # We might be able to get away without fetching all the hits in
            # some cases, but that would increase code complexity, so let's not
            # do it unless profiling shows it to be necessary.
            hits = self.get_hits(queryid, 0, startrank) + hits
            startrank = 0
            startchunk = 0

            for rank, docid in unmatched:
                rank = min(rank, len(hits) - 1)
                while rank >= 0:
                    if hits[rank] == docid:
                        del hits[rank]
                        break
                    rank -= 1
                assert rank != -1

        self.set_hits_internal(queryid, hits, startchunk)
