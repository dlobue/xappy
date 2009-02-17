#!/usr/bin/env python
#
# Copyright (C) 2007 Lemur Consulting Ltd
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
r"""errors.py: Exceptions for the search engine core.

"""
__docformat__ = "restructuredtext en"

class SearchEngineError(Exception):
    r"""Base class for exceptions thrown by the search engine.

    Any errors generated by xappy itself, or by xapian, will be instances of
    this class or its subclasses.

    """

class IndexerError(SearchEngineError):
    r"""Class used to report errors relating to the indexing API.

    """

class DuplicatedIdError(IndexerError):
    r"""Class used to report errors when a new document is added to an index
    with an ID which already exists.
    
    """

class SearchError(SearchEngineError):
    r"""Class used to report errors relating to the search API.

    """


class XapianError(SearchEngineError):
    r"""Base class for exceptions thrown by the xapian.

    Any errors generated by xapian will be instances of this class or its
    subclasses.

    """

def _rebase_xapian_exceptions():
    """Add new base classes for all the xapian exceptions.

    """
    import xapian
    for name in (
                 'AssertionError',
                 'DatabaseCorruptError',
                 'DatabaseCreateError',
                 'DatabaseError',
                 'DatabaseLockError',
                 'DatabaseModifiedError',
                 'DatabaseOpeningError',
                 'DatabaseVersionError',
                 'DocNotFoundError',
                 # We skip 'Error' because it inherits directly from exception
                 # and this causes problems with method resolution order.
                 # However, we probably don't need it anyway, because it's
                 # just a base class, and shouldn't ever actually be raised.
                 # Users can catch xappy.XapianError instead.
                 'FeatureUnavailableError',
                 'InternalError',
                 'InvalidArgumentError',
                 'InvalidOperationError',
                 'LogicError',
                 'NetworkError',
                 'NetworkTimeoutError',
                 'QueryParserError',
                 'RangeError',
                 'RuntimeError',
                 'UnimplementedError',
                 ):
        xapian_exception = getattr(xapian, name, None)
        if xapian_exception is not None:
            xapian_exception.__bases__ += (XapianError, )
            globals()['Xapian' + name] = xapian_exception

_rebase_xapian_exceptions()
