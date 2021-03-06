# Copyright (C) 2007 Lemur Consulting Ltd
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
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
