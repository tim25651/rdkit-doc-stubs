"""
rdkit.utils.listutils module¶
utility functions for lists
"""

def CompactListRepr(self, lst):
    """
    >>> CompactListRepr([0,1,1,1,1,0])
    '[0]+[1]*4+[0]'
    >>> CompactListRepr([0,1,1,2,1,1])
    '[0]+[1]*2+[2]+[1]*2'
    >>> CompactListRepr([])
    '[]'
    >>> CompactListRepr((0,1,1,1,1))
    '[0]+[1]*4'
    >>> CompactListRepr('foo')
    "['f']+['o']*2" """
    ...
