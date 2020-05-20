"""Tests for :class:`CAT._setattr.SetAttr`."""

import copy
import reprlib
from assertionlib import assertion

from CAT.utils import SetAttr


class _Test:
    a = True


OBJ = SetAttr(_Test, 'a', False)


def test_setattr() -> None:
    """Test :class:`CAT._setattr.SetAttr`."""
    assertion.is_(OBJ.obj, _Test)
    assertion.eq(OBJ.name, 'a')
    assertion.is_(OBJ.value, False)

    assertion.is_(OBJ.attr, _Test.a)
    try:
        OBJ.attr = False
        assertion.is_(OBJ.attr, False)
    finally:
        OBJ.attr = True

    assertion.contains(repr(OBJ), SetAttr.__name__)
    assertion.contains(repr(OBJ), object.__repr__(OBJ.obj))
    assertion.contains(repr(OBJ), reprlib.repr(OBJ.value))
    assertion.contains(repr(OBJ), 'a')

    obj2 = SetAttr(_Test, 'a', False)
    obj3 = SetAttr(_Test, 'a', True)
    assertion.eq(OBJ, obj2)
    assertion.ne(OBJ, obj3)
    assertion.ne(OBJ, 0)

    assertion.assert_(OBJ.__reduce__, exception=TypeError)
    assertion.is_(copy.copy(OBJ), OBJ)
    assertion.is_(copy.deepcopy(OBJ), OBJ)
    assertion.truth(hash(OBJ))
    assertion.is_(hash(OBJ), OBJ._hash)

    with OBJ:
        assertion.is_(_Test.a, False)
    assertion.is_(_Test.a, True)
