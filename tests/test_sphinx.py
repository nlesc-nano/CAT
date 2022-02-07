"""Test the :mod:`sphinx` documentation generation."""

import warnings
from os.path import join

import pytest
from nanoutils import delete_finally

try:
    from sphinx.application import Sphinx
    from sphinx.errors import SphinxWarning
except ImportError:
    HAS_SPHINX = False
else:
    HAS_SPHINX = True

SRCDIR = CONFDIR = 'docs'
OUTDIR = join('tests', 'test_files', 'build')
DOCTREEDIR = join('tests', 'test_files', 'build', 'doctrees')


@delete_finally(OUTDIR)
@pytest.mark.skipif(not HAS_SPHINX, reason="Requires Sphinx")
def test_sphinx_build() -> None:
    """Test :meth:`sphinx.application.Sphinx.build`."""
    try:
        app = Sphinx(SRCDIR, CONFDIR, OUTDIR, DOCTREEDIR, buildername='html', warningiserror=True)
        app.build(force_all=True)
    except SphinxWarning as ex:
        str_ex = str(ex)
        if "requests.exceptions.HTTPError" in str_ex:
            warnings.warn(str_ex)
            pytest.xfail(str_ex)
        else:
            raise
