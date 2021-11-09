"""Test the :mod:`sphinx` documentation generation."""

from os.path import join

import pytest
from nanoutils import delete_finally

try:
    from sphinx.application import Sphinx
except ImportError:
    HAS_SPHINX = False
else:
    HAS_SPHINX = True

SRCDIR = CONFDIR = 'docs'
OUTDIR = join('tests', 'test_files', 'build')
DOCTREEDIR = join('tests', 'test_files', 'build', 'doctrees')


@delete_finally(OUTDIR)
@pytest.mark.xfail  # TODO: Sort out the recent sphinx failure
@pytest.mark.skipif(not HAS_SPHINX, reason="Requires Sphinx")
def test_sphinx_build() -> None:
    """Test :meth:`sphinx.application.Sphinx.build`."""
    app = Sphinx(SRCDIR, CONFDIR, OUTDIR, DOCTREEDIR, buildername='html', warningiserror=True)
    app.build(force_all=True)
