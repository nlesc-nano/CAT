"""Test the :mod:`sphinx` documentation generation."""

import sys
from os.path import join

import pytest
from nanoutils import delete_finally

if sys.version_info >= (3, 7):
    from sphinx.application import Sphinx
    GE37 = True
else:
    GE37 = False

SRCDIR = CONFDIR = 'docs'
OUTDIR = join('tests', 'test_files', 'build')
DOCTREEDIR = join('tests', 'test_files', 'build', 'doctrees')


@delete_finally(OUTDIR)
@pytest.mark.skipif(not GE37, reason="Requires python >=3.7")
@pytest.mark.xfail  # TODO: Sort out the recent sphinx failure
def test_sphinx_build() -> None:
    """Test :meth:`sphinx.application.Sphinx.build`."""
    app = Sphinx(SRCDIR, CONFDIR, OUTDIR, DOCTREEDIR, buildername='html', warningiserror=True)
    app.build(force_all=True)
