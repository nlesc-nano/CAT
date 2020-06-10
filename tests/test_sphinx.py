"""Test the :mod:`sphinx` documentation generation."""

from os.path import join
from sphinx.application import Sphinx

from nanoutils import delete_finally

SRCDIR = CONFDIR = 'docs'
OUTDIR = join('tests', 'test_files', 'build')
DOCTREEDIR = join('tests', 'test_files', 'build', 'doctrees')


@delete_finally(OUTDIR)
def test_sphinx_build() -> None:
    """Test :meth:`sphinx.application.Sphinx.build`."""
    app = Sphinx(SRCDIR, CONFDIR, OUTDIR, DOCTREEDIR, buildername='html', warningiserror=True)
    app.build(force_all=True)
