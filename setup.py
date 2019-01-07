import sys
if sys.version_info < (3,):
    sys.exit('anndata requires Python >= 3.5')
from pathlib import Path

from setuptools import setup, find_packages
import versioneer

try:
    from anndata import __author__, __email__
except ImportError:  # deps not yet installed
    __author__ = __email__ = ''

setup(
    name='anndata',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='Annotated Data.',
    long_description=Path('README.rst').read_text('utf-8'),
    url='http://github.com/theislab/anndata',
    author=__author__,
    author_email=__email__,
    license='BSD-3-Clause',
    install_requires=[
        l.strip() for l in
        Path('requirements.txt').read_text('utf-8').splitlines()
    ],
    python_requires='>=3.5',
    packages=find_packages(),
    zip_safe=False,
    classifiers=[
        'Environment :: Console',
        'Framework :: Jupyter',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
)
