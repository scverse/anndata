from setuptools import setup, find_packages
from io import open
import versioneer

with open('requirements.txt', encoding='utf-8') as requirements:
    requires = [l.strip() for l in requirements]

with open('README.rst', encoding='utf-8') as readme_f:
    readme = readme_f.read()

setup(
    name='anndata',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='Annotated Data.',
    long_description=readme,
    url='http://github.com/theislab/anndata',
    author='Alex Wolf, Philipp Angerer, Sergei Rybakov',
    author_email='alex.wolf@helmholtz-muenchen.de',
    license='BSD-3-Clause',
    install_requires=requires,
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
