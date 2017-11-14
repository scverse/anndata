from setuptools import setup, find_packages
from pathlib import Path
import versioneer

package_name = 'anndata'

req_path = Path('requires.txt')
if not req_path.is_file():
    req_path = Path(package_name + '.egg-info') / req_path
with req_path.open() as requirements:
    requires = [l.strip() for l in requirements]

with open('README.rst') as readme_f:
    readme = readme_f.read()

setup(
    name=package_name,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='An annotated data matrix.',
    long_description=readme,
    url='http://github.com/theislab/anndata',
    author='Philipp Angerer, Alex Wolf',
    author_email='alex.wolf@helmholtz-muenchen.de',
    license='BSD-3-Clause',
    install_requires=requires,
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
