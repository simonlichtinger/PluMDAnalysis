#
# PluMDAnalysis setuptools script
#
from setuptools import setup, find_packages


def get_version():
    import os
    import sys

    sys.path.append(os.path.abspath('PluMDAnalysis'))
    from PluMDAnalysis.version_info import VERSION as version
    sys.path.pop()

    return version


def get_readme():
    """
    Load README.md text for use as description.
    """
    with open('README.md') as f:
        return f.read()


# Go!
setup(
    # Module name (lowercase)
    name='plumdanalysis',

    # Version
    version=get_version(),

    description='A PLUMED input file generator',

    long_description=get_readme(),

    license='MIT license',

    # author='',

    # author_email='',

    maintainer='Simon Lichtinger',

    maintainer_email='simon.lichtinger@sjc.ox.ac.uk',

    url='https://github.com/simonlichtinger/PluMDAnalysis',

    # Packages to include
    packages=find_packages(include=('PluMDAnalysis', 'PluMDAnalysis.*')),

    # List of dependencies
    install_requires=[
        # Dependencies go here!
        'numpy',
        'MDAnalysis'
    ],
    extras_require={
        'docs': [
            # Sphinx for doc generation. Version 1.7.3 has a bug:
            'sphinx>=1.5, !=1.7.3',
            # Nice theme for docs
            'sphinx_rtd_theme',
        ],
        'dev': [
            # Flake8 for code style checking
            #'flake8>=3',
            'pytest',
        ],
    },
)