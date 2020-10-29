from setuptools import setup, find_packages

VERSION = "0.2.3"
DESCRIPTION = "a Python library to run highlevel functions of DGGRIDv7"
LONG_DESCRIPTION = "a set of python modules for creating and manipulating Discrete Global Grids with DGGRID version 7.0 software which was created and maintained by Kevin Sahr.  - more info to come!"

# Setting up
setup(
    # the name must match the folder name 'dggrid4py'
    name="dggrid4py",
    version=VERSION,
    author="Alexander Kmoch",
    author_email="<alexander.kmoch@ut.ee>",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    url='https://github.com/allixender/dggrid4py',
    license='AGPLv3',
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "fiona",
        "shapely",
        "geopandas"
    ],  # add any additional packages that
    # needs to be installed along with your package. Eg: 'caer'
    keywords=["python", "GIS", "DGGS", "DGGRID", "hexagons", "grids", "spatial statistics"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
        "Programming Language :: Python :: 3",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: GIS",
    ],
)
