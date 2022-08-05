from setuptools import setup
import json

with open("./supy/supy_version.json") as f:
    dict_ver = json.load(f)

print(dict_ver)
__version__ = f"{dict_ver['ver_milestone']}.{dict_ver['ver_major']}.{dict_ver['ver_minor']}{dict_ver['ver_remark']}"


def readme():
    try:
        with open("../README.md", encoding="utf-8") as f:
            return f.read()
    except:
        return f"SuPy package"


setup(
    name="supy",
    version=__version__,
    description="the SUEWS model that speaks python",
    long_description=readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/UMEP-Dev/SuPy",
    author=", ".join(
        [
            "Dr Ting Sun",
            "Dr Hamidreza Omidvar",
            "Prof Sue Grimmond",
        ]
    ),
    author_email=", ".join(
        [
            "ting.sun@reading.ac.uk",
            "h.omidvar@reading.ac.uk",
            "c.s.grimmond@reading.ac.uk",
        ]
    ),
    license="GPL-V3.0",
    packages=["supy"],
    package_data={
        "supy": [
            "sample_run/*",
            "sample_run/Input/*",
            "*.json",
            "util/*",
            "cmd/*",
        ]
    },
    # distclass=BinaryDistribution,
    ext_modules=[],
    install_requires=[
        "pandas>=1.3",
        "matplotlib",
        "scipy",
        "dask",  # needs dask for parallel tasks
        "f90nml",  # utility for namelist files
        "seaborn",  # stat plotting
        "atmosp",  # my own `atmosp` module forked from `atmos-python`
        "cdsapi",  # ERA5 data
        "xarray",  # utility for high-dimensional datasets
        "multiprocess",  # a better multiprocessing library
        "click",  # cmd tool
        "lmfit",  # optimiser
        "numdifftools",  # required by `lmfit` for uncertainty estimation
        "pvlib",  # TMY-related solar radiation calculations
        "platypus-opt==1.0.4",  # a multi-objective optimiser
        "supy_driver==2021a9",  # a separate f2py-based driver
    ],
    extras_require={
        "hdf": [
            "tables",  # for dumping in hdf5
        ]
    },
    entry_points={
        #   command line tools
        "console_scripts": [
            "suews-run=supy.cmd.SUEWS:SUEWS",
            "suews-convert=supy.cmd.table_converter:convert_table_cmd",
        ]
    },
    include_package_data=True,
    python_requires="~=3.7",
    classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
    ],
    zip_safe=False,
)
