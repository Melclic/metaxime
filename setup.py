import os
from setuptools import setup, find_packages

major_version = 0
minor_version = 0
patch_version = 1

if os.environ.get("PATCH_NUMBER"):
    patch_version = os.environ.get("PATCH_NUMBER")

version = str(major_version) + "." + str(minor_version) + "." + str(patch_version)

here = os.path.abspath(os.path.dirname(__file__))
README = open(os.path.join(here, "README.md")).read()

third_party_require = [
    "rdkit>=2025.03.3",
    "pandas>=2.3.0",
    "numpy>=2.2.6",
    "compress_json>=1.1.1",
    "networkx>=3.4.2",
    "biopathopt"
]

tests_require = [
    "nosexcover",
    "coverage",
    "nose-timer",
    "nose-xunitmp",
    "pylint",
    "nose",
]

require = third_party_require + tests_require


setup(
    name="metaxime",
    version=version,
    description="A Python package for parsing and completing retrosynthesis outputs, performing FBA simulations, and exploring metabolic pathway features to predict the production potential of target compounds.",
    long_description=README,
    classifiers=[
        "Programming Language :: Python :: 3.10",
    ],  # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
    keywords="",
    author="Melchior du Lac",
    author_email="",
    license="",
    packages=find_packages(exclude=["ez_setup"]),
    package_dir={},
    package_data={},
    include_package_data=True,
    zip_safe=False,
    test_suite="nose.collector",
    install_requires=require,
    tests_require=tests_require,
)
