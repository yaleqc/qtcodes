"""
Setup for qtcodes
"""
import os

from setuptools import setup, find_namespace_packages

REQUIREMENTS = [
    "numpy",
    "matplotlib>=3.3.0",
    "retworkx",
    "qiskit",
    "tqdm",
    "pydot",
    "pylatexenc",
    "IPython",
]

# Read long description from README.
README_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), "README.md")
with open(README_PATH) as readme_file:
    README = readme_file.read()

version_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "qtcodes", "VERSION.txt")
)

with open(version_path, "r") as fd:
    version_str = fd.read().rstrip()

setup(
    name="qtcodes",
    version=version_str,
    description="Qiskit Topological Codes",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/yaleqc/qtcodes",
    author="Shantanu Jha",
    author_email="shantanu.rajesh.jha@gmail.com",
    license="Apache 2.0",
    packages=find_namespace_packages(exclude=["test*"]),
    install_requires=REQUIREMENTS,
    classifiers=[
        "Environment :: Console",
        "License :: OSI Approved :: Apache Software License",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering",
    ],
    keywords="qiskit topological surface codes qec quantum error correction",
    python_requires=">=3.6",
    project_urls={
        "Bug Tracker": "https://github.com/yaleqc/qtcodes/issues",
        "Documentation": "https://github.com/yaleqc/qtcodes",
        "Source Code": "https://github.com/yaleqc/qtcodes",
        "Tutorials": "https://github.com/yaleqc/qtcodes/tutorials",
        "Tests": "https://github.com/yaleqc/qtcodes/test",
    },
    include_package_data=True,
)
