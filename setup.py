#!/usr/bin/env python
# encoding: utf-8

from setuptools import setup, find_packages

with open("README.md", "r") as readme_file:
    readme = readme_file.read()

requirements = ["astropy", "numpy"]

setup(
    name="pysrcor",
    version="1.0",
    author="Guang Yang",
    author_email="gyang206265@gmail.com",
    description="A quick and easy function to correlate source positions from two catalogs",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/Guang91/pysrcor/",
    packages=find_packages(),
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
    ],
)
