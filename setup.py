import os
from setuptools import find_packages, setup

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='chem-spectra-app',
    version='0.9.0',
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    install_requires=required,
)
