import os
from setuptools import find_packages, setup

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='chem-spectra-app',
    version='0.10.17-patch.1',
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    install_requires=required,
)
