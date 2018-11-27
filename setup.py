from setuptools import find_packages, setup

setup(
    name='chem_spectra',
    version='0.1.0',
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    install_requires=[ # 'nmrglue',
        'flask',
        'numpy',
        'scipy',
        'waitress',
    ],
)
