from setuptools import find_packages, setup

setup(
    name='chem_spectra',
    version='0.1.0',
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    install_requires=[  # 'nmrglue',
        'flask==1.0.2',
        'numpy',
        'scipy==1.2.0',
        'pymzml==2.0.6',
        'matplotlib==2.0.2',
        'gunicorn',
        'pytest==4.0.0',
        'coverage==4.5.3',
        'flake8==3.7.7'
    ],
)
