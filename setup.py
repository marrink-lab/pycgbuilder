#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

from pkg_resources import parse_requirements

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('requirements/prod.txt') as req:
    requirements = [str(r) for r in parse_requirements(req)]

setup(
    name='pycgbuilder',
    version='0.1.0',
    description="A Qt5 based GUI for making mapping files",
    long_description=readme,
    author="Peter C Kroon",
    author_email='p.c.kroon@rug.nl',
    url='https://github.com/pckroon/pycgbuilder',
    packages=[
        'pycgbuilder',
    ],
    package_dir={'pycgbuilder':
                 'pycgbuilder'},
    entry_points={
        "gui_scripts": [
            "pycgbuilder = pycgbuilder.__main__:main"
        ],
    },
    include_package_data=True,
    install_requires=requirements,
    license="Apache Software License",
    zip_safe=False,
    keywords='pycgbuilder',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
    ],
)
