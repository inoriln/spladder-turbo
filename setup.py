#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md', encoding='utf-8') as readme_file:
    readme = readme_file.read()

with open('requirements.txt') as f:
    requirements = list(f.readlines())

setup(
    author="Yihao Zhang",
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    description="Unofficial, performance-optimized fork of SplAdder: detection and quantification of alternative splicing events from RNA-Seq data.",
    entry_points = {
        'console_scripts': [
            'spladder=spladder.spladder:main',
            'spladder-turbo=spladder.spladder:main',
        ],

    },
    install_requires=requirements,
    license="BSD-3-Clause",
    license_files=['LICENSE', 'COPYRIGHT', 'AUTHORS', 'ATTRIBUTION.md'],
    long_description=readme,
    long_description_content_type='text/markdown',
    include_package_data=True,
    keywords='spladder alternative-splicing rna-seq splicing-graph',
    name='spladder-turbo',
    packages=find_packages(exclude=['tests', 'tests.*']),
    project_urls={
        'Source': 'https://github.com/inoriln/spladder-turbo',
        'Issues': 'https://github.com/inoriln/spladder-turbo/issues',
        'Upstream (SplAdder)': 'https://github.com/ratschlab/spladder',
        'SplAdder documentation': 'https://spladder.readthedocs.io',
    },
    python_requires='>=3.9',
    url='https://github.com/inoriln/spladder-turbo',
    version='0.1.0',
    zip_safe=False,
)
