import os
from setuptools import find_packages, setup

# determining the directory containing setup.py
setup_path = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(setup_path, 'README.rst'), encoding='utf-8') as f:
    readme = f.read()

setup(
    # package information
    name = 'baghera_tool',
    packages = find_packages(),
    version = '1.0.0',
    description = 'Bayesian Gene Heritability Analysis tool',
    long_description = readme,
    license = 'MIT',
    url='git@git.ecdf.ed.ac.uk:stracquadaniolab/projects/baghera_tool.git',
    keywords='',

    # author information
    author = 'Viola Fanfani',
    author_email = 'v.fanfani@sms.ed.ac.uk',

    # installation info and requirements
    install_requires=['argh',
                    'pandas',
                    'numpy==1.15',
                    'htseq',
                    'gseapy',
                    'pymc3>=3.6',
                    ],
    setup_requires=[
                    ],


    # test info and requirements
    test_suite='tests',
    tests_require=[],

    # package deployment info
    include_package_data=True,
    zip_safe=False,

    # all tools have cli interface
    entry_points={
        'console_scripts': [
            'baghera_tool=baghera_tool.cli:main',
        ],
    },
)
