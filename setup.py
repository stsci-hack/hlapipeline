#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""
import os
import pkgutil
import shutil
import sys
import importlib
from setuptools import setup, find_packages

requirements = ['Click>=6.0', 'astropy', 'scipy', 'matplotlib', \
                'photutils', 'pyyaml', 'astroquery',\
                'scipy', 'sphinx','sphinx_rtd_theme', \
                'stsci_rtd_theme','stsci.tools',\
                'stwcs','setuptools']

setup_requirements = ['pytest-runner', 'sphinx']

test_requirements = ['pytest', 'requests_mock', 'ci_watson']


if not pkgutil.find_loader('relic'):
    relic_local = os.path.exists('relic')
    relic_submodule = (relic_local and
                       os.path.exists('.gitmodules') and
                       not os.listdir('relic'))
    try:
        if relic_submodule:
            check_call(['git', 'submodule', 'update', '--init', '--recursive'])
        elif not relic_local:
            check_call(['git', 'clone', 'https://github.com/spacetelescope/relic.git'])

        sys.path.insert(1, 'relic')
    except CalledProcessError as e:
        print(e)
        exit(1)

import relic.release

PACKAGENAME = 'hlapipeline'

# Due to overriding `install` and `build_sphinx` we need to download
# setup_requires dependencies before reaching `setup()`. This allows
# `sphinx` to exist before the `BuildSphinx` class is injected.
_install_setup_requires(dict(setup_requires=setup_requirements))

for dep_pkg in setup_requirements:
    try:
        importlib.import_module(dep_pkg)
    except ImportError:
        print("{0} is required in order to install '{1}'.\n"
              "Please install {0} first.".format(dep_pkg, PACKAGENAME),
              file=sys.stderr)
        exit(1)
from sphinx.cmd.build import build_main
from sphinx.setup_command import BuildDoc

version = relic.release.get_info()
relic.release.write_template(version, PACKAGENAME)

class BuildSphinx(BuildDoc):
    """Build Sphinx documentation after compiling C extensions"""

    description = 'Build Sphinx documentation'

    def initialize_options(self):
        BuildDoc.initialize_options(self)

    def finalize_options(self):
        BuildDoc.finalize_options(self)

    def run(self):
        build_cmd = self.reinitialize_command('build_ext')
        build_cmd.inplace = 1
        self.run_command('build_ext')
        build_main(['-b', 'html', 'doc/source', 'build/sphinx/html'])

        # Bundle documentation inside of drizzlepac
        if os.path.exists(docs_compiled_src):
            if os.path.exists(docs_compiled_dest):
                shutil.rmtree(docs_compiled_dest)

            shutil.copytree(docs_compiled_src, docs_compiled_dest)

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

setup(
    author="Warren J. Hack",
    author_email='hack@stsci.edu',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    cmdclass={
        'build_sphinx': BuildSphinx,
    },
    description="Code for implementing HLA-type processing in the HST Pipeline",
    entry_points={
        'console_scripts': [
            'hlapipeline=hlapipeline.cli:main',
        ],
    },
    install_requires=requirements,
    license="BSD license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='hlapipeline',
    name='hlapipeline',
    packages=find_packages(),
    project_urls={
        'Bug Reports': 'https://github.com/spacetelescope/hlapipeline/issues/',
        'Source': 'https://github.com/spacetelescope/hlapipeline/',
        'Help': 'https://hsthelp.stsci.edu/',
    },
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/spacetelescope/hlapipeline',
    version='0.1.0',
    zip_safe=False,
)
