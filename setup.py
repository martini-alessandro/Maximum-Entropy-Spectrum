from setuptools import setup
import os
import sys

try:
	from sphinx.setup_command import BuildDoc
	cmdclass = {'build_sphinx': BuildDoc} #to build with sphinx
except ImportError:
	if sys.argv[1] == 'build_sphinx': warnings.warn("sphinx module not found: impossibile to build the documents")
	pass



def readme():
    with open('README.md') as f:
        return f.read()

#setuptools.find_packages()

setup(
    name='memspectrum',
    version='1.3.0',
    author='Alessandro Martini',
    author_email='martini.alessandr@gmail.com',
    packages=['memspectrum'],
    package_dir = {'memspectrum':'memspectrum'},
    url="https://github.com/martini-alessandro/Maximum-Entropy-Spectrum/",
    license='CC by 4.0',
    description='maximum entropy spectral analysis',
    long_description=readme(),
    long_description_content_type='text/markdown',
	#long_description_content_type = 'text/x-md',
    include_package_data = True,
#    package_data={'memspectrum': get_package_data()},
    install_requires=[
        "numpy",
        "scipy",
    ],
	command_options={
        'build_sphinx': {
            'source_dir': ('setup.py', 'docs'),
            'build_dir': ('setup.py', 'docs/__build'),
            }},
)

