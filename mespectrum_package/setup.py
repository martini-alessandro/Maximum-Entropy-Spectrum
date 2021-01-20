from setuptools import setup
import os

def readme():
    with open('README.rst') as f:
        return f.read()

#setuptools.find_packages()

setup(
    name='mespectrum',
    version='1.0.0',
    author='Alessandro Martini, Walter Del Pozzo',
    author_email='martini.alessandr@gmail.com,walter.delpozzo@ligo.org',
    packages=['mespectrum'],
    package_dir = {'mespectrum':'./mespectrum'},
    url="https://github.com/martini-alessandro/Maximum-Entropy-Spectrum/",
    license='CC by 4.0',
    description='maximum entropy spectral analysis',
    long_description=readme(),
    include_package_data = True,
#    package_data={'mespectrum': get_package_data()},
    install_requires=[
        "numpy >= 1.16.4",
    ],
	long_description_content_type = 'text/x-rst'
)

