from setuptools import setup
import os

def readme():
    with open('README.rst') as f:
        return f.read()

#setuptools.find_packages()

setup(
    name='mesa',
    version='1.0.0',
    author='Alessandro Martini, Walter Del Pozzo',
    author_email='martini.alessandr@gmail.com,walter.delpozzo@ligo.org',
    packages=['mesa'],
    package_dir = {'mesa':'./mesa'},
    url="https://github.com/martini-alessandro/Maximum-Entropy-Spectrum/",
    license='CC by 4.0',
    description='maximum entropy spectral analysis',
    long_description=readme(),
    include_package_data = True,
#    package_data={'mesa': get_package_data()},
    install_requires=[
        "numpy >= 1.16.4",
		"scipy >= 1.3.1",
		#"lalsuite >= 6.62"
    ],
	long_description_content_type = 'text/x-rst'
)
