from setuptools import setup
import os

def readme():
    with open('README.rst') as f:
        return f.read()

#setuptools.find_packages()

setup(
    name='memspectrum',
    version='1.1.0.post',
    author='Alessandro Martini',
    author_email='martini.alessandr@gmail.com',
    packages=['memspectrum'],
    package_dir = {'memspectrum':'./memspectrum'},
    url="https://github.com/martini-alessandro/Maximum-Entropy-Spectrum/",
    license='CC by 4.0',
    description='maximum entropy spectral analysis',
    long_description=readme(),
    include_package_data = True,
#    package_data={'memspectrum': get_package_data()},
    install_requires=[
        "numpy >= 1.16.4",
    ],
	long_description_content_type = 'text/x-rst'
)

