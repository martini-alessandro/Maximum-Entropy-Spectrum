Installation
============

You have several options to install `memspectrum`!

## From the latest distributed version

The easiest way is to get the latest distributed version from the [PyPI](https://pypi.org/project/memspectrum/) archives.
To do this, you just need to type:

```
pip install memspectrum
```

This will install the package and its dependencies in your current working python environment.

## From source

If you want to get the latest developments, you can install the package from source.
To do this you need to clone the git [repository](https://github.com/martini-alessandro/Maximum-Entropy-Spectrum) and to build and install the code manually. A handy makefile will help you to do that.

These are the steps:


```Bash
git clone git@github.com:martini-alessandro/Maximum-Entropy-Spectrum.git
cd Maximum-Entropy-Spectrum
make install
```

This will build the package and will install it: `pip` will keep track of the dependencies.

## Build the docs

To get a local copy of this documentation:

```Bash
cd Maximum-Entropy-Spectrum
make docs
```

The documentation will be generated on `docs/__build`.



