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

## From the repo

You can get the latest changes (not yet released with pip) with:

```Bash
pip install git+https://github.com/martini-alessandro/Maximum-Entropy-Spectrum
```

This will install the package as is after the latest commit on the repo.

## From source

If you want to do some developments, you may want to install the package from source after cloning the git [repository](https://github.com/martini-alessandro/Maximum-Entropy-Spectrum).
A handy makefile will help you to install the code.

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



