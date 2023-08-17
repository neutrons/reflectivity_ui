=====================
How to Make a Release
=====================

.. contents::
    :local:


Overview
--------

The release of ``reflectivity_ui`` is configured to be done by an automated pipeline via GitHub action.
However, developers might have to create manual releases in case the automated system fails, the process of which requires a working local development.
For ``reflectivity_ui`` there are three official release channels, and this document will provide necessary information on how to release ``reflectivity_ui`` to all three.


Release to Conda
----------------

Release candidate and stable versions of ``reflectivity_ui`` are automatically released to the project channel, `neutrons`_
whenever a new tag is pushed to the repository.

To manually build a conda package (for testing puruposes only), the following steps are required:

* Make sure the local development environment has both ``anaconda-client`` and ``conda-build`` installed.
* Checkout the desired feature branch for testing.
* Build the package with ``conda build .``.
* Use ``conda build . --output`` to locate the built package.

.. _neutrons: https://anaconda.org/neutrons
