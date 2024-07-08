===========
Gctf plugin
===========

This plugin provides a wrapper for `Gctf <https://www2.mrc-lmb.cam.ac.uk/research/locally-developed-software/zhang-software/>`_ program.

.. image:: https://img.shields.io/pypi/v/scipion-em-gctf.svg
        :target: https://pypi.python.org/pypi/scipion-em-gctf
        :alt: PyPI release

.. image:: https://img.shields.io/pypi/l/scipion-em-gctf.svg
        :target: https://pypi.python.org/pypi/scipion-em-gctf
        :alt: License

.. image:: https://img.shields.io/pypi/pyversions/scipion-em-gctf.svg
        :target: https://pypi.python.org/pypi/scipion-em-gctf
        :alt: Supported Python versions

.. image:: https://img.shields.io/sonar/quality_gate/scipion-em_scipion-em-gctf?server=https%3A%2F%2Fsonarcloud.io
        :target: https://sonarcloud.io/dashboard?id=scipion-em_scipion-em-gctf
        :alt: SonarCloud quality gate

.. image:: https://img.shields.io/pypi/dm/scipion-em-gctf
        :target: https://pypi.python.org/pypi/scipion-em-gctf
        :alt: Downloads

Installation
------------

You will need to use 3.0+ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

   .. code-block::

      scipion installp -p scipion-em-gctf

b) Developer's version

   * download repository

   .. code-block::

      git clone -b devel https://github.com/scipion-em/scipion-em-gctf.git

   * install

   .. code-block::

      scipion installp -p /path/to/scipion-em-gctf --devel

Configuration variables
-----------------------
*CONDA_ACTIVATION_CMD*: If undefined, it will rely on conda command being in the
PATH (not recommended), which can lead to execution problems mixing scipion
python with conda ones. One example of this could can be seen below but
depending on your conda version and shell you will need something different:
CONDA_ACTIVATION_CMD = eval "$(/extra/miniconda3/bin/conda shell.bash hook)"

*GCTF_ENV_ACTIVATION* (default = conda activate gctf):
Command to activate the Gctf environment. It will have cudatoolkit=10.1 installed.

*GCTF_HOME* (default = software/em/gctf-1.18):
Path to Gctf installation folder.

*GCTF_BIN* (default = Gctf_v1.18_sm30-75_cu10.1):
Binary to use.

Verifying
---------
To check the installation, simply run one of the following Scipion tests: 

.. code-block::

   scipion test gctf.tests.test_protocols_gctf.TestGctf
   scipion tests gctf.tests.test_protocols_gctf_ts.TestGctfTs

Supported versions
------------------

1.18

Protocols
---------

* ctf estimation
* tilt-series gctf

References
----------

1. Zhang K. (2016). Gctf: Real-time CTF determination and correction. JSB, 193: 1-12.
