===========
Gctf plugin
===========

This plugin provide wrappers around `Gctf <https://www2.mrc-lmb.cam.ac.uk/research/locally-developed-software/zhang-software/>`_ program.

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


+--------------+----------------+--------------------+
| prod: |prod| | devel: |devel| | support: |support| |
+--------------+----------------+--------------------+

.. |prod| image:: http://scipion-test.cnb.csic.es:9980/badges/gctf_prod.svg
.. |devel| image:: http://scipion-test.cnb.csic.es:9980/badges/gctf_devel.svg
.. |support| image:: http://scipion-test.cnb.csic.es:9980/badges/gctf_support.svg


Installation
------------

You will need to use `3.0 <https://github.com/I2PC/scipion/releases/tag/V3.0.0>`_ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

   .. code-block::

      scipion installp -p scipion-em-gctf

b) Developer's version

   * download repository

   .. code-block::

      git clone https://github.com/scipion-em/scipion-em-gctf.git

   * install

   .. code-block::

      scipion installp -p path_to_scipion-em-gctf --devel

Gctf binaries will be installed automatically with the plugin, but you can also link an existing installation. 
Default installation path assumed is ``software/em/gctf-1.18``, if you want to change it, set *GCTF_HOME* in ``scipion.conf`` file to the folder where the Gctf is installed. Depending on your CUDA version and GPU card compute capability you might want to change the default binary from ``Gctf_v1.18_sm30-75_cu10.1`` to a different one by explicitly setting *GCTF* variable. If you need to use CUDA different from the one used during Scipion installation (defined by CUDA_LIB), you can add *GCTF_CUDA_LIB* variable to the config file. Various binaries can be downloaded from the official Gctf website.

To check the installation, simply run one of the following Scipion tests: 

.. code-block::

   scipion test gctf.tests.test_protocols_gctf.TestGctfRefine
   scipion test gctf.tests.test_protocols_gctf.TestGctf

Supported versions
------------------

1.06 and 1.18.

1.18 is a special version designed for VPP data, it does not support local/movie CTF refinement or validation options. 

Protocols
---------

* ctf estimation
* ctf refinement

References
----------

1. Zhang K. (2016). Gctf: Real-time CTF determination and correction. JSB, 193: 1-12.
