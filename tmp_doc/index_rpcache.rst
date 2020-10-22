rpCache's Documentation
=======================

Indices and tables
##################

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Introduction
############

.. _MetaNetX: https://www.metanetx.org/
.. _RetroRules: https://retrorules.org/
.. _Component-Contribution: https://github.com/eladnoor/component-contribution
.. _rpBase: https://github.com/Galaxy-SynBioCAD/rpBase

Welcome to rpCache's documentation. This project generates the cache that is used throughout many different projects. The files are automatically downloaded, from three sources MetaNetX_, RetroRules_ and Component-Contribution_ and are parsed to generate the cache. 

Before building the docker, you must build the rpBase_ docker before building it using the local docker:

.. code-block:: bash

   docker build -t brsynth/rpcache:v2 .

API
###

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. currentmodule:: rpCache

.. autoclass:: rpCache
    :show-inheritance:
    :members:
    :inherited-members:
