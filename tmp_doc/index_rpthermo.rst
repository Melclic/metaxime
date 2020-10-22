rpThermo's Documentation
========================

Indices and tables
##################

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Introduction
############

.. _rpBase: https://github.com/Galaxy-SynBioCAD/rpBase
.. _Equilibrator: http://equilibrator.weizmann.ac.il/
.. _Marvin: https://chemaxon.com/products/marvin/download
.. _MDF: http://equilibrator.weizmann.ac.il/static/classic_rxns/pathway.html#mdf
.. _SBtab: http://equilibrator.weizmann.ac.il/pathway/

Welcome to the documentation of the rpThermodynamics project. The different scripts extend the Equilibrator_ project to perform thermodynamics analyis of heterologous pathways. Furthermore, one can calculate the MDF of pathways and convert the rpSBML files to SBtab that is Equilibrator_ friendly.

Usage
#####

First build the rpBase_ docker before building the local docker using the follwing command:

.. code-block:: bash

   docker build -t brsynth/rpthermo-standalone:v2 -f Dockerfile .

rpThermo depends on the Marvin library and thus, you need to have a valid license.cxl that is in the same folder as the Dockerfile. Furthermore, you need to download the deb Marvin_ file (tested version is 20.9) and also place it in the same folder as the Dockerfile.

Then you can call the docker using the command line interface using the following script:

.. code-block:: bash

   python run_rpThermo.py -input /path/to/file_rpsbml.tar.gz -input_format tar -output /path/to/outfile.tar.gz

The same type of script may be called locally for calculating MDF_:

.. code-block:: bash

   python run_MDF.py -input /path/to/file_rpsbml.tar.gz -input_format tar -output /path/to/outfile.tar.gz
 
Lastly, the rpSBML pathway files can be converted SBtab_'s that are Equilibrator friendly:

.. code-block:: bash

   python run_eqSBtab.py -input /path/to/file_rpsbml.tar.gz -input_format tar -output /path/to/outfile.tar.gz

API
###

.. toctree::
   :maxdepth: 1
   :caption: Contents:

.. currentmodule:: rpToolServe

.. currentmodule:: rpEquilibrator

.. autoclass:: rpEquilibrator
    :show-inheritance:
    :members:
    :inherited-members:

.. currentmodule:: run_rpThermo

.. autoclass:: main
    :show-inheritance:
    :members:
    :inherited-members:

.. currentmodule:: run_MDF

.. autoclass:: main
    :show-inheritance:
    :members:
    :inherited-members:

.. currentmodule:: run_eqSBtab

.. autoclass:: main
    :show-inheritance:
    :members:
    :inherited-members:

