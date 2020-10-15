rpSBML's Documentation
======================

Indices and tables
##################

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Introduction
############

.. _libSBML: http://model.caltech.edu/software/libsbml/5.18.0/docs/formatted/python-api/
.. _RetroPath2.0: https://www.myexperiment.org/workflows/4987.html

Welcome to the documentation of rpSBML, a metabolic engeinnering extension of libSBML_ project to handle heterologous metabolic pathways. The project was created as a means of converting the results of RetroPath2.0_ to SBML files for easier handling and interoperability.

Usage
#####

Basic usage include opening a rpSBML file, extracting some of its information:

.. note::
   You can open any SBML file using rpSBML.py and manipulate it.

To open a SBML file:

.. code-block:: python

   rpsbml = rpSBML.rpSBML('rp_1_1', path='/path/to/your/file_rpsbml.xml')

To extract all the information of a pathway in the form of a dictionnary you can use:

.. code-block:: python

   json_friendly_pathway = rpsbml.genJSON()

To extract the species that are qualified as being central you can use:

.. code-block:: python

   groups = rpsbml.model.getPlugin('groups')
   central_species_group = groups.getGroup('central_species')
   central_species = [member.getIdRef() for member in central_species_group.getListOfMembers()]

API
###

.. toctree::
   :maxdepth: 1
   :caption: Contents:

.. currentmodule:: rpSBML

.. autoclass:: rpSBML
    :show-inheritance:
    :members:
    :inherited-members:
