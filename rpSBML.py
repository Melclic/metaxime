"""rpSBML
.. moduleauthor:: Melchior du Lac
"""


import libsbml
from hashlib import md5
import os
import logging
import copy


"""
TODO

1) fonction qui retourne l’identifiant du target “target_id"
2) fonction qui retourne directement la liste des identifiants des précurseurs “precursor_ids”.

"""

## @package RetroPath SBML writer
# Documentation for SBML representation of the different model
#
# To exchange between the different workflow nodes, the SBML (XML) format is used. This
# implies using the libSBML library to create the standard definitions of species, reactions, etc...
# Here we also define our own annotations that are used internally in that we call BRSYNTH nodes.
# The object holds an SBML object and a series of methods to write and access BRSYNTH related annotations

##################################################################
############################### rpSBML ###########################
##################################################################


class rpSBML:
    """This class uses the libSBML object and handles it by adding BRSynth annotation
    """
    def __init__(self, modelName, document=None, path=None):
        """Constructor for the rpSBML class

        Note that the user can pass either a document libSBML object or a path to a SBML file. If a path is passed it overwrite the passed document object.

        :param modelName: The Name of the model
        :param document: The libSBML document class (Default: None)
        :param path: The path of a SBML file (Default: None)

        :type modelName: str
        :type path: str
        :type document: libsbml.SBMLDocument
        """
        self.logger = logging.getLogger(__name__)
        #WARNING: change this to reflect the different debugging levels
        self.logger.debug('Started instance of rpSBML')
        #self.logger.setLevel(logging.INFO)
        self.modelName = modelName
        self.document = document
        if self.document==None:
            self.model = None
        else:
            self.model = self.document.getModel()
        if not path==None:
            self.path = path
            self.readSBML(path)
        else:
            self.path = None
        self.miriam_header = {'compartment': {'mnx': 'metanetx.compartment/', 'bigg': 'bigg.compartment/', 'seed': 'seed/', 'name': 'name/'}, 'reaction': {'mnx': 'metanetx.reaction/', 'rhea': 'rhea/', 'reactome': 'reactome/', 'bigg': 'bigg.reaction/', 'sabiork': 'sabiork.reaction/', 'ec': 'ec-code/', 'biocyc': 'biocyc/', 'lipidmaps': 'lipidmaps/', 'uniprot': 'uniprot/'}, 'species': {'inchikey': 'inchikey/', 'pubchem': 'pubchem.compound/','mnx': 'metanetx.chemical/', 'chebi': 'chebi/CHEBI:', 'bigg': 'bigg.metabolite/', 'hmdb': 'hmdb/', 'kegg_c': 'kegg.compound/', 'kegg_d': 'kegg.drug/', 'biocyc': 'biocyc/META:', 'seed': 'seed.compound/', 'metacyc': 'metacyc.compound/', 'sabiork': 'sabiork.compound/', 'reactome': 'reactome/R-ALL-'}}
        self.header_miriam = {'compartment': {'metanetx.compartment': 'mnx', 'bigg.compartment': 'bigg', 'seed': 'seed', 'name': 'name'}, 'reaction': {'metanetx.reaction': 'mnx', 'rhea': 'rhea', 'reactome': 'reactome', 'bigg.reaction': 'bigg', 'sabiork.reaction': 'sabiork', 'ec-code': 'ec', 'biocyc': 'biocyc', 'lipidmaps': 'lipidmaps', 'uniprot': 'uniprot'}, 'species': {'inchikey': 'inchikey', 'pubchem.compound': 'pubchem', 'metanetx.chemical': 'mnx', 'chebi': 'chebi', 'bigg.metabolite': 'bigg', 'hmdb': 'hmdb', 'kegg.compound': 'kegg_c', 'kegg.drug': 'kegg_d', 'biocyc': 'biocyc', 'seed.compound': 'seed', 'metacyc.compound': 'metacyc', 'sabiork.compound': 'sabiork', 'reactome': 'reactome'}}


    #######################################################################
    ############################# PRIVATE FUNCTIONS ####################### 
    #######################################################################


    def _checklibSBML(self, value, message):
        """Private function that checks the libSBML calls.

        Check that the libSBML python calls do not return error INT and if so, display the error. Taken from: http://sbml.org/Software/libSBML/docs/python-api/create_simple_model_8py-example.html

        :param value: The libSBML command returned int
        :param message: The string that describes the call
        :type value: int
        :type message: str

        :raises AttributeError: If the libSBML command encounters an error or the input value is None

        :return: None
        :rtype: None
        """
        if value is None:
            self.logger.error('LibSBML returned a null value trying to ' + message + '.')
            raise AttributeError
        elif type(value) is int:
            if value==libsbml.LIBSBML_OPERATION_SUCCESS:
                return
            else:
                err_msg = 'Error encountered trying to ' + message + '.' \
                        + 'LibSBML returned error code ' + str(value) + ': "' \
                        + libsbml.OperationReturnValue_toString(value).strip() + '"'
                self.logger.error(err_msg)
                raise AttributeError
        else:
            #self.logger.debug(message)
            return None


    def _nameToSbmlId(self, name):
        """String to SBML id's
        
        Convert any String to one that is compatible with the SBML meta_id formatting requirements
        
        :param name: The input string
        :type name: str

        :return: SBML valid string
        :rtype: str
        """
        IdStream = []
        count = 0
        end = len(name)
        if '0' <= name[count] and name[count] <= '9':
            IdStream.append('_')
        for count in range(0, end):
            if (('0' <= name[count] and name[count] <= '9') or
                    ('a' <= name[count] and name[count] <= 'z') or
                    ('A' <= name[count] and name[count] <= 'Z')):
                IdStream.append(name[count])
            else:
                IdStream.append('_')
        Id = ''.join(IdStream)
        if Id[len(Id) - 1] != '_':
            return Id
        return Id[:-1]


    def _genMetaID(self, name):
        """String to hashed id

        Hash an input string and then pass it to _nameToSbmlId()

        :param name: Input string
        :type name: str

        :return: Hashed string id
        :rtype: str
        """
        return self._nameToSbmlId(md5(str(name).encode('utf-8')).hexdigest())


    def _compareXref(self, current, toadd):
        """Compare two dictionaries of lists that describe the cross-reference and return the difference

        :param current: The source cross-reference dictionary
        :param toadd: The target cross-reference dictionary
        :type current: dict
        :type toadd: dict

        :return: Difference between the two cross-reference dictionaries
        :rtype: dict
        """
        toadd = copy.deepcopy(toadd)
        for database_id in current:
            try:
                list_diff = [i for i in toadd[database_id] if i not in current[database_id]]
                if not list_diff:
                    toadd.pop(database_id)
                else:
                    toadd[database_id] = list_diff
            except KeyError:
                pass
        return toadd


    ######################################################################
    ####################### Annotations ##################################
    ######################################################################


    def _defaultBothAnnot(self, meta_id):
        """Returns a default annotation string that include MIRIAM and BRSynth annotation

        :param meta_id: The meta ID to be added to the default annotation
        :type meta_id: str

        :return: The default annotation string
        :rtype: str
        """
        return '''<annotation>
  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
    <rdf:Description rdf:about="#'''+str(meta_id or '')+'''">
      <bqbiol:is>
        <rdf:Bag>
        </rdf:Bag>
      </bqbiol:is>
    </rdf:Description>
    <rdf:BRSynth rdf:about="#'''+str(meta_id or '')+'''">
      <brsynth:brsynth xmlns:brsynth="http://brsynth.eu">
      </brsynth:brsynth>
    </rdf:BRSynth>
  </rdf:RDF>
</annotation>'''


    def _defaultBRSynthAnnot(self, meta_id):
        """Returns BRSynth default annotation string

        :param meta_id: The meta ID to be added to the annotation string
        :type meta_id: str

        :return: The default annotation string
        :rtype: str
        """
        return '''<annotation>
  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
    <rdf:BRSynth rdf:about="#'''+str(meta_id or '')+'''">
      <brsynth:brsynth xmlns:brsynth="http://brsynth.eu">
      </brsynth:brsynth>
    </rdf:BRSynth>
  </rdf:RDF>
</annotation>'''


    def _defaultMIRIAMAnnot(self, meta_id):
        """Returns MIRIAM default annotation string

        :param meta_id: The meta ID to be added to the annotation string
        :type meta_id: str

        :return: The default annotation string
        :rtype: str
        """
        return '''<annotation>
  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
    <rdf:Description rdf:about="#'''+str(meta_id or '')+'''">
      <bqbiol:is>
        <rdf:Bag>
        </rdf:Bag>
      </bqbiol:is>
    </rdf:Description>
  </rdf:RDF>
</annotation>'''



    def addUpdateBRSynth(self, sbase_obj, annot_header, value, units=None, isAlone=False, isList=False, isSort=True, meta_id=None):
        """Append or update an entry to the BRSynth annotation of the passed libsbml.SBase object.
        
        If the annot_header isn't contained in the annotation it is created. If it already exists it overwrites it

        :param sbase_obj: The libSBML object to add the different 
        :param annot_header: The annotation header that defines the type of entry
        :param value: The value(s) to add
        :param units: Add a values unit to the entry  
        :param isAlone: Add the entry without any units or defined within a value child (Setting this to True will ignore any units)
        :param isList: Define if the value entry is a list or not
        :param isSort: Sort the list that is passed (Only if the isList is True)
        :param meta_id: The meta ID to be added to the annotation string
        
        :type units: str
        :type isAlone: bool
        :type isList: bool
        :type isSort: bool

        :type sbase_obj: libsbml.SBase
        :type annot_header: str
        :type value: Union[str, int, float, list]
        :type meta_id: str

        :rtype: bool
        :return: Sucess or failure of the function
        """
        self.logger.debug('############### '+str(annot_header)+' ################')
        if isList:
            annotation = '''<annotation>
      <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">
        <rdf:BRSynth rdf:about="#adding">
          <brsynth:brsynth xmlns:brsynth="http://brsynth.eu">
            <brsynth:'''+str(annot_header)+'''>'''
            if isSort:
                for name in sorted(value, key=value.get, reverse=True):
                    if isAlone:
                        annotation += '<brsynth:'+str(name)+'>'+str(value[name])+'</brsynth:'+str(name)+'>'
                    else:
                        if units:
                            annotation += '<brsynth:'+str(name)+' units="'+str(units)+'" value="'+str(value[name])+'" />'
                        else:
                            annotation += '<brsynth:'+str(name)+' value="'+str(value[name])+'" />'
            else:
                for name in value:
                    if isAlone:
                        annotation += '<brsynth:'+str(name)+'>'+str(value[name])+'</brsynth:'+str(name)+'>'
                    else:
                        if units:
                            annotation += '<brsynth:'+str(name)+' units="'+str(units)+'" value="'+str(value[name])+'" />'
                        else:
                            annotation += '<brsynth:'+str(name)+' value="'+str(value[name])+'" />'
            annotation += '''
            </brsynth:'''+str(annot_header)+'''>
          </brsynth:brsynth>
        </rdf:BRSynth>
      </rdf:RDF>
    </annotation>'''
        else:
            #### create the string
            annotation = '''<annotation>
      <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">
        <rdf:BRSynth rdf:about="#adding">
          <brsynth:brsynth xmlns:brsynth="http://brsynth.eu">'''
            if isAlone:
                annotation += '<brsynth:'+str(annot_header)+'>'+str(value)+'</brsynth:'+str(annot_header)+'>'
            else:
                if units:
                    annotation += '<brsynth:'+str(annot_header)+' units="'+str(units)+'" value="'+str(value)+'" />'
                else:
                    annotation += '<brsynth:'+str(annot_header)+' value="'+str(value)+'" />'
            annotation += '''
          </brsynth:brsynth>
        </rdf:BRSynth>
      </rdf:RDF>
    </annotation>'''
        annot_obj = libsbml.XMLNode.convertStringToXMLNode(annotation)
        if annot_obj==None:
            self.logger.error('Cannot conver this string to annotation object: '+str(annotation))
            return False
        #### retreive the annotation object
        brsynth_annot = None
        obj_annot = sbase_obj.getAnnotation()
        if not obj_annot:
            sbase_obj.setAnnotation(libsbml.XMLNode.convertStringToXMLNode(self._defaultBRSynthAnnot(meta_id)))
            obj_annot = sbase_obj.getAnnotation()
            if not obj_annot:
                self.logger.error('Cannot update BRSynth annotation')
                return False
        brsynth_annot = obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
        if not brsynth_annot:
             self.logger.error('Cannot find the BRSynth annotation')
             return False
        #add the annotation and replace if it exists
        isfound_target = False
        #self.logger.debug(brsynth_annot.toXMLString())
        for i in range(brsynth_annot.getNumChildren()):
            self.logger.debug(annot_header+' -- '+str(brsynth_annot.getChild(i).getName()))
            if annot_header==brsynth_annot.getChild(i).getName():
                isfound_target = True
                '''
                self._checklibSBML(brsynth_annot.removeChild(brsynth_annot.getIndex(i)),
                    'Removing annotation '+str(annot_header))
                '''
                self._checklibSBML(brsynth_annot.removeChild(i), 'Removing annotation '+str(annot_header))
                isfound_source = False
                source_brsynth_annot = annot_obj.getChild('RDF').getChild('BRSynth').getChild('brsynth')
                for y in range(source_brsynth_annot.getNumChildren()):
                    self.logger.debug('\t'+annot_header+' -- '+str(source_brsynth_annot.getChild(y).getName()))
                    if str(annot_header)==str(source_brsynth_annot.getChild(y).getName()):
                        isfound_source = True
                        self.logger.debug('Adding annotation to the brsynth annotation: '+str(source_brsynth_annot.getChild(y).toXMLString()))
                        towrite_annot = source_brsynth_annot.getChild(y)
                        self.logger.debug(brsynth_annot.toXMLString())
                        self._checklibSBML(brsynth_annot.addChild(towrite_annot), ' 1 - Adding annotation to the brsynth annotation')
                        self.logger.debug(brsynth_annot.toXMLString())
                        break
                if not isfound_source:
                    self.logger.error('Cannot find '+str(annot_header)+' in source annotation')
        if not isfound_target:
            self.logger.debug('Cannot find '+str(annot_header)+' in target annotation') 
            isfound_source = False
            source_brsynth_annot = annot_obj.getChild('RDF').getChild('BRSynth').getChild('brsynth')
            for y in range(source_brsynth_annot.getNumChildren()):
                self.logger.debug('\t'+annot_header+' -- '+str(source_brsynth_annot.getChild(y).getName()))
                if str(annot_header)==str(source_brsynth_annot.getChild(y).getName()):
                    isfound_source = True
                    self.logger.debug('Adding annotation to the brsynth annotation: '+str(source_brsynth_annot.getChild(y).toXMLString()))
                    towrite_annot = source_brsynth_annot.getChild(y)
                    self.logger.debug(brsynth_annot.toXMLString())
                    self._checklibSBML(brsynth_annot.addChild(towrite_annot), '2 - Adding annotation to the brsynth annotation')
                    self.logger.debug(brsynth_annot.toXMLString())
                    break
            if not isfound_source:
                self.logger.error('Cannot find '+str(annot_header)+' in source annotation')
            #toWrite_annot = annot_obj.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(annot_header)
            #self._checklibSBML(brsynth_annot.addChild(toWrite_annot), 'Adding annotation to the brsynth annotation')
                return False
        '''
        if brsynth_annot.getChild(annot_header).toXMLString()=='':
            toWrite_annot = annot_obj.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(annot_header)
            self._checklibSBML(brsynth_annot.addChild(toWrite_annot), 'Adding annotation to the brsynth annotation')
        else:
            #try:
            self.logger.debug('==============================')
            found_child = False
            for i in range(brsynth_annot.getNumChildren()):
                if annot_header==brsynth_annot.getChild(i).getName():
                    self.logger.debug('Found the same name to remove: '+str(annot_header))
                    self._checklibSBML(brsynth_annot.removeChild(brsynth_annot.getIndex(i)),
                        'Removing annotation '+str(annot_header))
                    toWrite_annot = annot_obj.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(annot_header)
                    self._checklibSBML(brsynth_annot.addChild(toWrite_annot), 'Adding annotation to the brsynth annotation')
                    found_child = True
                    break
            #cause by a bbug with string lookup
            if not found_child:
                self.logger.warning('Bug with lookup adding it now: '+str(annot_header))
                toWrite_annot = annot_obj.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(annot_header)
                self._checklibSBML(brsynth_annot.addChild(toWrite_annot), 'Adding annotation to the brsynth annotation')
            #except OverflowError:
            #    self.logger.warning('TODO: Overflow error that must be dealt with')
            #    self.logger.warning(brsynth_annot.getChild(annot_header).toXMLString())
            #    return False
        '''
        return True


    def addUpdateMIRIAM(self, sbase_obj, type_param, xref, meta_id=None):
        """Append or update an entry to the MIRIAM annotation of the passed libsbml.SBase object.
        
        If the annot_header isn't contained in the annotation it is created. If it already exists it overwrites it

        :param sbase_obj: libsbml.SBase, The libSBML object to add the different 
        :param type_param: The type of parameter entered. Valid include ['compartment', 'reaction', 'species']
        :param xref: Dictionnary of the cross reference
        :param meta_id: The meta ID to be added to the annotation string
        
        :type type_param: str
        :type xref: dict

        :type meta_id: str

        :rtype: bool
        :return: Sucess or failure of the function
        """
        if not type_param in ['compartment', 'reaction', 'species']:
            self.logger.error('type_param must be '+str(['compartment', 'reaction', 'species'])+' not '+str(type_param))
            return False
        miriam_annot = None
        isReplace = False
        try:
            miriam_annot = sbase_obj.getAnnotation().getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
            miriam_elements = self.readMIRIAMAnnotation(sbase_obj.getAnnotation())
            if not miriam_elements:
                isReplace = True
                if not meta_id:
                    meta_id = self._genMetaID('tmp_addUpdateMIRIAM')
                miriam_annot_1 = libsbml.XMLNode.convertStringToXMLNode(self._defaultBothAnnot(meta_id))
                miriam_annot = miriam_annot_1.getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
            else:
                miriam_elements = None
        except AttributeError:
            try:
                #Cannot find MIRIAM annotation, create it
                isReplace = True
                if not meta_id:
                    meta_id = self._genMetaID('tmp_addUpdateMIRIAM')
                miriam_annot = libsbml.XMLNode.convertStringToXMLNode(self._defaultMIRIAMAnnot(meta_id))
                miriam_annot = miriam_annot.getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
            except AttributeError:
                self.logger.error('Fatal error fetching the annotation')
                return False
        #compile the list of current species
        inside = {}
        for i in range(miriam_annot.getNumChildren()):
            single_miriam = miriam_annot.getChild(i)
            if single_miriam.getAttributes().getLength()>1:
                self.logger.error('MIRIAM annotations should never have more than 1: '+str(single_miriam.toXMLString()))
                continue
            single_miriam_attr = single_miriam.getAttributes()
            if not single_miriam_attr.isEmpty():
                try:
                    db = single_miriam_attr.getValue(0).split('/')[-2]
                    v = single_miriam_attr.getValue(0).split('/')[-1]
                    inside[self.header_miriam[type_param][db]].append(v)
                except KeyError:
                    try:
                        db = single_miriam_attr.getValue(0).split('/')[-2]
                        v = single_miriam_attr.getValue(0).split('/')[-1]
                        inside[self.header_miriam[type_param][db]] = [v]
                    except KeyError:
                        self.logger.warning('Cannot find the self.header_miriram entry '+str(db))
                        continue
            else:
                self.logger.warning('Cannot return MIRIAM attribute')
                pass
        #add or ignore
        toadd = self._compareXref(inside, xref)
        for database_id in toadd:
            for species_id in toadd[database_id]:
                #not sure how to avoid having it that way
                if database_id in self.miriam_header[type_param]:
                    try:
                        #determine if the dictionnaries
                        annotation = '''<annotation>
    <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">
    <rdf:Description rdf:about="#tmp">
      <bqbiol:is>
        <rdf:Bag>'''
                        if type_param=='species':
                            if database_id=='kegg' and species_id[0]=='C':
                                annotation += '''
              <rdf:li rdf:resource="http://identifiers.org/'''+self.miriam_header[type_param]['kegg_c']+str(species_id)+'''"/>'''
                            elif database_id=='kegg' and species_id[0]=='D':
                                annotation += '''
              <rdf:li rdf:resource="http://identifiers.org/'''+self.miriam_header[type_param]['kegg_d']+str(species_id)+'''"/>'''
                            else:
                                annotation += '''
              <rdf:li rdf:resource="http://identifiers.org/'''+self.miriam_header[type_param][database_id]+str(species_id)+'''"/>'''
                        else:
                            annotation += '''
              <rdf:li rdf:resource="http://identifiers.org/'''+self.miriam_header[type_param][database_id]+str(species_id)+'''"/>'''
                        annotation += '''
        </rdf:Bag>
      </bqbiol:is>
    </rdf:Description>
    </rdf:RDF>
    </annotation>'''
                        toPass_annot = libsbml.XMLNode.convertStringToXMLNode(annotation)
                        toWrite_annot = toPass_annot.getChild('RDF').getChild('Description').getChild('is').getChild('Bag').getChild(0)
                        miriam_annot.insertChild(0, toWrite_annot)
                    except KeyError:
                        #WARNING need to check this
                        self.logger.warning('Cannot find '+str(database_id)+' in self.miriam_header for '+str(type_param))
                        continue
        if isReplace:
            ori_miriam_annot = sbase_obj.getAnnotation()
            if ori_miriam_annot==None:
                sbase_obj.unsetAnnotation()
                sbase_obj.setAnnotation(miriam_annot)
            else:
                self._checklibSBML(ori_miriam_annot.getChild('RDF').getChild('Description').getChild('is').removeChild(0), 'Removing annotation "is"')
                self._checklibSBML(ori_miriam_annot.getChild('RDF').getChild('Description').getChild('is').addChild(miriam_annot), 'Adding annotation to the brsynth annotation')
        return True


    def genJSON(self, pathway_id='rp_pathway'):
        """Generate the dictionnary of all the annotations of a pathway species, reaction and pathway annotations

        :param pathway_id: The pathway ID (Default: rp_pathway)
        
        :type pathway_id: str

        :rtype: dict
        :return: Dictionnary of the pathway annotation
        """
        groups = self.model.getPlugin('groups')
        rp_pathway = groups.getGroup(pathway_id)
        reactions = rp_pathway.getListOfMembers()
        #pathway
        rpsbml_json = {}
        rpsbml_json['pathway'] = {}
        rpsbml_json['pathway']['brsynth'] = self.readBRSYNTHAnnotation(rp_pathway.getAnnotation())
        #reactions
        rpsbml_json['reactions'] = {}
        for member in reactions:
            reaction = self.model.getReaction(member.getIdRef())
            annot = reaction.getAnnotation()
            rpsbml_json['reactions'][member.getIdRef()] = {}
            rpsbml_json['reactions'][member.getIdRef()]['brsynth'] = self.readBRSYNTHAnnotation(annot)
            rpsbml_json['reactions'][member.getIdRef()]['miriam'] = self.readMIRIAMAnnotation(annot)
        #loop though all the species
        rpsbml_json['species'] = {}
        for spe_id in self.readUniqueRPspecies(pathway_id):
            species = self.model.getSpecies(spe_id)
            annot = species.getAnnotation()
            rpsbml_json['species'][spe_id] = {}
            rpsbml_json['species'][spe_id]['brsynth'] = self.readBRSYNTHAnnotation(annot)
            rpsbml_json['species'][spe_id]['miriam'] = self.readMIRIAMAnnotation(annot)
        return rpsbml_json


    #####################################################################
    ########################## INPUT/OUTPUT #############################
    #####################################################################


    def readSBML(self, inFile):
        """Open an SBML file to the object

        :param inFile: Path to the input SBML file
        
        :type inFile: str

        :raises FileNotFoundError: If the file cannot be found
        :raises AttributeError: If the libSBML command encounters an error or the input value is None

        :rtype: None
        :return: Dictionnary of the pathway annotation
        """
        if not os.path.isfile(inFile):
            self.logger.error('Invalid input file')
            raise FileNotFoundError
        document = libsbml.readSBMLFromFile(inFile)
        self._checklibSBML(document, 'reading input file')
        errors = document.getNumErrors()
        #display the errors in the log accordning to the severity
        for err in [document.getError(i) for i in range(document.getNumErrors())]:
            #TODO if the error is related to packages not enabled (like groups or fbc) activate them
            if err.isFatal:
                self.logger.error('libSBML reading error: '+str(err.getShortMessage()))
                raise FileNotFoundError
            else:
                self.logger.warning('libSBML reading warning: '+str(err.getShortMessage()))
        model = document.getModel()
        if not model:
            self.logger.error('Either the file was not read correctly or the SBML is empty')
            raise FileNotFoundError
        self.document = document
        self.model = model
        #enabling the extra packages if they do not exists when reading a model
        if not self.model.isPackageEnabled('groups'):
            self._checklibSBML(self.model.enablePackage(
                'http://www.sbml.org/sbml/level3/version1/groups/version1',
                'groups',
                True),
                    'Enabling the GROUPS package')
            self._checklibSBML(self.document.setPackageRequired('groups', False), 'enabling groups package')
        if not self.model.isPackageEnabled('fbc'):
            self._checklibSBML(self.model.enablePackage(
                'http://www.sbml.org/sbml/level3/version1/fbc/version2',
                'fbc',
                True),
                    'Enabling the FBC package')
            self._checklibSBML(self.document.setPackageRequired('fbc', False), 'enabling FBC package')


    def writeSBML(self, path):
        """Export the metabolic network to a SBML file

        :param path: Path to the output SBML file
        
        :type path: str

        :raises FileNotFoundError: If the file cannot be found
        :raises AttributeError: If the libSBML command encounters an error or the input value is None

        :rtype: bool
        :return: Success or failure of the command
        """
        ####### check the path #########
        #need to determine where are the path id's coming from
        p = None
        if path:
            if path[-1:]=='/':
                path = path[:-1]
            if not os.path.isdir(path):
                if self.path:
                    p = self.path
                else:
                    self.logger.error('The output path is not a directory: '+str(path))
                    return False
            else:
                p = path
        else:
            p = self.path
        ########## check and create folder #####
        if not os.path.exists(p):
            os.makedirs(p)
        libsbml.writeSBMLToFile(self.document, p+'/'+str(self.modelName)+'.sbml')
        return True


    #####################################################################
    ########################## FindCreate ###############################
    #####################################################################


    #TODO: seperate the find and the create
    def findCreateObjective(self, reactions, coefficients, isMax=True, objective_id=None):
        """Find the objective (with only one reaction associated) based on the reaction ID and if not found create it

        :param reactions: List of the reactions id's to set as objectives
        :param coefficients: List of the coefficients about the objectives
        :param isMax: Maximise or minimise the objective
        :param objective_id: overwite the default id if created (from obj_[reactions])
        
        :type reactions: list
        :type coefficients: list
        :type isMax: bool
        :type objective_id: str

        :raises FileNotFoundError: If the file cannot be found
        :raises AttributeError: If the libSBML command encounters an error or the input value is None

        :rtype: str
        :return: Objective ID
        """
        fbc_plugin = self.model.getPlugin('fbc')
        self._checklibSBML(fbc_plugin, 'Getting FBC package')
        if not objective_id:
            objective_id = 'obj_'+'_'.join(reactions)
            self.logger.debug('Setting objective as '+str(objective_id))
        for objective in fbc_plugin.getListOfObjectives():
            if objective.getId()==objective_id:
                self.logger.warning('The specified objective id ('+str(objective_id)+') already exists')
                return objective_id
            if not set([i.getReaction() for i in objective.getListOfFluxObjectives()])-set(reactions):
                #TODO: consider setting changing the name of the objective
                self.logger.warning('The specified objective id ('+str(objective_id)+') has another objective with the same reactions: '+str(objective.getId()))
                return objective.getId()
        #If cannot find a valid objective create it
        self.createMultiFluxObj(objective_id,
                                reactions,
                                coefficients,
                                isMax)
        return objective_id


    #####################################################################
    ########################## READ #####################################
    #####################################################################


    #TODO: rename this function to readGroupsMembers
    #TODO: add error handling if the groups does not exist
    #TODO: change the pathway_id to groups_id
    def readRPpathwayIDs(self, pathway_id='rp_pathway'):
        """Return the members of a groups entry

        :param pathway_id: The pathway ID (Default: rp_pathway)
        
        :type pathway_id: str

        :rtype: list
        :return: List of member id's of a particular group
        """
        groups = self.model.getPlugin('groups')
        rp_pathway = groups.getGroup(pathway_id)
        self._checklibSBML(rp_pathway, 'retreiving groups rp_pathway')
        toRet = []
        for member in rp_pathway.getListOfMembers():
            toRet.append(member.getIdRef())
        return toRet


    def readRPrules(self, pathway_id='rp_pathway'):
        """Return the list of reaction rules contained within a pathway

        :param pathway_id: The pathway ID (Default: rp_pathway)
        
        :type pathway_id: str

        :rtype: dict
        :return: Dictionnary of reaction rules (rule_id as key)
        """
        toRet = {}
        for reacId in self.readRPpathwayIDs(pathway_id):
            reac = self.model.getReaction(reacId)
            brsynth_annot = self.readBRSYNTHAnnotation(reac.getAnnotation())
            if not brsynth_annot['rule_id']=='' and not brsynth_annot['smiles']=='':
                toRet[brsynth_annot['rule_id']] = brsynth_annot['smiles'].replace('&gt;', '>')
        return toRet


    #TODO: merge with unique species
    #TODO: change the name of the function to read
    def readRPspecies(self, pathway_id='rp_pathway'):
        """Return the species stoichiometry of a pathway

        :param pathway_id: The pathway ID (Default: rp_pathway)
        
        :type pathway_id: str

        :rtype: dict
        :return: Dictionary of the pathway species and reactions
        """
        reacMembers = {}
        for reacId in self.readRPpathwayIDs(pathway_id):
            reacMembers[reacId] = {}
            reacMembers[reacId]['products'] = {}
            reacMembers[reacId]['reactants'] = {}
            reac = self.model.getReaction(reacId)
            for pro in reac.getListOfProducts():
                reacMembers[reacId]['products'][pro.getSpecies()] = pro.getStoichiometry()
            for rea in reac.getListOfReactants():
                reacMembers[reacId]['reactants'][rea.getSpecies()] = rea.getStoichiometry()
        return reacMembers


    def readUniqueRPspecies(self, pathway_id='rp_pathway'):
        """Return the unique species of a pathway

        :param pathway_id: The pathway ID (Default: rp_pathway)
        
        :type pathway_id: str

        :rtype: list
        :return: List of unique species
        """
        rpSpecies = self.readRPspecies()
        toRet = []
        for i in rpSpecies:
            for y in rpSpecies[i]:
                for z in rpSpecies[i][y]:
                    if not z in toRet:
                        toRet.append(z)
        return toRet
        #reacMembers = self.readRPspecies(pathway_id)
        #return set(set(ori_rp_path['products'].keys())|set(ori_rp_path['reactants'].keys()))


    def readTaxonAnnotation(self, annot):
        """Return he taxonomy ID from an annotation

        :param annot: libsbml.XMLNode, The annotation object of libSBML

        :rtype: dict
        :return: Dictionary of all taxonomy id's
        """
        try:
            toRet = {}
            bag = annot.getChild('RDF').getChild('Description').getChild('hasTaxon').getChild('Bag')
            for i in range(bag.getNumChildren()):
                str_annot = bag.getChild(i).getAttrValue(0)
                if str_annot=='':
                    self.logger.warning('This contains no attributes: '+str(bag.getChild(i).toXMLString()))
                    continue
                dbid = str_annot.split('/')[-2].split('.')[0]
                if len(str_annot.split('/')[-1].split(':'))==2:
                    cid = str_annot.split('/')[-1].split(':')[1]
                else:
                    cid = str_annot.split('/')[-1]
                if not dbid in toRet:
                    toRet[dbid] = []
                toRet[dbid].append(cid)
            return toRet
        except AttributeError:
            return {}


    def readMIRIAMAnnotation(self, annot):
        """Return the MIRIAM annotations of species

        :param annot: libsbml.XMLNode, The annotation object of libSBML

        :rtype: dict
        :return: Dictionary of all the annotation of species
        """
        try:
            toRet = {}
            bag = annot.getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
            for i in range(bag.getNumChildren()):
                str_annot = bag.getChild(i).getAttrValue(0)
                if str_annot=='':
                    self.logger.warning('This contains no attributes: '+str(bag.getChild(i).toXMLString()))
                    continue
                dbid = str_annot.split('/')[-2].split('.')[0]
                if len(str_annot.split('/')[-1].split(':'))==2:
                    cid = str_annot.split('/')[-1].split(':')[1]
                else:
                    cid = str_annot.split('/')[-1]
                if not dbid in toRet:
                    toRet[dbid] = []
                toRet[dbid].append(cid)
            return toRet
        except AttributeError:
            return {}


    def readBRSYNTHAnnotation(self, annot):
        """Return a dictionnary of all the information in a BRSynth annotations

        :param annot: libsbml.XMLNode, The annotation object of libSBML

        :rtype: dict
        :return: Dictionary of all the BRSynth annotations
        """
        toRet = {'dfG_prime_m': None,
                 'dfG_uncert': None,
                 'dfG_prime_o': None,
                 'path_id': None,
                 'step_id': None,
                 'sub_step_id': None,
                 'rule_score': None,
                 'smiles': None,
                 'inchi': None,
                 'inchikey': None,
                 'selenzyme': None,
                 'rule_id': None,
                 'rule_ori_reac': None,
                 'rule_score': None,
                 'global_score': None}
        if annot==None:
            self.logger.warning('The passed annotation is None')
            return {}
        bag = annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
        for i in range(bag.getNumChildren()):
            ann = bag.getChild(i)
            if ann=='':
                self.logger.warning('This contains no attributes: '+str(ann.toXMLString()))
                continue
            if ann.getName()=='dfG_prime_m' or ann.getName()=='dfG_uncert' or ann.getName()=='dfG_prime_o' or ann.getName()[0:4]=='fba_' or ann.getName()=='flux_value':
                try:
                    toRet[ann.getName()] = {
                            'units': ann.getAttrValue('units'),
                            'value': float(ann.getAttrValue('value'))}
                except ValueError:
                    self.logger.warning('Cannot interpret '+str(ann.getName())+': '+str(ann.getAttrValue('value')+' - '+str(ann.getAttrValue('units'))))
                    toRet[ann.getName()] = {
                            'units': None,
                            'value': None}
            elif ann.getName()=='path_id' or ann.getName()=='step_id' or ann.getName()=='sub_step_id':
                try:
                    #toRet[ann.getName()] = int(ann.getAttrValue('value'))
                    toRet[ann.getName()] = {'value': int(ann.getAttrValue('value'))}
                except ValueError:
                    toRet[ann.getName()] = None
            elif ann.getName()=='rule_score' or ann.getName()=='global_score' or ann.getName()[:5]=='norm_':
                try:
                    #toRet[ann.getName()] = float(ann.getAttrValue('value'))
                    toRet[ann.getName()] = {'value': float(ann.getAttrValue('value'))}
                except ValueError:
                    toRet[ann.getName()] = None
            elif ann.getName()=='smiles':
                toRet[ann.getName()] = ann.getChild(0).toXMLString().replace('&gt;', '>')
            #lists in the annotation
            #elif ann.getName()=='selenzyme' or ann.getName()=='rule_ori_reac':
            elif ann.getName()=='selenzyme':
                toRet[ann.getName()] = {}
                for y in range(ann.getNumChildren()):
                    selAnn = ann.getChild(y)
                    try:
                        toRet[ann.getName()][selAnn.getName()] = float(selAnn.getAttrValue('value'))
                    except ValueError:
                        toRet[ann.getName()][selAnn.getName()] = selAnn.getAttrValue('value')
            else:
                toRet[ann.getName()] = ann.getChild(0).toXMLString()
        #to delete empty
        return {k: v for k, v in toRet.items() if v is not None}
        #return toRet


    # TODO: delete
    def readReactionSpecies_old(self, reaction, isID=False):
        """Return the products and the species associated with a reaction
        
        WARNING: DEPRECATED
        
        :param reaction: libsbml.Reaction, Reaction object of libSBML

        :rtype: dict
        :return: Dictionary of the reactions stoichiometry
        """
        #TODO: check that reaction is either an sbml species; if not check that its a string and that
        # it exists in the rpsbml model
        toRet = {'left': {}, 'right': {}}
        #reactants
        for i in range(reaction.getNumReactants()):
            reactant_ref = reaction.getReactant(i)
            reactant = self.model.getSpecies(reactant_ref.getSpecies())
            if isID:
                toRet['left'][reactant.getId()] = int(reactant_ref.getStoichiometry())
            else:
                toRet['left'][reactant.getName()] = int(reactant_ref.getStoichiometry())
        #products
        for i in range(reaction.getNumProducts()):
            product_ref = reaction.getProduct(i)
            product = self.model.getSpecies(product_ref.getSpecies())
            if isID:
                toRet['right'][product.getId()] = int(product_ref.getStoichiometry())
            else:
                toRet['right'][product.getName()] = int(product_ref.getStoichiometry())
            toRet['reversible'] = reaction.getReversible()
        return toRet


    def readReactionSpecies(self, reaction):
        """Return the products and the species associated with a reaction
        
        :param reaction: libsbml.Reaction, Reaction object of libSBML

        :rtype: dict
        :return: Dictionary of the reaction stoichiometry
        """
        #TODO: check that reaction is either an sbml species; if not check that its a string and that
        # it exists in the rpsbml model
        toRet = {'left': {}, 'right': {}}
        #reactants
        for i in range(reaction.getNumReactants()):
            reactant_ref = reaction.getReactant(i)
            toRet['left'][reactant_ref.getSpecies()] = int(reactant_ref.getStoichiometry())
        #products
        for i in range(reaction.getNumProducts()):
            product_ref = reaction.getProduct(i)
            toRet['right'][product_ref.getSpecies()] = int(product_ref.getStoichiometry())
        return toRet


    #####################################################################
    ######################### INQUIRE ###################################
    #####################################################################


    def speciesExists(self, speciesName, compartment_id='MNXC3'):
        """Determine if the model already contains a species according to its ID
        
        :param reaction: libsbml.Reaction, Reaction object of libSBML

        :rtype: bool
        :return: True if exists and False if not
        """
        if speciesName in [i.getName() for i in self.model.getListOfSpecies()] or speciesName+'__64__'+compartment_id in [i.getId() for i in self.model.getListOfSpecies()]:
            return True
        return False


    def isSpeciesProduct(self, species_id, ignoreReactions=[]):
        """Function to determine if a species can be a product of any reaction.
        
        :param species_id: ID of the species to find
        
        :type species_id: str
        :paran ignoreReactions: list, List of all the reaction id's to ignore

        :rtype: bool
        :return: True if its a product of a reaction False if not
        """
        #return all the parameters values
        param_dict = {i.getId(): i.getValue() for i in self.model.parameters}
        for reaction in self.model.getListOfReactions():
            if reaction.getId() not in ignoreReactions:
                #check that the function is reversible by reversibility and FBC bounds
                if reaction.reversible:
                    reaction_fbc = reaction.getPlugin('fbc')
                    #strict left to right
                    if param_dict[reaction_fbc.getLowerFluxBound()]>=0 and param_dict[reaction_fbc.getUpperFluxBound()]>0:
                        if species_id in [i.getSpecies() for i in reaction.getListOfProducts()]:
                            return True
                    #can go both ways
                    elif param_dict[reaction_fbc.getLowerFluxBound()]<0 and param_dict[reaction_fbc.getUpperFluxBound()]>0:
                        if species_id in [i.getSpecies() for i in reaction.getListOfProducts()]:
                            return True
                        elif species_id in [i.getSpecies() for i in reaction.getListOfReactants()]:
                            return True
                    #strict right to left
                    elif param_dict[reaction_fbc.getLowerFluxBound()]<0 and param_dict[reaction_fbc.getUpperFluxBound()]<=0 and param_dict[reaction_fbc.getLowerFluxBound()]<param_dict[reaction_fbc.getUpperFluxBound()]:
                        if species_id in [i.getSpecies() for i in reaction.getListOfReactants()]:
                            return True
                    else:
                        self.logger.warning('isSpeciesProduct does not find the directionailty of the reaction for reaction: '+str(species_id))
                        return True
                else:
                    #if the reaction is not reversible then product are the only way to create it
                    if species_id in [i.getSpecies() for i in reaction.getListOfProducts()]:
                        return True
        return False


    #########################################################################
    ################### CONVERT BETWEEEN FORMATS ############################
    #########################################################################


    def outPathsDict(self, pathway_id='rp_pathway'):
        """Function to return in a dictionary in the same format as the out_paths rp2paths file dictionary object
        
        Example format returned: {'rule_id': 'RR-01-503dbb54cf91-49-F', 'right': {'TARGET_0000000001': 1}, 'left': {'MNXM2': 1, 'MNXM376': 1}, 'pathway_id': 1, 'step': 1, 'sub_step': 1, 'transformation_id': 'TRS_0_0_17'}. Really used to complete the monocomponent reactions

        :param pathway_id: The pathway ID (Default: rp_pathway)
        
        :type pathway_id: str

        :rtype: dict
        :return: Dictionary of the pathway
        """
        pathway = {}
        for member in self.readRPpathwayIDs(pathway_id):
            #TODO: need to find a better way
            reaction = self.model.getReaction(member)
            brsynthAnnot = self.readBRSYNTHAnnotation(reaction.getAnnotation())
            speciesReac = self.readReactionSpecies(reaction)
            self.logger.debug('brsynthAnnot:'+str(brsynthAnnot))
            step = {'reaction_id': member,
                    'reaction_rule': brsynthAnnot['smiles'],
                    'rule_score': brsynthAnnot['rule_score'],
                    'rule_id': brsynthAnnot['rule_id'],
                    'rule_ori_reac': brsynthAnnot['rule_ori_reac'],
                    'right': speciesReac['right'],
                    'left': speciesReac['left'],
                    'path_id': brsynthAnnot['path_id'],
                    'step': brsynthAnnot['step_id'],
                    'sub_step': brsynthAnnot['sub_step_id']}
            self.logger.debug('Step: '+str(step))
            pathway[brsynthAnnot['step_id']['value']] = step
        return pathway


    #########################################################################
    ############################# COMPARE MODELS ############################
    #########################################################################


    def compareBRSYNTHAnnotations(self, source_annot, target_annot):
        """Determine if two libsbml species or reactions have members in common in BRSynth annotation
        
        Compare two dictionnaries and if any of the values of any of the same keys are the same then the function return True, and if none are found then return False

        :param source_annot: libsbml.Reaction, Source object of libSBML
        :param target_annot: libsbml.Reaction, Target object of libSBML

        :rtype: bool
        :return: True if there is at least one similar and False if none
        """
        source_dict = self.readBRSYNTHAnnotation(source_annot)
        target_dict = self.readBRSYNTHAnnotation(target_annot)
        #ignore thse when comparing reactions
        for i in ['path_id', 'step', 'sub_step', 'rule_score', 'rule_ori_reac']:
            try:
                del source_dict[i]
            except KeyError:
                pass
            try:
                del target_dict[i]
            except KeyError:
                pass
        #list the common keys between the two
        for same_key in list(set(list(source_dict.keys())).intersection(list(target_dict.keys()))):
            if source_dict[same_key]==target_dict[same_key]:
                return True
        return False


    def compareMIRIAMAnnotations(self, source_annot, target_annot):
        """Determine if two libsbml species or reactions have members in common in MIRIAM annotation
        
        Compare two dictionnaries and if any of the values of any of the same keys are the same then the function return True, and if none are found then return False

        :param source_annot: libsbml.Reaction, Source object of libSBML
        :param target_annot: libsbml.Reaction, Target object of libSBML

        :rtype: bool
        :return: True if there is at least one similar and False if none
        """
        source_dict = self.readMIRIAMAnnotation(source_annot)
        target_dict = self.readMIRIAMAnnotation(target_annot)
        #list the common keys between the two
        for com_key in set(list(source_dict.keys()))-(set(list(source_dict.keys()))-set(list(target_dict.keys()))):
            #compare the keys and if same is non-empty means that there
            #are at least one instance of the key that is the same
            if bool(set(source_dict[com_key]) & set(target_dict[com_key])):
                return True
        return False


    def compareAnnotations_annot_dict(self, source_annot, target_dict):
        """Compare an annotation object and annotation dictionary
        
        :param source_annot: libsbml.Reaction, Source object of libSBML
        :param target_annot: Target dictionary
        
        :type target_annot: dict

        :rtype: bool
        :return: True if there is at least one similar and False if none
        """
        source_dict = self.readMIRIAMAnnotation(source_annot)
        #list the common keys between the two
        for com_key in set(list(source_dict.keys()))-(set(list(source_dict.keys()))-set(list(target_dict.keys()))):
            #compare the keys and if same is non-empty means that there
            #are at least one instance of the key that is the same
            if bool(set(source_dict[com_key]) & set(target_dict[com_key])):
                return True
        return False


    def compareAnnotations_dict_dict(self, source_dict, target_dict):
        """Compare an annotation as dictionaries
        
        :param source_annot: Source dictionary
        :param target_annot: Target dictionary
        
        :type source_annot: dict
        :type target_annot: dict

        :rtype: bool
        :return: True if there is at least one similar and False if none
        """
        #list the common keys between the two
        for com_key in set(list(source_dict.keys()))-(set(list(source_dict.keys()))-set(list(target_dict.keys()))):
            #compare the keys and if same is non-empty means that there
            #are at least one instance of the key that is the same
            if bool(set(source_dict[com_key]) & set(target_dict[com_key])):
                return True
        return False


    def compareRPpathways(self, measured_sbml):
        """Function to compare two SBML's RP pathways
        
        Function that compares the annotations of reactions and if not found, the annotations of all
        species in that reaction to try to recover the correct ones. Because we are working with
        intermediate cofactors for the RP generated pathways, the annotation crossreference will
        not work. Best is to use the cross-reference to the original reaction

        :param measured_sbml: rpSBML object
        
        :type measured_sbml: rpsbml

        :rtype: bool, dict
        :return: True if there is at least one similar and return the dict of similarities and False if none with empty dictionary 
        """
        #return all the species annotations of the RP pathways
        try:
            meas_rp_species = measured_sbml.readRPspecies()
            found_meas_rp_species = measured_sbml.readRPspecies()
            for meas_step_id in meas_rp_species:
                meas_rp_species[meas_step_id]['annotation'] = measured_sbml.model.getReaction(meas_step_id).getAnnotation()
                found_meas_rp_species[meas_step_id]['found'] = False
                for spe_name in meas_rp_species[meas_step_id]['reactants']:
                    meas_rp_species[meas_step_id]['reactants'][spe_name] = measured_sbml.model.getSpecies(spe_name).getAnnotation()
                    found_meas_rp_species[meas_step_id]['reactants'][spe_name] = False
                for spe_name in meas_rp_species[meas_step_id]['products']:
                    meas_rp_species[meas_step_id]['products'][spe_name] = measured_sbml.model.getSpecies(spe_name).getAnnotation()
                    found_meas_rp_species[meas_step_id]['products'][spe_name] = False
            rp_rp_species = self.readRPspecies()
            for rp_step_id in rp_rp_species:
                rp_rp_species[rp_step_id]['annotation'] = self.model.getReaction(rp_step_id).getAnnotation()
                for spe_name in rp_rp_species[rp_step_id]['reactants']:
                    rp_rp_species[rp_step_id]['reactants'][spe_name] = self.model.getSpecies(spe_name).getAnnotation()
                for spe_name in rp_rp_species[rp_step_id]['products']:
                    rp_rp_species[rp_step_id]['products'][spe_name] = self.model.getSpecies(spe_name).getAnnotation()
        except AttributeError:
            self.logger.error('TODO: debug, for some reason some are passed as None here')
            return False, {}
        #compare the number of steps in the pathway
        if not len(meas_rp_species)==len(rp_rp_species):
            self.logger.warning('The pathways are not of the same length')
            return False, {}
        ############## compare using the reactions ###################
        for meas_step_id in measured_sbml.readRPpathwayIDs():
            for rp_step_id in rp_rp_species:
                if self.compareMIRIAMAnnotations(rp_rp_species[rp_step_id]['annotation'], meas_rp_species[meas_step_id]['annotation']):
                    found_meas_rp_species[meas_step_id]['found'] = True
                    found_meas_rp_species[meas_step_id]['rp_step_id'] = rp_step_id
                    break
        ############## compare using the species ###################
        for meas_step_id in measured_sbml.readRPpathwayIDs():
            #if not found_meas_rp_species[meas_step_id]['found']:
            for rp_step_id in rp_rp_species:
                # We test to see if the meas reaction elements all exist in rp reaction and not the opposite
                #because the measured pathways may not contain all the elements
                ########## reactants ##########
                for meas_spe_id in meas_rp_species[meas_step_id]['reactants']:
                    for rp_spe_id in rp_rp_species[rp_step_id]['reactants']:
                        if self.compareMIRIAMAnnotations(meas_rp_species[meas_step_id]['reactants'][meas_spe_id], rp_rp_species[rp_step_id]['reactants'][rp_spe_id]):
                            found_meas_rp_species[meas_step_id]['reactants'][meas_spe_id] = True
                            break
                        else:
                            if self.compareBRSYNTHAnnotations(meas_rp_species[meas_step_id]['reactants'][meas_spe_id], rp_rp_species[rp_step_id]['reactants'][rp_spe_id]):
                                found_meas_rp_species[meas_step_id]['reactants'][meas_spe_id] = True
                                break
                ########### products ###########
                for meas_spe_id in meas_rp_species[meas_step_id]['products']:
                    for rp_spe_id in rp_rp_species[rp_step_id]['products']:
                        if self.compareMIRIAMAnnotations(meas_rp_species[meas_step_id]['products'][meas_spe_id], rp_rp_species[rp_step_id]['products'][rp_spe_id]):
                            found_meas_rp_species[meas_step_id]['products'][meas_spe_id] = True
                            break
                        else:
                            if self.compareBRSYNTHAnnotations(meas_rp_species[meas_step_id]['products'][meas_spe_id], rp_rp_species[rp_step_id]['products'][rp_spe_id]):
                                found_meas_rp_species[meas_step_id]['products'][meas_spe_id] = True
                                break
                ######### test to see the difference
                pro_found = [found_meas_rp_species[meas_step_id]['products'][i] for i in found_meas_rp_species[meas_step_id]['products']]
                rea_found = [found_meas_rp_species[meas_step_id]['reactants'][i] for i in found_meas_rp_species[meas_step_id]['reactants']]
                if pro_found and rea_found:
                    if all(pro_found) and all(rea_found):
                        found_meas_rp_species[meas_step_id]['found'] = True
                        found_meas_rp_species[meas_step_id]['rp_step_id'] = rp_step_id
                        break
        ################# Now see if all steps have been found ############
        if all(found_meas_rp_species[i]['found'] for i in found_meas_rp_species):
            found_meas_rp_species['measured_model_id'] = measured_sbml.model.getId()
            found_meas_rp_species['rp_model_id'] = self.model.getId()
            return True, found_meas_rp_species
        else:
            return False, {}


    #########################################################################
    ############################# MODEL APPEND ##############################
    #########################################################################


    def setReactionConstraints(self,
                               reaction_id,
                               upper_bound,
                               lower_bound,
                               unit='mmol_per_gDW_per_hr',
                               is_constant=True):
        """Set a given reaction's upper and lower bounds

        Sets the upper and lower bounds of a reaction. Note that if the numerical values passed
        are not recognised, new parameters are created for each of them
        
        :param reaction_id: The id of the reaction
        :param upper_bound: Reaction upper bound
        :param lower_bound: Reaction lower bound
        :param unit: Unit to the bounds (Default: mmol_per_gDW_per_hr)
        :param is_constant: Set if the reaction is constant
        
        :type reaction_id: str
        :type upper_bound: float
        :type lower_bound: float
        :type unit: str
        :type is_constant: bool

        :rtype: tuple or bool
        :return: bool if there is an error and tuple of the lower and upper bound
        """
        reaction = self.model.getReaction(reaction_id)
        if not reaction:
            self.logger.error('Cannot find the reaction: '+str(reaction_id))
            return False
        reac_fbc = reaction.getPlugin('fbc')
        self._checklibSBML(reac_fbc, 'extending reaction for FBC')
        ########## upper bound #############
        old_upper_value = self.model.getParameter(reac_fbc.getUpperFluxBound()).value
        upper_param = self.createReturnFluxParameter(upper_bound, unit, is_constant)
        self._checklibSBML(reac_fbc.setUpperFluxBound(upper_param.getId()),
            'setting '+str(reaction_id)+' upper flux bound')
        ######### lower bound #############
        old_lower_value = self.model.getParameter(reac_fbc.getLowerFluxBound()).value
        lower_param = self.createReturnFluxParameter(lower_bound, unit, is_constant)
        self._checklibSBML(reac_fbc.setLowerFluxBound(lower_param.getId()),
            'setting '+str(reaction_id)+' lower flux bound')
        return old_upper_value, old_lower_value


    ##### ADD SOURCE FROM ORPHAN #####
    #if the heterologous pathway from the self.model contains a sink molecule that is not included in the 
    # original model (we call orhpan species) then add another reaction that creates it
    #TODO: that transports the reactions that creates the species in the
    # extracellular matrix and another reaction that transports it from the extracellular matrix to the cytoplasm
    #TODO: does not work
    def fillOrphan(self,
            rpsbml=None,
            pathway_id='rp_pathway',
            compartment_id='MNXC3',
            upper_flux_bound=999999,
            lower_flux_bound=10):
        """Fill the orgpan

        WARNING: in progress

        :rtype: tuple or bool
        :return: bool if there is an error and tuple of the lower and upper bound
        """
        if rpsbml==None:
            model = self.model
        else:
            model = rpsbml.model
        self.logger.debug('Adding the orphan species to the GEM model')
        #only for rp species
        groups = model.getPlugin('groups')
        rp_pathway = groups.getGroup(pathway_id)
        reaction_id = sorted([(int(''.join(x for x in i.id_ref if x.isdigit())), i.id_ref) for i in rp_pathway.getListOfMembers()], key=lambda tup: tup[0], reverse=True)[0][1]
        #for reaction_id in [i.getId() for i in self.model.getListOfReactions()]:
        for species_id in set([i.getSpecies() for i in model.getReaction(reaction_id).getListOfReactants()]+[i.getSpecies() for i in model.getReaction(reaction_id).getListOfProducts()]):
            if rpsbml==None:
                isSpePro = self.isSpeciesProduct(species_id, [reaction_id])
            else:
                isSpePro = rpsbml.isSpeciesProduct(species_id, [reaction_id])
            if not isSpePro:
                #create the step
                createStep = {'rule_id': None,
                              'left': {species_id.split('__')[0]: 1},
                              'right': {},
                              'step': None,
                              'sub_step': None,
                              'path_id': None,
                              'transformation_id': None,
                              'rule_score': None,
                              'rule_ori_reac': None}
                #create the model in the
                if rpsbml==None:
                    self.createReaction('create_'+species_id,
                                        upper_flux_bound,
                                        lower_flux_bound,
                                        createStep,
                                        compartment_id)
                else:
                    rpsbml.createReaction('create_'+species_id,
                                        upper_flux_bound,
                                        lower_flux_bound,
                                        createStep,
                                        compartment_id)


    #########################################################################
    ############################# MODEL CREATION FUNCTIONS ##################
    #########################################################################


    def createModel(self, name, model_id, meta_id=None):
        """Create libSBML model instance

        Function that creates a new libSBML model instance and initiates it with the appropriate packages. Creates a cytosol compartment

        :param name: The name of the of the model
        :param model_id: The id of the model
        :param meta_id: Meta ID of the model
        
        :type name: str
        :type model_id: str
        :type meta_id: str

        :rtype: None
        :return: None
        """
        ## sbmldoc
        self.sbmlns = libsbml.SBMLNamespaces(3,1)
        self._checklibSBML(self.sbmlns, 'generating model namespace')
        self._checklibSBML(self.sbmlns.addPkgNamespace('groups',1), 'Add groups package')
        self._checklibSBML(self.sbmlns.addPkgNamespace('fbc',2), 'Add FBC package')
        #sbmlns = libsbml.SBMLNamespaces(3,1,'groups',1)
        self.document = libsbml.SBMLDocument(self.sbmlns)
        self._checklibSBML(self.document, 'generating model doc')
        #!!!! must be set to false for no apparent reason
        self._checklibSBML(self.document.setPackageRequired('fbc', False), 'enabling FBC package')
        #!!!! must be set to false for no apparent reason
        self._checklibSBML(self.document.setPackageRequired('groups', False), 'enabling groups package')
        ## sbml model
        self.model = self.document.createModel()
        self._checklibSBML(self.model, 'generating the model')
        self._checklibSBML(self.model.setId(model_id), 'setting the model ID')
        model_fbc = self.model.getPlugin('fbc')
        model_fbc.setStrict(True)
        if meta_id==None:
            meta_id = self._genMetaID(model_id)
        self._checklibSBML(self.model.setMetaId(meta_id), 'setting model meta_id')
        self._checklibSBML(self.model.setName(name), 'setting model name')
        self._checklibSBML(self.model.setTimeUnits('second'), 'setting model time unit')
        self._checklibSBML(self.model.setExtentUnits('mole'), 'setting model compartment unit')
        self._checklibSBML(self.model.setSubstanceUnits('mole'), 'setting model substance unit')


    #TODO: set the compName as None by default. To do that you need to regenerate the compXref to
    #TODO: consider seperating it in another function if another compartment is to be created
    #TODO: use MNX ids as keys instead of the string names
    def createCompartment(self, size, compId, compName, compXref, meta_id=None):
        """Create libSBML compartment
        
        :param size: Size of the compartment
        :param compId: Compartment id
        :param compName: Compartment Name
        :param compXref: Cross reference dictionary of the compartment
        :param meta_id: Meta id (Default: None)
        
        :type size: float
        :type compId: str
        :type compName: str
        :type compXref: dict
        :type meta_id: str

        :rtype: None
        :return: None
        """
        comp = self.model.createCompartment()
        self._checklibSBML(comp, 'create compartment')
        self._checklibSBML(comp.setId(compId), 'set compartment id')
        if compName:
            self._checklibSBML(comp.setName(compName), 'set the name for the cytoplam')
        self._checklibSBML(comp.setConstant(True), 'set compartment "constant"')
        self._checklibSBML(comp.setSize(size), 'set compartment "size"')
        self._checklibSBML(comp.setSBOTerm(290), 'set SBO term for the cytoplasm compartment')
        if meta_id==None:
            meta_id = self._genMetaID(compId)
        self._checklibSBML(comp.setMetaId(meta_id), 'set the meta_id for the compartment')
        ############################ MIRIAM ############################
        comp.setAnnotation(libsbml.XMLNode.convertStringToXMLNode(self._defaultMIRIAMAnnot(meta_id)))
        self.addUpdateMIRIAM(comp, 'compartment', compXref, meta_id)


    def createUnitDefinition(self, unit_id, meta_id=None):
        """Create libSBML unit definition

        Function that creates a unit definition (composed of one or more units)

        :param unit_id: Unit id definition
        :param meta_id: Meta id (Default: None)
        
        :type unit_id: str
        :type meta_id: str

        :rtype: libsbml.UnitDefinition
        :return: Unit definition object created
        """
        unitDef = self.model.createUnitDefinition()
        self._checklibSBML(unitDef, 'creating unit definition')
        self._checklibSBML(unitDef.setId(unit_id), 'setting id')
        if meta_id==None:
            meta_id = self._genMetaID(unit_id)
        self._checklibSBML(unitDef.setMetaId(meta_id), 'setting meta_id')
        #self.unitDefinitions.append(unit_id)
        return unitDef


    def createUnit(self, unitDef, libsbmlunit, exponent, scale, multiplier):
        """Set or update the parameters of a libSBML unit definition

        :param unitDef: libsbml.Unit, libSBML Unit 
        :param libsbmlunit: String unit
        :param exponent: Exponent unit
        :param sale: Scale of the unit 
        :param multiplier: Multiplier of the unit 
        
        :type libsbmlunit: str
        :type exponent: int
        :type sale: int
        :type multiplier: int

        :rtype: None
        :return: None
        """
        unit = unitDef.createUnit()
        self._checklibSBML(unit, 'creating unit')
        self._checklibSBML(unit.setKind(libsbmlunit), 'setting the kind of unit')
        self._checklibSBML(unit.setExponent(exponent), 'setting the exponenent of the unit')
        self._checklibSBML(unit.setScale(scale), 'setting the scale of the unit')
        self._checklibSBML(unit.setMultiplier(multiplier), 'setting the multiplier of the unit')


    def createReturnFluxParameter(self,
            value,
            unit='mmol_per_gDW_per_hr',
            is_constant=True,
            parameter_id=None,
            meta_id=None):
        """Create libSBML flux parameters

        Parameters are used for the bounds for FBA analysis. Unit parameter must be an instance of unitDefinition.
        If the parameter id exists, then the function returns the libsbml.Parameter object

        :param value: Value set for the parameter
        :param unit: The unit id of the parameter
        :param is_constant: Define if the parameter is constant
        :param parameter_id: Overwrite the default naming convention (Default: None)
        :param meta_id: Meta id (Default: None)
        
        :type value: float
        :type unit: str
        :type is_constant: bool
        :type parameter_id: str
        :type meta_id: str

        :rtype: libsbml.Parameter
        :return: The newly created libsbml.Parameter
        """
        if parameter_id:
            param_id = parameter_id
        else:
            if value>=0:
                param_id = 'B_'+str(round(abs(value), 4)).replace('.', '_')
            else:
                param_id = 'B__'+str(round(abs(value), 4)).replace('.', '_')
        if param_id in [i.getId() for i in self.model.getListOfParameters()]:
            return self.model.getParameter(param_id)
        else:
            newParam = self.model.createParameter()
            self._checklibSBML(newParam, 'Creating a new parameter object')
            self._checklibSBML(newParam.setConstant(is_constant), 'setting as constant')
            self._checklibSBML(newParam.setId(param_id), 'setting ID')
            self._checklibSBML(newParam.setValue(value), 'setting value')
            self._checklibSBML(newParam.setUnits(unit), 'setting units')
            self._checklibSBML(newParam.setSBOTerm(625), 'setting SBO term')
            if meta_id==None:
                meta_id = self._genMetaID(parameter_id)
            self._checklibSBML(newParam.setMetaId(meta_id), 'setting meta ID')
            #self.parameters.append(parameter_id)
            return newParam


    #TODO as of now not generic, works when creating a new SBML file, but no checks if modifying existing SBML file
    def createReaction(self,
                       reac_id,
                       fluxUpperBound,
                       fluxLowerBound,
                       step,
                       compartment_id,
                       reaction_smiles=None,
                       reacXref={},
                       pathway_id=None,
                       meta_id=None):
        """Create libSBML reaction

        Create a reaction that is added to the self.model in the input compartment id. fluxBounds is a list of libSBML.UnitDefinition, length of exactly 2 with the first position that is the upper bound and the second is the lower bound. reactants_dict and reactants_dict are dictionnaries that hold the following parameters: name, compartment, stoichiometry

        :param name: Name of the reaction
        :param fluxUpperBound: The reaction fbc upper bound
        :param fluxLowerBound: The reaction fbc lower bound
        :param step: The id's of the reactant and products of the reactions. Example: {'left': [], 'right': []}
        :param compartment_id: The id of the compartment to add the reaction
        :param reaction_smiles: The reaction rule to add to the BRSynth annotation of the reaction (Default: None)
        :param reacXref: The dict containing the MIRIAM annotation (Default: {})
        :param pathway_id: The Groups id of the reaction to which the reacion id will be added (Default: None)
        :param meta_id: Meta id (Default: None)
        
        :type name: str
        :type fluxUpperBound: float
        :type fluxLowerBound: float
        :type step: dict
        :type compartment_id: str
        :type reaction_smiles: str
        :type reacXref: dict
        :type pathway_id: str
        :type meta_id: str

        :rtype: None
        :return: None
        """
        reac = self.model.createReaction()
        self._checklibSBML(reac, 'create reaction')
        ################ FBC ####################
        reac_fbc = reac.getPlugin('fbc')
        self._checklibSBML(reac_fbc, 'extending reaction for FBC')
        #bounds
        upper_bound = self.createReturnFluxParameter(fluxUpperBound)
        lower_bound = self.createReturnFluxParameter(fluxLowerBound)
        self._checklibSBML(reac_fbc.setUpperFluxBound(upper_bound.getId()), 'setting '+str(reac_id)+' upper flux bound')
        self._checklibSBML(reac_fbc.setLowerFluxBound(lower_bound.getId()), 'setting '+str(reac_id)+' lower flux bound')
        #########################################
        #reactions
        self._checklibSBML(reac.setId(reac_id), 'set reaction id') #same convention as cobrapy
        self._checklibSBML(reac.setSBOTerm(176), 'setting the system biology ontology (SBO)') #set as process
        #TODO: consider having the two parameters as input to the function
        self._checklibSBML(reac.setReversible(True), 'set reaction reversibility flag')
        self._checklibSBML(reac.setFast(False), 'set reaction "fast" attribute')
        if meta_id==None:
            meta_id = self._genMetaID(reac_id)
        self._checklibSBML(reac.setMetaId(meta_id), 'setting species meta_id')
        #TODO: check that the species exist
        #reactants_dict
        for reactant in step['left']:
            spe = reac.createReactant()
            self._checklibSBML(spe, 'create reactant')
            #use the same writing convention as CobraPy
            self._checklibSBML(spe.setSpecies(str(reactant)+'__64__'+str(compartment_id)), 'assign reactant species')
            #TODO: check to see the consequences of heterologous parameters not being constant
            self._checklibSBML(spe.setConstant(True), 'set "constant" on species '+str(reactant))
            self._checklibSBML(spe.setStoichiometry(float(step['left'][reactant])),
                'set stoichiometry ('+str(float(step['left'][reactant]))+')')
        #TODO: check that the species exist
        #products_dict
        for product in step['right']:
            pro = reac.createProduct()
            self._checklibSBML(pro, 'create product')
            self._checklibSBML(pro.setSpecies(str(product)+'__64__'+str(compartment_id)), 'assign product species')
            #TODO: check to see the consequences of heterologous parameters not being constant
            self._checklibSBML(pro.setConstant(True), 'set "constant" on species '+str(product))
            self._checklibSBML(pro.setStoichiometry(float(step['right'][product])),
                'set the stoichiometry ('+str(float(step['right'][product]))+')')
        ############################ MIRIAM ############################
        self._checklibSBML(reac.setAnnotation(self._defaultBothAnnot(meta_id)), 'creating annotation')
        self.addUpdateMIRIAM(reac, 'reaction', reacXref, meta_id)
        ###### BRSYNTH additional information ########
        if reaction_smiles:
            self.addUpdateBRSynth(reac, 'smiles', reaction_smiles, None, True, False, False, meta_id)
        if step['rule_id']:
            self.addUpdateBRSynth(reac, 'rule_id', step['rule_id'], None, True, False, False, meta_id)
        #TODO: need to change the name and content (to dict) upstream
        if step['rule_ori_reac']:
            self.addUpdateBRSynth(reac, 'rule_ori_reac', step['rule_ori_reac'], None, True, False, False, meta_id)
        if step['rule_score']:
            self.addUpdateBRSynth(reac, 'rule_score', step['rule_score'], None, False, False, False, meta_id)
        if step['path_id']:
            self.addUpdateBRSynth(reac, 'path_id', step['path_id'], None, False, False, False, meta_id)
        if step['step']:
            self.addUpdateBRSynth(reac, 'step_id', step['step'], None, False, False, False, meta_id)
        if step['sub_step']:
            self.addUpdateBRSynth(reac, 'sub_step_id', step['sub_step'], None, False, False, False, meta_id)
        #### GROUPS #####
        if not pathway_id==None:
            groups_plugin = self.model.getPlugin('groups')
            hetero_group = groups_plugin.getGroup(pathway_id)
            if not hetero_group:
                self.logger.warning('The pathway_id '+str(pathway_id)+' does not exist in the model')
            else:
                newM = hetero_group.createMember()
                self._checklibSBML(newM, 'Creating a new groups member')
                self._checklibSBML(newM.setIdRef(reac_id), 'Setting name to the groups member')


    def createSpecies(self,
                      species_id,
                      compartment_id,
                      species_name=None,
                      chemXref={},
                      inchi=None,
                      inchikey=None,
                      smiles=None,
                      species_group_id=None,
                      in_sink_group_id=None,
                      meta_id=None):
                      #TODO: add these at some point -- not very important
                      #charge=0,
                      #chemForm=''):
        """Create libSBML species

        Create a species that is added to self.model

        :param species_id: The id of the created species
        :param compartment_id: The id of the compartment to add the reaction
        :param species_name: Overwrite the default name of the created species (Default: None)
        :param chemXref: The dict containing the MIRIAM annotation (Default: {})
        :param inchi: The InChI string to be added to BRSynth annotation (Default: None)
        :param inchikey: The InChIkey string to be added to BRSynth annotation (Default: None)
        :param smiles: The SMLIES string to be added to BRSynth annotation (Default: None)
        :param species_group_id: The Groups id to add the species (Default: None)
        :param in_sink_group_id: The Groups id sink species to add the species (Default: None)
        :param meta_id: Meta id (Default: None)
        
        :type species_id: str
        :type compartment_id: str
        :type species_name: str
        :type chemXref: dict
        :type inchi: str
        :type inchikey: str
        :type smiles: str
        :type species_group_id: str
        :type in_sink_group_id: str
        :type meta_id: str

        :rtype: None
        :return: None
        """
        spe = self.model.createSpecies()
        self._checklibSBML(spe, 'create species')
        ##### FBC #####
        spe_fbc = spe.getPlugin('fbc')
        self._checklibSBML(spe_fbc, 'creating this species as an instance of FBC')
        #spe_fbc.setCharge(charge) #### These are not required for FBA
        #spe_fbc.setChemicalFormula(chemForm) #### These are not required for FBA
        #if compartment_id:
        self._checklibSBML(spe.setCompartment(compartment_id), 'set species spe compartment')
        #else:
        #    #removing this could lead to errors with xref
        #    self._checklibSBML(spe.setCompartment(self.compartment_id), 'set species spe compartment')
        #ID same structure as cobrapy
        #TODO: determine if this is always the case or it will change
        self._checklibSBML(spe.setHasOnlySubstanceUnits(False), 'set substance units')
        self._checklibSBML(spe.setBoundaryCondition(False), 'set boundary conditions')
        self._checklibSBML(spe.setConstant(False), 'set constant')
        #useless for FBA (usefull for ODE) but makes Copasi stop complaining
        self._checklibSBML(spe.setInitialConcentration(1.0), 'set an initial concentration')
        #same writting convention as COBRApy
        self._checklibSBML(spe.setId(str(species_id)+'__64__'+str(compartment_id)), 'set species id')
        self.logger.debug('Setting species id as: '+str(species_id)+'__64__'+str(compartment_id))
        if meta_id==None:
            meta_id = self._genMetaID(species_id)
        self._checklibSBML(spe.setMetaId(meta_id), 'setting reaction meta_id')
        if species_name==None:
            self._checklibSBML(spe.setName(species_id), 'setting name for the metabolite '+str(species_id))
        else:
            self._checklibSBML(spe.setName(species_name), 'setting name for the metabolite '+str(species_name))
        #this is setting MNX id as the name
        #this is setting the name as the input name
        #self._checklibSBML(spe.setAnnotation(self._defaultBRSynthAnnot(meta_id)), 'creating annotation')
        self._checklibSBML(spe.setAnnotation(self._defaultBothAnnot(meta_id)), 'creating annotation')
        ###### annotation ###
        self.addUpdateMIRIAM(spe, 'species', chemXref, meta_id)
        ###### BRSYNTH additional information ########
        if smiles:
            self.addUpdateBRSynth(spe, 'smiles', smiles, None, True, False, False, meta_id)
            #                   sbase_obj, annot_header, value, units=None, isAlone=False, isList=False, isSort=True, meta_id=None)
        if inchi:
            self.addUpdateBRSynth(spe, 'inchi', inchi, None, True, False, False, meta_id)
        if inchikey:
            self.addUpdateBRSynth(spe, 'inchikey', inchikey, None, True, False, False, meta_id)
            self.addUpdateMIRIAM(spe, 'species', {'inchikey': [inchikey]})
        #### GROUPS #####
        #TODO: check that it actually exists
        self.logger.debug('species_group_id: '+str(species_group_id))
        if not species_group_id==None:
            groups_plugin = self.model.getPlugin('groups')
            hetero_group = groups_plugin.getGroup(species_group_id)
            if not hetero_group:
                self.logger.warning('The species_group_id '+str(species_group_id)+' does not exist in the model')
                #TODO: consider creating it if
            else:
                newM = hetero_group.createMember()
                self._checklibSBML(newM, 'Creating a new groups member')
                self._checklibSBML(newM.setIdRef(str(species_id)+'__64__'+str(compartment_id)), 'Setting name to the groups member')
        #TODO: check that it actually exists
        #add the species to the sink species
        self.logger.debug('in_sink_group_id: '+str(in_sink_group_id))
        if not in_sink_group_id==None:
            groups_plugin = self.model.getPlugin('groups')
            sink_group = groups_plugin.getGroup(in_sink_group_id)
            if not sink_group:
                self.logger.warning('The species_group_id '+str(in_sink_group_id)+' does not exist in the model')
                #TODO: consider creating it if
            else:
                newM = sink_group.createMember()
                self._checklibSBML(newM, 'Creating a new groups member')
                self._checklibSBML(newM.setIdRef(str(species_id)+'__64__'+str(compartment_id)), 'Setting name to the groups member') 


    #TODO: change the name of this function to createGroup
    def createPathway(self, pathway_id, meta_id=None):
        """Create libSBML pathway

        Create a pathway that is added to self.model

        :param pathway_id: The Groups id of the pathway id (Default: None)
        :param meta_id: Meta id (Default: None)
        
        :type pathway_id: str
        :type meta_id: str

        :rtype: None
        :return: None
        """
        groups_plugin = self.model.getPlugin('groups')
        new_group = groups_plugin.createGroup()
        self.logger.debug('setting new group id: '+str(pathway_id))
        new_group.setId(pathway_id)
        if meta_id==None:
            meta_id = self._genMetaID(pathway_id)
        new_group.setMetaId(meta_id)
        new_group.setKind(libsbml.GROUP_KIND_COLLECTION)
        new_group.setAnnotation(self._defaultBRSynthAnnot(meta_id))


    def createGene(self, reac, step_id, meta_id=None):
        """Create libSBML gene

        Create a gene that is associated with a reaction

        :param reac: The id of the reaction that is associated with the gene
        :param step_id: The id of the reaction to name the gene
        :param meta_id: Meta id (Default: None)
        
        :type reac: str
        :type step_id: str
        :type meta_id: str

        :rtype: None
        :return: None
        """
        #TODO: pass this function to Pablo for him to fill with parameters that are appropriate for his needs
        geneName = 'RP'+str(step_id)+'_gene'
        fbc_plugin = self.model.getPlugin('fbc')
        #fbc_plugin = reac.getPlugin("fbc")
        gp = fbc_plugin.createGeneProduct()
        gp.setId(geneName)
        if meta_id==None:
            meta_id = self._genMetaID(str(geneName))
        gp.setMetaId(meta_id)
        gp.setLabel('gene_'+str(step_id))
        gp.setAssociatedSpecies('RP'+str(step_id))
        ##### NOTE: The parameters here require the input from Pablo to determine what he needs
        #gp.setAnnotation(self._defaultBothAnnot(meta_id))


    def createFluxObj(self, fluxobj_id, reactionName, coefficient, isMax=True, meta_id=None):
        """Create libSBML flux objective

        WARNING DEPRECATED -- use the createMultiFluxObj() with lists of size one to define an objective function
        with a single reaction
        Using the FBC package one can add the FBA flux objective directly to the model. This function sets a particular reaction as objective with maximization or minimization objectives
        
        :param fluxobj_id: The id of the flux objective 
        :param reactionName: The id of the reaction that is associated with the reaction
        :param coefficient: The coefficient of the flux objective
        :param isMax: Define if the objective is coefficient
        :param meta_id: Meta id (Default: None)
        
        :type fluxobj_id: str
        :type reactionName: str
        :type coefficient: int
        :type isMax: bool
        :type meta_id: str

        :rtype: None
        :return: None
        """
        fbc_plugin = self.model.getPlugin('fbc')
        target_obj = fbc_plugin.createObjective()
        #TODO: need to define inpiut metaID
        target_obj.setAnnotation(self._defaultBRSynthAnnot(meta_id))
        target_obj.setId(fluxobj_id)
        if isMax:
            target_obj.setType('maximize')
        else:
            target_obj.setType('minimize')
        fbc_plugin.setActiveObjectiveId(fluxobj_id) # this ensures that we are using this objective when multiple
        target_flux_obj = target_obj.createFluxObjective()
        target_flux_obj.setReaction(reactionName)
        target_flux_obj.setCoefficient(coefficient)
        if meta_id==None:
            meta_id = self._genMetaID(str(fluxobj_id))
        target_flux_obj.setMetaId(meta_id)
        target_flux_obj.setAnnotation(self._defaultBRSynthAnnot(meta_id))


    def createMultiFluxObj(self, fluxobj_id, reactionNames, coefficients, isMax=True, meta_id=None):
        """Create libSBML flux objective

        Using the FBC package one can add the FBA flux objective directly to the model. Can add multiple reactions. This function sets a particular reaction as objective with maximization or minimization objectives
        
        :param fluxobj_id: The id of the flux objective 
        :param reactionNames: The list of string id's of the reaction that is associated with the reaction
        :param coefficients: The list of int defining the coefficients of the flux objective
        :param isMax: Define if the objective is coefficient
        :param meta_id: Meta id (Default: None)
        
        :type fluxobj_id: str
        :type reactionNames: list
        :type coefficients: list
        :type isMax: bool
        :type meta_id: str

        :rtype: None
        :return: None
        """
        if not len(reactionNames)==len(coefficients):
            self.logger.error('The size of reactionNames is not the same as coefficients')
            return False
        fbc_plugin = self.model.getPlugin('fbc')
        target_obj = fbc_plugin.createObjective()
        target_obj.setAnnotation(self._defaultBRSynthAnnot(meta_id))
        target_obj.setId(fluxobj_id)
        if isMax:
            target_obj.setType('maximize')
        else:
            target_obj.setType('minimize')
        fbc_plugin.setActiveObjectiveId(fluxobj_id) # this ensures that we are using this objective when multiple
        for reac, coef in zip(reactionNames, coefficients):
            target_flux_obj = target_obj.createFluxObjective()
            target_flux_obj.setReaction(reac)
            target_flux_obj.setCoefficient(coef)
            if meta_id==None:
                meta_id = self._genMetaID(str(fluxobj_id))
            target_flux_obj.setMetaId(meta_id)
            target_flux_obj.setAnnotation(self._defaultBRSynthAnnot(meta_id))


    ##############################################################################################
    ############################### Generic Model ################################################
    ##############################################################################################


    def genericModel(self,
                     modelName,
                     model_id,
                     compXref,
                     compartment_id,
                     upper_flux_bound=999999,
                     lower_flux_bound=0):
        """Generate a generic model

        Since we will be using the same type of parameters for the RetroPath model, this function
        generates a libSBML model with parameters that will be mostly used
        
        :param modelName: The given name of the model
        :param model_id: The id of the model
        :param compXref: The model MIRIAM annotation
        :param upper_flux_bound: The upper flux bounds unit definitions default when adding new reaction (Default: 999999.0)
        :param lower_flux_bound: The lower flux bounds unit definitions default when adding new reaction (Defaul: 0.0)
        
        :type modelName: str
        :type model_id: str
        :type compXref: dict
        :type upper_flux_bound: float
        :type lower_flux_bound: float

        :rtype: None
        :return: None
        """
        self.createModel(modelName, model_id)
        # mmol_per_gDW_per_hr -- flux
        unitDef = self.createUnitDefinition('mmol_per_gDW_per_hr')
        self.createUnit(unitDef, libsbml.UNIT_KIND_MOLE, 1, -3, 1)
        self.createUnit(unitDef, libsbml.UNIT_KIND_GRAM, 1, 0, 1)
        self.createUnit(unitDef, libsbml.UNIT_KIND_SECOND, 1, 0, 3600)
        # kj_per_mol -- thermodynamics
        gibbsDef = self.createUnitDefinition('kj_per_mol')
        self.createUnit(gibbsDef, libsbml.UNIT_KIND_JOULE, 1, 3, 1)
        self.createUnit(gibbsDef, libsbml.UNIT_KIND_MOLE, -1, 1, 1)
        ### set the bounds
        upBound = self.createReturnFluxParameter(upper_flux_bound)
        lowBound = self.createReturnFluxParameter(lower_flux_bound)
        #compartment
        #TODO: try to recover the name from the Xref
        try:
            name = compXref['name'][0]
        except KeyError:
            name = compartment_id+'_name'
        self.createCompartment(1, compartment_id, name, compXref)

