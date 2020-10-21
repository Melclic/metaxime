from rdkit.Chem import MolFromInchi
from rdkit.Chem import Draw
import copy
import json
import re
import logging
import drawSvg as draw
import svgutils.transform as sg


## Class that contains a collection to draw a rpSBML file
#
#
class rpDraw:
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.mnx_cofactors = json.load(open('data/mnx_cofactors.json', 'r'))
        #some drawing constants
        self.arrowhead = draw.Marker(-0.1, -0.5, 0.9, 0.5, scale=4, orient='auto', id='normal_arrow')
        self.arrowhead.append(draw.Lines(-0.1, -0.5, -0.1, 0.5, 0.9, 0, fill='black', close=True))
        self.arrowhead_flat = draw.Marker(-0.1, -0.5, 0.9, 0.5, scale=4, orient=0, id='flat_arrow')
        self.arrowhead_flat.append(draw.Lines(-0.1, -0.5, -0.1, 0.5, 0.9, 0, fill='black', close=True))
        self.rev_arrowhead_flat = draw.Marker(-0.1, -0.5, 0.9, 0.5, scale=4, orient=180, id='rev_flat_arrow')
        self.rev_arrowhead_flat.append(draw.Lines(-0.1, -0.5, -0.1, 0.5, 0.9, 0, fill='black', close=True))
        '''
        self.rev_arrowhead = draw.Marker(-0.1, -0.5, 0.9, 0.5, scale=4, orient=0, id='rev_flat_arrow')
        self.rev_arrowhead.append(draw.Lines(-0.1, -0.5, -0.1, 0.5, 0.9, 0, fill='black', close=True))
        '''
        self.arrowhead_comp_x = 7.0
        self.arrowhead_comp_y = 7.0


    ######################################################################################################
    ######################################### Private Function ###########################################
    ######################################################################################################


    ## Organise the heterologous pathway in a tree like structure
    #
    # This method iterates from the TARGET and organises the heterologous in a tree like manner. The 
    # TODO: add the global filter of cofactors to remove the species that are shared among many (ex: H+, O2, etc...)
    #BUG: /Users/melchior/Downloads/Galaxy111/rp_50_28.sbml.xml not having the  ink entry added
    #TODO: move this to rpGraph
    def _hierarchy_pos(self,
                       G,
                       root,
                       width=1.0,
                       y_gap=0.2,
                       xcenter=0.5,
                       plot_only_central=False,
                       filter_cofactors=True,
                       filter_sink_species=False,
                       central_species_group_id='central_species',
                       sink_species_group_id='sink_species'):
        #### find all the nodes that are in a layer ###########
        def _make_layer(parent_nodes):
            layer = []
            for current_node in parent_nodes:
                self.logger.debug('\t'+str(current_node)+' --> '+str(list(G.predecessors(current_node))+list(G.successors(current_node))))
                for nei in list(G.predecessors(current_node))+list(G.successors(current_node)):
                #self.logger.debug('\t'+str(current_node)+' --> '+str(list(G.predecessors(current_node))))
                #for nei in list(G.predecessors(current_node)):
                    if nei in toadd_nodes and nei not in layer:
                        self.logger.debug('\tAdding node: '+str(nei))
                        layer.append(nei)
                        #now check that the predeceessor of that node is a reaction and if yes add its successors to the current layer
                        self.logger.debug('\t'+str(nei)+' predecessors: '+str(list(G.predecessors(nei))))
                        for nei_pre in list(G.predecessors(nei)):
                            self.logger.debug('\t\t'+str(nei_pre)+' type: '+str(G.node.get(nei_pre)['type']))
                            if G.node.get(nei_pre)['type']=='reaction':
                                self.logger.debug('\t\t'+str(nei_pre)+' successors: '+str(list(G.successors(nei_pre))))
                                for reac_nei_suc in list(G.successors(nei_pre)):
                                    if reac_nei_suc in toadd_nodes and reac_nei_suc not in layer:
                                        self.logger.debug('\tAdding node: '+str(reac_nei_suc))
                                        layer.append(reac_nei_suc)
            return layer
        ##### filter the nodes that will not be used #######
        filtered_species = []
        reac_cofactors_id = {}
        toadd_nodes = list(set(list(G.nodes)))
        for node in list(set(list(G.nodes))):
            self.logger.debug('--------- '+str(node)+' ---------')
            node_obj = G.node.get(node)
            self.logger.debug(node_obj)
            if node_obj['type']=='reaction':
                reac_cofactors_id[node] = {'substrates': [], 'products': []}
            elif node_obj['type']=='species':
                #above all else, if there are no InChI structuures then filter
                if not 'inchi' in node_obj['brsynth']:
                    self.logger.warning(str(node)+': has no InChI structure associated with it, adding as cofactor')
                    toadd_nodes.remove(node)
                    filtered_species.append(node)
                    continue
                if filter_cofactors:
                    if 'metanetx' in node_obj['miriam']:
                        if any([i in list(self.mnx_cofactors.keys()) for i in node_obj['miriam']['metanetx']]):
                            self.logger.debug('filter_cofactors: '+str(filter_cofactors))

                            self.logger.warning(str(node)+' is a list cofactor and is filtered')
                            toadd_nodes.remove(node)
                            filtered_species.append(node)
                            continue
                if not filter_sink_species and node_obj[sink_species_group_id]:
                    self.logger.debug('filter_sink_species: '+str(filter_sink_species))
                    self.logger.debug('node_obj['+str(sink_species_group_id)+']: '+str(node_obj[sink_species_group_id]))
                    self.logger.warning(str(node)+' is a sink species and will not be filtered')
                    continue
                if plot_only_central and not node_obj[central_species_group_id]:
                    self.logger.debug('plot_only_central: '+str(plot_only_central))
                    self.logger.debug('node_obj['+str(central_species_group_id)+']: '+str(node_obj[central_species_group_id]))
                    self.logger.warning(str(node)+' is not a central species and is filtered')
                    toadd_nodes.remove(node)
                    filtered_species.append(node)
                    continue
        ###### create the cofactor list based on the ones you are removing
        for edge in list(G.edges):
            node_obj_source = G.node.get(edge[0])
            node_obj_target = G.node.get(edge[1])
            if node_obj_source['type']=='reaction':
                if edge[1] in filtered_species:
                    reac_cofactors_id[edge[0]]['products'].append(edge[1])
            elif node_obj_target['type']=='reaction':
                if edge[0] in filtered_species:
                    reac_cofactors_id[edge[1]]['substrates'].append(edge[0])
        #############  Add the parent nodes first
        pos = {}
        parent_layer = [root]
        y_layer = 0.0
        self.logger.debug('\t'+str(root)+' --> '+str(list(G.predecessors(root))+list(G.successors(root))))
        for nei in list(G.predecessors(root))+list(G.successors(root)):
            self.logger.debug('\t'+str(nei)+' predecessors: '+str(list(G.predecessors(nei))))
            if G.node.get(nei)['type']=='reaction':
                for reac_nei_suc in list(G.successors(nei)):
                    if reac_nei_suc in toadd_nodes and reac_nei_suc not in parent_layer:
                        parent_layer.append(reac_nei_suc)
        self.logger.debug('parent_layer: '+str(parent_layer))
        if parent_layer==[]:
            self.logger.warning('parent_layer is empty')
            return False
        dx = width/len(parent_layer)
        nextx = xcenter-width/2-dx/2
        for l in parent_layer:
            nextx += dx
            pos[l] = (nextx, y_layer)
            toadd_nodes.remove(l)
        #toadd_nodes.remove(root)
        y_layer = -y_gap
        layer_num = 1
        while not len(toadd_nodes)==0:
            self.logger.debug('==================================')
            self.logger.debug('toadd_nodes: '+str(toadd_nodes))
            self.logger.debug('parent_layer: '+str(parent_layer))
            layer = _make_layer(parent_layer)
            #layer = list(set(layer))
            self.logger.debug('layer: '+str(layer))
            if layer==[]:
                self.logger.warning('layer is empty')
                break
            dx = width/len(layer)
            nextx = xcenter-width/2-dx/2
            for l in layer:
                nextx += dx
                pos[l] = (nextx, y_layer)
                toadd_nodes.remove(l)
            #set for next loop
            self.logger.debug('pos: '+str(pos))
            layer_num += 1
            y_layer -= y_gap
            parent_layer = layer
        plotG = G.copy()
        for node_id in list(plotG.nodes):
            if node_id not in list(pos.keys()):
                plotG.remove_node(node_id)
        #normalise the values to between 0 and 1
        #round the pos since we can have division errors
        self.logger.debug('----------------------------')
        self.logger.debug('pos: '+str(pos))
        all_x = []
        all_y = []
        for node in pos:
            all_x.append(pos[node][0])
            all_y.append(pos[node][1])
        for node in pos:
            #convert from up/down to right/left
            #reverse the x and y locations
            try:
                x = round((pos[node][1]-min(all_y))/(max(all_y)-min(all_y)), 5)
            except ZeroDivisionError:
                x = 0
            try:
                y = round((pos[node][0]-min(all_x))/(max(all_x)-min(all_x)), 5)
            except ZeroDivisionError:
                y = 0
            pos[node] = (x, y)
        self.logger.debug('pos: '+str(pos))
        ##### TODO: need to adjust the position layers if they are not equidistant
        return plotG, pos, reac_cofactors_id


    ######################################################################################################
    ########################################## Public Functions ##########################################
    ######################################################################################################


    ## Draw the chemicals using their InChI
    #
    # @param 
    #
    #
    def drawChemicalList(self, id_inchi, subplot_size=[200, 200]):
        toRet = {}
        inchi_list = list(set([id_inchi[i] for i in id_inchi]))
        list_mol = [MolFromInchi(inchi) for inchi in inchi_list]
        for i in range(len(list_mol)):
            cp_list_mol = copy.deepcopy(list_mol)
            cp_list_mol.pop(i)
            tmp_list_mol = [list_mol[i]]+cp_list_mol
            img = Draw.MolsToGridImage(tmp_list_mol, molsPerRow=1, subImgSize=(subplot_size[0], subplot_size[1]), useSVG=True)
            #add the groups tag with the id's of the reactions -- should have be size width=subplot_size[0] height=subplot_size[1]*len(list_mol)
            bond_0_count = 0
            svg_str = ''
            for line in img.splitlines():
                add_line = True
                m0 = re.findall("(\d+\.\d+)", line)
                if m0:
                    for y in m0:
                        if float(y)>subplot_size[1]:
                            add_line = False
                m1 = re.findall("height=\'\d+", line)
                if m1:
                    line = re.sub(r"height=\'\d+", "height=\'"+str(subplot_size[1]), line)
                    #line.replace(str(subplot_size[i]*len(list_mol)), str(subplot_size[1]))
                if add_line:
                    svg_str += line+'\n'
            for y in id_inchi:
                if id_inchi[y]==inchi_list[i]:
                    toRet[y] = svg_str
        return toRet



    def _iscross(self, edge1, edge2):
        """Determine if the line
        """
        #Starting points or end points are the same
        if edge1['source'][0]==edge2['source'][0] and edge1['source'][1]==edge2['source'][1]:
            return True
        if edge1['target'][0]==edge2['target'][0] and edge1['target'][1]==edge2['target'][1]:
            return True
        #if the breakpoints are the same
        if edge1['L1'][0]==edge2['L1'][0] and edge1['L1'][1]==edge2['L1'][1]:
            return True
        if edge1['L2'][0]==edge2['L2'][0] and edge1['L2'][1]==edge2['L2'][1]:
            return True
        if edge1['L2'][1]==edge2['L2'][1]:
            return True


    def line_intersection(line1, line2):
        """ Taken from https://stackoverflow.com/questions/20677795/how-do-i-compute-the-intersection-point-of-two-lines
            Determine the location of an intersection between two lines
        """
        xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
        ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])
        def det(a, b):
            return a[0] * b[1] - a[1] * b[0]
        div = det(xdiff, ydiff)
        if div == 0:
           #raise Exception('lines do not intersect')
           return None, None
        d = (det(*line1), det(*line2))
        x = det(d, xdiff) / div
        y = det(d, ydiff) / div
        return x, y


    def drawsvg(self, rpgraph,
                target,
                subplot_size=[200,200],
                #reac_size=[20,60],
                reac_fill_color='#ddd',
                reac_stroke_color='black',
                reac_stroke_width=2,
                arrow_stroke_color='black',
                arrow_stroke_width=2,
                font_family='sans-serif',
                font_size=10,
                font_color='black',
                plot_only_central=True,
                filter_cofactors=True,
                filter_sink_species=False):
        """Generate a reaction SVG image from the rpgraph object.

        :param rpgraph: rpGraph object to draw the SVG from
        :param target: source node to calculate the hierarchy tree organisation of the reaction (should be TARGET)
        :param suboplot_size: The size in pixels of the subplot boxes used to draw the SVG  (default: [200, 200])
        :param reac_fill_color: Hex (or name) color code to fill of the reaction box (default: '#ddd')
        :param reac_stroke_color: Hex (or name) color code of the reaction box stroke (default: 'black')
        :param reac_stroke_width: Size of the reaction rectangle stroke width (default: 2)
        :param arrow_stroke_color: Hex (or name) color code of the reaction arrows (default: 'black')
        :param arrow_stroke_width: Size of the reaction arrows (default: 2)
        :param font_family: The font of the cofactors (default: 'sans-serif'
        :param font_size: The font size of the cofactors (default: 10)
        :param font_color: The font color of the cofactors (default: 'black')
        :param plot_only_central: Do not draw the chemical structure of the non-central species (default: True)
        :param filter_cofactors: Do not draw the chemical structire of the identified cofactors (see: data/mnx_cofactors.json) (default: True)
        :param filter_sink_species: Do not draw the chemical structure of sink species (default: False)

        :type rpgraph: rpGraph
        :type target: str
        :type suboplot_size: list
        :type reac_fill_color: str
        :type reac_stroke_color: str
        :type reac_stroke_width: int
        :type arrow_stroke_color: str
        :type arrow_stroke_width: int
        :type font_family: str
        :type font_size: int
        :type font_color: str
        :type plot_only_central: bool
        :type filter_cofactors: bool
        :type filter_sink_species: bool

        :returns: tuple (svg, resG, mod_pos, reac_cofactors_name)
            - svg - SVG as string
            - regG - Result networkx object with the cofactors removed
            - mod_pos - The calculates positions for the objects
            - reac_cofactors_name - Dictionnary of reactions with the subtrates and products ID's

        :rtype: tuple
        """
        #TODO: Check this one: /Users/melchior/Downloads/rpglobalscore_77/rp_109_2.sbml.xml
        reac_size = [subplot_size[0]/8, subplot_size[1]/2]
        #gather all the inchis and convert to svg
        resG, pos, reac_cofactors_id = self._hierarchy_pos(rpgraph.G,
                                                           target,
                                                           plot_only_central=plot_only_central,
                                                           filter_cofactors=filter_cofactors,
                                                           filter_sink_species=filter_sink_species)
        self.logger.debug('+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=')
        ######## convert the id's to string name #####
        reac_cofactors_name = {}
        for reac_id in reac_cofactors_id:
            if not reac_id in reac_cofactors_name:
                reac_cofactors_name[reac_id] = {'substrates': [], 'products': []}
            for sub in reac_cofactors_id[reac_id]['substrates']:
                try:
                    name = rpgraph.G.node.get(sub)['name']
                except KeyError:
                    name = sub
                if name=='':
                    name = sub
                reac_cofactors_name[reac_id]['substrates'].append(name)
            for pro in reac_cofactors_id[reac_id]['products']:
                try:
                    name = rpgraph.G.node.get(pro)['name']
                except KeyError:
                    name = pro
                if name=='':
                    name = pro
                reac_cofactors_name[reac_id]['products'].append(name)
        ############# Calculate the size of the image ##### 
        id_inchi = {}
        self.logger.debug('positions: '+str(pos))
        for node in list(resG.nodes):
            if resG.node.get(node)['type']=='species':
                id_inchi[node] = resG.node.get(node)['brsynth']['inchi']
        id_svg = self.drawChemicalList(id_inchi, subplot_size)
        self.logger.debug('============================')
        a = {}
        for n in pos:
            if not pos[n][0] in a:
                a[pos[n][0]] = []
            a[pos[n][0]].append(pos[n][1])
        self.logger.debug('a: '+str(a))
        largest_y = 0
        for i in a:
            if len(a[i])>largest_y:
                largest_y = len(a[i])
        if largest_y==0:
            largest_y = 1
        self.logger.debug('largest_y: '+str(largest_y))
        u_x_layers = sorted(list(set([pos[i][0] for i in pos])))
        u_y_layers = sorted(list(set([pos[i][1] for i in pos])))
        self.logger.debug('u_x_layers: '+str(u_x_layers))
        self.logger.debug('u_y_layers: '+str(u_y_layers))
        background_len_x = subplot_size[0]*len(u_x_layers)
        background_len_y = subplot_size[1]*largest_y
        self.logger.debug('background_len_x: '+str(background_len_x))
        self.logger.debug('background_len_y: '+str(background_len_y))
        mod_pos = {}
        #adjust the x axis lications for the boxes to be close to each other
        #TODO: readjust to equidistant on the y-axis nodes that are not
        for node in pos:
            mod_pos[node] = (pos[node][0]*(subplot_size[0]*(len(u_x_layers)-1)), #not sure why I have to correct that
                             (pos[node][1]*background_len_y))
        self.logger.debug('============================')
        self.logger.debug('mod_pos: '+str(mod_pos))
        ########### draw the background #############
        len_fig_x = background_len_x
        len_fig_y = background_len_y
        self.logger.debug('len_fig_x: '+str(len_fig_x))
        self.logger.debug('len_fig_y: '+str(len_fig_y))
        fig = sg.SVGFigure(str(len_fig_x),
                           str(len_fig_y))
        background_white = draw.Drawing(len_fig_x, len_fig_y, origin=(0,0))
        background_white.append(draw.Rectangle(0, 0, len_fig_x, len_fig_y, fill='#ffffff'))
        a = sg.fromstring(background_white.asSvg())
        b_w = a.getroot()
        b_w.moveto(0, background_len_y)#WARNING: not sure why I have to + subpot
        fig.append(b_w)
        nodes_attach_locs = {}
        for node_id in mod_pos:
            node = rpgraph.G.node.get(node_id)
            self.logger.debug('\tSpecies: '+str(node_id))
            if node['type']=='species':
                self.logger.debug('\tNode pos: '+str(mod_pos[node_id]))
                self.logger.debug('\tx: '+str(mod_pos[node_id][0]))
                self.logger.debug('\ty: '+str(mod_pos[node_id][1]))
                move_x = mod_pos[node_id][0]
                # because of the nature of the x y locations, need to reverse them here
                if mod_pos[node_id][1]==0.0:
                    y_coord = mod_pos[node_id][1]+subplot_size[1]/2
                    move_y = len_fig_y-mod_pos[node_id][1]-subplot_size[1]
                elif mod_pos[node_id][1]==len_fig_y:
                    y_coord = mod_pos[node_id][1]-subplot_size[1]/2
                    move_y = len_fig_y-mod_pos[node_id][1]
                else:
                    y_coord = mod_pos[node_id][1]
                    move_y = len_fig_y-mod_pos[node_id][1]-subplot_size[1]/2
                self.logger.debug('\tmove_x: '+str(move_x))
                self.logger.debug('\tmove_y: '+str(move_y))
                f = sg.fromstring(id_svg[node_id])
                p = f.getroot()
                #Rememeber that you are moving the object on the x and y axis while the coordinates are coordinates so its reversed
                p.moveto(move_x, move_y)
                fig.append(p)
                nodes_attach_locs[node_id] = {'left': (mod_pos[node_id][0], y_coord),
                                              'right': (mod_pos[node_id][0]+subplot_size[0], y_coord)}
            elif node['type']=='reaction':
                d = draw.Drawing(subplot_size[0], subplot_size[1], origin=(0,0))
                d.append(draw.Rectangle(0, 0, subplot_size[0], subplot_size[1], fill='#FFFFFF'))
                edge_x = subplot_size[0]/2-reac_size[1]/2
                edge_y = subplot_size[1]/2-reac_size[0]/2
                self.logger.debug('\tedge_x: '+str(edge_x))
                self.logger.debug('\tedge_y: '+str(edge_y))
                self.logger.debug('\tx: '+str(mod_pos[node_id][0]))
                self.logger.debug('\ty: '+str(mod_pos[node_id][1]+subplot_size[1]))
                if reac_cofactors_name[node_id]['substrates']:
                    d.append(draw.Line(subplot_size[0]/2, edge_y-reac_size[0],
                                       subplot_size[0]/2, edge_y-self.arrowhead_comp_x,
                                       stroke=arrow_stroke_color, stroke_width=arrow_stroke_width, fill='none', marker_end=self.arrowhead))
                if reac_cofactors_name[node_id]['products']:
                    d.append(draw.Line(subplot_size[0]/2, edge_y+reac_size[0],
                                       subplot_size[0]/2, (subplot_size[1]/2)+reac_size[1]/3,
                                       stroke=arrow_stroke_color, stroke_width=arrow_stroke_width, fill='none', marker_end=self.arrowhead))
                y_shift = 0.0
                for sub in reac_cofactors_name[node_id]['substrates']:
                    self.logger.debug(sub)
                    d.append(draw.Text(sub, font_size,
                                       subplot_size[0]/2, ((subplot_size[1]/2)-reac_size[1]/3)-y_shift-font_size,
                                       font_family=font_family, center=True, fill=font_color))
                    y_shift += font_size
                y_shift = 0.0
                for pro in reac_cofactors_name[node_id]['products']:
                    self.logger.debug(pro)
                    d.append(draw.Text(pro, font_size,
                                       subplot_size[0]/2, ((subplot_size[1]/2)+reac_size[1]/3)+y_shift+font_size+self.arrowhead_comp_x,
                                       font_family=font_family, center=True, fill=font_color))
                    y_shift += font_size
                d.append(draw.Rectangle(edge_x, edge_y,
                                        reac_size[1], reac_size[0],
                                        fill=reac_fill_color,
                                        stroke_width=reac_stroke_width,
                                        stroke=reac_stroke_color))
                a = sg.fromstring(d.asSvg())
                a_r = a.getroot()
                move_x = mod_pos[node_id][0]
                if mod_pos[node_id][1]==0.0:
                    move_y = len_fig_y-mod_pos[node_id][1]
                elif mod_pos[node_id][1]==len_fig_y:
                    move_y = len_fig_y-mod_pos[node_id][1]+subplot_size[1]
                else:
                    move_y = len_fig_y-mod_pos[node_id][1]-subplot_size[1]/2+subplot_size[1]
                self.logger.debug('\tmove_x: '+str(move_x))
                self.logger.debug('\tmove_y: '+str(move_y))
                a_r.moveto(move_x, move_y)
                fig.append(a_r)
                nodes_attach_locs[node_id] = {'left': (move_x+edge_x,
                                                       move_y-subplot_size[1]/2),
                                              'right': (move_x+subplot_size[0]-edge_x+reac_stroke_width/2,
                                                        move_y-subplot_size[1]/2)}
            self.logger.debug('\t-------------------------------')
        self.logger.debug('nodes_attach_locs: '+str(nodes_attach_locs))
        self.logger.debug(list(resG.edges))
        arrow_box = draw.Drawing(len_fig_x, len_fig_y, origin=(0,0))
        ################ Add the arrowhead depending on edge direction #######
        #depending on the directions of the node, switch the source and the target
        #and calculate the center pocitions of the arrows
        edge_pos = {}
        strict_edge_pos = {}
        for edge in list(resG.edges):
            #left to right
            if pos[edge[0]][0]>pos[edge[1]][0]:
                source_x = nodes_attach_locs[edge[0]]['left'][0]
                source_y = nodes_attach_locs[edge[0]]['left'][1]
                target_x = nodes_attach_locs[edge[1]]['right'][0]
                target_y = nodes_attach_locs[edge[1]]['right'][1]
                edge_pos[edge] = {'source': (round(source_x, 2), round(source_y, 2)),
                                  'L1': (round(source_x+(target_x-source_x)/2, 2), round(source_y, 2)),
                                  'L2': (round(source_x+(target_x-source_x)/2, 2), round(target_y, 2)),
                                  'target': (round(target_x+self.arrowhead_comp_x, 2), round(target_y, 2)),
                                  'arrow_direction': self.rev_arrowhead_flat}
                strict_edge_pos[edge] = {'source': (int(round(source_x, 0)), int(round(source_y, 0))),
                                  'L1': (int(round(source_x+(target_x-source_x)/2, 0)), int(round(source_y, 0))),
                                  'L2': (int(round(source_x+(target_x-source_x)/2, 0)), int(round(target_y, 0))),
                                  'target': (int(round(target_x, 0)), int(round(target_y, 0))),
                                  'arrow_direction': self.rev_arrowhead_flat}
                '''
                edge_pos[edge] = {'source': [source_x, source_y],
                                  'L1': [source_x+(target_x-source_x)/2+self.arrowhead_comp_y/2, source_y],
                                  'L2': [source_x+(target_x-source_x)/2+self.arrowhead_comp_y/2, target_y],
                                  'target': [target_x+self.arrowhead_comp_x, target_y],
                                  'arrow_direction': self.rev_arrowhead_flat}
                '''
            #right to left
            elif pos[edge[0]][0]<pos[edge[1]][0]:
                source_x = nodes_attach_locs[edge[0]]['right'][0]
                source_y = nodes_attach_locs[edge[0]]['right'][1]
                target_x = nodes_attach_locs[edge[1]]['left'][0]
                target_y = nodes_attach_locs[edge[1]]['left'][1]
                edge_pos[edge] = {'source': (round(source_x, 2), round(source_y, 2)),
                                  'L1': (round(source_x+(target_x-source_x)/2, 2), round(source_y, 2)),
                                  'L2': (round(source_x+(target_x-source_x)/2, 2), round(target_y, 2)),
                                  'target': (round(target_x-self.arrowhead_comp_x, 2), round(target_y, 2)),
                                  'arrow_direction': self.arrowhead_flat}
                strict_edge_pos[edge] = {'source': (int(round(source_x, 0)), int(round(source_y, 0))),
                                  'L1': (int(round(source_x+(target_x-source_x)/2, 0)), int(round(source_y, 0))),
                                  'L2': (int(round(source_x+(target_x-source_x)/2, 0)), int(round(target_y, 0))),
                                  'target': (int(round(target_x, 0)), int(round(target_y, 0))),
                                  'arrow_direction': self.arrowhead_flat}
                '''
                edge_pos[edge] = {'source': [source_x, source_y],
                                  'L1': [source_x+(target_x-source_x)/2-self.arrowhead_comp_y/2, source_y],
                                  'L2': [source_x+(target_x-source_x)/2-self.arrowhead_comp_y/2, target_y],
                                  'target': [target_x-self.arrowhead_comp_x, target_y],
                                  'arrow_direction': self.arrowhead_flat}
                '''
            else:
                self.logger.error('Cannot connect same y-axi')
        #calculate the center positions of the arrows
        self.logger.debug('edge_pos: '+str(edge_pos))
        self.logger.debug('strict_edge_pos: '+str(strict_edge_pos))
        ############# Calculate the overlaps ##########
        #find the edges that have the same source/target locations - if more than that do not go in the same direction
        # then create a input/output location and update the positions
        #NOTE: have to make strings from the edge locations to be able to be
        overlaps_arrow = {'node': {}, 'L': {}}
        #overlaps_edge = {'node': {}, 'L': {}}
        for edge in strict_edge_pos:
            source_id = str(strict_edge_pos[edge]['source'][0])+'-'+str(strict_edge_pos[edge]['source'][1])
            target_id = str(strict_edge_pos[edge]['target'][0])+'-'+str(strict_edge_pos[edge]['target'][1])
            l1_id = str(strict_edge_pos[edge]['L1'][0])+'-'+str(strict_edge_pos[edge]['L1'][1])
            l2_id = str(strict_edge_pos[edge]['L2'][0])+'-'+str(strict_edge_pos[edge]['L1'][1])
            ####### make the nodes and arrow overlaps_arrow #######
            if not source_id in overlaps_arrow['node']:
                overlaps_arrow['node'][source_id] = [strict_edge_pos[edge]['arrow_direction']]
            else:
                overlaps_arrow['node'][source_id].append(strict_edge_pos[edge]['arrow_direction'])
            if not target_id in overlaps_arrow['node']:
                overlaps_arrow['node'][target_id] = [strict_edge_pos[edge]['arrow_direction']]
            else:
                overlaps_arrow['node'][target_id].append(strict_edge_pos[edge]['arrow_direction'])
            if not l1_id in overlaps_arrow['L']:
                overlaps_arrow['L'][l1_id] = [strict_edge_pos[edge]['arrow_direction']]
            else:
                overlaps_arrow['L'][l1_id].append(strict_edge_pos[edge]['arrow_direction'])
            if not l2_id in overlaps_arrow['L']:
                overlaps_arrow['L'][l2_id] = [strict_edge_pos[edge]['arrow_direction']]
            else:
                overlaps_arrow['L'][l2_id].append(strict_edge_pos[edge]['arrow_direction'])
            ##### make the overlap edge ####
            """
            if not source_id in overlaps_edge['node']:
                overlaps_edge['node'][source_id] = [edge]
            else:
                overlaps_edge['node'][source_id].append(edge)
            if not target_id in overlaps_edge['node']:
                overlaps_edge['node'][target_id] = [edge]
            else:
                overlaps_edge['node'][target_id].append(edge)
            if not l1_id in overlaps_edge['L']:
                overlaps_edge['L'][l1_id] = [edge]
            else:
                overlaps_edge['L'][l1_id].append(edge)
            if not l2_id in overlaps_edge['L']:
                overlaps_edge['L'][l2_id] = [edge]
            else:
                overlaps_edge['L'][l2_id].append(edge)
            """
        ########## Add entry/exit of node if same side node has multiple types #######
        #adjust the perpendecular arrows if there is overlap with reversed directions
        #TODO: adjust arrows that overlap in the same direction but do not have the same destination
        self.logger.debug('overlaps_arrow: '+str(overlaps_arrow))
        self.logger.debug('overlaps_edge: '+str(overlaps_edge))
        for pos_id in overlaps_arrow['node']:
            #if the direction of the node locations are not the same then you need seperate the input/output
            if len(overlaps_arrow['node'][pos_id])>1:
                if not overlaps_arrow['node'][pos_id].count(self.arrowhead_flat)==len(overlaps_arrow['node'][pos_id]) or overlaps_arrow['node'][pos_id].count(self.rev_arrowhead_flat)==len(overlaps_arrow['node'][pos_id]):
                    for edge in strict_edge_pos:
                        source_id = str(strict_edge_pos[edge]['source'][0])+'-'+str(strict_edge_pos[edge]['source'][1])
                        target_id = str(strict_edge_pos[edge]['target'][0])+'-'+str(strict_edge_pos[edge]['target'][1])
                        if source_id==pos_id and strict_edge_pos[edge]['arrow_direction']==self.arrowhead_flat:
                            edge_pos[edge]['source'] = (edge_pos[edge]['source'][0], edge_pos[edge]['source'][1]-self.arrowhead_comp_y/2)
                            edge_pos[edge]['L1'] = (edge_pos[edge]['L1'][0], edge_pos[edge]['L1'][1]-self.arrowhead_comp_y/2)
                        elif source_id==pos_id and strict_edge_pos[edge]['arrow_direction']==self.rev_arrowhead_flat:
                            edge_pos[edge]['source'] = (edge_pos[edge]['source'][0], edge_pos[edge]['source'][1]+self.arrowhead_comp_y/2)
                            edge_pos[edge]['L1'] = (edge_pos[edge]['L1'][0], edge_pos[edge]['L1'][1]+self.arrowhead_comp_y/2)
                        if target_id==pos_id and strict_edge_pos[edge]['arrow_direction']==self.arrowhead_flat:
                            edge_pos[edge]['target'] = (edge_pos[edge]['target'][0], edge_pos[edge]['target'][1]-self.arrowhead_comp_y/2)
                            edge_pos[edge]['L2'] = (edge_pos[edge]['L2'][0], edge_pos[edge]['L2'][1]-self.arrowhead_comp_y/2)
                        elif target_id==pos_id and strict_edge_pos[edge]['arrow_direction']==self.rev_arrowhead_flat:
                            edge_pos[edge]['target'] = (edge_pos[edge]['target'][0], edge_pos[edge]['target'][1]+self.arrowhead_comp_y/2)
                            edge_pos[edge]['L2'] = (edge_pos[edge]['L2'][0], edge_pos[edge]['L2'][1]+self.arrowhead_comp_y/2)
        #TODO: problem of overlap of arrows
        #1) loop through all the arrows in the same layer
        #2) for each node and each entry/exit, record the perpendicular locations
        #3) if two overlap when they should not then add y_shift
        #           - if they are entry/exit
        #           - if they are entry or exit that overlap with other entry/exit from another reaction
        #do the same for the the perpendicular
        '''
        for pos_id in overlaps_arrow['edge']:
        for edge in strict_edge_pos:
            source_id = str(strict_edge_pos[edge]['source'][0])+'-'+str(strict_edge_pos[edge]['source'][1])
            target_id = str(strict_edge_pos[edge]['target'][0])+'-'+str(strict_edge_pos[edge]['target'][1])
        '''
        #calculate the perpendicular overlaps
        for edge in edge_pos:
            perpendicular_layers = []
            for comp_edge in edge_pos:
                #if edge_pos[comp_edge]['L1']==edge_pos[edge]['L1']:
                x, y = line_intersection(edge_pos[comp_edge]['L1'][0], edge_pos[comp_edge]['L1'][1], edge_pos[edge]['L1'][0], edge_pos[edge]['L1'][1])
                if not x==None and y==None:
                    perpendicular_layers.append(comp_edge)
            if len(perpendicular_layers)>1:
                #cases to ignore:
                #   - when they overlap but go in the same direction and the same target or source
                #case when overlap but they go to the same direction and not the same target
                #case when 


        for pos_id in overlaps_arrow['L']:
            if len(overlaps_arrow['L'][pos_id])>1:
                #TODO: need a better overlap algo that takes into consideration if they do not overlap
                #example /Users/melchior/Downloads/rpglobalscore_101/rp_1_1.sbml.xml
                #BUG: /Users/melchior/Downloads/rpglobalscore_101/rp_3_2.sbml.xml --> line not detected to be moved
                #BUG: three perpendicular lines drawn: /Users/melchior/Downloads/rpglobalscore_101/rp_2_2.sbml.xml
                #BUG: HO being shown, should be cofactor: /Users/melchior/Downloads/rpglobalscore_91/rp_2_1.sbml.xml
                #BUG: this one should not shift perpendivular lines: /Users/melchior/Downloads/rpglobalscore_101/rp_12_1.sbml.xml
                if not overlaps_arrow['L'][pos_id].count(self.arrowhead_flat)==len(overlaps_arrow['L'][pos_id]) or overlaps_arrow['L'][pos_id].count(self.rev_arrowhead_flat)==len(overlaps_arrow['L'][pos_id]):
                    '''Close but removes some that should have seperated
                    #ignore cases where there are only 2 that go in opposite direction and actually never meet
                    if len(overlaps_arrow['L'][pos_id])==2:
                        if not strict_edge_pos[overlaps_edge['L'][pos_id][0]]['L1'][1]-strict_edge_pos[overlaps_edge['L'][pos_id][0]]['L2'][1]>0 and not strict_edge_pos[overlaps_edge['L'][pos_id][1]]['L1'][1]-strict_edge_pos[overlaps_edge['L'][pos_id][1]]['L2'][1]>0:
                            continue
                        if not strict_edge_pos[overlaps_edge['L'][pos_id][0]]['L1'][1]-strict_edge_pos[overlaps_edge['L'][pos_id][0]]['L2'][1]<0 and not strict_edge_pos[overlaps_edge['L'][pos_id][1]]['L1'][1]-strict_edge_pos[overlaps_edge['L'][pos_id][1]]['L2'][1]<0:
                            continue
                    '''
                    #TODO: need to detect if there are any criss-cross of perpendecular arrows at their target of source to seperate in either direction
                    for edge in strict_edge_pos:
                        l1_id = str(strict_edge_pos[edge]['L1'][0])+'-'+str(strict_edge_pos[edge]['L1'][1])
                        l2_id = str(strict_edge_pos[edge]['L2'][0])+'-'+str(strict_edge_pos[edge]['L1'][1])
                        if l1_id==pos_id and strict_edge_pos[edge]['arrow_direction']==self.arrowhead_flat:
                            '''
                            edge_pos[edge]['L1'] = (edge_pos[edge]['L1'][0]-self.arrowhead_comp_y/2, edge_pos[edge]['L1'][1])
                            edge_pos[edge]['L2'] = (edge_pos[edge]['L2'][0]-self.arrowhead_comp_y/2, edge_pos[edge]['L2'][1])
                            '''
                            edge_pos[edge]['L1'] = (edge_pos[edge]['L1'][0]+self.arrowhead_comp_y/2, edge_pos[edge]['L1'][1])
                            edge_pos[edge]['L2'] = (edge_pos[edge]['L2'][0]+self.arrowhead_comp_y/2, edge_pos[edge]['L2'][1])
                        elif l1_id==pos_id and strict_edge_pos[edge]['arrow_direction']==self.rev_arrowhead_flat:
                            '''
                            edge_pos[edge]['L1'] = (edge_pos[edge]['L1'][0]+self.arrowhead_comp_y/2, edge_pos[edge]['L1'][1])
                            edge_pos[edge]['L2'] = (edge_pos[edge]['L2'][0]+self.arrowhead_comp_y/2, edge_pos[edge]['L2'][1])
                            '''
                            edge_pos[edge]['L1'] = (edge_pos[edge]['L1'][0]-self.arrowhead_comp_y/2, edge_pos[edge]['L1'][1])
                            edge_pos[edge]['L2'] = (edge_pos[edge]['L2'][0]-self.arrowhead_comp_y/2, edge_pos[edge]['L2'][1])
                        if l2_id==pos_id and strict_edge_pos[edge]['arrow_direction']==self.arrowhead_flat:
                            '''
                            edge_pos[edge]['L1'] = (edge_pos[edge]['L1'][0]-self.arrowhead_comp_y/2, edge_pos[edge]['L1'][1])
                            edge_pos[edge]['L2'] = (edge_pos[edge]['L2'][0]-self.arrowhead_comp_y/2, edge_pos[edge]['L2'][1])
                            '''
                            edge_pos[edge]['L1'] = (edge_pos[edge]['L1'][0]+self.arrowhead_comp_y/2, edge_pos[edge]['L1'][1])
                            edge_pos[edge]['L2'] = (edge_pos[edge]['L2'][0]+self.arrowhead_comp_y/2, edge_pos[edge]['L2'][1])
                        elif l2_id==pos_id and strict_edge_pos[edge]['arrow_direction']==self.rev_arrowhead_flat:
                            '''
                            edge_pos[edge]['L1'] = (edge_pos[edge]['L1'][0]+self.arrowhead_comp_y/2, edge_pos[edge]['L1'][1])
                            edge_pos[edge]['L2'] = (edge_pos[edge]['L2'][0]+self.arrowhead_comp_y/2, edge_pos[edge]['L2'][1])
                            '''
                            edge_pos[edge]['L1'] = (edge_pos[edge]['L1'][0]-self.arrowhead_comp_y/2, edge_pos[edge]['L1'][1])
                            edge_pos[edge]['L2'] = (edge_pos[edge]['L2'][0]-self.arrowhead_comp_y/2, edge_pos[edge]['L2'][1])
        ############## Finally draw ##################
        for edge in edge_pos:
            p = draw.Path(stroke=arrow_stroke_color,
                          stroke_width=arrow_stroke_width,
                          fill='none',
                          marker_end=edge_pos[edge]['arrow_direction'])
            p.M(edge_pos[edge]['source'][0], edge_pos[edge]['source'][1]).L(edge_pos[edge]['L1'][0], edge_pos[edge]['L1'][1]).L(edge_pos[edge]['L2'][0], edge_pos[edge]['L2'][1]).L(edge_pos[edge]['target'][0], edge_pos[edge]['target'][1])
            arrow_box.append(p)
        a = sg.fromstring(arrow_box.asSvg())
        a_b = a.getroot()
        a_b.moveto(0, background_len_y)#WARNING: not sure why I have to + subpot
        fig.append(a_b)
        svg = fig.to_str().decode("utf-8")
        return svg, resG, mod_pos, reac_cofactors_name, nodes_attach_locs


