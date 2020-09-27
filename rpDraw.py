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


    ##
    #
    # We do this to keep the proportions between the different molecules
    #
    def _drawChemicalList(self, inchi_list, subplot_size=[200, 200]):
        toRet = {}
        inchi_list = list(set(inchi_list))
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
            toRet[inchi_list[i]] = svg_str
        return toRet


    ##
    #
    # Tree diagram for the reactions
    #
    def _drawReactionArrows(self, mol_list, left_to_right=True, gap_size=100, subplot_size=[200, 200], stroke_color='black', stroke_width=2):
        x_len = gap_size
        y_len = len(mol_list)*subplot_size[1]
        d = draw.Drawing(x_len, y_len, origin=(0,0))
        self.logger.debug('############ _drawReactionArrows ##########')
        self.logger.debug('x_len: '+str(x_len))
        self.logger.debug('y_len: '+str(y_len))
        #white rectangle
        d.append(draw.Rectangle(0, 0, gap_size+subplot_size[0], len(mol_list)*subplot_size[1], fill='#FFFFFF'))
        self.logger.debug('rectangle: '+str(0)+','+str(0))
        self.logger.debug('rectangle: '+str(gap_size+subplot_size[0])+','+str(len(mol_list)*subplot_size[1]))
        #calculate the y starts of the lines
        if left_to_right:
            x_species = 0
            self.logger.debug('x_species: '+str(x_species))
            y_species = [y_len-(i*subplot_size[1])-(subplot_size[1]/2) for i in range(len(mol_list)) if mol_list[i]]
            #calculate the y end of the lines
            x_reaction = gap_size
            self.logger.debug('x_reaction: '+str(x_reaction))
            y_reaction = subplot_size[1]*len(mol_list)/2
            self.logger.debug('y_reaction: '+str(y_reaction))
            for y_spe in y_species:
                self.logger.debug('y_spe: '+str(y_spe))
                p = draw.Path(stroke=stroke_color, stroke_width=stroke_width, fill='transparent', marker_end=self.arrowhead_flat)
                p.M(x_species, y_spe).C(gap_size+x_species, y_spe,
                        x_species, y_reaction,
                        x_reaction-self.arrowhead_comp_x, y_reaction)
                d.append(p)
        else:
            x_species = gap_size
            self.logger.debug('x_species: '+str(x_species))
            y_species = [y_len-(i*subplot_size[1])-(subplot_size[1]/2) for i in range(len(mol_list)) if mol_list[i]]
            #calculate the y end of the lines
            x_reaction = 0
            self.logger.debug('x_reaction: '+str(x_reaction))
            y_reaction = subplot_size[1]*len(mol_list)/2
            self.logger.debug('y_reaction: '+str(y_reaction))
            for y_spe in y_species:
                self.logger.debug('y_spe: '+str(y_spe))
                p = draw.Path(stroke=stroke_color, stroke_width=stroke_width, fill='transparent', marker_start=self.rev_arrowhead)
                p.M(x_species-self.arrowhead_comp_x, y_spe).C(gap_size-x_species, y_spe,
                        x_species, y_reaction,
                        x_reaction, y_reaction)
                d.append(p)
        #d.saveSvg('test_arrow.svg')
        return d.asSvg(), x_len, y_len


    # make sure that inchi_list is ordered in such a way that [1,2,3] -> [1,
    #                                                                     2,
    #                                                                     3]
    #
    #
    def _drawPartReaction(self,
                          mol_list,
                          arrow_list,
                          left_to_right=True,
                          gap_size=100,
                          subplot_size=[200, 200],
                          is_inchi=True,
                          draw_mol=True,
                          stroke_color='black',
                          stroke_width=2):
        self.logger.debug('---- _drawPartReaction ----')
        assert len(mol_list)==len(arrow_list)
        if draw_mol:
            x_len = subplot_size[0]+gap_size
        else:
            x_len = gap_size
        self.logger.debug('x_len: '+str(x_len))
        y_len = len(mol_list)*subplot_size[1]
        self.logger.debug('y_len: '+str(y_len))
        fig = sg.SVGFigure(str(x_len), str(y_len))
        if is_inchi:
            mol_dict = self._drawChemicalList(mol_list, subplot_size)
            mol_list = [mol_dict[i] for i in mol_list]
        if left_to_right:
            count = 0
            for svg in mol_list:
                if svg:
                    f = sg.fromstring(svg)
                    p = f.getroot()
                    p.moveto(0, subplot_size[1]*count)
                    if draw_mol:
                        fig.append(p)
                count += 1
            lines, arrow_len_x, arrow_len_y = self._drawReactionArrows(arrow_list, left_to_right, gap_size, subplot_size, stroke_color=stroke_color, stroke_width=stroke_width)
            line = sg.fromstring(lines)
            line_p = line.getroot()
            self.logger.debug('Move line x: '+str(subplot_size[0]))
            self.logger.debug('Move line y: '+str(subplot_size[1]*len(mol_list)))
            line_p.moveto(x_len-gap_size, subplot_size[1]*len(mol_list))
            fig.append(line_p)
        else:
            lines, arrow_len_x, arrow_len_y = self._drawReactionArrows(arrow_list, left_to_right, gap_size, subplot_size, stroke_color=stroke_color, stroke_width=stroke_width)
            line = sg.fromstring(lines)
            line_p = line.getroot()
            self.logger.debug('Move line y: '+str(subplot_size[1]*len(mol_list)))
            line_p.moveto(0, subplot_size[1]*len(mol_list))
            fig.append(line_p)
            count = 0
            for svg in mol_list:
                f = sg.fromstring(svg)
                p = f.getroot()
                p.moveto(gap_size, subplot_size[1]*count)
                count += 1
                if draw_mol:
                    fig.append(p)
        return fig.to_str().decode("utf-8"), x_len, y_len


    ## draw the cofactors with a box and a line going through
    #
    #
    def _drawCofactors(self,
                       cofactor_reactants,
                       cofactor_products,
                       x_len=100,
                       y_len=200,
                       rec_size_y=20,
                       font_family='sans-serif',
                       font_size=10,
                       font_color='black',
                       stroke_color='black',
                       fill_color='#ddd',
                       stroke_width=2):
        d = draw.Drawing(x_len, y_len, origin=(0,0))
        d.append(draw.Rectangle(0, 0, x_len, y_len, fill='#FFFFFF'))
        #draw the cofactors arrow
        start_line_x = x_len/2
        start_line_y = (y_len/2-rec_size_y)-rec_size_y
        end_line_x = x_len/2
        end_line_y = (y_len/2+rec_size_y)+rec_size_y
        self.logger.debug('start_line_x: '+str(start_line_x))
        self.logger.debug('start_line_y: '+str(start_line_y))
        self.logger.debug('end_line_x: '+str(end_line_x))
        self.logger.debug('end_line_y: '+str(end_line_y))
        if cofactor_reactants or cofactor_products:
            #d.append(draw.Line(start_line_x, start_line_y, end_line_x, end_line_y,
            #                    stroke=stroke_color, stroke_width=stroke_width, fill='none'))
            ### arrowhead ##
            p = draw.Path(stroke=stroke_color, stroke_width=stroke_width, fill='transparent', marker_end=self.arrowhead)
            p.M(start_line_x, start_line_y).Q(start_line_x-rec_size_y/1.5, y_len/2,
                    end_line_x, end_line_y-self.arrowhead_comp_y)
            d.append(p)
        #Draw the reaction rectangle
        self.logger.debug('react_x1: '+str(0))
        self.logger.debug('react_y1: '+str(y_len/2+rec_size_y))
        react_width = x_len
        react_height = rec_size_y
        self.logger.debug('react_width: '+str(react_width))
        self.logger.debug('react_height: '+str(react_height))
        d.append(draw.Rectangle(0, y_len/2-rec_size_y/2, react_width, react_height, fill=fill_color, stroke_width=stroke_width, stroke=stroke_color))
        #draw the text
        #reactants
        y_shift = 0
        for react in cofactor_reactants:
            x = start_line_x
            y = start_line_y-font_size-y_shift-5
            self.logger.debug('x: '+str(x))
            self.logger.debug('y: '+str(y))
            d.append(draw.Text(react, font_size, x, y, font_family=font_family, center=0, fill=font_color))
            y_shift += font_size
        #product
        y_shift = 0
        for pro in cofactor_products:
            x = end_line_x
            y = end_line_y+y_shift+5
            self.logger.debug('x: '+str(x))
            self.logger.debug('y: '+str(y))
            d.append(draw.Text(pro, font_size, x, y, font_family=font_family, center=0, fill=font_color))
            y_shift += font_size
        return d.asSvg()


    ## Draw a reaction 
    #
    #
    def _drawReaction(self,
                     reactants,
                     products,
                     react_arrow=[],
                     pro_arrow=[],
                     cofactor_reactants=[],
                     cofactor_products=[],
                     is_inchi=True,
                     react_arrow_size=100,
                     arrow_gap_size=100,
                     subplot_size=[200, 200],
                     draw_left_mol=True,
                     stroke_color='black',
                     stroke_width=2):
        self.logger.debug('========= _drawReaction =========')
        if not react_arrow:
            react_arrow = [True for i in reactants]
        if not pro_arrow:
            pro_arrow = [True for i in products]
        if is_inchi:
            #Need to combine them into single request to keep the same proportions
            inchi_svg_dict = self._drawChemicalList([i for i in reactants if i]+[i for i in products if i], subplot_size)
            svg_left, x_len_left, y_len_left = self._drawPartReaction([inchi_svg_dict[i] for i in reactants],
                                                                     react_arrow,
                                                                     left_to_right=True,
                                                                     gap_size=arrow_gap_size,
                                                                     subplot_size=subplot_size,
                                                                     is_inchi=False,
                                                                     draw_mol=draw_left_mol,
                                                                     stroke_color=stroke_color,
                                                                     stroke_width=stroke_width)
            svg_right, x_len_right, y_len_right = self._drawPartReaction([inchi_svg_dict[i] for i in products],
                                                                        pro_arrow,
                                                                        left_to_right=False,
                                                                        gap_size=arrow_gap_size,
                                                                        subplot_size=subplot_size,
                                                                        is_inchi=False,
                                                                        draw_mol=True,
                                                                        stroke_color=stroke_color,
                                                                        stroke_width=stroke_width)
        else:
            svg_left, x_len_left, y_len_left = self._drawPartReaction(reactants,
                                                                     react_arrow,
                                                                     left_to_right=True,
                                                                     gap_size=arrow_gap_size,
                                                                     subplot_size=subplot_size,
                                                                     is_inchi=False,
                                                                     draw_mol=draw_left_mol,
                                                                     stroke_color=stroke_color,
                                                                     stroke_width=stroke_width)
            svg_right, x_len_right, y_len_right = self._drawPartReaction(products,
                                                                        pro_arrow,
                                                                        left_to_right=False,
                                                                        gap_size=arrow_gap_size,
                                                                        subplot_size=subplot_size,
                                                                        is_inchi=False,
                                                                        draw_mol=True,
                                                                        stroke_color=stroke_color,
                                                                        stroke_width=stroke_width)
        self.logger.debug('x_len_left: '+str(x_len_left))
        self.logger.debug('y_len_left: '+str(y_len_left))
        self.logger.debug('x_len_right: '+str(x_len_right))
        self.logger.debug('y_len_right: '+str(y_len_right))
        y_len = 0
        if y_len_left>y_len_right:
            y_len = y_len_left
        else:
            y_len = y_len_right
        x_len = x_len_left+x_len_right+react_arrow_size
        self.logger.debug('x_len: '+str(x_len))
        self.logger.debug('y_len: '+str(y_len))
        #create the final svg
        fig = sg.SVGFigure(str(x_len), str(y_len))
        f_left = sg.fromstring(svg_left)
        p_left = f_left.getroot()
        p_left.moveto(0, (y_len-y_len_left)/2)
        fig.append(p_left)
        f_right = sg.fromstring(svg_right)
        p_right = f_right.getroot()
        p_right.moveto(x_len_left+react_arrow_size, (y_len-y_len_right)/2)
        fig.append(p_right)
        co_reac_svg = self._drawCofactors(cofactor_reactants,
                                        cofactor_products,
                                        x_len=react_arrow_size,
                                        y_len=y_len)
        f_rect = sg.fromstring(co_reac_svg)
        p_rect = f_rect.getroot()
        #300, 500
        p_rect.moveto(x_len_left, y_len)
        fig.append(p_rect)
        return fig.to_str().decode("utf-8"), x_len, y_len


    ##
    #
    # pathway_list: List of dict where each entry contains a list of reaction inchi, with None for positions and
    def _drawPathway(self,
                    pathway_list,
                    react_arrow_size=100,
                    arrow_gap_size=100,
                    subplot_size=[200, 200],
                    stroke_color='black',
                    stroke_width=2):
        self.logger.debug('________ _drawPathway ________')
        #Need to make single command to keep all in the same proportions
        all_inchis = []
        for reaction in pathway_list:
            for inchi in reaction['reactants_inchi']+reaction['products_inchi']:
                if inchi:
                    all_inchis.append(inchi)
        inchi_svg_dict = self._drawChemicalList(all_inchis, subplot_size)
        ###### if there is only one sinply return the reaction 
        if len(pathway_list)==1:
            svg, x_len, y_len = self._drawReaction([inchi_svg_dict[i] for i in pathway_list[0]['reactants_inchi']],
                                                   [inchi_svg_dict[i] for i in pathway_list[0]['products_inchi']],
                                                   react_arrow=[],
                                                   pro_arrow=[],
                                                   is_inchi=False,
                                                   cofactor_reactants=reaction['cofactor_reactants'],
                                                   cofactor_products=reaction['cofactor_products'],
                                                   react_arrow_size=react_arrow_size,
                                                   arrow_gap_size=arrow_gap_size,
                                                   subplot_size=subplot_size,
                                                   draw_left_mol=True,
                                                   stroke_color=stroke_color,
                                                   stroke_width=stroke_width)
            return svg, x_len, y_len
        is_first = True
        pathway_svgs = []
        count = 0
        prev_inchi = []
        for reaction in pathway_list:
            self.logger.debug('count: '+str(count))
            if count==0: #first reaction
                self.logger.debug('This is the first one: '+str(count))
                react_inchi = [i for i in reaction['reactants_inchi']]
                react_arrow = [True for i in reaction['reactants_inchi']]
                #combine the current products with the next reactants
                prev_inchi = list(set([i for i in reaction['products_inchi']]+pathway_list[count+1]['reactants_inchi']))
                pro_inchi = [i for i in prev_inchi]
                pro_arrow = [True if i in reaction['products_inchi'] else False for i in pro_inchi]
            elif count==len(pathway_list)-1: #last reaction
                self.logger.debug('This is the last one: '+str(count)+' - '+str(len(pathway_list)-1))
                react_inchi = [i for i in prev_inchi]
                react_arrow = [True if i in reaction['reactants_inchi'] else False for i in react_inchi]
                pro_inchi = reaction['products_inchi']
                pro_arrow = [True for i in reaction['products_inchi']]
            else: #middle
                self.logger.debug('This is a middle one: '+str(count)+' - '+str(len(pathway_list)-1))
                react_inchi = [i for i in prev_inchi]
                ''' this should be unessecary since the previous reaction had all the inchi mashed together
                for inchi in reaction['reactants_inchi']:
                    if not inchi in react_inchi:
                        react_inchi.append(inchi)
                '''
                react_arrow = [True if i in reaction['reactants_inchi'] else False for i in react_inchi]
                #combine the current products with the next reactants
                prev_inchi = list(set([i for i in reaction['products_inchi']]+pathway_list[count+1]['reactants_inchi']))
                pro_inchi = prev_inchi
                pro_arrow = [True if i in reaction['products_inchi'] else False for i in pro_inchi]
            self.logger.debug('prev_inchi: '+str(prev_inchi))
            self.logger.debug('react_inchi: '+str(react_inchi))
            self.logger.debug('react_arrow: '+str(react_arrow))
            self.logger.debug('pro_inchi: '+str(react_inchi))
            self.logger.debug('pro_arrow: '+str(react_arrow))
            svg, x_len, y_len = self._drawReaction([inchi_svg_dict[i] for i in react_inchi],
                                                   [inchi_svg_dict[i] for i in pro_inchi],
                                                   react_arrow=react_arrow,
                                                   pro_arrow=pro_arrow,
                                                   is_inchi=False,
                                                   cofactor_reactants=reaction['cofactor_reactants'],
                                                   cofactor_products=reaction['cofactor_products'],
                                                   react_arrow_size=react_arrow_size,
                                                   arrow_gap_size=arrow_gap_size,
                                                   subplot_size=subplot_size,
                                                   draw_left_mol=is_first,
                                                   stroke_color=stroke_color,
                                                   stroke_width=stroke_width)
            pathway_svgs.append({'svg': svg, 'x_len': x_len, 'y_len': y_len})
            if is_first:
                is_first = False
            count += 1
        y_len = 0
        x_len = 0
        for reaction_svg in pathway_svgs:
            x_len += reaction_svg['x_len']
            if reaction_svg['y_len']>y_len:
                y_len = reaction_svg['y_len']
        self.logger.debug('x_len: '+str(x_len))
        self.logger.debug('y_len: '+str(y_len))
        #create the final svg
        fig = sg.SVGFigure(str(x_len), str(y_len))
        cum_x_len = 0
        for reaction_svg in pathway_svgs:
            f = sg.fromstring(reaction_svg['svg'])
            p = f.getroot()
            p.moveto(cum_x_len, (y_len-reaction_svg['y_len'])/2)
            self.logger.debug('cum_x_len: '+str(cum_x_len))
            cum_x_len += reaction_svg['x_len']
            fig.append(p)
        return fig.to_str().decode("utf-8"), x_len, y_len






    """
    ## Organise the heterologous pathway in a tree like structure using recursive function
    # #DEPRECATED? this function works well but the other works better for drawing
    # This method iterates from the TARGET and organises the heterologous in a tree like manner. The 
    # TODO: add the global filter of cofactors to remove the species that are shared among many (ex: H+, O2, etc...)
    def _hierarchy_pos_recursive(self, G, root, width=1.0, vert_gap=0.2, vert_loc=0, xcenter=0.5, plot_only_central=False, filter_cofactors=True):
        global_xcenter = xcenter
        def h_recur(G, root, width=1.0, vert_gap=0.2, vert_loc=0, xcenter=0.5,
                      pos=None, parent=None, parsed=[], saw_first=[], parent_neighbors=[]):
            self.logger.debug('####################### '+str(root)+' #######################')
            self.logger.debug('parsed:\t\t'+str(parsed))
            self.logger.debug('saw_first:\t'+str(saw_first))#G.node.get('MNXM4__64__MNXC3')['central_species']
            if root not in parsed and root not in saw_first:
                parsed.append(root)
                if pos==None:
                    self.logger.debug('x --> '+str(xcenter))
                    self.logger.debug('y --> '+str(vert_loc))
                    pos = {root:(xcenter,vert_loc)}
                else:
                    self.logger.debug('x --> '+str(xcenter))
                    self.logger.debug('y --> '+str(vert_loc))
                    pos[root] = (xcenter, vert_loc)
                #if you plot_only_central then remove non central
                neighbors = []
                for nei in list(G.predecessors(root))+list(G.successors(root)):
                    if not nei in saw_first:
                        node_obj = G.node.get(nei)
                        #filters only apply to species
                        if node_obj['type']=='species':
                            if plot_only_central and not node_obj['central_species']:
                                self.logger.debug('\t'+str(nei)+' is not a central species')
                                continue
                            if filter_cofactors:
                                if 'metanetx' in node_obj['miriam']:
                                    if any([i in list(self.mnx_cofactors.keys()) for i in node_obj['miriam']['metanetx']]):
                                        self.logger.debug('\t'+str(nei)+' is a list cofactor')
                                        continue
                            neighbors.append(nei)
                        else:
                            neighbors.append(nei)
                #neighbors = [i for i in list(G.predecessors(root))+list(G.successors(root)) if i not in saw_first]
                self.logger.debug('neighbors:\t\t'+str(neighbors))
                if parent!=None:
                    try:
                        neighbors.remove(parent)
                    except ValueError:
                        pass
                if len(neighbors)!=0:
                    #remove the species that have been already added
                    layer_neighbors = [i for i in neighbors]
                    for nh in list(set([i for i in saw_first]))+[root]:
                        try:
                            layer_neighbors.remove(nh)
                        except ValueError:
                            pass
                    #if you want to plot only the central species
                    self.logger.debug('layer_neighbors:\t'+str(layer_neighbors))
                    self.logger.debug('parent_neighbors:\t'+str(parent_neighbors))
                    #dx = width/len(neighbors)
                    try:
                        dx = width/len(layer_neighbors)
                    except ZeroDivisionError:
                        dx = width
                    self.logger.debug('dx: '+str(dx))
                    self.logger.debug(G.node.get(root))
                    '''
                    if len(layer_neighbors)==1:
                        try:
                            xcenter = width/len(parent_neighbors)
                        except ZeroDivisionError:
                            pass
                    '''
                    xcenter = 0.5
                    #nextx = xcenter - width/2 - dx/2
                    self.logger.debug('xcenter: '+str(xcenter))
                    self.logger.debug('width: '+str(width))
                    nextx = xcenter - width/2 - dx/2
                    self.logger.debug('nextx: '+str(nextx))
                    for neighbor in neighbors:
                        self.logger.debug('\tneighbor -> '+str(neighbor))
                        pass_saw_first = list(set([i for i in neighbors]+saw_first))
                        #self.logger.debug('\tpass_saw_first: '+str(pass_saw_first))
                        pass_saw_first.remove(neighbor)
                        nextx += dx
                        #self.logger.debug('\tpass_saw_first: '+str(pass_saw_first))
                        #overwrite to have the reactions center
                        if len(neighbors)==1:
                            nextx = global_xcenter
                        self.logger.debug('\twidth: '+str(dx))
                        self.logger.debug('\tvert_gap: '+str(vert_gap))
                        self.logger.debug('\txcenter: '+str(nextx))
                        self.logger.debug('\tvert_loc: '+str(vert_loc-vert_gap))
                        self.logger.debug('==============================================')
                        pos = h_recur(G, neighbor, width=width, vert_gap=vert_gap,
                                            vert_loc=vert_loc-vert_gap, xcenter=nextx, pos=pos,
                                            parent=root, parsed=parsed, saw_first=pass_saw_first,
                                            parent_neighbors=layer_neighbors)
            return pos
        #if we plot only the center species then we must remove the non-central species
        plotG = G.copy()
        pos = h_recur(plotG, root, width=width, vert_gap=vert_gap, vert_loc=vert_loc, xcenter=xcenter)
        for node_id in list(plotG.nodes):
            if node_id not in list(pos.keys()):
                plotG.remove_node(node_id)
        return plotG, pos
        #return G, h_recur(G, root, width=1., vert_gap = 0.2, vert_loc = 0, xcenter = 0.5)
    """



    ## draw the cofactors with a box and a line going through
    #
    #
    def _YYdrawReaction(self,
                      reaction_name,
                      cofactor_reactants,
                      cofactor_products,
                      x_len=100,
                      y_len=200,
                      rec_size_y=20,
                      font_family='sans-serif',
                      font_size=10,
                      font_color='black',
                      stroke_color='black',
                      fill_color='#ddd',
                      stroke_width=2):
        d = draw.Drawing(x_len, y_len, origin=(0,0))
        d.append(draw.Rectangle(0, 0, x_len, y_len, fill='#FFFFFF'))
        d.append(draw.Text(reaction_name, font_size, x_len/2, 0, font_family=font_family, center=0, fill=font_color))
        #draw the cofactors arrow
        start_line_x = x_len/2
        start_line_y = (y_len/2-rec_size_y)-rec_size_y
        end_line_x = x_len/2
        end_line_y = (y_len/2+rec_size_y)+rec_size_y
        self.logger.debug('start_line_x: '+str(start_line_x))
        self.logger.debug('start_line_y: '+str(start_line_y))
        self.logger.debug('end_line_x: '+str(end_line_x))
        self.logger.debug('end_line_y: '+str(end_line_y))
        if cofactor_reactants or cofactor_products:
            ### arrowhead ##
            p = draw.Path(stroke=stroke_color, stroke_width=stroke_width, fill='transparent', marker_end=self.arrowhead)
            p.M(start_line_x, start_line_y).Q(start_line_x-rec_size_y/1.5, y_len/2,
                    end_line_x, end_line_y-self.arrowhead_comp_y)
            d.append(p)
        #Draw the reaction rectangle
        self.logger.debug('react_x1: '+str(0))
        self.logger.debug('react_y1: '+str(y_len/2+rec_size_y))
        react_width = x_len
        react_height = rec_size_y
        self.logger.debug('react_width: '+str(react_width))
        self.logger.debug('react_height: '+str(react_height))
        d.append(draw.Rectangle(0, y_len/2-rec_size_y/2, react_width, react_height, fill=fill_color, stroke_width=stroke_width, stroke=stroke_color))
        #draw the text
        #reactants
        y_shift = 0
        for react in cofactor_reactants:
            x = start_line_x
            y = start_line_y-font_size-y_shift-5
            self.logger.debug('x: '+str(x))
            self.logger.debug('y: '+str(y))
            d.append(draw.Text(react, font_size, x, y, font_family=font_family, center=0, fill=font_color))
            y_shift += font_size
        #product
        y_shift = 0
        for pro in cofactor_products:
            x = end_line_x
            y = end_line_y+y_shift+5
            self.logger.debug('x: '+str(x))
            self.logger.debug('y: '+str(y))
            d.append(draw.Text(pro, font_size, x, y, font_family=font_family, center=0, fill=font_color))
            y_shift += font_size
        return d.asSvg()














    ## Organise the heterologous pathway in a tree like structure
    #
    # This method iterates from the TARGET and organises the heterologous in a tree like manner. The 
    # TODO: add the global filter of cofactors to remove the species that are shared among many (ex: H+, O2, etc...)
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
        return plotG, pos, reac_cofactors_id



    ##
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



    ##
    #
    #
    def drawsvg(self, G,
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
        reac_size = [subplot_size[0]/8, subplot_size[1]/2]
        #gather all the inchis and convert to svg
        resG, pos, reac_cofactors_id = self._hierarchy_pos(G,
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
                    name = G.node.get(sub)['name']
                except KeyError:
                    name = sub
                if name=='':
                    name = sub
                reac_cofactors_name[reac_id]['substrates'].append(name)
            for pro in reac_cofactors_id[reac_id]['products']:
                try:
                    name = G.node.get(pro)['name']
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
        '''
        a = {}
        for n in pos:
            if not pos[n][1] in a:
                a[pos[n][1]] = []
            a[pos[n][1]].append(pos[n][0])
        self.logger.debug('a: '+str(a))
        largest_x = 0
        for i in a:
            if len(a[i])>largest_x:
                largest_x = len(a[i])
        if largest_x==0:
            largest_x = 1
        '''
        u_x_layers = sorted(list(set([pos[i][0] for i in pos])))
        u_y_layers = sorted(list(set([pos[i][1] for i in pos])))



        '''
        layer_x_conv = {}
        count = 1
        for i in u_x_layers:
            if count==1:
                layer_x_conv[i] = 0.0
            else:
                layer_x_conv[i] = count/len(u_x_layers)
            count += 1
        layer_y_conv = {}
        count = 1
        for i in u_y_layers:
            layer_y_conv[i] = count
            count += 1
        '''
        self.logger.debug('u_x_layers: '+str(u_x_layers))
        self.logger.debug('u_y_layers: '+str(u_y_layers))
        #background_len_x = subplot_size[0]*largest_x
        background_len_x = subplot_size[0]*len(u_x_layers)
        #background_len_y = subplot_size[1]*len(u_y_layers)
        background_len_y = subplot_size[1]*largest_y
        self.logger.debug('background_len_x: '+str(background_len_x))
        self.logger.debug('background_len_y: '+str(background_len_y))
        mod_pos = {}
        #adjust the x axis lications for the boxes to be close to each other
        '''
        u_x_layers.remove(1.0)
        adjust_x = min(u_x_layers)
        '''
        for node in pos:
            mod_pos[node] = (pos[node][0]*(subplot_size[0]*(len(u_x_layers)-1)), #not sure why I have to correct that
                             (pos[node][1]*background_len_y))
            '''
            mod_pos[node] = (pos[node][0]*background_len_x*adjust_x, #not sure why I have to correct that
                             (pos[node][1]*background_len_y))
            '''
            '''
            mod_pos[node] = (pos[node][0]/(1/len(u_x_layers))*background_len_x, #not sure why I have to correct that
                             (pos[node][1]/(1/len(u_y_layers))*background_len_y))
            '''
            '''
            mod_pos[node] = (pos[node][0]*background_len_x*adjust_x, #not sure why I have to correct that
                             (pos[node][1]*background_len_y))
            '''
        self.logger.debug('============================')
        self.logger.debug('mod_pos: '+str(mod_pos))
        ########### draw the background #############
        #len_fig_x = (background_len_x+subplot_size[0])*adjust_x #not sure why I have to correct that
        #len_fig_x = background_len_x+subplot_size[0]
        len_fig_x = background_len_x
        #len_fig_y = background_len_y+subplot_size[1]
        len_fig_y = background_len_y
        self.logger.debug('len_fig_x: '+str(len_fig_x))
        self.logger.debug('len_fig_y: '+str(len_fig_y))
        fig = sg.SVGFigure(str(len_fig_x),
                           str(len_fig_y))
        background_white = draw.Drawing(len_fig_x, len_fig_y, origin=(0,0))
        '''
        background_white.append(draw.Rectangle(0, 0, len_fig_x, len_fig_y, fill='#FFFFFF'))
        '''
        background_white.append(draw.Rectangle(0, 0, len_fig_x, len_fig_y, fill='#e0d7d7'))
        a = sg.fromstring(background_white.asSvg())
        b_w = a.getroot()
        #b_w.moveto(0, background_len_y+subplot_size[1])#WARNING: not sure why I have to + subpot
        b_w.moveto(0, background_len_y)#WARNING: not sure why I have to + subpot
        fig.append(b_w)
        nodes_attach_locs = {}
        for node_id in mod_pos:
            node = G.node.get(node_id)
            self.logger.debug('\tSpecies: '+str(node_id))
            if node['type']=='species':
                self.logger.debug('\tNode pos: '+str(mod_pos[node_id]))
                self.logger.debug('\tx: '+str(mod_pos[node_id][0]))
                self.logger.debug('\ty: '+str(mod_pos[node_id][1]))
                move_x = mod_pos[node_id][0]
                #move_y = len_fig_y-mod_pos[node_id][1]-subplot_size[1]/2
                # because of the nature of the x y locations, need to reverse them here
                if mod_pos[node_id][1]==0.0:
                    y_coord = mod_pos[node_id][1]+subplot_size[1]/2
                    move_y = len_fig_y-mod_pos[node_id][1]-subplot_size[1]
                    #move_y = len_fig_y-mod_pos[node_id][1]
                elif mod_pos[node_id][1]==len_fig_y:
                    y_coord = mod_pos[node_id][1]-subplot_size[1]/2
                    move_y = len_fig_y-mod_pos[node_id][1]
                    #move_y = len_fig_y-mod_pos[node_id][1]-subplot_size[1]
                else:
                    y_coord = mod_pos[node_id][1]
                    move_y = len_fig_y-mod_pos[node_id][1]-subplot_size[1]/2
                self.logger.debug('\tmove_x: '+str(move_x))
                self.logger.debug('\tmove_y: '+str(move_y))
                f = sg.fromstring(id_svg[node_id])
                p = f.getroot()
                #Rememeber that you are moving the object on the x and y axis while the coordinates are coordinates so its reversed
                '''
                p.moveto(mod_pos[node_id][0],
                         mod_pos[node_id][1])
                '''
                p.moveto(move_x, move_y)
                fig.append(p)
                #if mod_pos[node_id][1]==0:
                '''
                if mod_pos[node_id][1]-subplot_size[1]/2==-100:
                    newy = 100
                else:
                    newy = mod_pos[node_id][1]-subplot_size[1]/2
                nodes_attach_locs[node_id] = {'left': (mod_pos[node_id][0],
                                                       newy),
                                              'right': (mod_pos[node_id][0]+subplot_size[0],
                                                        newy)}
                '''
                '''
                nodes_attach_locs[node_id] = {'left': (mod_pos[node_id][0],
                                                       mod_pos[node_id][1]-subplot_size[1]/2),
                                              'right': (mod_pos[node_id][0]+subplot_size[0],
                                                        mod_pos[node_id][1]-subplot_size[1]/2)}
                '''
                nodes_attach_locs[node_id] = {'left': (mod_pos[node_id][0], y_coord),
                                              'right': (mod_pos[node_id][0]+subplot_size[0], y_coord)}
                '''
                nodes_attach_locs[node_id] = {'left': (move_x,
                                                       move_y+subplot_size[1]-subplot_size[1]/2),
                                              'right': (move_x+subplot_size[0],
                                                        move_y+subplot_size[1]-subplot_size[1]/2)}
                '''
            elif node['type']=='reaction':
                d = draw.Drawing(subplot_size[0], subplot_size[1], origin=(0,0))
                d.append(draw.Rectangle(0, 0, subplot_size[0], subplot_size[1], fill='#FFFFFF'))
                edge_x = subplot_size[0]/2-reac_size[1]/2
                edge_y = subplot_size[1]/2-reac_size[0]/2
                self.logger.debug('\tedge_x: '+str(edge_x))
                self.logger.debug('\tedge_y: '+str(edge_y))
                self.logger.debug('\tx: '+str(mod_pos[node_id][0]))
                self.logger.debug('\ty: '+str(mod_pos[node_id][1]+subplot_size[1]))
                '''
                d.append(draw.Line(subplot_size[0]/2, (subplot_size[1]/2)-reac_size[1]/3,
                                   subplot_size[0]/2, (subplot_size[1]/2)+reac_size[1]/3,
                                   stroke=arrow_stroke_color, stroke_width=arrow_stroke_width, fill='none', marker_end=self.arrowhead))
                '''
                if reac_cofactors_name[node_id]['substrates']:
                    d.append(draw.Line(subplot_size[0]/2, edge_y-reac_size[0],
                                       subplot_size[0]/2, edge_y-self.arrowhead_comp_x,
                                       stroke=arrow_stroke_color, stroke_width=arrow_stroke_width, fill='none', marker_end=self.arrowhead))
                if reac_cofactors_name[node_id]['products']:
                    d.append(draw.Line(subplot_size[0]/2, edge_y+reac_size[0],
                                       subplot_size[0]/2, (subplot_size[1]/2)+reac_size[1]/3,
                                       stroke=arrow_stroke_color, stroke_width=arrow_stroke_width, fill='none', marker_end=self.arrowhead))
                #y_shift = (subplot_size[1]/2)-reac_size[1]
                y_shift = 0.0
                for sub in reac_cofactors_name[node_id]['substrates']:
                    self.logger.debug(sub)
                    d.append(draw.Text(sub, font_size,
                                       subplot_size[0]/2, ((subplot_size[1]/2)-reac_size[1]/3)-y_shift-font_size,
                                       font_family=font_family, center=True, fill=font_color))
                    y_shift += font_size
                #y_shift = (subplot_size[1]/2)+reac_size[1]
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
                #move_y = len_fig_y-mod_pos[node_id][1]-subplot_size[1]/2
                if mod_pos[node_id][1]==0.0:
                    move_y = len_fig_y-mod_pos[node_id][1]
                elif mod_pos[node_id][1]==len_fig_y:
                    move_y = len_fig_y-mod_pos[node_id][1]+subplot_size[1]
                else:
                    move_y = len_fig_y-mod_pos[node_id][1]-subplot_size[1]/2+subplot_size[1]
                self.logger.debug('\tmove_x: '+str(move_x))
                self.logger.debug('\tmove_y: '+str(move_y))
                a_r.moveto(move_x, move_y)
                '''
                a_r.moveto(mod_pos[node_id][0],
                           mod_pos[node_id][1]+subplot_size[1]) #WARNING: not sure why I have to + subpot
                '''
                '''
                a_r.moveto(mod_pos[node_id][0],
                           mod_pos[node_id][1]) #WARNING: not sure why I have to + subpot
                '''
                fig.append(a_r)
                #if mod_pos[node_id][1]==0:
                '''
                if mod_pos[node_id][1]-subplot_size[1]/2==-100:
                    newy = 100
                else:
                    newy = mod_pos[node_id][1]-subplot_size[1]/2
                nodes_attach_locs[node_id] = {'left': (mod_pos[node_id][0]+edge_x,
                                                       newy),
                                              'right': (mod_pos[node_id][0]+subplot_size[0]-edge_x+reac_stroke_width/2,
                                                        newy)}
                '''
                '''
                nodes_attach_locs[node_id] = {'left': (mod_pos[node_id][0]+edge_x,
                                                       mod_pos[node_id][1]-subplot_size[1]/2),
                                              'right': (mod_pos[node_id][0]+subplot_size[0]-edge_x+reac_stroke_width/2,
                                                        mod_pos[node_id][1]-subplot_size[1]/2)}
                '''
                '''
                nodes_attach_locs[node_id] = {'left': (move_x+edge_x,
                                                       move_y-subplot_size[1]/2),
                                              'right': (move_x+subplot_size[0]-edge_x+reac_stroke_width/2,
                                                        move_y-subplot_size[1]/2)}
                '''
                nodes_attach_locs[node_id] = {'left': (move_x+edge_x,
                                                       move_y-subplot_size[1]/2),
                                              'right': (move_x+subplot_size[0]-edge_x+reac_stroke_width/2,
                                                        move_y-subplot_size[1]/2)}
            self.logger.debug('\t-------------------------------')
        self.logger.debug('nodes_attach_locs: '+str(nodes_attach_locs))
        self.logger.debug(list(resG.edges))
        arrow_box = draw.Drawing(len_fig_x, len_fig_y, origin=(0,0))
        for edge in list(resG.edges):
            self.logger.debug('\t---------- edge: '+str(edge)+' -----------')
            node_obj_source = G.node.get(edge[0])
            node_obj_target = G.node.get(edge[1])
            #if node_obj_source['type']=='reaction' and node_obj_target['type']=='species':
            self.logger.debug('source_x: '+str(pos[edge[0]][0]))
            self.logger.debug('target_x: '+str(pos[edge[1]][0]))
            if pos[edge[0]][0]>pos[edge[1]][0]:
                source_x = nodes_attach_locs[edge[0]]['left'][0]
                source_y = nodes_attach_locs[edge[0]]['left'][1]
                #source_y -= len_fig_y-subplot_size[1]/2
                target_x = nodes_attach_locs[edge[1]]['right'][0]
                target_y = nodes_attach_locs[edge[1]]['right'][1]
                #target_y -= len_fig_y-subplot_size[1]/2
                self.logger.debug('\tsource_x: '+str(source_x))
                self.logger.debug('\tsource_y: '+str(source_y))
                self.logger.debug('\ttarget_x: '+str(target_x))
                self.logger.debug('\ttarget_y: '+str(target_y))
                '''
                len_x = source_x-target_x
                len_y = source_y-target_y
                self.logger.debug('\tlen_x: '+str(len_x))
                self.logger.debug('\tlen_y: '+str(len_y))
                d = draw.Drawing(len_x+10, len_y+10, origin=(0,0))
                '''
                '''
                if source_y==-100.0:
                    source_y = 100.0
                if target_y==-100.0:
                    target_y = 100.0
                if source_y==400.0:
                    source_y = 100.0
                if target_y==400.0:
                    target_y = 100.0
                '''
                self.logger.debug('\tsource_x: '+str(source_x))
                self.logger.debug('\tsource_y: '+str(source_y))
                self.logger.debug('\ttarget_x: '+str(target_x))
                self.logger.debug('\ttarget_y: '+str(target_y))
                p = draw.Path(stroke=arrow_stroke_color,
                              stroke_width=arrow_stroke_width,
                              fill='none',
                              marker_end=self.rev_arrowhead_flat)
                '''
                p.M(source_x, source_y).C(target_x+self.arrowhead_comp_x-1.5,
                                          source_y,
                                          source_x,
                                          target_y,
                                          target_x+self.arrowhead_comp_x-1.5,
                                          target_y)
                '''
                '''
                p.M(0.0, 0.0).L((target_x-source_x)/2, 0.0).L((target_x-source_x)/2, len_y).L(len_x, len_y)
                d.append(p)
                a = sg.fromstring(d.asSvg())
                a_r = a.getroot()
                a_r.moveto(source_x, source_y)
                fig.append(a_r)
                '''
                self.logger.debug((source_x, source_y))
                self.logger.debug((source_x+(target_x-source_x)/2, source_y))
                self.logger.debug((source_x+(target_x-source_x)/2, target_y))
                self.logger.debug((target_x-self.arrowhead_comp_x, target_y))
                p.M(source_x, source_y).L(source_x+(target_x-source_x)/2, source_y).L(source_x+(target_x-source_x)/2, target_y).L(target_x+self.arrowhead_comp_x, target_y)
                #background.append(p)
                arrow_box.append(p)
            #elif node_obj_source['type']=='species' and node_obj_target['type']=='reaction':
            elif pos[edge[0]][0]<pos[edge[1]][0]:
                source_x = nodes_attach_locs[edge[0]]['right'][0]
                source_y = nodes_attach_locs[edge[0]]['right'][1]
                target_x = nodes_attach_locs[edge[1]]['left'][0]
                target_y = nodes_attach_locs[edge[1]]['left'][1]
                self.logger.debug('\tsource_x: '+str(source_x))
                self.logger.debug('\tsource_y: '+str(source_y))
                self.logger.debug('\ttarget_x: '+str(target_x))
                self.logger.debug('\ttarget_y: '+str(target_y))
                '''
                len_x = source_x-target_x
                len_y = source_y-target_y
                self.logger.debug('\tlen_x: '+str(len_x))
                self.logger.debug('\tlen_y: '+str(len_y))
                d = draw.Drawing(len_x+10, len_y+10, origin=(0,0))
                '''
                '''
                if source_y==-100.0:
                    source_y = 100.0
                if target_y==-100.0:
                    target_y = 100.0
                if source_y==300.0:
                    source_y = 100.0
                if target_y==300.0:
                    target_y = 100.0
                self.logger.debug('\tsource_x: '+str(source_x))
                self.logger.debug('\tsource_y: '+str(source_y))
                self.logger.debug('\ttarget_x: '+str(target_x))
                self.logger.debug('\ttarget_y: '+str(target_y))
                '''
                p = draw.Path(stroke=arrow_stroke_color,
                              stroke_width=arrow_stroke_width,
                              fill='none',
                              marker_end=self.arrowhead_flat)
                '''
                p.M(source_x, source_y).C(subplot_size[0]/2-self.arrowhead_comp_x+source_x,
                                          source_y,
                                          source_x,
                                          target_y,
                                          target_x-self.arrowhead_comp_x,
                                          target_y)
                '''
                '''
                p.M(0.0, 0.0).L((target_x-source_x)/2, 0.0).L((target_x-source_x)/2, len_y).L(len_x, len_y)
                d.append(p)
                a = sg.fromstring(d.asSvg())
                a_r = a.getroot()
                a_r.moveto(source_x, source_y)
                fig.append(a_r)
                '''
                self.logger.debug((source_x, source_y))
                self.logger.debug((source_x+(target_x-source_x)/2, source_y))
                self.logger.debug((source_x+(target_x-source_x)/2, target_y))
                self.logger.debug((target_x-self.arrowhead_comp_x, target_y))
                p.M(source_x, source_y).L(source_x+(target_x-source_x)/2, source_y).L(source_x+(target_x-source_x)/2, target_y).L(target_x-self.arrowhead_comp_x, target_y)
                #background.append(p)
                arrow_box.append(p)
            else:
                self.logger.error('Cannot connect same y-axi')
        a = sg.fromstring(arrow_box.asSvg())
        a_b = a.getroot()
        a_b.moveto(0, background_len_y)#WARNING: not sure why I have to + subpot
        fig.append(a_b)
        svg = fig.to_str().decode("utf-8")
        return svg, resG, pos, reac_cofactors_name






    ######################################################################################################
    ########################################## Public Function ###########################################
    ######################################################################################################

    def drawsvg1(self,
                G,
                path=None,
                plot_only_central=True,
                filter_cofactors=True,
                react_arrow_size=100,
                arrow_gap_size=100,
                subplot_size=[200,200],
                stroke_color='black',
                stroke_width=2):
        #TODO: need to find a better way to find the root node without name
        root_node = [i for i in list(rpgraph.G.nodes) if 'TARGET' in i]
        if not len(root_node)==1:
            self.logger.debug('There are multiple nodes with TARGET')
            return False
        pathway_list = []
        drawG, pos = self._hierarchy_pos(G, root_node[0], plot_only_central=plot_only_central, filter_cofactors=filter_cofactors)
        #order the layers on the y-axis and then order the layers on the x-axis

        svg, len_x, len_y = self._drawPathway(pathway_list,
                                              react_arrow_size=react_arrow_size,
                                              arrow_gap_size=arrow_gap_size,
                                              subplot_size=subplot_size,
                                              stroke_color=stroke_color,
                                              stroke_width=stroke_width)
        if path:
            open(path, 'w').write(svg)
        return svg

