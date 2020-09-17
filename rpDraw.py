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
        self.rev_arrowhead = draw.Marker(-0.1, -0.5, 0.9, 0.5, scale=4, orient=0, id='rev_flat_arrow')
        self.rev_arrowhead.append(draw.Lines(-0.1, -0.5, -0.1, 0.5, 0.9, 0, fill='black', close=True))
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


    ######################################################################################################
    ########################################## Public Function ###########################################
    ######################################################################################################

    
    ##
    #
    # TODO: add the conparison by inchikey to determine the cofactors
    def drawsvg(self,
                G,
                path=None,
                react_arrow_size=100,
                arrow_gap_size=100,
                subplot_size=[200,200],
                stroke_color='black',
                stroke_width=2):
        ordered_reactions = self.orderedRetroReactions()
        if not ordered_reactions:
            self.logger.error('Ordered reaction returned emtpy results')
            return ''
        self.logger.debug('ordered_reactions: '+str(ordered_reactions))
        flat_ordered_reactions = []
        for i in ordered_reactions:
            if len(i)==1:
                flat_ordered_reactions.append(i[0])
            else:
                self.logger.error('This pathway contains steps with multiple reactions')
                return ''
        self.logger.debug('flat_ordered_reactions: '+str(flat_ordered_reactions))
        pathway_list = []
        for reaction_id in flat_ordered_reactions:
            #edges
            to_add = {'reactants_inchi': [],
                      'products_inchi': [],
                      'cofactor_reactants': [],
                      'cofactor_products': []}
            for pred_id in G.predecessors(reaction_id):
                pred = G.node.get(pred_id)
                if pred['type']=='species':
                    if pred['central_species']:
                        if 'metanetx' in pred['miriam']:
                            if any([i in list(self.mnx_cofactors.keys()) for i in pred['miriam']['metanetx']]):
                                #double check but lazy
                                for mnx in pred['miriam']['metanetx']:
                                    if mnx in self.mnx_cofactors.keys():
                                        to_add['cofactor_reactants'].append(self.mnx_cofactors[mnx]['name'])
                                        break
                                #to_add['cofactor_reactants'].append(pred['name'])
                            else:
                                to_add['reactants_inchi'].append(pred['brsynth']['inchi'])
                        else:
                            to_add['reactants_inchi'].append(pred['brsynth']['inchi'])
                    else:
                        is_found = False
                        if 'metanetx' in pred['miriam']:
                            for mnx in pred['miriam']['metanetx']:
                                if mnx in self.mnx_cofactors.keys():
                                    to_add['cofactor_reactants'].append(self.mnx_cofactors[mnx]['name'])
                                    is_found = True
                                    break
                        if not is_found:
                            to_add['cofactor_reactants'].append(pred['name'])
                else:
                    self.logger.warning('The current node is not a species: '+str(pred))
            prev_products = [] #reset the list
            for succ_id in G.successors(reaction_id):
                succ = G.node.get(succ_id)
                if succ['type']=='species':
                    if succ['central_species']:
                        if 'metanetx' in succ['miriam']:
                            if any([i in list(self.mnx_cofactors.keys()) for i in succ['miriam']['metanetx']]):
                                #double check i know but lazy
                                for mnx in succ['miriam']['metanetx']:
                                    if mnx in self.mnx_cofactors.keys():
                                        to_add['cofactor_products'].append(self.mnx_cofactors[mnx]['name'])
                                        break
                                #to_add['cofactor_products'].append(succ['name'])
                            else:
                                to_add['products_inchi'].append(succ['brsynth']['inchi'])
                        else:
                            to_add['products_inchi'].append(succ['brsynth']['inchi'])
                    else:
                        is_found = False
                        if 'metanetx' in succ['miriam']:
                            for mnx in succ['miriam']['metanetx']:
                                if mnx in self.mnx_cofactors.keys():
                                    to_add['cofactor_products'].append(self.mnx_cofactors[mnx]['name'])
                                    is_found = True
                                    break
                        if not is_found:
                            to_add['cofactor_products'].append(succ['name'])
                else:
                    self.logger.warning('The current node is not a species: '+str(succ))
            pathway_list.append(to_add)
        self.logger.debug('------------------------------------------')
        self.logger.debug(pathway_list)
        self.logger.debug('------------------------------------------')
        svg, len_x, len_y = self._drawPathway(pathway_list,
                                              react_arrow_size=react_arrow_size,
                                              arrow_gap_size=arrow_gap_size,
                                              subplot_size=subplot_size,
                                              stroke_color=stroke_color,
                                              stroke_width=stroke_width)
        if path:
            open(path, 'w').write(svg)
        return svg
