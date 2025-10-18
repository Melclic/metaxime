from biopathopt import Data
import compress_json
import logging
import requests

from .utils import read_compressed_tsv

class RR_Data(Data):
    def __init__(
            self, 
            use_progressbar=False, 
            low_memory_mode=False,
        ):
        """Class that inherits Data used to build a cobra model
        """
        super().__init__(low_memory_mode=low_memory_mode, use_progressbar=use_progressbar)

    def rr_recipes(self):
        """Return the chemical properties of molecules

        This function will download the following file https://www.metanetx.org/cgi-bin/mnxget/mnxref/mnxr_prop.tsv that
        describes the chemical structure, etc...

        Args:
        Returns:
            dict: Reaction properties from MetaNetX
        """
        if not self._rr_recipes:
            logging.debug("------ rr_recipes -----")
            logging.debug("\t-> Populating...")
            rr_recipes_path = os.path.join(self.base_dir, "flatfiles/rr_recipes.json.gz")
            if not os.path.exists(rr_recipes_path):
                rr_recipes = read_compressed_tsv('data/rxn_recipes.tsv.tar.gz')
                rr_recipes.columns = [
                    "Reaction_ID", 
                    "Equation", 
                    "Description", 
                    "Direction", 
                    "EC_number", 
                    "Name", 
                    "Type", 
                    "UniProt_IDs",
                    "Additional_RIDs", 
                    "Main_left", 
                    "Main_right", 
                    "Secondary_left", 
                    "Secondary_right",
                ]
                rr_recipes = rr_recipes.set_index("Reaction_ID")
                rr_recipes = rr_recipes.transpose().to_dict()
                for i in rr_recipes:
                    try:
                        rr_recipes[i]["EC_number"] = rr_recipes[i]["EC_number"].split(
                            ","
                        )
                    except AttributeError:
                        rr_recipes[i]["EC_number"] = []
                ### generate the reaction based in inchikeys
                for i in rr_recipes:
                    try:
                        rr_recipes[i]['inchikey2_equation'] = \
                            self.convert_mnxr_equation(rr_recipes[i]['Equation'], inchikey_levels=2)
                    except (ValueError, KeyError) as e:
                        rr_recipes[i]['inchikey2_equation'] = ''
                    try:
                        rr_recipes[i]['inchikey_equation'] = \
                            self.convert_mnxr_equation(rr_recipes[i]['Equation'], inchikey_levels=3)
                    except (ValueError, KeyError) as e:
                        rr_recipes[i]['inchikey_equation'] = ''
                ### generate main left and main right
                for i in rr_recipes:
                    reactants, products = self.parse_mnxr_equation(rr_recipes[i]['Equation'])
                    rr_recipes[i]['main_reactants'] = {y[1]: y[0] for y in reactants if y[1] not in self.mnxm_cofactors}
                    rr_recipes[i]['secondary_reactants'] = {y[1]: y[0] for y in reactants if y[1] in self.mnxm_cofactors}
                    rr_recipes[i]['main_products'] = {y[1]: y[0] for y in products if y[1] not in self.mnxm_cofactors}
                    rr_recipes[i]['secondary_products'] = {y[1]: y[0] for y in products if y[1]  in self.mnxm_cofactors}
                #save it
                compress_json.dump(rr_recipes, rr_recipes_path)
            self._rr_recipes = compress_json.load(rr_recipes_path)
        return self._rr_recipes
