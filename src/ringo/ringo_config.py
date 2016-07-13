import ConfigParser
import os


class RingoConfig(object):
    # Singleton pattern, to avoid reading the config file multiple times.
    __instance = None
    def __new__(cls):
        if RingoConfig.__instance is None:
            instance = object.__new__(cls)
            instance.config = ConfigParser.ConfigParser()
            instance.config.readfp(open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'ringo.cfg')))
            RingoConfig.__instance = instance
        return RingoConfig.__instance

    # RINGO
    def ringo_output_tree(self):
        return self.config.get("Filenames", "ringo_output_tree")

    def ringo_output_genomes(self):
        return self.config.get("Filenames", "ringo_output_genomes")

    def ringo_output_parameters(self):
        return self.config.get("Filenames", "ringo_output_parameters")

    # DeClone
    def declone_path(self):
        return self.config.get("Paths", "declone_path")

    def declone_output_extant_weight(self, kT):
        return self.config.get("Filenames", "declone_extant_weight") + '{0:.6g}'.format(kT)

    def declone_output_single_leaf(self, kT):
        return self.config.get("Filenames", "declone_single_leaf") + '{0:.6g}'.format(kT)

    def declone_output_internal_weight(self, kT):
        return self.config.get("Filenames", "declone_internal_weight") + '{0:.6g}'.format(kT)

    def declone_nhx_tree(self):
        return self.config.get("Filenames", "declone_nhx_tree")

    # MGRA
    def mgra_path(self):
        return self.config.get("Paths", "mgra_path")

    # PhySca
    def physca_path(self):
        return self.config.get("Paths", "physca_path")

    def physca_reconstructed_adjacencies(self):
        return self.config.get("Filenames", "physca_reconstructed_adjacencies")

    # Blossom5
    def blossom5_path(self):
        return self.config.get("Paths", "blossom5_path")

    # Simulations
    def sim_leaf_genomes(self):
        return self.config.get("Filenames", "sim_leaf_genomes")

    def sim_ancestral_genomes(self):
        return self.config.get("Filenames", "sim_ancestral_genomes")

    def sim_tree(self):
        return self.config.get("Filenames", "sim_tree")

    def sim_tree_no_lengths(self):
        return self.config.get("Filenames", "sim_tree_no_lengths")

    def sim_logfile(self):
        return self.config.get("Filenames", "sim_logfile")

    def sim_mgra_config(self):
        return self.config.get("Filenames", "sim_mgra_config")

    def sim_paramfile(self):
        return self.config.get("Filenames", "sim_paramfile")

    def scj_genomes(self):
        return self.config.get("Filenames", "scj_genomes")

    def mgra_output_folder(self):
        return self.config.get("Filenames", "mgra_output_folder")
