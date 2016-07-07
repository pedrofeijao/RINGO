import ConfigParser
import os


class RingoConfig():

    def __init__(self):
        self.config = ConfigParser.ConfigParser()
        self.config.readfp(
            open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'ringo.cfg')))

    def declone_path(self):
        return self.config.get("Paths", "declone_path")

    def declone_output_extant_weight(self, kT):
        return self.config.get("Filenames", "declone_extant_weight") + '{0:.3g}'.format(kT)

    def declone_output_single_leaf(self, kT):
        return self.config.get("Filenames", "declone_single_leaf") + '{0:.3g}'.format(kT)

    def declone_output_internal_weight(self, kT):
        return self.config.get("Filenames", "declone_internal_weight") + '{0:.3g}'.format(kT)

    def declone_nhx_tree(self):
        return self.config.get("Filenames", "declone_nhx_tree")

    def mgra_path(self):
        return self.config.get("Paths", "mgra_path")

    def sim_leaf_genomes(self):
        return self.config.get("Filenames", "sim_leaf_genomes")

    def sim_ancestral_genomes(self):
        return self.config.get("Filenames", "sim_ancestral_genomes")

    def sim_tree(self):
        return self.config.get("Filenames", "sim_tree")

    def sim_basic_tree(self):
        return self.config.get("Filenames", "sim_basic_tree")

    def sim_logfile(self):
        return self.config.get("Filenames", "sim_logfile")

    def sim_mgra_config(self):
        return self.config.get("Filenames", "sim_mgra_config")

    def sim_paramfile(self):
        return self.config.get("Filenames", "sim_paramfile")
