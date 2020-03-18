from sequal.base_block import BaseBlock
from sequal.modification import Modification
from sequal.resources import AA_mass


class AminoAcid(BaseBlock):
    def __init__(self, value, position=None, mass=None):
        super().__init__(value, position, branch=False, mass=mass)
        self.mods = []
        if not self.mass:
            if value in AA_mass:
                self.mass = AA_mass[value]

    def set_modification(self, i: Modification):
        self.mods.append(i)

    def __repr__(self):
        s = self.value
        for i in self.mods:
            s += "[{}]".format(i.value)
        return s

    def __str__(self):
        s = self.value
        for i in self.mods:
            s += "[{}]".format(i.value)
        return s


