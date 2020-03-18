from sequal.amino_acid import AminoAcid
from sequal.sequence import Sequence
from sequal import resources


def calculate_mass(seq, mass_dict=None, N_terminus=0, O_terminus=0, with_water=True):
    mass = 0
    if with_water:
        mass += resources.H*2 + resources.O
    for i in seq:
        if not i.mass:
            if mass_dict:
                if i.value in mass_dict:
                    mass += mass_dict[i.value]
                else:
                    raise ValueError('Block {} not found in mass_dict'.format(i.value))
            else:
                raise ValueError('Block {} mass is not available in mass attribute and no additional mass_dict was supplied'.format(i.value))
        else:
            mass += i.mass
        if i.mods:
            for m in i.mods:
                if m.mass != 0 and not m.mass:
                    if mass_dict:
                        if m.value in mass_dict:
                            mass += mass_dict[m.value]
                        else:
                            raise ValueError('Block {} not found in mass_dict'.format(m.value))
                    else:
                        raise ValueError(
                            'Block {} mass is not available in mass attribute and no additional mass_dict was supplied'.format(
                                m.value))
                else:
                    mass += m.mass
    return mass + N_terminus + O_terminus






