import chemlib
from rdkit import Chem
import pubchempy as pcp


def get_formula_and_smiles() -> list:
    """ Return formula and SMILES representation of compound entered by user"""
    compound = input("Enter compound name: ")
    result = pcp.get_compounds(compound.strip(), 'name')  # Get compound info from PubChem database

    while len(result) == 0:
        print("PubChem did not recognize compound. Try again")
        compound = input("Enter compound: ")
        result = pcp.get_compounds(compound.strip(), 'name')  # Can do name and formula

    compound = result[0]

    formula = compound.to_dict(properties=['molecular_formula'])['molecular_formula']
    smiles = compound.isomeric_smiles

    return [formula, smiles]


def _get_bond_numbers(smiles: str) -> list:
    """ Return list of bond numbers based on compound SMILES representation"""
    moles = Chem.MolFromSmiles(smiles)

    return [bond.GetBondTypeAsDouble() for bond in moles.GetBonds()]


def _create_vsepr_dict(keys: list, values: list) -> dict:
    """ Return dictionary containing vsepr info based on keys and values of compound """
    vsepr_dict = dict(map(lambda i, j: (i, j), keys, values))
    return vsepr_dict


def return_vsepr_info() -> dict:
    """ Return VSEPR info"""
    formula, smiles = get_formula_and_smiles()

    bond_list = _get_bond_numbers(smiles)

    compound = chemlib.Compound(formula)

    element_dict = compound.occurences

    electronegativity_dict = dict()

    for key in element_dict.keys():
        element = chemlib.Element(key)
        electronegativity_dict[key] = element.properties['Electronegativity']

    elect_dict_keys = list(electronegativity_dict.keys())
    elect_dict_values = list(electronegativity_dict.values())

    position = elect_dict_values.index(min(elect_dict_values))

    central_element = chemlib.Element(elect_dict_keys[position])

    valence_elec_central = int(central_element.properties['Valence'])

    valence_electrons_needed = 0

    if len(bond_list) > 0:
        valence_electrons_needed = int(sum(bond_list))
    else:
        for key in elect_dict_keys:
            if key != central_element.properties['Symbol']:
                element = chemlib.Element(key)
                electrons_needed = (8 - int(element.properties['Valence'])) * element_dict[key]
                valence_electrons_needed += electrons_needed
            else:
                continue

    num_lone_pairs = (valence_elec_central - valence_electrons_needed) // 2

    if num_lone_pairs < 0:
        num_lone_pairs = 0

    num_bonding_pairs = sum(element_dict[key] for key in element_dict if key !=
                            central_element.properties['Symbol'])

    steric_num = num_lone_pairs + num_bonding_pairs

    key_list = ['Shape', 'Molecular Geometry', 'Hybridization', 'Ideal Bond Angle', 'Bonded Pairs (X)',
                'Lone Pairs (E)']

    match steric_num:
        case 2:
            value_list = ['Linear', 'Linear', 'sp', '180', num_bonding_pairs, num_lone_pairs]
            return _create_vsepr_dict(key_list, value_list)
        case 3:
            match num_lone_pairs:
                case 0:
                    value_list = ['Trigonal Planar', 'Trigonal Planar', 'sp2', '180', num_bonding_pairs,
                                  num_lone_pairs]
                    return _create_vsepr_dict(key_list, value_list)
                case 1:
                    value_list = ['Trigonal Planar', 'Bent', 'sp2', '< 120', num_bonding_pairs, num_lone_pairs]
                    return _create_vsepr_dict(key_list, value_list)

                case _:
                    value_list = ['Trigonal Planar', 'Trigonal Planar', 'sp2', '180', num_bonding_pairs,
                                  num_lone_pairs]
                    return _create_vsepr_dict(key_list, value_list)

        case 4:
            match num_lone_pairs:
                case 0:
                    value_list = ['Tetrahedral', 'Tetrahedral', 'sp3', '109.5', num_bonding_pairs,
                                  num_lone_pairs]
                    return _create_vsepr_dict(key_list, value_list)
                case 1:
                    value_list = ['Tetrahedral', 'Trigonal Pyramidal', 'sp3', '< 109.5', num_bonding_pairs,
                                  num_lone_pairs]
                    return _create_vsepr_dict(key_list, value_list)
                case 2:
                    value_list = ['Tetrahedral', 'Bent', 'sp3', '<< 109.5', num_bonding_pairs, num_lone_pairs]
                    return _create_vsepr_dict(key_list, value_list)

                case _:
                    value_list = ['Tetrahedral', 'Tetrahedral', 'sp3', '109.5', num_bonding_pairs,
                                  num_lone_pairs]
                    return _create_vsepr_dict(key_list, value_list)

        case 5:
            match num_lone_pairs:
                case 0:
                    value_list = ['Trigonal Bipyramidal', 'Trigonal Bipyramidal', 'sp3',
                                  '90 (Axial) 120 (Equatorial)', num_bonding_pairs, num_lone_pairs]
                    return _create_vsepr_dict(key_list, value_list)
                case 1:
                    value_list = ['Trigonal Bipyramidal', 'Seesaw', 'sp3', '90 120 180', num_bonding_pairs,
                                  num_lone_pairs]
                    return _create_vsepr_dict(key_list, value_list)
                case 2:
                    value_list = ['Trigonal Bipyramidal', 'T-shaped', 'sp3', '90 180', num_bonding_pairs,
                                  num_lone_pairs]
                    return _create_vsepr_dict(key_list, value_list)
                case 3:
                    value_list = ['Trigonal Bipyramidal', 'Linear', 'sp3', '180', num_bonding_pairs,
                                  num_lone_pairs]
                    return _create_vsepr_dict(key_list, value_list)
                case _:
                    value_list = ['Trigonal Bipyramidal', 'Trigonal Bipyramidal', 'sp3',
                                  '90 (Axial) 120 (Equatorial)', num_bonding_pairs, num_lone_pairs]
                    return _create_vsepr_dict(key_list, value_list)

        case 6:
            match num_lone_pairs:
                case 0:
                    value_list = ['Octahedral', 'Octahedral', 'sp3', '90', num_bonding_pairs, num_lone_pairs]
                    return _create_vsepr_dict(key_list, value_list)
                case 1:
                    value_list = ['Octahedral', 'Square Pyramidal', 'sp3', '90 180', num_bonding_pairs,
                                  num_lone_pairs]
                    return _create_vsepr_dict(key_list, value_list)
                case 2:
                    value_list = ['Octahedral', 'Square Planar', 'sp3', '90 180', num_bonding_pairs,
                                  num_lone_pairs]
                    return _create_vsepr_dict(key_list, value_list)
                case 3:
                    value_list = ['Octahedral', 'T-shaped', 'sp3', '90 180', num_bonding_pairs,
                                  num_lone_pairs]
                    return _create_vsepr_dict(key_list, value_list)
                case 4:
                    value_list = ['Octahedral', 'Linear', 'sp3', '180', num_bonding_pairs, num_lone_pairs]
                    return _create_vsepr_dict(key_list, value_list)
                case _:
                    value_list = ['Octahedral', 'Octahedral', 'sp3', '90', num_bonding_pairs, num_lone_pairs]
                    return _create_vsepr_dict(key_list, value_list)
