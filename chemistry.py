# Libraries
import chemlib
from pprint import pprint
import pandas as pd

# Files
import vsepr
import draw_molecule
import titration
import element_properties


def _balance_equation(reactants_list: list, products_list: list) -> chemlib.Reaction:
    """ Balance and return reaction based on list of reactants and products"""
    reaction = chemlib.Reaction(reactants=reactants_list, products=products_list)
    reaction.balance()

    return reaction


def _find_limiting_reagent(reaction: chemlib.Reaction, amounts_list: list, unit: str) -> chemlib.Compound:
    """ Calculate and return limiting reagent based on reaction
        Return"""
    match len(amounts_list):
        case 1:
            return reaction.limiting_reagent(float(amounts_list[0]), mode=unit)
        case 2:
            return reaction.limiting_reagent(float(amounts_list[0]), float(amounts_list[1]), mode=unit)
        case 3:
            return reaction.limiting_reagent(float(amounts_list[0]), float(amounts_list[1]),
                                             float(amounts_list[-1]), mode=unit)
        case 4:
            return reaction.limiting_reagent(float(amounts_list[0]), float(amounts_list[1]),
                                             float(amounts_list[2]), float(amounts_list[-1]), mode=unit)
        case 5:
            return reaction.limiting_reagent(float(amounts_list[0]), float(amounts_list[1]),
                                             float(amounts_list[2]), float(amounts_list[3]),
                                             float(amounts_list[-1]), mode=unit)


def _convert_dict_to_dataframe(info_dict: dict) -> pd.DataFrame:
    """ Convert dictionary to pandas dataframe """
    return pd.DataFrame([info_dict])


def main() -> None:
    """ Run Chemistry Calculator"""
    while True:
        command = input("Enter command (h for help): ")

        match command.strip().lower():
            case 'b' | 'balance':
                reactants_input = input("Enter reactants (separate by commas and spaces): ")
                products_input = input("Enter products (separate by commas and spaces): ")

                reactants_list = reactants_input.split(', ')  # List of reactants represented by str
                reactants_list = [chemlib.Compound(i) for i in reactants_list]

                products_list = products_input.split(', ')  # List of reactants represented by str
                products_list = [chemlib.Compound(i) for i in products_list]

                balanced_reaction = _balance_equation(reactants_list, products_list)

                print(f"Balanced Reaction: {balanced_reaction.formula}")  # Balanced Reaction

            case 'l' | 'limiting reagent' | 'limiting':
                reactants_input = input("Enter reactants (separate by commas and spaces): ")
                products_input = input("Enter products (separate by commas and spaces): ")

                reactants_list = reactants_input.split(', ')  # List of reactants represented by str
                reactants_list = [chemlib.Compound(i) for i in reactants_list]

                products_list = products_input.split(', ')  # List of reactants represented by str
                products_list = [chemlib.Compound(i) for i in products_list]

                balanced_reaction = _balance_equation(reactants_list, products_list)

                amounts = input("Enter amounts of reactants (separated by space): ")
                amounts_list = amounts.split(' ')

                unit = input("Enter unit (grams, moles, molecules): ")

                match unit.strip().lower():
                    case 'g' | 'gram' | 'grams':
                        limiting_reagent = _find_limiting_reagent(balanced_reaction, amounts_list, 'grams')
                        print(limiting_reagent.formula)
                    case 'moles':
                        limiting_reagent = _find_limiting_reagent(balanced_reaction, amounts_list, 'moles')
                        print(limiting_reagent.formula)
                    case 'molecules':
                        limiting_reagent = _find_limiting_reagent(balanced_reaction,
                                                                  amounts_list, 'molecules')
                        print(limiting_reagent.formula)

            case 'c' | 'combustion':
                carbon_dioxide_amount = float(input("Enter grams of carbon dioxide produced: "))
                water_amount = float(input("Enter grams of water produced: "))

                print(chemlib.thermochemistry.combustion_analysis(carbon_dioxide_amount, water_amount))

            case 'a' | 'acidity' | 'acid':  # Find pH, pOH, [H+], [OH-], acidity status
                molarity_type = input("Enter type (m for molarity and p for pH): ")

                match molarity_type.strip().lower():
                    case 'm':
                        molarity = float(input("Enter molarity: "))
                        h_status = input('Enter h for H and oh for OH: ')

                        match h_status.strip().lower():
                            case 'h':
                                info_dict = chemlib.pH(H=molarity)

                                print(_convert_dict_to_dataframe(info_dict))

                            case 'oh':
                                info_dict = chemlib.pH(OH=molarity)

                                print(_convert_dict_to_dataframe(info_dict))

                    case 'p':
                        acidity_num = float(input("Enter pH: "))

                        info_dict = chemlib.pH(pH=acidity_num)

                        print(_convert_dict_to_dataframe(info_dict))

            case 'w' | 'wave' | 'waves':  # Find wavelength, frequency, energy
                value_type = input("Enter type of value (wavelength, frequency, or energy): ")
                value = float(input("Enter value (use e- instead of 10^-): "))

                match value_type.strip().lower():
                    case 'w' | 'wavelength':
                        wave = chemlib.Wave(wavelength=value)
                        print(pd.DataFrame([wave.properties]))
                    case 'f' | 'frequency':
                        wave = chemlib.Wave(frequency=value)
                        print(pd.DataFrame([wave.properties]))
                    case 'e' | 'energy':
                        wave = chemlib.Wave(energy=value)
                        print(pd.DataFrame([wave.properties]))

            case 'v' | 'vespr' | 'vsepr':
                print(vsepr.return_vsepr_info())

            case 'd' | 'draw':
                draw_molecule.draw_2d_molecule()

            case 'e' | 'electrolysis' | 'electrochem':  # Find cell, cathode, anode, cell potential
                electrodes = input("Enter electrodes (separate by commas and spaces): ")

                electrode_list = electrodes.split(', ')

                electrode_1 = electrode_list[0]
                electrode_2 = electrode_list[-1]

                galvanic_cell = chemlib.Galvanic_Cell(electrode_1, electrode_2)

                print(_convert_dict_to_dataframe(galvanic_cell.properties).set_index('Cell'))

            case 't' | 'titration':
                starting = input("Acid or base: ")
                acid_volume = float(input("Enter volume (mL): "))
                acid_conc = float(input("Enter first concentration (M): "))

                init_k = float(input("Enter initial k value (use e- instead of 10^-): "))
                other_conc = float(input("Enter other concentration (M): "))

                t = titration.Titration(acid_volume, acid_conc, init_k, other_conc, starting.strip().lower())

                pprint(t.titration_info())

                diagram_input = input("Diagram (Y or N): ")

                if diagram_input.strip().lower() == 'y':
                    t.plot_titration_curve()

            case 'p' | 'property' | 'properties':
                elements_input = str(input('Enter elements separated by commas and spaces: '))

                elements_list = elements_input.split(', ')

                prop_df = element_properties.return_prop_df(elements_list)

                print(prop_df.set_index('Symbol'))

                diagram_input = input("Diagram (Y or N): ")

                if diagram_input.strip().lower() == 'y':
                    element_properties.draw_comparison(prop_df)

            case 's' | 'stoich' | 'stoichiometry':
                while True:
                    try:
                        compound_input = str(input('Enter compound: ')).strip()
                        compound = chemlib.Compound(compound_input)
                    except IndexError:
                        print("Invalid compound. Try again.")
                        continue
                    else:
                        break

                amount = float(input("Enter amount of compound: "))

                unit = input("Enter unit (grams, moles, molecules): ")

                def _add_compound_name_to_output(stoich_dict: dict) -> dict:
                    stoich_dict['compound'] = compound_input
                    stoich_dict['molar mass'] = f'{compound.molar_mass()} g/mol'
                    return stoich_dict

                match unit.strip().lower():
                    case 'g' | 'gram' | 'grams':
                        stoich_dict = _add_compound_name_to_output(compound.get_amounts(grams=amount))
                        print(_convert_dict_to_dataframe(stoich_dict).set_index('compound'))
                    case 'moles':
                        stoich_dict = _add_compound_name_to_output(compound.get_amounts(moles=amount))
                        print(_convert_dict_to_dataframe(stoich_dict).set_index('compound'))
                    case 'molecules':
                        stoich_dict = _add_compound_name_to_output(compound.get_amounts(molecules=amount))
                        print(_convert_dict_to_dataframe(stoich_dict).set_index('compound'))

            case 'h' | 'help':
                with open('README.md') as f:
                    lines: list = f.readlines()
                    command = lines.index('Commands:\n')
                    command_lines = lines[command+2:-1]  # Excludes help command from being printed
                    removed_space_lines = [line[:-1] for line in command_lines]

                    for line in removed_space_lines:
                        print(f'\t\t\t\t{line}')


            case 'q' | 'quit':
                break

            case _:
                print("Feature not supported")


if __name__ == '__main__':
    main()
