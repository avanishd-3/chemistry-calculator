import math
import sympy
from sympy.abc import x

import matplotlib.pyplot as plt


class Titration:
    """
    Parameters

      initial_volume: (float) The volume of the acid/base solution in milliliters (mL).
      initial_concentration: (float) The concentration of the acid/base solution in molarity (M).
      initial_ka: (float) The Ka/Kb value of the acid/base
      other_concentration: (float) The concentration of the base/acid solution in molarity (M).
      other_k: (float) The Kb/Ka value of the base/acid
      starting: (str) Whether the reaction starts with an acid or a base

    """

    def __init__(self, acid_volume, acid_conc, k_a, base_conc, starting):

        # Data values
        self._initial_volume = acid_volume  # In milliliters
        self._initial_concentration = acid_conc  # In M

        self._initial_k = k_a

        self._other_concentration = base_conc  # In M

        self._other_k = 1e-14 / self._initial_k

        self._starting = starting

    def _initial_ph(self) -> dict:
        """ Calculate and return [H+], pH, [OH-], pOh in initial state"""
        if self._initial_k < 1:
            ion_concentration = [i for i in sympy.solve(((x ** 2) /
                                 (self._initial_concentration) - self._initial_k), x) if i > 0][0]
        else:
            ion_concentration = self._initial_concentration

        if self._starting == 'a' or self._starting == 'acid':
            ph_val = -math.log10(ion_concentration)

            return {"Initial [H+]": ion_concentration, "Initial pH": ph_val,
                    "Initial [OH-]": 1e-14 / ion_concentration, "Initial pOH": 14 - ph_val}
        else:
            poh_val = -math.log10(ion_concentration)

            return {"Initial [H+]": 1e-14 / ion_concentration, "Initial pH": 14 - poh_val,
                    "Initial [OH-]": ion_concentration, "Initial pOH": poh_val}

    def _equivalence_point(self) -> dict:
        """ Calculate and return equivalence point, [H+], pH, [OH-], pOh at equivalence point"""
        mol_acid = (self._initial_volume / 1000) * self._initial_concentration

        eq_pt = (mol_acid / self._other_concentration) * 1000

        conc = mol_acid / (self._initial_volume / 1000 + eq_pt / 1000)  # In M

        ion_concentration = [i for i in sympy.solve(((x ** 2) / (conc) - self._other_k), x) if i > 0][0]

        if self._starting == 'a' or self._starting == 'acid':
            if self._initial_k > 1 and self._other_k > 1:
                ph_val = 7
            else:
                ph_val = 14 + math.log10(ion_concentration)  # 14 - pOH

            return {"Equivalence point": f"{eq_pt} mL", "Equivalence [H+]": 1e-14 / ion_concentration,
                    "Equivalence pH": ph_val, "Equivalence [OH-]": ion_concentration,
                    "Equivalence pOH": 14 - ph_val}
        else:
            if self._initial_k > 1 and self._other_k > 1:
                poh_val = 7
            else:
                poh_val = 14 + math.log10(ion_concentration)  # 14 - pH

            return {"Equivalence point": f"{eq_pt} mL", "Equivalence [H+]": ion_concentration,
                    "Equivalence pH": 14 - poh_val, "Equivalence [OH-]": 1e-14 / ion_concentration,
                    "Equivalence pOH": poh_val}

    def _halfway_point(self) -> dict:
        """ Calculate and return halfway point, [H+], pH, [OH-], pOh in initial state"""
        eq_pt_info = self._equivalence_point()['Equivalence point']

        eq_pt = float(eq_pt_info.split(' ')[0])

        halfway_pt = eq_pt / 2

        if self._starting == 'a' or self._starting == 'acid':

            ph_val = -math.log10(self._initial_k)
            h_ion_concentration = 10 ** (-ph_val)

            return {"Halfway point": f"{halfway_pt} mL", "Halfway [H+]": h_ion_concentration,
                    "Halfway pH": ph_val, "Halfway [OH-]": 1e-14 / h_ion_concentration,
                    "Halfway pOH": 14 - ph_val}
        else:

            poh_val = -math.log10(self._initial_k)
            oh_ion_concentration = 10 ** (-poh_val)

            return {"Halfway point": f"{halfway_pt} mL", "Halfway [H+]": 1e-14 / oh_ion_concentration,
                    "Halfway pH": 14 - poh_val, "Halfway [OH-]": oh_ion_concentration,
                    "Halfway pOH": poh_val}

    def _before_and_after_equivalence_point(self) -> list[tuple]:
        """ Calculate and return pH at 20 points before reaching equivalence point"""
        eq_pt_info = self._equivalence_point()['Equivalence point']

        eq_pt = float(eq_pt_info.split(' ')[0])

        amount_of_pts = 20

        iterate_amount = eq_pt / amount_of_pts
        init = 0

        mol_init = (self._initial_volume / 1000) * self._initial_concentration

        data_list = []

        # Before equivalence point
        for i in range(amount_of_pts - 1):  # Exclude equivalence point
            init += iterate_amount

            mol_other = (init / 1000) * self._other_concentration

            molarity_init = (mol_init - mol_other) / ((self._initial_volume / 1000) + (init / 1000))

            molarity_other = mol_other / ((self._initial_volume / 1000) + (init / 1000))

            if self._starting == 'a' or self._starting == 'acid':
                ph_val = -math.log10(self._initial_k) + math.log10(molarity_other / molarity_init)

            else:
                ph_val = -math.log10(self._other_k) + math.log10(molarity_init / molarity_other)

            data_list.append((init, ph_val))

        # After equivalence point
        init = eq_pt

        eq_mol = (eq_pt * self._other_concentration) / 1000

        for i in range(amount_of_pts - 1):  # Keep equal to before equivalence point
            init += iterate_amount

            mol_other = (init / 1000) * self._other_concentration

            if self._starting == 'a' or self._starting == 'acid':
                oh_concentration = (mol_other - eq_mol) / ((self._initial_volume / 1000) + (init / 1000))

                ph_val = 14 + math.log10(oh_concentration)
            else:
                h_concentration = (mol_other - eq_mol) / ((self._initial_volume / 1000) + (init / 1000))

                ph_val = -math.log10(h_concentration)

            data_list.append((init, ph_val))

        return data_list

    def titration_info(self) -> dict:
        titration_info_dict = dict()

        initial = self._initial_ph()

        equivalence = self._equivalence_point()

        halfway = self._halfway_point()

        titration_info_dict.update(initial)
        titration_info_dict.update(halfway)
        titration_info_dict.update(equivalence)

        return titration_info_dict

    def plot_titration_curve(self) -> None:
        initial_ph_val = self._initial_ph()['Initial pH']

        initial_pt = (0, initial_ph_val)

        eq_pt_info = self._equivalence_point()['Equivalence point']

        eq_volume = float(eq_pt_info.split(' ')[0])

        eq_ph_val = self._equivalence_point()['Equivalence pH']

        eq_pt = (eq_volume, eq_ph_val)

        halfway_pt_info = self._halfway_point()['Halfway point']

        halfway_volume = float(halfway_pt_info.split(' ')[0])

        halfway_ph_val = self._halfway_point()['Halfway pH']

        halfway_pt = (halfway_volume, halfway_ph_val)

        data_list = self._before_and_after_equivalence_point()

        # Add initial point, halfway point, equivalence point

        data_list.append(initial_pt)
        data_list.append(eq_pt)
        data_list.append(halfway_pt)

        x_points = [i[0] for i in data_list]
        y_points = [i[-1] for i in data_list]

        plt.plot(x_points, y_points, 'o')

        plt.xlabel('mL added')
        plt.ylabel('pH')
        plt.title('Titration Curve')

        plt.show()


if __name__ == "__main__":  # For testing
    t = Titration(25, 0.15, 1.8e-5, 0.10, 'a')

    t.plot_titration_curve()
