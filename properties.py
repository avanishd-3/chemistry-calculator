import pandas as pd
import chemlib


def return_prop_df(input_list: list) -> pd.DataFrame:
    """ Return dataframe containing Atomic Radius, Electronegativity, FirstIonization
        of all elements in input list"""

    elements_list = []

    for obj in input_list:
        try:
            element = chemlib.Element(obj)
        except IndexError:
            continue
        else:
            elements_list.append(element.properties)

    info_df = pd.DataFrame(elements_list)[['AtomicRadius', 'Electronegativity', 'FirstIonization', 'Symbol']]

    return info_df.drop_duplicates()
