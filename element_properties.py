import chemlib
import pandas as pd
import matplotlib.pyplot as plt


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

    return info_df.drop_duplicates().fillna(0)


def draw_comparison(comp_df: pd.DataFrame) -> None:
    """ Draw 3D scatterplot of comparison b/n molecules selected"""

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    ax.scatter(comp_df['AtomicRadius'], comp_df['Electronegativity'], comp_df['FirstIonization'])

    for i, label in enumerate(comp_df['Symbol']):
        ax.text(x=comp_df['AtomicRadius'][i],
                y=comp_df['Electronegativity'][i],
                z=comp_df['FirstIonization'][i],
                s=comp_df['Symbol'][i])  # Text

    ax.set_xlabel('Atomic Radius')
    ax.set_ylabel('Electronegativity')
    ax.set_zlabel('FirstIonization')

    plt.show()
