import assemblycalculator as ac  #imported on Cradle
import json
import pandas as pd


def link_smiles(labels):
    """Links KEGG ids to smiles strings in the KEGG database
    Assumes a csv file in KeggData directory with "C" and "S" columns

    Args:
        labels (list): List of KEGG IDs

    Returns:
        list: all smiles strings associated with a list of KEGG IDs
    """
    #Open compound label file as a pandas dataframe
    kegg_df = pd.read_csv("KeggData\\KEGG_chiral_molweight_formula_labels.csv")

    smiles = []
    for l in labels:
        #Find the associated smiles string for each label
        smiles.append(kegg_df[kegg_df["C"] == l]["S"].item())

    return smiles


def find_KEGG_smiles(classes):
    """ Finds the smiles strings for a given set of KEGG BRITE classes.
    Assumes a brite.json file in KeggData directory.

    Args:
        classes (list): List of KEGG BRITE classes (e.g., "organic acids" or "lipids")

    Returns:
        list: list of all smiles strings found in KEGG associated with a
        given class
    """
    with open("KeggData\\brite.json") as json_file:
        data = json.load(json_file)

    cpd_smiles = []
    cpd_labels = []
    for c in classes:
        #Go into each cpd label (goes down 5 levels)
        for group in data["children"]:
            if group["name"] == c:
                for l in group["children"]:
                    for l2 in l["children"]:
                        for l3 in l2["children"]:
                            cpd_labels.append(l3["name"].split()[0])

        print("Size of", c, "=", len(cpd_labels))

    cpd_smiles = link_smiles(cpd_labels)

    return cpd_smiles


def main():
    #Reads in opioid data
    df = pd.read_csv("DrugData\opioid_structures.csv")
    print("Size of opioids =", len(df))
    opioid_smiles = df["Smiles"].tolist()

    print(opioid_smiles)


if __name__ == "__main__":
    main()
