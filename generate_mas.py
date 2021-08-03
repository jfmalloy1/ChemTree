import assemblycalculator as ac  #imported on Cradle
import json
import pandas as pd
from tqdm import tqdm


def link_smiles(labels):
    """Links KEGG ids to smiles strings in the KEGG database
    Assumes a csv file in KeggData directory with "C" and "S" columns

    Args:
        labels (list): List of KEGG IDs

    Returns:
        list: all smiles strings associated with a list of KEGG IDs
    """
    #Open compound label file as a pandas dataframe
    kegg_df = pd.read_csv("KeggData/KEGG_chiral_molweight_formula_labels.csv")

    smiles = []
    for l in labels:
        #Find the associated smiles string for each label
        smiles.append(kegg_df[kegg_df["C"] == l]["S"].item())

    return smiles


def find_KEGG_smiles(classes):
    """ Finds the smiles strings for a given set of KEGG BRITE classes.
    Assumes a brite.json file in KeggData directory.

    Args:
        classes (list): List of KEGG BRITE classes (e.g., "Organic acids" or "Lipids")

    Returns:
        list: list of all smiles strings found in KEGG associated with a
        given class
    """
    with open("KeggData/brite.json") as json_file:
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


def find_min_pathway_frags(cpd):
    """Generates all fragments from the minimal assembly theory pathway

    Args:
        cpd (string): smiles description of a chemical compound

    Returns:
        None (so far)
    """
    calc = ac.MACalculation(cpd, 30, "exact")

    calc.execute()

    print("Compound:", calc.compound)
    print("Valid mol:", calc.valid_mol)
    print("Original Compound:", calc.original_compound)
    print("Timeout:", calc.timeout)
    print("Method:", calc.method)
    #print("ncpus:", calc.ncpus)
    print("MA value from Calculation: ", calc.result)
    print("Number of minimal Pathways: ", len(calc.pathway))
    print("Full pathways:", calc.pathway)
    print("Time used:", calc.time_used)
    print("Calculation completed: ", calc.completed)
    print("Failed:", calc.failed)
    print("Modified:", calc.modified)
    # print("Fragments:", calc.fragments)
    # print("Histogram:", calc.num_frags_hist)
    # print("Path samples:", calc.path_samples)
    # print("Max atoms:", calc.max_atoms)
    print("Error message:", calc.error_message)
    #print("Full output:", calc.full_output)

    # for min_path in calc.pathway:
    #     print("Possible Pathway: ")
    #     for (frag, count) in min_path:
    #         print("Fragment ", frag, " used ", count, " times")


def main():
    #Reads in opioid data
    df = pd.read_csv("DrugData/opioid_structures.csv")
    opioid_smiles = df["Smiles"].tolist()
    # find_min_pathway_frags(opioid_smiles[0])

    #Find organic acid compounds
    kegg_smiles = find_KEGG_smiles(["Organic acids"])
    print(kegg_smiles)
    find_min_pathway_frags("CC(=O)OC1=CC=CC=C1C(=O)O")

    # for kegg_cpd in kegg_smiles:
    #     find_min_pathway_frags(kegg_cpd)

    #     print()


if __name__ == "__main__":
    main()
