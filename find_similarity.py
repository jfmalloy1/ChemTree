from rdkit import Chem
from rdkit import DataStructs
import json
import pandas as pd
import numpy as np
import itertools
import pickle

""" Link KEGG cpd labels with smiles strings
    Input: list of compound labels
    Output: list of smiles strings
"""
def link_smiles(labels):
    #Open compound label file as a pandas dataframe
    kegg_df = pd.read_csv("KeggData\\KEGG_chiral_molweight_formula_labels.csv")

    smiles = []
    for l in labels:
        #Find the associated smiles string for each label
        smiles.append(kegg_df[kegg_df["C"] == l]["S"].item())

    return smiles

""" Find the smiles strings associated with specific classes of KEGG
    Input: cs, list of KEGG classes (lipids, organic acids, etc...)
    Output: smiles strings belonging to those classes, ordered by placement in input list
"""
def find_KEGG_smiles(classes):
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

""" Find smiles strings associated with opioids (obtained via Drugbank)
    Input: None (assumes data is in 'DrugData' folder)
    Output: list of smiles strings
"""
def find_opioid_smiles():
    df = pd.read_csv("DrugData\opioid_structures.csv")
    print("Size of opioids =", len(df))
    return df["Smiles"].tolist()

def main():
    # #Test from rdkit tutorial
    # ms = [Chem.MolFromSmiles('CCOC'), Chem.MolFromSmiles('CCO'), Chem.MolFromSmiles('COC')]
    # fps = [Chem.RDKFingerprint(x) for x in ms]
    # print(DataStructs.FingerprintSimilarity(fps[0],fps[1]))
    # print(DataStructs.FingerprintSimilarity(fps[0],fps[2]))
    # print(DataStructs.FingerprintSimilarity(fps[1],fps[2]))
    # print(DataStructs.FingerprintSimilarity(fps[2],fps[1]))

    #Goal: build a tanimoto similarity matrix between lipids & organic acids

    #Get KEGG cpds (input is a list of BRITE classes)
    cpds = []
    cpds.append(find_KEGG_smiles(["Organic acids"]))
    #Get opioids
    cpds.append(find_opioid_smiles())
    #Flatten list
    cpds = sum(cpds, [])

    #Convert to RDKFingerprints
    fps = []
    fps = [Chem.RDKFingerprint(Chem.MolFromSmiles(c)) for c in cpds]

    #Build & save tanimoto matrix
    matrix = np.zeros((len(fps), len(fps)))

    combinations = itertools.combinations(np.arange(0,len(fps),1), 2)
    for c in combinations:
        matrix[c[0]][c[1]] = DataStructs.FingerprintSimilarity(fps[c[0]],fps[c[1]])

    print(matrix)
    pickle.dump(matrix, open("tanimoto_matrixOpioidsOrgAcids.p", "wb"))

    #Data viz things (probably will use jupyter for this)


if __name__ == "__main__":
    main()
