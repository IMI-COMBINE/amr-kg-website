# Description: This file contains the code to display the TMAP of the chemicals in AMR-KG.

import datetime
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from rdkit.Chem import CanonSmiles, MolFromSmiles, MolToInchiKey


st.set_page_config(
    layout="wide",
    page_title="AMR-KG",
    page_icon=":microscope:",
    initial_sidebar_state="auto",
)

st.markdown(
    "<h1 style='text-align: center; color: #006c8b;'>AMR-KG Chemical Space Exploration</h1>",
    unsafe_allow_html=True,
)

st.header(
    "🔍 TMAP of the chemicals in AMR-KG Database",
    divider="orange",
    help="TMAP distribution of compounds in AMR-KG.",
)


# Add some styling with CSS selectors
st.markdown(
    """
    <style>

    a[href] {
        color: #78bc1e;
    }

    container {
        background-color: grey;
        color: white;
        border-radius: 20px;
    }
    </style>
    """,
    unsafe_allow_html=True,
)


st.markdown(
    """Explore the chemical space of the AMR-KG database using the interactive [TMAP](https://tmap.gdb.tools/) below. \
    The TMAP is a 2D representation of the chemical space of the compounds in the database. \
    Each point represents a compound, and the distance between points indicates the similarity between compounds. \
    You can zoom in, pan, and hover over the points to view the compound structure."""
)

HtmlFile = open("amrkg_chemspace.html", "r", encoding="utf-8")
source_code = HtmlFile.read()
components.html(source_code, height=500, scrolling=True)

st.markdown(
    "Want to view the interactive version locally? \
    [Download this file](https://raw.githubusercontent.com/IMI-COMBINE/broad_spectrum_prediction/main/figures/amrkg_chemspace.html)"
)


st.header(
    "👀 Search substructure in database",
    divider="orange",
    help="Look for sub-structures in the database based in InChI keys.",
)

df = pd.read_csv("data/combined_bioassay_data.tsv", sep="\t")
df["scaffold_inchikey"] = df["compound_inchikey"].str.split("-").str[0]

user_smiles = st.text_input(
    "## Search for compounds with the same scaffold as a given compound in the database. \
    Enter the SMILES of the compound you want to search for below:",
    "N[C@H]1CN(c2c(F)cc3c(=O)c(C(=O)O)cn(C4CC4F)c3c2Cl)CC12CC2",
)
can_smiles = CanonSmiles(user_smiles)
mol = MolFromSmiles(can_smiles)
inchi_key = MolToInchiKey(mol)

scaffold_inchikey = inchi_key.split("-")[0]  # remove the stereochemistry information

st.markdown(
    f"**Input: SMILES - :green[{can_smiles}] \
    \n Generated Scaffold - :green[{scaffold_inchikey}]**"
)

df_scaffold = df[df["scaffold_inchikey"] == scaffold_inchikey]
if df_scaffold.empty:
    st.write("**No compounds with this scaffold in the database.**")
else:
    st.dataframe(df_scaffold)
    st.write(f"Found :green[{df_scaffold.shape[0]}] compounds with the same scaffold.")

# last updated
date = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
st.markdown(
    f"<p style='text-align: center;'>Last updated: {date}</p>",
    unsafe_allow_html=True,
)

# footer with text and green background
st.markdown(
    "<footer style='background-color: #006c8b; padding: 10px; border-radius: 10px;'>"
    "<p style='color: white; text-align: center;'>Fraunhofer ITMP © 2021</p>"
    "<p style='color: white;'>This project has received funding from the Innovative Medicines Initiative 2 Joint Undertaking under Grant Agreement No 853967. This Joint Undertaking receives support from the European Union’s Horizon 2020 research and innovation programme and EFPIA companies’ in kind contribution.</p>"
    "</footer>",
    unsafe_allow_html=True,
)
