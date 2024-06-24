# Description: This file contains the code for model prediction.

import logging
import pandas as pd
import streamlit as st
import torch
import pickle

from rdkit.Chem import (
    CanonSmiles,
    MolFromSmiles,
    rdFingerprintGenerator,
    MACCSkeys,
    Descriptors,
    rdMolDescriptors,
    rdReducedGraphs,
    GraphDescriptors,
)
from mhfp.encoder import MHFPEncoder
import seaborn as sns
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

st.set_page_config(
    layout="wide",
    page_title="AMR-KG",
    page_icon=":microscope:",
    initial_sidebar_state="auto",
)

st.markdown(
    "<h1 style='text-align: center; color: #78bc1e;'>Chemical - Strain activity prediction</h1>",
    unsafe_allow_html=True,
)

st.header(
    "🔍 Compound broad spectrum activity prediction",
    divider="orange",
    help="Prediction of chemical strain activity by best model.",
)

# Data input
cols = st.columns(2)

with cols[0]:
    uploaded_file = st.file_uploader(
        "Upload a file containing SMILES",
        help="Upload a file with one SMILES per row or a comma-separated list (CSV or TSV)",
        type=["csv"],
    )
with cols[1]:
    text_input = st.text_area(
        "Or enter SMILES (one per line)",
        "",
        help="Paste yours SMILES in each line",
    )


st.write(
    "**The best model combination is Random Forest with MHFP6 fingerprint. \
    By default, this is used for prediction.**"
)

st.header(
    "🧠 Prediction results",
    divider="orange",
    help="Results of the model prediction.",
)

if st.button("Predict"):
    if uploaded_file is not None:
        smiles_df = pd.read_csv(uploaded_file, header=None)
    else:
        smiles_df = pd.DataFrame(text_input.split("\n"), columns=["smiles"])

    model_name = "rf"
    fingerprint_name = "mhfp6"

    model = torch.load(f"./models/{fingerprint_name}_{model_name}.pkl")

    # warnings.simplefilter("error", InconsistentVersionWarning)

    mhfp_encoder = MHFPEncoder(n_permutations=2048, seed=42)  # MHFP6 fingerprint

    mhfp6_fingerprints = []

    if smiles_df.empty:
        st.write("No SMILES provided.")
        st.stop()

    for smiles in smiles_df["smiles"].values:
        # Canonicalize the smiles
        try:
            can_smiles = CanonSmiles(smiles)
        except Exception as e:
            can_smiles = mhfp6_fingerprints.append(None)
            continue

        # Generate the mol object
        mol = MolFromSmiles(can_smiles)

        if not mol:
            mhfp6_fingerprints.append(None)
            continue

        mhfp6_fingerprints.append(mhfp_encoder.encode(can_smiles, radius=3))

        vals = Descriptors.CalcMolDescriptors(mol)

    smiles_df["mhfp6"] = mhfp6_fingerprints

    smiles_df_subset = smiles_df.dropna(subset=[fingerprint_name])[
        ["smiles", fingerprint_name]
    ]

    if smiles_df_subset.empty:
        st.write("No valid SMILES provided.")
        st.stop()

    predictions = model.predict(smiles_df_subset[fingerprint_name].tolist())
    prediction_proba = model.predict_proba(smiles_df_subset[fingerprint_name].tolist())
    label_classes = model.classes_.tolist()

    smiles_df_subset["Prediction"] = predictions
    probs = []
    for idx, probability in enumerate(prediction_proba):
        predicted_class = predictions[idx]
        probs.append(probability[label_classes.index(predicted_class)])
    smiles_df_subset["Probability"] = probs

    st.dataframe(smiles_df_subset[["smiles", "Prediction", "Probability"]])
    st.write("Note: The compounds that could not generate fingerprints are not shown.")
