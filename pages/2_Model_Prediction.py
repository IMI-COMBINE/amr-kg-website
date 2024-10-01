# Description: This file contains the code for model prediction.

import datetime
import logging
import pandas as pd
import streamlit as st
import torch

from rdkit.Chem import (
    CanonSmiles,
    MolFromSmiles,
    rdFingerprintGenerator,
    MACCSkeys,
    rdReducedGraphs,
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
    "<h1 style='text-align: center; color: #006c8b;'>Chemical - Strain activity prediction</h1>",
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


# Model selection
st.markdown("### Select the model to use for prediction:")

cols = st.columns(2)
with cols[0]:
    st.write("#### Select the fingerprint:")
    fingerprint = st.radio(
        "Fingerprint",
        ("MHFP6", "ECFP4", "RDKIT", "MACCS", "ErG"),
        index=0,
        help="Select the fingerprint representation to use for prediction",
    )

    st.write("**The best model combination is Random Forest with MHFP6 fingerprint.**")

with cols[1]:

    metric_df = pd.read_csv("data/test_metrics.tsv", sep="\t").reset_index(drop=True)
    metric_df.rename(columns={"Unnamed: 0": "model_name"}, inplace=True)
    metric_df["model"] = metric_df["model_name"].apply(
        lambda x: (
            x.split("_")[1].upper()
            if len(x.split("_")) < 3
            else x.split("_")[-1].upper()
        )
    )
    metric_df["fingerprints"] = metric_df["model_name"].apply(
        lambda x: (
            x.split("_")[0].upper()
            if len(x.split("_")) < 3
            else x.split("_")[0].upper() + "_" + x.split("_")[1].upper()
        )
    )
    metric_df["fingerprints"] = metric_df["fingerprints"].replace(
        {
            "ERG": "ErG",
            "CHEM_PHYS": "ChemPhys",
        }
    )

    colors = {
        "MHFP6": "#3a2c20",
        "ECFP4": "#b65c11",
        "RDKIT": "#e7a504",
        "MACCS": "#719842",
        "ErG": "#3d8ebf",
        "ChemPhys": "#901b1b",
        "RF": "#3a2c20",
        "XGBOOST": "#719842",
    }

    metric_df["accuracy"] = metric_df["accuracy"] * 100
    plt.figure(figsize=(5, 5))
    sns.violinplot(x="model", y="accuracy", data=metric_df, palette=colors, hue="model")
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel("Model", fontsize=15)
    plt.ylabel("Accuracy", fontsize=15)
    st.pyplot(plt)

if st.button("Predict"):
    if uploaded_file is not None:
        smiles_df = pd.read_csv(uploaded_file, header=None)
    else:
        smiles_df = pd.DataFrame(text_input.split("\n"), columns=["smiles"])

    if smiles_df.empty:
        st.write("No SMILES provided.")
        st.stop()

    if fingerprint == "MHFP6":
        fingerprint_name = "mhfp6"
    elif fingerprint == "ECFP4":
        fingerprint_name = "ecfp4"
    elif fingerprint == "RDKIT":
        fingerprint_name = "rdkit"
    elif fingerprint == "MACCS":
        fingerprint_name = "maccs"
    else:
        fingerprint_name = "erg"

    logger.info("⏳ Loading models")

    model_name = "rf"
    model = torch.load(f"models/{fingerprint_name}_{model_name}.pkl")

    logger.info("🔮 Processing SMILES to fingerprints")

    mfpgen = rdFingerprintGenerator.GetMorganGenerator(
        radius=4, fpSize=1024
    )  # ECFP4 fingerprint
    rdkgen = rdFingerprintGenerator.GetRDKitFPGenerator(
        fpSize=1024
    )  # RDKit fingerprint
    mhfp_encoder = MHFPEncoder(n_permutations=2048, seed=42)  # MHFP6 fingerprint

    ecfp_fingerprints = []
    rdkit_fingerprints = []
    maccs_fingerprints = []
    mhfp6_fingerprints = []
    erg_fingerprints = []

    if smiles_df.empty:
        st.write("No SMILES provided.")
        st.stop()

    for smiles in smiles_df["smiles"].values:
        # Canonicalize the smiles
        try:
            can_smiles = CanonSmiles(smiles)
        except Exception:
            can_smiles = smiles

        # Generate the mol object
        mol = MolFromSmiles(can_smiles)

        if not mol:
            ecfp_fingerprints.append(None)
            rdkit_fingerprints.append(None)
            maccs_fingerprints.append(None)
            mhfp6_fingerprints.append(None)
            erg_fingerprints.append(None)
            continue

        ecfp_fingerprints.append(mfpgen.GetFingerprint(mol))
        rdkit_fingerprints.append(rdkgen.GetFingerprint(mol))
        maccs_fingerprints.append(MACCSkeys.GenMACCSKeys(mol))
        mhfp6_fingerprints.append(mhfp_encoder.encode(can_smiles, radius=3))
        erg_fingerprints.append(rdReducedGraphs.GetErGFingerprint(mol))

    smiles_df["ecfp4"] = ecfp_fingerprints
    smiles_df["rdkit"] = rdkit_fingerprints
    smiles_df["maccs"] = maccs_fingerprints
    smiles_df["mhfp6"] = mhfp6_fingerprints
    smiles_df["erg"] = erg_fingerprints

    smiles_df["mhfp6"] = mhfp6_fingerprints

    logger.info("🏃 Running model")

    smiles_df_subset = smiles_df.dropna(subset=[fingerprint_name])[
        ["smiles", fingerprint_name]
    ]

    if smiles_df.empty:
        st.write("No SMILES provided.")
        st.stop()

    predictions = model.predict(smiles_df_subset[fingerprint_name].tolist())
    prediction_proba = model.predict_proba(smiles_df_subset[fingerprint_name].tolist())
    label_classes = model.classes_.tolist()

    logger.info("✅ Finished task")

    st.write("### Predictions")
    smiles_df_subset["Prediction"] = predictions
    probs = []
    for idx, probability in enumerate(prediction_proba):
        predicted_class = predictions[idx]
        probs.append(probability[label_classes.index(predicted_class)])
    smiles_df_subset["Probability"] = probs

    st.dataframe(smiles_df_subset[["smiles", "Prediction", "Probability"]])

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
