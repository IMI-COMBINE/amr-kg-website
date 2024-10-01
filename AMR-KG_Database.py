import datetime
import pandas as pd
from collections import defaultdict
import streamlit as st
import matplotlib.pyplot as plt

st.set_page_config(
    layout="wide",
    page_title="AMR-KG",
    page_icon=":microscope:",
    initial_sidebar_state="auto",
)

# Customize sidebar
markdown = """
**Info**: This is AMR-KG repository and the broad-spectrum prediction models trained on the datasets. 
This work was done under the IMI project COMBINE.

**Developers**:
* [Yojana Gadiya](https://orcid.org/0000-0002-7683-0452)
* [Andrea Zaliani](https://orcid.org/0000-0002-1740-8390)
* [Philip Gribbon](https://orcid.org/0000-0001-7655-2459)

**External links**:
* [GitHub](https://github.com/IMI-COMBINE/broad_spectrum_prediction)
* [Website](https://www.itmp.fraunhofer.de/)
"""

st.sidebar.title("About")
st.sidebar.markdown(markdown)
st.sidebar.image("data/COMBINE_logo.jpg")

st.markdown(
    "<h1 style='text-align: center; color: #006c8b;'>AMR-KG Database</h1>",
    unsafe_allow_html=True,
)
st.markdown(
    "<h4 style='text-align: center;'>An exhaustive data warehouse of experimentally \
    validated antibacterial chemicals</h4>",
    unsafe_allow_html=True,
)

# Add some styling with CSS selectors
st.markdown(
    """
    <style>

    a[href] {
        color: #1e85bc;
    }

    sidebar .sidebar-content {
        background-color: #111 !important;
    }

    [data-testid="column"]:nth-child(1){background-color: #78bc1e;}
    [data-testid="column"]:nth-child(2){background-color: #78bc1e;}
    [data-testid="column"]:nth-child(3){background-color: #78bc1e;}

    </style>
    """,
    unsafe_allow_html=True,
)

# AMR-KG Description
st.header(
    "ℹ️ About the resources",
    divider="orange",
    help="Information on the data in AMR-KG.",
)

st.markdown(
    "Antimicrobial Resistant Knowledge Graph (AMR-KG) is an exhaustive data warehouse of experimentally validated antibacterial chemicals \
    covering Gram-positive, Gram-negative, acid-fast bacteria and fungi. The construction of the AMR-KG involved collecting \
    [minimum inhibitory concentration (MIC)](http://purl.obolibrary.org/obo/ARO_3004370) data from three different public data resources:"
)


col = st.columns(3)

with col[0]:
    container = st.container(border=True, height=280)
    container.write("### [ChEMBL](https://www.ebi.ac.uk/chembl/)")

    container.write(
        "ChEMBL (v.34) is a manually curated database of bioactive molecules with drug-like properties. \
        It brings together chemical and bioactivity to aid the translation of experimental information into effective new drugs.",
    )

with col[1]:
    container = st.container(border=True, height=280)
    container.write("### [CO-ADD](https://co-add.org/)")
    container.write(
        "Community for Open Antimicrobial Drug Discovery (CO-ADD) is a not-for-profit initiative led by academics at The University of Queensland. \
        It provides free antimicrobial screening for researchers worldwide."
    )

with col[2]:
    container = st.container(border=True, height=280)
    container.write("### [SPARK](spark.co-add.org)")
    container.write(
        "Shared Platform for Antibiotic Research (SPARK), now integrated and maintained by the CO-ADD community, was initially created by the\
        Pew Charitable Trusts to expand research around antibiotics targeting Gram-negative bacteria."
    )


# """Stats about the data"""
st.header(
    "📊 Data overview",
    divider="orange",
    help="Stats on the underlying data.",
)
df = pd.read_csv("data/combined_bioassay_data.tsv", sep="\t")


def get_base_stats():

    chembl_cmpds = set(
        df[df["compound_source"] == "chembl_34"]["compound_inchikey"].unique()
    )
    coadd_cmpds = set(
        df[df["compound_source"] == "coadd_03_01-02-2020"]["compound_inchikey"].unique()
    )
    spark_cmpds = set(
        df[df["compound_source"] == "spark"]["compound_inchikey"].unique()
    )

    pchem_dist_dict = defaultdict(list)

    for idx, row in df.iterrows():
        (
            inchikey,
            smiles,
            source,
            gram_pos,
            gram_neg,
            fungi,
            acid_fast,
            _,
            _,
            _,
            gp_class,
            gn_class,
            fungi_class,
            af_class,
            best_class,
        ) = row

        if pd.notna(gp_class):
            pchem_dist_dict["gram-positive"].append(gram_pos)
        if pd.notna(gn_class):
            pchem_dist_dict["gram-negative"].append(gram_neg)
        if pd.notna(fungi_class):
            pchem_dist_dict["fungi"].append(fungi)
        if pd.notna(af_class):
            pchem_dist_dict["acid-fast"].append(acid_fast)

    return (chembl_cmpds, coadd_cmpds, spark_cmpds), pchem_dist_dict


(chembl_cmpds, coadd_cmpds, spark_cmpds), pchem_dist_dict = get_base_stats()

fig = plt.figure(figsize=(15, 10))

plt.subplot(2, 2, 1)
plt.bar(
    [
        "ChEMBL",
        "SPARK",
        "CO-ADD",
    ],
    [len(chembl_cmpds), len(spark_cmpds), len(coadd_cmpds)],
)
plt.title("Compound distribution across resources", fontsize=15, fontweight="bold")
plt.yscale("log")
# display the number of compounds on the bars
for i, v in enumerate([len(chembl_cmpds), len(spark_cmpds), len(coadd_cmpds)]):
    plt.text(i, v + 10, str(v), ha="center", va="bottom")
plt.ylabel("Number of deduplicated compounds", fontsize=12)
plt.xlabel("Data source", fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.subplot(2, 2, 2)
plt.hist(
    pchem_dist_dict["gram-positive"], alpha=0.5, color="red", label="gram-positive"
)
plt.hist(
    pchem_dist_dict["gram-negative"], alpha=0.5, color="blue", label="gram-negative"
)
plt.hist(pchem_dist_dict["fungi"], alpha=0.5, color="green", label="fungi")
plt.hist(pchem_dist_dict["acid-fast"], alpha=0.5, color="orange", label="acid-fast")
plt.axvline(5, color="red", linestyle="--", label="pChEMBL threshold")
plt.legend(title="Pathogen class", fontsize=12)
plt.title("Distribution of pChEMBL values", fontsize=15, fontweight="bold")
plt.ylabel("Number of deduplicated compounds", fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel("pChEMBL value", fontsize=12)

st.pyplot(fig)


st.header(
    ":arrow_down: Data Download",
    divider="orange",
    help="Downloading the compound-pathogen data in AMR-KG.",
)


@st.cache_data
def convert_df(df):
    return df.to_csv(index=False, sep="\t").encode("utf-8")


csv = convert_df(df)
st.write("The data files contains the following columns in the tab-seperated manner:")
st.markdown(
    """
    * `compound_inchikey` - InChI key of the compound \n
    * `compound_smiles` - SMILES representation of the compound \n
    * `compound_source` - Source database of the compound \n
    * `gram_positive` - pMIC value for compound activity against Gram-positive pathogen \n
    * `gram_negative` - pMIC value for compound activity against Gram-negative pathogen \n
    * `fungi` - pMIC value for compound activity against fungi \n
    * `acid_fast` - pMIC value for compound activity against acid-fast pathogen \n
    * `chemical_class` - Chemical class disitribution of the compounds based on NPClassifier \n
    * `compound_superclass` - Superclass of the compound based on NPClassifier \n
    * `compound_pathway` - Pathway of the compound based on NPClassifier \n
    * `gram_positive_label` - Indicating whether a compound is active or inactive for a Gram-positive pathogen based on pChEMBL threshold (i.e. 5) \n
    * `gram_negative_label` - Indicating whether a compound is active or inactive for a Gram-negative pathogen based on pChEMBL threshold (i.e. 5) \n
    * `fungi_label` - Indicating whether a compound is active or inactive for a Fungi pathogen based on pChEMBL threshold (i.e. 5) \n
    * `acid_fast_label` - Indicating whether a compound is active or inactive for a Acid-fast pathogen based on pChEMBL threshold (i.e. 5) \n
    * `best_class` - Pathogen class of the compound based on the highest pChEMBL values \n
    """
)
st.dataframe(df.head(3))
st.download_button(
    "Press to Download", csv, "amrkg_data_dump.tsv", "text/tsv", key="download-tsv"
)

# Publucation note
with st.expander(
    "If you have found our resource or model useful in your work, please consider citing us: "
):
    st.write("""Manuscript in preparation. Please check back soon for more details.""")


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
