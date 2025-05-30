import streamlit as st
import GEOparse
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import seaborn as sns
import os
import requests
from geo_public import extract_sample_groups


# st.title("DEG Analysis from GEO Dataset")

# geo_id = st.text_input("Enter GEO Series ID (e.g., GSE7305)", value="GSE7306")

# # Load GEO dataset
# @st.cache_data
# def load_geo(geo_id):
#     try:
#         gse = GEOparse.get_GEO(geo=geo_id, destdir="./")
    
#     except Exception as e:
#         st.error(f"Failed to load GEO dataset: {e}")
#         st.stop()

#     platform = list(gse.gpls.values())[0]
#     annotation = platform.table
#     samples = list(gse.gsms.values())
#     return gse, annotation, samples

# gse, annotation, samples = load_geo(geo_id)

@st.cache_data
def load_geo(geo_id):
    def download_soft_file(geo_id, dest_dir="./data"):
        geo_prefix = geo_id[:6] + "nnn"
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_prefix}/{geo_id}/soft/{geo_id}_family.soft.gz"
        dest_path = os.path.join(dest_dir, f"{geo_id}_family.soft.gz")

        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)

        r = requests.get(url)
        if r.status_code == 200:
            with open(dest_path, "wb") as f:
                f.write(r.content)
            return dest_path
        else:
            return None

    try:
        gse = GEOparse.get_GEO(geo=geo_id, destdir="./data", how="full")
    except Exception as e:
        st.warning(f"FTP failed: {e}. Trying manual download...")
        filepath = download_soft_file(geo_id)
        if filepath:
            try:
                gse = GEOparse.get_GEO(filepath=filepath)
            except Exception as ex:
                st.error(f"Failed to parse manually downloaded file: {ex}")
                st.stop()
        else:
            st.error("Could not download GEO dataset. It may not be public.")
            st.stop()

    platform = list(gse.gpls.values())[0]
    annotation = platform.table
    samples = list(gse.gsms.values())
    return gse, annotation, samples


st.title("DEG Analysis from GEO Dataset")

geo_id = st.text_input("Enter GEO Series ID (e.g., GSE7305)", value="GSE7306")

if geo_id:
    gse, annotation, samples = load_geo(geo_id)
    st.success(f"Loaded {len(samples)} samples from {geo_id}")
    group_labels = extract_sample_groups(gse)
    # for gsm_id, gsm in gse.gsms.items():
    #     st.write(f"Sample ID: {gsm_id}")
    #     st.write(gsm.metadata)
    #     st.write("---")
    display_samples = [f"{gsm} ({label})" for gsm, label in group_labels.items()]

    # st.write("Samples with dynamic groups:")
    # st.write(display_samples)
    sample_titles = {}
    for gsm in samples:
        title = gsm.metadata.get("title", ["No title"])[0]  # safer access
        sample_titles[gsm.name] = title

    # Build markdown table with clickable sample IDs linking to GEO page
    md_table = "| Sample ID | Title |\n|---|---|\n"
    for gsm_id, title in sample_titles.items():
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm_id}"
        link = f"[{gsm_id}]({url})"
        md_table += f"| {link} | {title} |\n"

    st.markdown("### Sample Metadata")
    st.markdown(md_table, unsafe_allow_html=True)
    st.write("Annotation preview:")
    st.dataframe(annotation.head())

    # Add group selection widgets
    group_a_display = st.multiselect("Select samples for group A (e.g., Healthy)", display_samples)
    group_b_display = st.multiselect("Select samples for group B (e.g., Disease)", display_samples)

    def extract_gsm(selected_list):
        return [s.split(" (")[0] for s in selected_list]

    group_a = extract_gsm(group_a_display)
    group_b = extract_gsm(group_b_display)

    if group_a and group_b:
        st.write("Extracting expression data...")
        expression_data = gse.pivot_samples("VALUE")
        expression_data = expression_data.dropna()

        # Map gene symbols if available
        if "Symbol" in annotation.columns and "ID" in annotation.columns:
            symbol_map = annotation.set_index("ID")["Symbol"]
            expression_data["Symbol"] = expression_data.index.map(symbol_map)
        else:
            expression_data["Symbol"] = "NA"

        st.write("Performing t-test...")
        logfc = expression_data[group_b].mean(axis=1) - expression_data[group_a].mean(axis=1)
        pvals = ttest_ind(expression_data[group_b], expression_data[group_a], axis=1, equal_var=False).pvalue

        results = pd.DataFrame({
            "logFC": logfc,
            "pvalue": pvals,
            "symbol": expression_data["Symbol"]
        })
        results["-log10(pval)"] = -np.log10(results["pvalue"])
        results = results.dropna()

        st.subheader("Top Differentially Expressed Genes")
        st.dataframe(results.sort_values("pvalue").head(20))

        st.subheader("Volcano Plot")
        fig, ax = plt.subplots()
        sns.scatterplot(data=results, x="logFC", y="-log10(pval)", hue=results["pvalue"] < 0.05, ax=ax)
        plt.axhline(y=-np.log10(0.05), color='red', linestyle='--')
        st.pyplot(fig)

        selected_samples = st.multiselect("Select samples for boxplot", display_samples)

        if selected_samples:
            selected_samples_clean = [s.split(" (")[0] for s in selected_samples]
            data_for_plot = expression_data[selected_samples_clean]
            st.write(f"Data shape: {data_for_plot.shape}")
            st.write(f"Number of labels: {len(selected_samples)}")
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.boxplot(data_for_plot.values.T, labels=selected_samples, vert=True)
            ax.set_title("Expression Value Distribution per Sample")
            ax.set_ylabel("Expression Value")
            ax.set_xlabel("Samples")
            plt.xticks(rotation=45)
            plt.tight_layout()
            st.pyplot(fig)
        else:
            st.write("Select one or more samples above to see boxplot.")
    else:
        st.warning("Please select at least one sample for each group.")
