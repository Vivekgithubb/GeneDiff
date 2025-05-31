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
from sklearn.decomposition import PCA
import io
from google.cloud import vision


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

def plot_pca(expr_df, sample_labels):
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(expr_df.T)
    pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])
    pca_df['Sample'] = expr_df.columns
    # Extract group ("Healthy", "Disease", "Unknown") from sample_labels
    def extract_group(label):
        if "Healthy" in label:
            return "Healthy"
        elif "Disease" in label:
            return "Disease"
        else:
            return "Unknown"
    pca_df['Group'] = pca_df['Sample'].map(lambda x: extract_group(sample_labels.get(x, "")))
    palette = {"Healthy": "blue", "Disease": "red", "Unknown": "gray"}
    fig, ax = plt.subplots()
    sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='Group',palette=palette, s=50, ax=ax)
    ax.set_title('PCA of Samples')
    ax.legend(fontsize=8, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    st.pyplot(fig)

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
geo_options = {
"Breast Cancer vs. Normal (GSE15852)": "GSE15852",
"Cervical Cancer vs. Normal (GSE6791)": "GSE6791",
"Multiple Cancers vs. Normal (GSE5364)": "GSE5364",
"Ductal Carcinoma in Situ (GSE10950)": "GSE10950",
"Lung Cancer vs. Normal (GSE19188)": "GSE19188",
"Bladder Cancer (GSE109169)": "GSE109169",
"Prostate Cancer (GSE3325)": "GSE3325"
}
selected_label = st.selectbox("Select a GEO Dataset for Analysis", list(geo_options.keys()))
geo_id = geo_options[selected_label]

st.write(f"Selected GEO ID: {geo_id}")

# geo_id = st.text_input("Enter GEO Series ID (e.g., GSE7305)", value="GSE10950")

if geo_id:
    gse, annotation, samples = load_geo(geo_id)
    st.success(f"Loaded {len(samples)} samples from {geo_id}")
    sample_labels,  group_labels = extract_sample_groups(gse)
    # for gsm_id, gsm in gse.gsms.items():
    #     st.write(f"Sample ID: {gsm_id}")
    #     st.write(gsm.metadata)
    #     st.write("---")
    display_samples = [f"{gsm} ({label})" for gsm, label in group_labels.items()]

    # st.write("Samples with dynamic groups:")
    # st.write(display_samples)
    sample_titles = {}
    # for gsm in samples:
    #     title = gsm.metadata.get("title", ["No title"])[0]  # safer access
    #     sample_titles[gsm.name] = title

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
    st.subheader("Assign Groups")
    # Correct unpacking
    sample_labels, group_dict = extract_sample_groups(gse)

# sample_labels is dict: {gsm_id: "gsm_id (label)"}
# group_dict is dict: {label: [gsm_id list]}

    display_samples = list(sample_labels.values())

    group_a_display = st.multiselect(
        "Select samples for Group A (e.g., Healthy)",
        options=[sample_labels[gsm] for gsm in group_dict.get("Healthy", [])]
    )

    group_b_display = st.multiselect(
        "Select samples for Group B (e.g., Disease)",
        options=[sample_labels[gsm] for gsm in group_dict.get("Disease", [])]
    )

    def extract_gsm(selected_list):
        return [s.split(" ")[0] for s in selected_list]

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
        
                # --- Heatmap of Differential Expression ---
        st.subheader("Heatmap of Top Differentially Expressed Genes")
        num_genes = st.slider("Number of top genes to show in heatmap", min_value=10, max_value=50, value=20, step=1)

        # Get top N genes by p-value
# Get top N genes by p-value (use index if symbol is missing)
        top_gene_indices = results.sort_values("pvalue").head(num_genes).index

        # Get expression data for these genes (use all samples)
        heatmap_data = expression_data.loc[top_gene_indices]

        # Remove the Symbol column for plotting if it exists
        if "Symbol" in heatmap_data.columns:
            heatmap_matrix = heatmap_data.drop(columns=["Symbol"])
        else:
            heatmap_matrix = heatmap_data.copy()

        # Z-score normalize each gene (row) for better visualization
        if not heatmap_matrix.empty:
            heatmap_matrix = heatmap_matrix.sub(heatmap_matrix.mean(axis=1), axis=0)
            heatmap_matrix = heatmap_matrix.div(heatmap_matrix.std(axis=1), axis=0)
            # Set probe IDs as row labels
            heatmap_matrix.index = top_gene_indices

            fig, ax = plt.subplots(figsize=(min(1.5 + 0.3 * len(heatmap_matrix.columns), 20), min(0.5 * num_genes, 20)))
            sns.heatmap(heatmap_matrix, cmap="vlag", ax=ax, cbar_kws={'label': 'Z-score'})
            ax.set_title(f"Top {num_genes} Differentially Expressed Probes (Z-score)")
            ax.set_xlabel("Sample")
            ax.set_ylabel("Probe ID")
            st.pyplot(fig)
        else:
            st.warning("No valid genes found for heatmap. Try increasing the number of top genes or check your data.")


        selected_samples = st.multiselect("Select samples for boxplot", display_samples)

        # if selected_samples:
        #     selected_samples_clean = [s.split(" (")[0] for s in selected_samples]
        #     data_for_plot = expression_data[selected_samples_clean]
        #     # Ensure labels match columns
        #     labels = [sample_labels[gsm] for gsm in selected_samples_clean]
        #     st.write(f"Data shape: {data_for_plot.shape}")
        #     st.write(f"Number of labels: {len(labels)}")
        #     fig, ax = plt.subplots(figsize=(10, 6))
        #     ax.boxplot(data_for_plot.values, labels=labels, vert=True)
        #     ax.set_title("Expression Value Distribution per Sample")
        #     ax.set_ylabel("Expression Value")
        #     ax.set_xlabel("Samples")
        #     plt.xticks(rotation=45)
        #     plt.tight_layout()
        #     st.pyplot(fig)

        if selected_samples:
            selected_samples_clean = [s.split(" (")[0] for s in selected_samples]
            data_for_plot = expression_data[selected_samples_clean]
            labels = [sample_labels[gsm] for gsm in selected_samples_clean]
            st.write(f"Data shape: {data_for_plot.shape}")
            st.write(f"Number of labels: {len(labels)}")
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.boxplot(data_for_plot.values, labels=labels, vert=True)
            ax.set_title("Expression Value Distribution per Sample")
            ax.set_ylabel("Expression Value")
            ax.set_xlabel("Samples")
            plt.xticks(rotation=45)
            plt.tight_layout()
            st.pyplot(fig)
        else:
            st.write("Select one or more samples above to see boxplot.")
        if st.button("Show PCA Plot"):
                 # Use only numeric columns (samples)
                st.markdown("### PCA Plot: Sample Clustering")
                st.write("This plot shows how samples cluster based on their gene expression profiles using Principal Component Analysis (PCA).")
                expr_df = expression_data.drop(columns=["Symbol"], errors="ignore")
                plot_pca(expr_df, sample_labels)
    else:
        st.warning("Please select at least one sample for each group.")
