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

# st.markdown(
#     """
#     <style>
#     body {
#         background-color: #ffffff;
#     }
#     .main .block-container {
#         background-color: #ffffff;
#         border-radius: 12px;
#         padding: 2rem 2rem 2rem 2rem;
#         box-shadow: 0 2px 8px rgba(0,0,0,0.05);
#     }
#     header, .css-18e3th9, .css-1d391kg { 
#         background-color: #e6f2ff !important;
#     }
#     /* Optional: Style the sidebar*/ 
#     .css-1d391kg {
#         background-color: #cce0ff !important;
#     }
#     </style>
#     """,
#     unsafe_allow_html=True
# )
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

def plot_ma(results, expression_data, group_a, group_b):
    # Calculate mean expression for each group
    mean_group_a = expression_data[group_a].mean(axis=1)
    mean_group_b = expression_data[group_b].mean(axis=1)
    mean_expr = (mean_group_a + mean_group_b) / 2

    ma_df = results.copy()
    ma_df['mean_expr'] = mean_expr

    fig, ax = plt.subplots()
    sns.scatterplot(
        data=ma_df,
        x='mean_expr',
        y='logFC',
        hue=ma_df['pvalue'] < 0.05,
        palette={True: 'red', False: 'gray'},
        ax=ax
    )
    ax.axhline(0, color='black', linestyle='--')
    ax.set_xlabel('Mean Expression')
    ax.set_ylabel('log2 Fold Change')
    ax.set_title('MA Plot')
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
# Sidebar for dataset selection
st.sidebar.title("Navigation")
selected_label = st.sidebar.selectbox("Select a GEO Dataset for Analysis", list(geo_options.keys()))
geo_id = geo_options[selected_label]

st.write(f"Selected GEO ID: {geo_id}")

# geo_id = st.text_input("Enter GEO Series ID (e.g., GSE7305)", value="GSE10950")

if geo_id:
    gse, annotation, samples = load_geo(geo_id)
    st.success(f"Loaded {len(samples)} samples from {geo_id}")
    sample_labels, group_labels = extract_sample_groups(gse)
    display_samples = [f"{gsm} ({label})" for gsm, label in group_labels.items()]
    sample_titles = {gsm.name: gsm.metadata.get("title", ["No title"])[0] for gsm in samples}

    # Prepare for group selection
    sample_labels, group_dict = extract_sample_groups(gse)
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

    # Only proceed if both groups are selected
    if group_a and group_b:
        expression_data = gse.pivot_samples("VALUE").dropna()
        if "Symbol" in annotation.columns and "ID" in annotation.columns:
            symbol_map = annotation.set_index("ID")["Symbol"]
            expression_data["Symbol"] = expression_data.index.map(symbol_map)
        else:
            expression_data["Symbol"] = "NA"
        logfc = expression_data[group_b].mean(axis=1) - expression_data[group_a].mean(axis=1)
        pvals = ttest_ind(expression_data[group_b], expression_data[group_a], axis=1, equal_var=False).pvalue
        results = pd.DataFrame({
            "logFC": logfc,
            "pvalue": pvals,
            "symbol": expression_data["Symbol"]
        })
        results["-log10(pval)"] = -np.log10(results["pvalue"])
        results = results.dropna()

        # --- TABS ---
        tabs = st.tabs([
            "Sample Metadata", "Sample Groups Overview", "Annotation Table", "DEG Table", "Volcano Plot",
            "Heatmap", "Boxplot", "PCA Plot", "MA Plot"
        ])

        # 1. Sample Metadata
        with tabs[0]:
            st.subheader("Sample Metadata Table (with clickable links)")
            st.markdown("""
**Summary:**  
This section provides a comprehensive overview of all samples in your selected GEO dataset. Each sample is listed with its unique GEO accession ID (which links directly to the GEO database) and its descriptive title.

**Table Columns:**  
- **Sample ID:** The unique identifier for each sample in GEO. Click to view the sample’s full metadata on the GEO website.
- **Title:** The descriptive title provided by the original study authors.

**Assumptions:**  
- Sample titles and IDs are accurate as provided by GEO.
- Clicking a Sample ID will open the corresponding GEO page in your browser.
            """)

            # Clickable links table
            md_table = "| Sample ID | Title |\n|---|---|\n"
            for gsm_id, title in list(sample_titles.items()):
                url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm_id}"
                link = f"[{gsm_id}]({url})"
                md_table += f"| {link} | {title} |\n"
            st.markdown(md_table, unsafe_allow_html=True)

        # 2. Sample Groups Overview
        with tabs[1]:
            st.subheader("Sample Groups Overview")
            st.markdown("""
**Summary:**  
This table groups all sample IDs by their inferred biological status.

**Groups:**  
- **Healthy:** Samples labeled as healthy, normal, or control.
- **Disease/Unhealthy:** Samples labeled as disease, cancer, tumor, or unknown.

**How grouping is determined:**  
- If the sample label contains "Healthy", it is classified as Healthy.
- All other samples are classified as Unhealthy.

**Assumptions:**  
- The group assignment is based on keyword matching in sample metadata.
- Some samples may be misclassified if their metadata is ambiguous.
            """)
            # Grouped sample IDs table
            group_dict_display = {"Healthy": [], "Unhealthy": []}
            for gsm_id, label in sample_labels.items():
                if "Healthy" in label:
                    group_dict_display["Healthy"].append(gsm_id)
                else:
                    group_dict_display["Unhealthy"].append(gsm_id)

            group_md = "| Group | Sample IDs |\n|---|---|\n"
            for group, ids in group_dict_display.items():
                if ids:
                    id_links = ", ".join([f"[{gsm_id}](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm_id})" for gsm_id in ids])
                    group_md += f"| {group} | {id_links} |\n"
            st.markdown(group_md, unsafe_allow_html=True)

        # 3. Annotation Table
        with tabs[2]:
            st.subheader("Annotation Table Preview")
            st.markdown("""
**Summary:**  
The annotation table links each probe (measurement spot on the microarray) to biological information such as gene names and genomic locations.

**Table Columns (may vary by platform):**  
- **ID:** The probe identifier (unique for each spot on the array).
- **Symbol:** The gene symbol (e.g., TP53) that the probe measures.
- **Other columns:** May include gene title, chromosome, or database cross-references.

**Assumptions:**  
- The mapping between probe IDs and gene symbols is correct as provided by the platform annotation.
- Some probes may not map to a known gene (symbol may be NA).

**Why this matters:**  
This table is essential for interpreting your results in terms of real genes and their biological functions.
            """)
            st.dataframe(annotation.head())

        # 4. DEG Table
        with tabs[3]:
            st.subheader("Top Differentially Expressed Genes")
            st.markdown("""
**Summary:**  
This table lists the top 20 genes (or probes) with the most statistically significant differences in expression between your selected groups.

**Table Columns:**  
- **logFC (Log2 Fold Change):**  
  - The log2 ratio of average expression between Group B (e.g., Disease) and Group A (e.g., Healthy).
  - Positive values: Higher in Group B; Negative: Higher in Group A.
- **pvalue:**  
  - The probability that the observed difference is due to random chance (from a t-test).
  - Lower values indicate more significant differences.
- **symbol:**  
  - The gene symbol for each probe.
- **-log10(pval):**  
  - The negative log10 of the p-value, used for easier visualization (higher = more significant).

**Assumptions:**  
- Expression values are approximately normally distributed for each gene.
- The t-test is appropriate for comparing the two groups.
- Multiple testing correction is not shown here (raw p-values).

**How to use:**  
Focus on genes with large absolute logFC and low p-values for biological interpretation.
            """)
            st.dataframe(results.sort_values("pvalue").head(20))

        # 5. Volcano Plot
        with tabs[4]:
            st.subheader("Volcano Plot")
            st.markdown("""
**Summary:**  
A volcano plot is a type of scatter plot that helps you quickly identify genes that are both statistically significant and have large fold changes between two groups.

**Axes:**  
- **X-axis (logFC):**  
  - Log2 fold change between groups.
  - Positive values: Higher expression in Group B (e.g., Disease).
  - Negative values: Higher expression in Group A (e.g., Healthy).
- **Y-axis (-log10(p-value)):**  
  - The negative log10 of the p-value for each gene.
  - Higher values mean more statistically significant differences.

**Colors:**  
- **Red points:** Genes with p-value < 0.05 (statistically significant).
- **Gray points:** Genes not meeting the significance threshold.

**What does this plot show?**  
- Genes far from zero on the x-axis have large differences in expression between groups.
- Genes high on the y-axis are more statistically significant.
- Genes that are both far from zero (large logFC) and high (low p-value) are the most interesting biologically.

**Assumptions:**  
- Each gene is tested independently.
- The p-value threshold (0.05) is used for significance, but multiple testing correction is not applied here.

**How to use:**  
- Look for red points far from zero on the x-axis and high on the y-axis—these are your top candidate genes.
- Use this plot to prioritize genes for further study or validation.
            """)
            fig, ax = plt.subplots()
            sns.scatterplot(data=results, x="logFC", y="-log10(pval)", hue=results["pvalue"] < 0.05, ax=ax)
            plt.axhline(y=-np.log10(0.05), color='red', linestyle='--')
            st.pyplot(fig)

        # 6. Heatmap
        with tabs[5]:
            st.subheader("Heatmap of Top Differentially Expressed Genes")
            st.markdown("""
**Summary:**  
A heatmap is a graphical representation of data where individual values are represented as colors. Here, it shows the expression patterns of the top N most differentially expressed genes across all samples in your dataset.

**Axes:**  
- **Rows:** Each row represents a gene or probe (specifically, the top N genes ranked by statistical significance).
- **Columns:** Each column represents a sample from your dataset.

**Colors:**  
- **Red (or lighter colors):** Higher-than-average expression for that gene in that sample (Z-score > 0).
- **Blue (or darker colors):** Lower-than-average expression for that gene in that sample (Z-score < 0).
- The color scale is based on Z-score normalization, which centers and scales each gene’s expression values so you can compare patterns across genes.

**What does this plot show?**  
- **Clusters of samples:** If columns (samples) group together with similar color patterns, it suggests those samples have similar gene expression profiles—often reflecting biological similarity (e.g., all Disease samples clustering together).
- **Clusters of genes:** If rows (genes) group together, it suggests those genes behave similarly across samples, which may indicate shared biological function.
- **Outliers:** Samples or genes with very different color patterns may be outliers or have unique biology.

**Assumptions:**  
- Only the top N genes (by p-value) are shown for clarity and interpretability.
- Expression values are Z-score normalized for each gene (row), so differences reflect relative, not absolute, expression.

**How to use:**  
- Look for blocks of color that indicate groups of samples or genes with similar expression.
- Check if your sample groups (e.g., Healthy vs. Disease) cluster together, which supports biological differences.
- Investigate outliers or unexpected patterns for possible technical issues or novel findings.
            """)
            num_genes = st.slider("Number of top genes to show in heatmap", min_value=10, max_value=50, value=20, step=1)
            top_gene_indices = results.sort_values("pvalue").head(num_genes).index
            heatmap_data = expression_data.loc[top_gene_indices]
            if "Symbol" in heatmap_data.columns:
                heatmap_matrix = heatmap_data.drop(columns=["Symbol"])
            else:
                heatmap_matrix = heatmap_data.copy()
            if not heatmap_matrix.empty:
                heatmap_matrix = heatmap_matrix.sub(heatmap_matrix.mean(axis=1), axis=0)
                heatmap_matrix = heatmap_matrix.div(heatmap_matrix.std(axis=1), axis=0)
                heatmap_matrix.index = top_gene_indices
                fig, ax = plt.subplots(figsize=(min(1.5 + 0.3 * len(heatmap_matrix.columns), 20), min(0.5 * num_genes, 20)))
                sns.heatmap(heatmap_matrix, cmap="vlag", ax=ax, cbar_kws={'label': 'Z-score'})
                ax.set_title(f"Top {num_genes} Differentially Expressed Probes (Z-score)")
                ax.set_xlabel("Sample")
                ax.set_ylabel("Probe ID")
                st.pyplot(fig)
            else:
                st.warning("No valid genes found for heatmap. Try increasing the number of top genes or check your data.")

        # 7. Boxplot
        with tabs[6]:
            st.subheader("Boxplot")
            st.markdown("""
**Summary:**  
The boxplot displays the distribution of expression values for each selected sample.

**Axes:**  
- **X-axis:** Selected samples (each box is one sample).
- **Y-axis:** Expression value for all genes in that sample.

**Boxplot elements:**  
- **Box:** Interquartile range (middle 50% of values).
- **Line in box:** Median expression value.
- **Whiskers:** Range of most values (excluding outliers).
- **Dots:** Outlier values.

**Assumptions:**  
- Expression values are comparable across samples.
- Outliers may indicate technical artifacts or true biological variation.

**How to use:**  
Compare the spread and median of each sample. Large differences may indicate batch effects or sample quality issues.
            """)
            selected_samples = st.multiselect("Select samples for boxplot", display_samples)
            if selected_samples:
                selected_samples_clean = [s.split(" (")[0] for s in selected_samples]
                data_for_plot = expression_data[selected_samples_clean]
                labels = [sample_labels[gsm] for gsm in selected_samples_clean]
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

        # 8. PCA Plot
        with tabs[7]:
            st.subheader("PCA Plot: Sample Clustering")
            st.markdown("""
**Summary:**  
Principal Component Analysis (PCA) is a statistical technique that reduces the complexity of high-dimensional data (like gene expression) to just a few dimensions, making it easier to visualize patterns and clusters among samples.

**Axes:**  
- **PC1 (Principal Component 1):** The direction (linear combination of genes) that captures the greatest variance in the data.
- **PC2 (Principal Component 2):** The direction orthogonal to PC1 that captures the next greatest variance.

**Points:**  
- **Each point:** Represents a single sample from your dataset.
- **Color:**  
    - **Blue:** Healthy samples  
    - **Red:** Disease samples  
    - **Gray:** Unknown or unclassified samples

**What does this plot show?**  
- Samples that are close together have similar overall gene expression profiles.
- If samples from the same group (e.g., all Healthy or all Disease) cluster together, it suggests strong biological differences between groups.
- Outliers (points far from others) may indicate technical artifacts, mislabeled samples, or unique biology.

**Assumptions:**  
- PCA is performed on the top N most significant genes (by p-value), as selected by the slider.
- Data is normalized so that differences reflect biology, not technical variation.

**How to use:**  
- Look for clear separation between groups (Healthy vs. Disease).
- Investigate outliers or unexpected clustering for possible data quality or biological insights.
            """)
            num_pca_genes = st.slider("Number of top genes to use for PCA", min_value=10, max_value=100, value=30, step=1)
            top_pca_indices = results.sort_values("pvalue").head(num_pca_genes).index
            expr_df = expression_data.loc[top_pca_indices].drop(columns=["Symbol"], errors="ignore")
            plot_pca(expr_df, sample_labels)

        # 9. MA Plot
        with tabs[8]:
            st.subheader("MA Plot (logFC vs. Mean Expression)")
            st.markdown("""
**Summary:**  
The MA plot visualizes the relationship between the average expression of each gene and its log fold change between groups.

**Axes:**  
- **X-axis (Mean Expression):** Average expression of each gene across all samples.
- **Y-axis (logFC):** Log2 fold change between groups.
-**Positive values**: Higher in Group B (e.g., Disease)
-**Negative values**: Higher in Group A (e.g., Healthy)

**Points:**  
- **True**: The gene is statistically significant (p-value < 0.05).
    These points are colored red
- **False**: The gene is not statistically significant (p-value ≥ 0.05).
    These points are colored gray.
- **How to interpret:**
    The p-value comes from a statistical test (t-test) comparing gene expression between your two groups.
    It represents the probability that the observed difference in expression happened by chance.
    Lower p-values (e.g., < 0.05) mean the difference is unlikely due to random chance and is considered statistically significant.

**Assumptions:**  
- Genes with low mean expression may have more variable fold changes.
- The t-test is used for significance.

**How to use:**  
Look for trends or biases (e.g., are highly expressed genes more likely to be differentially expressed?). Outliers may be of biological interest.
            """)
            mean_group_a = expression_data[group_a].mean(axis=1)
            mean_group_b = expression_data[group_b].mean(axis=1)
            mean_expr = (mean_group_a + mean_group_b) / 2

            ma_df = results.copy()
            ma_df['mean_expr'] = mean_expr

            fig, ax = plt.subplots()
            sns.scatterplot(
                data=ma_df,
                x='mean_expr',
                y='logFC',
                hue=ma_df['pvalue'] < 0.05,
                palette={True: 'red', False: 'gray'},
                ax=ax
            )
            ax.axhline(0, color='black', linestyle='--')
            ax.set_xlabel('Mean Expression')
            ax.set_ylabel('log2 Fold Change')
            ax.set_title('MA Plot')
            st.pyplot(fig)
    else:
        st.warning("Please select at least one sample for each group.")
