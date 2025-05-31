# 🧬 GeneDiff: Visual Gene Expression Analyzer

**GeneDiff** is a web-based application that allows users to visually compare gene expression between healthy and diseased samples using publicly available microarray datasets (GEO). It offers interactive plots and statistical filtering to aid researchers, students, and clinicians in understanding gene regulation patterns.

🔗 **Live App**: [https://genediff.streamlit.app](https://genediff.streamlit.app)

---

## 📌 Features

- 🧪 **Load Public Datasets** (e.g., GSE7308, GSE10950)
- 🔥 **Heatmap** of gene expression levels across conditions
- 🌋 **Volcano Plot** to identify significantly up/downregulated genes
- 📦 **Box Plot** to compare expression distribution across samples
- ⚖️ **MA Plot** for fold-change vs expression bias
- 🧭 **PCA Plot** for visualizing sample clustering
- 🔍 **Gene Search & Filter** by p-value, fold change, etc.

---

## 🛠️ Tech Stack

- **Frontend/App Framework**: Streamlit
- **Visualization**: Matplotlib, Seaborn, Plotly
- **Data Handling**: Pandas, NumPy
- **Datasets**: Publicly available GEO datasets (preprocessed)

---

## 🚀 Getting Started

### Prerequisites

Ensure you have Python 3.8+ installed.

### Installation

```bash
git clone https://github.com/Vivekgithubb/GeneDiff.git
cd GeneDiff
pip install -r requirements.txt
streamlit run app.py
```

---

## 📂 Project Structure

```
GeneDiff/
├── app.py               # Main application script
├── data/                # Dataset files
├── utils/               # Helper functions
├── screenshots/         # Visual examples (optional)
├── requirements.txt     # Python dependencies
└── README.md
```

---

## 🙋‍♂️ Team

**Team Byte-Me**
- Akshant A.A.
- D Vivek Pai
- Thrisha Santhosh
- Chaitanya Kamath

---

## 📈 Future Plans

- Upload custom user datasets
- Add clustering & enrichment analysis
- Download filtered results
- Improve performance on large datasets
---

## 📬 Contact

For questions or suggestions, please open an issue on the [GitHub repository](https://github.com/Vivekgithubb/GeneDiff/issues).
