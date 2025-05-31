import streamlit as st 

def extract_sample_groups(gse):
    sample_labels = {}
    group_dict = {"Healthy": [], "Disease": [], "Unknown": []}

    for gsm_id, gsm in gse.gsms.items():
        title = gsm.metadata.get("title", [""])[0].lower()
        characteristics = " ".join(gsm.metadata.get("characteristics_ch1", [])).lower()
        full_text = f"{title} {characteristics}"

        label = "Unknown"
        disease_type = None

        # Try to infer label
        if any(x in full_text for x in ["healthy", "normal", "control"]):
            label = "Healthy"
        elif any(x in full_text for x in ["disease", "cancer", "tumor", "carcinoma", "malignant", "sarcoma", "adenoma"]):
            label = "Disease"

            # Extract possible disease type dynamically from title or characteristics
            tokens = title.split() + characteristics.split()
            disease_keywords = [word.capitalize() for word in tokens if word not in ("disease", "tumor", "cancer", "sample") and len(word) > 3]

            if disease_keywords:
                disease_type = disease_keywords[0]  # pick the first valid word as hint

        if label == "Disease" and disease_type:
            display_label = f"Disease ({disease_type})"
        elif label != "Unknown":
            display_label = label
        else:
            display_label = "Unknown"

        # if label == "Healthy" and disease_type:
        #     display_label = f"Healthy ({disease_type})"
        # elif label != "Unknown":
        #     display_label = label
        # else:
        #     display_label = "Unknown"

        sample_labels[gsm_id] = f"{gsm_id} ({display_label})"
        group_dict[label].append(gsm_id)

    return sample_labels, group_dict
