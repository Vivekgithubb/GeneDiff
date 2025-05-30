import streamlit as st

def extract_sample_groups(gse):
    group_labels = {}
    for gsm_id, gsm in gse.gsms.items():
        # Try title metadata
        title = gsm.metadata.get('title', ['Unknown'])[0].lower()
        
        if 'healthy' in title or 'normal' in title:
            group_labels[gsm_id] = "Healthy"
        elif 'disease' in title or 'cancer' in title or 'tumor' in title:
            group_labels[gsm_id] = "Disease"
        else:
            # Try characteristics field if title not useful
            characteristics = gsm.metadata.get('characteristics_ch1', [])
            # This is a list, try to find group name
            label_found = False
            for c in characteristics:
                c_lower = c.lower()
                if 'healthy' in c_lower or 'normal' in c_lower:
                    group_labels[gsm_id] = "Healthy"
                    label_found = True
                    break
                elif 'disease' in c_lower or 'cancer' in c_lower or 'tumor' in c_lower:
                    group_labels[gsm_id] = "Disease"
                    label_found = True
                    break
            if not label_found:
                group_labels[gsm_id] = "Unknown"

    st.write("All sample IDs in GEO dataset:", list(gse.gsms.keys()))

    return group_labels
