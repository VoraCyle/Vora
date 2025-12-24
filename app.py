import streamlit as st
from openai import OpenAI

# --- 1. ACCESS ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except:
    st.error("🚨 API Key Missing.")
    st.stop()

# --- 2. THE SOURCE OF TRUTH (3-Pillar DNA) ---
dna_database = {
    "Multi-Layer Chicken Bags": {
        "bad": ["Nylon-6 Barrier Film", "Polyurethane Adhesives", "Carbon Black Pigment"],
        "fix": ["92% Mono-PE (The Body)", "5% Vora-C1 Catalyst (The Brain)", "3% Mineral-Anchor (The Skeleton)"]
    },
    "MAP Poultry Trays": {
        "bad": ["Rigid Polystyrene (PS)", "EVOH Oxygen Barrier", "Chemical Blowing Agents"],
        "fix": ["95% Mono-PP (The Body)", "3% Vora-C2 Catalyst (The Brain)", "2% Mineral-Anchor (The Skeleton)"]
    },
    "PVC Clamshells": {
        "bad": ["Polyvinyl Chloride (PVC)", "Phthalate Plasticizers", "Heavy Metal Stabilizers"],
        "fix": ["90% Cellulose-Fiber (The Body)", "7% Bio-Polymer Binder (The Brain)", "3% Mineral-Tracer (The Skeleton)"]
    }
}

# --- 3. UI SETUP ---
st.set_page_config(page_title="VoraCycle: Executive Command", layout="wide")
st.title("🛡️ VoraCycle: Strategic DNA Command Center")

dropdown_items = ["-- Select a Strategic Asset --"] + list(dna_database.keys())
choice = st.selectbox("Select Asset for Executive Audit:", dropdown_items)

if choice != "-- Select a Strategic Asset --":
    st.divider()
    st.subheader(f"🧬 DNA Forensic Transformation: {choice}")
    
    bad_dna = dna_database[choice]["bad"]
    vora_dna = dna_database[choice]["fix"]

    col_a, col_b = st.columns(2)
    
    with col_a:
        st.error("### 🔴 BEFORE: Status Quo DNA")
        for item in bad_dna:
            st.warning(f"❌ **{item}**")
        
    with col_b:
        st.success("### 🟢 AFTER: The 3 Pillars of Vora DNA")
        for item in vora_dna:
            st.info(f"🧬 **{item}**")

    st.divider()
    # Path Success Indicators
    p1, p2 = st.columns(2)
    p1.write("✅ **Path A (Waste):** Safe Bio-Mineralization")
    p2.write("✅ **Path B (Recycle):** 100% High-Value Circularity")

    with st.spinner("Generating Strategic Report..."):
        prompt = f"Explain how the 3-pillar Vora DNA fix for {choice} ensures safety in both landfill and recycling streams."
        report = client.chat.completions.create(model="gpt-4o", messages=[{"role": "user", "content": prompt}])
        st.markdown(report.choices[0].message.content)