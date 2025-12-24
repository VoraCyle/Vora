import streamlit as st
from openai import OpenAI

# --- 1. ACCESS ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except:
    st.error("🚨 API Key Missing.")
    st.stop()

# --- 2. THE SOURCE OF TRUTH (The 3-Pillar DNA) ---
dna_database = {
    "Multi-Layer Chicken Bags": {
        "bad": ["Nylon-6 Barrier Film", "Polyurethane Adhesives", "Carbon Black Pigment"],
        "fix": ["92% Mono-PE (THE BODY)", "5% Vora-C1 Catalyst (THE BRAIN)", "3% Mineral-Anchor (THE SKELETON)"]
    },
    "MAP Poultry Trays": {
        "bad": ["Rigid Polystyrene (PS)", "EVOH Oxygen Barrier", "Chemical Blowing Agents"],
        "fix": ["95% Mono-PP (THE BODY)", "3% Vora-C2 Catalyst (THE BRAIN)", "2% Mineral-Anchor (THE SKELETON)"]
    },
    "PVC Clamshells": {
        "bad": ["Polyvinyl Chloride (PVC)", "Phthalate Plasticizers", "Heavy Metal Stabilizers"],
        "fix": ["90% Cellulose-Fiber (THE BODY)", "7% Bio-Polymer Binder (THE BRAIN)", "3% Mineral-Tracer (THE SKELETON)"]
    },
    "Kirkland Bath Tissue Case-Wrap": {
        "bad": ["LDPE Low-Density Film", "High-Slip Chemical Additives", "Mixed Polymer Regrind"],
        "fix": ["97% High-Strength Vora-PE (THE BODY)", "3% Vora-C1 Catalyst (THE BRAIN)"]
    },
    "PFAS Wrappers": {
        "bad": ["Fluorinated Coatings (PFAS)", "Bleached Kraft Paper", "Synthetic Wax"],
        "fix": ["94% Natural Fiber (THE BODY)", "4% Aqueous Vora-Barrier (THE BRAIN)", "2% Mineral-Inert (THE SKELETON)"]
    }
}

# --- 3. UI SETUP ---
st.set_page_config(page_title="VoraCycle: Executive Command", layout="wide")
st.title("🛡️ VoraCycle: Strategic DNA Command Center")

dropdown_items = ["-- Select a Strategic Asset --"] + list(dna_database.keys())
choice = st.selectbox("Select Asset for Executive Audit:", dropdown_items)

if choice != "-- Select a Strategic Asset --":
    st.divider()
    
    # 🧬 DNA COMPARISON SECTION
    st.subheader(f"🧬 DNA Forensic Transformation: {choice}")
    col_a, col_b = st.columns(2)
    
    bad_dna = dna_database[choice]["bad"]
    vora_dna = dna_database[choice]["fix"]

    with col_a:
        st.error("### 🔴 BEFORE: Status Quo DNA")
        for item in bad_dna:
            st.warning(f"❌ **{item}**")
        st.markdown("---")
        st.write("🛑 **Path A (Waste):** 500+ Year Persistence")
        st.write("🛑 **Path B (Recycle):** Contaminates Stream")
        
    with col_b:
        st.success("### 🟢 AFTER: The 3 Pillars of Vora DNA")
        for item in vora_dna:
            st.info(f"🧬 **{item}**")
        st.markdown("---")
        st.write("✅ **Path A (Waste):** Safe Bio-Mineralization")
        st.write("✅ **Path B (Recycle):** 100% High-Value Circularity")

    st.divider()

    # 📊 THE BUSINESS CASE SUMMARY
    st.subheader("🏁 Executive Verdict & Company Benefits")
    
    with st.spinner("Finalizing Business Case..."):
        prompt = f"""
        Provide a 3-point business summary for the product '{choice}' using Vora DNA. 
        Focus on:
        1. FINANCIAL: How this stops EPR taxes and future plastic penalties.
        2. OPERATIONAL: How it simplifies the supply chain with mono-materials.
        3. BRAND: How it protects the stock price by removing 'PFAS' and 'Toxin' liabilities.
        Keep it sharp and executive-toned.
        """
        report = client.chat.completions.create(model="gpt-4o", messages=[{"role": "user", "content": prompt}])
        st.info(report.choices[0].message.content)

    st.caption("Strategic ROI: This transition creates an immediate EPR tax exemption and future-proofs the enterprise against 2026 bans.")

else:
    st.info("👆 Please select an asset to view the Business Case.")