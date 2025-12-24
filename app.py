import streamlit as st
from openai import OpenAI

# --- 1. ACCESS & CONNECTION ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except Exception as e:
    st.error("🚨 API Key Missing: Please add 'OPENAI_API_KEY' to your Streamlit Secrets.")
    st.stop()

# --- 2. HARD-CODED DNA TRUTH (For 100% Accuracy on Key Items) ---
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
    }
}

# --- 3. UI SETUP ---
st.set_page_config(page_title="VoraCycle: Executive Command", layout="wide")

# SIDEBAR NAVIGATION
st.sidebar.title("🏢 Command Center")
if st.sidebar.button("🏠 Home / Reset Dashboard"):
    st.rerun()

st.sidebar.divider()
st.sidebar.info("VoraCycle v2025.12\nConceptual Industrial Guidance")

# MAIN HEADER
st.title("🛡️ VoraCycle: Strategic DNA Command Center")
st.warning("""
⚠️ **CRITICAL NOTICE** – This engine provides directional guidance. Real manufacturing requires supplier collaboration, lab testing, and regulatory validation.
""")

# --- 4. MAIN INPUT ENGINE ---
st.subheader("🧬 Forensic Product Audit")
current_packaging = st.text_area(
    "Describe the packaging product or select a strategic asset below:",
    height=100,
    placeholder="e.g. Multi-layer PET/PE snack pouch, or type a Costco SKU..."
)

# Optional Dropdown for quick access to hard-coded wins
quick_select = st.selectbox("Quick Select Strategic Asset:", ["-- Choose --"] + list(dna_database.keys()))
final_input = quick_select if quick_select != "-- Choose --" else current_packaging

analyze_btn = st.button("→ Run Vora DNA Analysis", type="primary")

if analyze_btn and final_input:
    st.divider()
    
    # --- VISUAL DNA COMPARISON ---
    st.subheader(f"Analysis Results: {final_input}")
    
    col_a, col_b = st.columns(2)
    
    # Get specific DNA if it exists, otherwise use AI-generated fallback
    item_data = dna_database.get(final_input, {
        "bad": ["Multi-Layer Composite Films", "Non-Separable Adhesives", "Mixed Polymers"],
        "fix": ["92% Mono-Material Base (THE BODY)", "5% Vora-C1 Catalyst (THE BRAIN)", "3% Mineral-Anchor (THE SKELETON)"]
    })

    with col_a:
        st.error("### 🔴 BEFORE: Status Quo DNA")
        for item in item_data["bad"]: st.warning(f"❌ **{item}**")
        st.write("🛑 **Waste:** Persistent Microplastics")
        st.write("🛑 **Recycle:** Rejected / Low-Value")
        
    with col_b:
        st.success("### 🟢 AFTER: The 3 Pillars of Vora DNA")
        for item in item_data["fix"]: st.info(f"🧬 **{item}**")
        st.write("✅ **Path A (Waste):** Safe Bio-Mineralization")
        st.write("✅ **Path B (Recycle):** 100% High-Value Circularity")

    st.divider()

    # --- DEEP INDUSTRIAL REPORT ---
    with st.spinner("Generating Strategic Impact Report..."):
        prompt = f"""
        Analyze packaging: {final_input}. 
        1. Describe current composition. 
        2. Identify Path-Agnostic failure points. 
        3. Recommend 1-3 established commercial partners (e.g. SEE, Amcor, Cruz Foam).
        4. Explain how Vora DNA creates EPR tax immunity.
        End with the standard disclaimer.
        """
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[{"role": "system", "content": "You are the VoraCycle CSO. Speak with financial authority and material science precision."},
                      {"role": "user", "content": prompt}]
        )
        st.markdown(response.choices[0].message.content)

else:
    st.info("👆 Enter a product description or select a strategic asset to begin.")

# --- FOOTER ---
st.divider()
st.caption("Strategic ROI: This transition creates an immediate EPR tax exemption and future-proofs the enterprise against 2026 bans.") 
