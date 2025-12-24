import streamlit as st
from openai import OpenAI

# --- 1. ACCESS ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except:
    st.error("🚨 API Key Missing.")
    st.stop()

# --- 2. THE SOURCE OF TRUTH (Demo Data) ---
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
st.set_page_config(page_title="VoraCycle Industrial", layout="wide")

# --- SIDEBAR: NAVIGATION & RESET ---
st.sidebar.title("🏢 Command Center")

# THE HOME BUTTON
if st.sidebar.button("🏠 Home / Reset Dashboard"):
    st.rerun()

st.sidebar.divider()

st.sidebar.markdown("### 🏭 Industrial Intelligence")
industrial_query = st.sidebar.text_area("Input Chemical String or Raw Material Specs:", 
                                        placeholder="e.g. 70% LDPE, 20% Nylon-6, 10% Glue")
analyze_btn = st.sidebar.button("Run Forensic Analysis")

# --- MAIN PAGE ---
st.title("🛡️ VoraCycle: Strategic DNA Command Center")

# LOGIC FOR INDUSTRIAL SIDEBAR
if analyze_btn and industrial_query:
    st.divider()
    st.subheader("🔬 Live Forensic Analysis: Custom Compound")
    with st.spinner("Processing Molecular Signature..."):
        prompt = f"""
        Act as a Vora Material Scientist. Analyze this raw input: {industrial_query}.
        Identify 'Monster DNA', propose 'Vora 3-Pillar Fix', and rate success.
        """
        response = client.chat.completions.create(model="gpt-4o", messages=[{"role": "user", "content": prompt}])
        st.info(response.choices[0].message.content)
else:
    # STANDARD DROPDOWN VIEW
    dropdown_items = ["-- Select a Strategic Asset --"] + list(dna_database.keys())
    choice = st.selectbox("Select Asset for Executive Audit:", dropdown_items)

    if choice != "-- Select a Strategic Asset --":
        st.divider()
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
        st.subheader("🏁 Executive Verdict & Company Benefits")
        with st.spinner("Finalizing Business Case..."):
            report_prompt = f"Executive summary for {choice} transformation to Vora DNA."
            report = client.chat.completions.create(model="gpt-4o", messages=[{"role": "user", "content": report_prompt}])
            st.info(report.choices[0].message.content)
    else:
        # LANDING PAGE VIEW
        st.markdown("""
        ### Welcome to the VoraCycle Executive Command Center.
        Select a strategic asset from the menu above to begin the **DNA Forensic Audit**.
        
        **How it works:**
        1. **Forensic Audit:** We identify the "Monster DNA" currently causing tax liabilities.
        2. **3-Pillar Fix:** We engineer a new DNA (Body, Brain, Skeleton) to neutralize the threat.
        3. **Path-Agnostic Success:** We verify success in both Waste and Recycle streams.
        """)