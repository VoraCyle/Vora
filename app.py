

import streamlit as st
from openai import OpenAI

# --- 1. ACCESS ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except:
    st.error("🚨 API Key Missing.")
    st.stop()

# --- 2. THE DNA DATABASE ---
dna_database = {
    "Multi-Layer Chicken Bags": {
        "base_material": "Mono-PE (Polyethylene)",
        "catalyst_type": "Enzymatic/Pro-degradant (Vora-C1)",
        "bad": ["Nylon-6 Barrier", "Polyurethane Adhesives", "Carbon Black"],
        "fix": ["92% Mono-PE (BODY)", "5% Vora-C1 Catalyst (BRAIN)", "3% Mineral-Anchor (SKELETON)"]
    },
    "MAP Poultry Trays": {
        "base_material": "Mono-PP (Polypropylene)",
        "catalyst_type": "Bio-Mineral Catalyst (Vora-C2)",
        "bad": ["Rigid Polystyrene (PS)", "EVOH Oxygen Barrier"],
        "fix": ["95% Mono-PP (BODY)", "3% Vora-C2 Catalyst (BRAIN)", "2% Mineral-Anchor (SKELETON)"]
    },
    "PVC Clamshells": {
        "base_material": "Cellulose / Fiber Base",
        "catalyst_type": "Aqueous Bio-Binder",
        "bad": ["PVC", "Phthalates", "Lead Stabilizers"],
        "fix": ["90% Fiber (BODY)", "7% Bio-Binder (BRAIN)", "3% Mineral-Tracer (SKELETON)"]
    }
}

# --- 3. UI SETUP ---
st.set_page_config(page_title="VoraCycle Industrial", layout="wide")
st.title("🛡️ VoraCycle: Strategic DNA Command Center")

# --- SIDEBAR & HOME ---
st.sidebar.title("🏢 Command Center")
if st.sidebar.button("🏠 Home / Reset Dashboard"):
    st.rerun()

# --- MAIN SELECTION ---
dropdown_items = ["-- Select a Strategic Asset --"] + list(dna_database.keys())
choice = st.selectbox("Select Asset for Executive Audit:", dropdown_items)

if choice != "-- Select a Strategic Asset --":
    item_data = dna_database[choice]
    
    # 🧬 DNA COMPARISON
    st.divider()
    col_a, col_b = st.columns(2)
    with col_a:
        st.error("### 🔴 BEFORE: Status Quo DNA")
        for item in item_data["bad"]: st.warning(f"❌ **{item}**")
    with col_b:
        st.success("### 🟢 AFTER: Vora DNA Blueprint")
        for item in item_data["fix"]: st.info(f"🧬 **{item}**")

    # 🏭 THE SUPPLY CHAIN ENGINE (The "Real Life" Part)
    st.divider()
    st.subheader("🔗 Global Supplier Match (Tier 1 & 2)")
    
    if st.button("🔎 FIND CERTIFIED VORA-STYLE SUPPLIERS"):
        with st.spinner("Scanning Global Material Markets..."):
            # This prompt forces the AI to find REAL companies for the specific fix
            supplier_prompt = f"""
            Identify real-world global suppliers for the following transition:
            From: {item_data['bad']} 
            To: {item_data['fix']}
            
            List real companies for:
            1. The Base Mono-Material (e.g., Dow, LyondellBasell, Sabic).
            2. The 'Vora-style' Catalyst/Masterbatch (e.g., companies like Wells Plastics (Reverte), Symphony Environmental (d2w), or Evonik).
            3. Specify which ASTM/ISO standards these suppliers help meet.
            """
            response = client.chat.completions.create(model="gpt-4o", messages=[{"role": "user", "content": supplier_prompt}])
            st.write(response.choices[0].message.content)

    # 📝 WORK ORDER
    with st.expander("📝 VIEW MANUFACTURER WORK ORDER"):
        st.markdown(f"""
        **TO:** Costco Packaging Vendor
        **ACTION:** Transition {choice} to Path-Agnostic DNA.
        **RECIPE:** Mix {item_data['base_material']} with {item_data['catalyst_type']} additive at 5% load rate.
        """)

else:
    st.info("👆 Select a product to identify the DNA and the Suppliers.")

st.divider()
st.caption("Conceptual Guidance • Real production requires lab validation and certified supplier onboarding.")