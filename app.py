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
        "bad": ["Nylon-6 Barrier Film", "Polyurethane Adhesives", "Carbon Black Pigment"],
        "fix": ["92% Mono-PE (THE BODY)", "5% Vora-C1 Catalyst (THE BRAIN)", "3% Mineral-Anchor (THE SKELETON)"]
    },
    "MAP Poultry Trays": {
        "base_material": "Mono-PP (Polypropylene)",
        "catalyst_type": "Bio-Mineral Catalyst (Vora-C2)",
        "bad": ["Rigid Polystyrene (PS)", "EVOH Oxygen Barrier", "Chemical Blowing Agents"],
        "fix": ["95% Mono-PP (THE BODY)", "3% Vora-C2 Catalyst (THE BRAIN)", "2% Mineral-Anchor (THE SKELETON)"]
    },
    "PVC Clamshells": {
        "base_material": "Cellulose / Fiber Base",
        "catalyst_type": "Aqueous Bio-Binder",
        "bad": ["Polyvinyl Chloride (PVC)", "Phthalate Plasticizers", "Heavy Metal Stabilizers"],
        "fix": ["90% Cellulose-Fiber (THE BODY)", "7% Bio-Polymer Binder (THE BRAIN)", "3% Mineral-Tracer (THE SKELETON)"]
    }
}

# --- 3. UI SETUP ---
st.set_page_config(page_title="VoraCycle Industrial", layout="wide")
st.title("🛡️ VoraCycle: Strategic DNA Command Center")

st.sidebar.title("🏢 Command Center")
if st.sidebar.button("🏠 Home / Reset Dashboard"):
    st.rerun()

# --- MAIN SELECTION ---
dropdown_items = ["-- Select a Strategic Asset --"] + list(dna_database.keys())
choice = st.selectbox("Select Asset for Executive Audit:", dropdown_items)

if choice != "-- Select a Strategic Asset --":
    item_data = dna_database[choice]
    
    st.divider()
    col_a, col_b = st.columns(2)
    with col_a:
        st.error("### 🔴 BEFORE: Status Quo DNA")
        for item in item_data["bad"]: st.warning(f"❌ **{item}**")
    with col_b:
        st.success("### 🟢 AFTER: Vora DNA Blueprint")
        for item in item_data["fix"]: st.info(f"🧬 **{item}**")

    # --- 4. PATH-AGNOSTIC SCIENCE SUMMARY ---
    st.divider()
    st.subheader("🛡️ The Path-Agnostic Guarantee")
    col_s1, col_s2 = st.columns(2)
    with col_s1:
        st.markdown("""
        ### 🔄 Path A: Waste (Landfill)
        **The Science:** Enzymatic 'Molecular Scissors'
        - **Speed:** 90% breakdown in <5 years (vs. 450+ years).
        - **Outcome:** Full bio-mineralization; zero microplastics.
        - **Standard:** Meets **ASTM D5511**.
        """)
    with col_s2:
        st.markdown("""
        ### ♻️ Path B: Recycle (Circular)
        **The Science:** Mono-Material Purity
        - **Outcome:** High-value pellet for food-grade reuse.
        - **Benefit:** Zero EPR Plastic Tax liability.
        - **Standard:** Meets **APR Critical Guidance**.
        """)

    # --- 5. SUPPLIERS & WORK ORDER ---
    st.divider()
    st.subheader("🔗 Global Supplier Match")
    if st.button("🔎 FIND CERTIFIED VORA-STYLE SUPPLIERS"):
        with st.spinner("Scanning Global Markets..."):
            prompt = f"Identify real-world suppliers for {item_data['base_material']} and {item_data['catalyst_type']} masterbatch (e.g. Dow, Wells Plastics)."
            response = client.chat.completions.create(model="gpt-4o", messages=[{"role": "user", "content": prompt}])
            st.write(response.choices[0].message.content)

    with st.expander("📝 VIEW MANUFACTURER WORK ORDER"):
        st.markdown(f"**TO:** Packaging Vendor\n**ACTION:** Transition to {item_data['base_material']} + 5% {item_data['catalyst_type']}.")

else:
    st.info("👆 Select a product to begin the DNA audit.")