import streamlit as st
from openai import OpenAI

# --- 1. ACCESS ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except:
    st.error("🚨 API Key Missing.")
    st.stop()

# --- 2. THE EXPANDED DNA TRUTH (For 100% Accuracy) ---
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
    "PFAS Wrappers": {
        "bad": ["Fluorinated Coatings (PFAS)", "Bleached Kraft Paper", "Synthetic Wax"],
        "fix": ["94% Natural Fiber (THE BODY)", "4% Aqueous Vora-Barrier (THE BRAIN)", "2% Mineral-Inert (THE SKELETON)"]
    },
    "Kirkland Bath Tissue Case-Wrap": {
        "bad": ["LDPE Low-Density Film", "High-Slip Chemical Additives", "Mixed Polymer Regrind"],
        "fix": ["97% High-Strength Vora-PE (THE BODY)", "3% Vora-C1 Catalyst (THE BRAIN)"]
    },
    "Black Meat Trays": {
        "bad": ["Carbon Black Pigment (IR Invisible)", "Polystyrene Foam", "Mixed Resin Scrap"],
        "fix": ["96% Mono-PET (THE BODY)", "3% Vora-C2 Catalyst (THE BRAIN)", "1% IR-Visible Tracer (THE SKELETON)"]
    }
}

# --- 3. THE VORA 100 REGISTRY ---
vora_100 = {
    "🥩 POULTRY & FRESH MEATS": ["Multi-Layer Chicken Bags", "MAP Poultry Trays", "Absorbent Poultry Pads", "Black Meat Trays"],
    "🧻 PAPER & HYGIENE WRAPS": ["Kirkland Bath Tissue Case-Wrap", "Paper Towel Overwrap"],
    "🚩 HIGH-RISK LIABILITIES": ["PVC Clamshells", "PFAS Wrappers", "LLDPE Stretch Wrap"]
}

# --- 4. UI SETUP ---
st.set_page_config(page_title="VoraCycle: Executive Command", layout="wide")

st.sidebar.title("🏢 Command Center")
if st.sidebar.button("🏠 Home / Reset"):
    st.rerun()

st.title("🛡️ VoraCycle: Strategic DNA Command Center")

# Asset Selection
dropdown_items = ["-- Select a Strategic Asset --"]
for category, items in vora_100.items():
    dropdown_items.extend(items)

choice = st.selectbox("Select Asset for Executive Audit:", dropdown_items)

if choice != "-- Select a Strategic Asset --":
    st.divider()
    st.subheader(f"🧬 DNA Forensic Transformation: {choice}")
    
    # 5. THE DNA LISTS (BEFORE & AFTER)
    col_a, col_b = st.columns(2)
    
    # Logic to pull from database
    item_data = dna_database.get(choice, {
        "bad": ["Multi-Layer Composite", "Non-Separable Glues", "Mixed Polymers"],
        "fix": ["92% Vora-Base (BODY)", "5% Vora-Catalyst (BRAIN)", "3% Mineral-Anchor (SKELETON)"]
    })

    with col_a:
        st.error("### 🔴 BEFORE: Status Quo DNA")
        st.markdown("**Current Toxic Ingredients:**")
        for bad_item in item_data["bad"]:
            st.warning(f"❌ **{bad_item}**")
        st.markdown("---")
        st.write("🛑 **Path A (Waste):** 500+ Year Persistence")
        st.write("🛑 **Path B (Recycle):** Contaminates Stream")
        
    with col_b:
        st.success("### 🟢 AFTER: The 3 Pillars of Vora DNA")
        st.markdown("**New Engineering Manifest:**")
        for fix_item in item_data["fix"]:
            st.info(f"🧬 **{fix_item}**")
        st.markdown("---")
        st.write("✅ **Path A (Waste):** Safe Bio-Mineralization")
        st.write("✅ **Path B (Recycle):** 100% High-Value Circularity")

    st.divider()

    # 6. EXECUTIVE REPORT
    with st.spinner("Generating Strategic Business Case..."):
        prompt = f"Provide a CSO-level report for {choice}. Explain how the 3-pillar Vora DNA fix protects the company from EPR taxes and ensures success in both Waste and Recycle paths."
        report = client.chat.completions.create(model="gpt-4o", messages=[{"role": "user", "content": prompt}])
        st.markdown(report.choices[0].message.content)

else:
    st.info("👆 Please select an asset from the Vora 100 list to begin the DNA audit.")

