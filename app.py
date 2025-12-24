import streamlit as st
from openai import OpenAI
import json

# --- 1. ACCESS & CONNECTION ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except Exception as e:
    st.error("🚨 API Key Missing.")
    st.stop()

# --- 2. LOAD THE INGREDIENT DATABASE ---
try:
    with open('material_library.json', 'r') as f:
        monster_db = json.load(f)
except:
    monster_db = {}

# --- 3. THE VORA 100 REGISTRY ---
vora_100 = {
    "🥩 POULTRY & FRESH MEATS": ["MAP Poultry Trays", "Absorbent Poultry Pads", "Vacuum Wraps", "Black Meat Trays"],
    "🧻 PAPER & HYGIENE WRAPS": ["Kirkland Bath Tissue Case-Wrap", "Paper Towel Overwrap"],
    "🚩 HIGH-RISK LIABILITIES": ["PVC Clamshells", "PFAS Wrappers", "Lithium Battery Packs", "LLDPE Stretch Wrap"]
}

# --- 4. USER INTERFACE ---
st.set_page_config(page_title="VoraCycle: Executive Command", layout="wide")
st.title("🛡️ VoraCycle: Strategic DNA Command Center")

col1, col2 = st.columns(2)
with col1:
    dropdown_items = ["-- Select a Strategic Asset --"]
    for cat, items in vora_100.items(): dropdown_items.extend(items)
    choice = st.selectbox("Select Asset for Executive Audit:", dropdown_items)
with col2:
    search = st.text_input("Search Custom SKU:")

final_target = search if search else (choice if choice != "-- Select a Strategic Asset --" else None)

if final_target:
    st.divider()
    st.subheader(f"🧬 DNA Forensic Transformation: {final_target}")
    
    before_col, after_col = st.columns(2)
    
    with before_col:
        st.error("### 🔴 BEFORE: Status Quo DNA")
        st.markdown("**Core Liabilities:**")
        st.warning("⚠️ Multi-layer Composite Films")
        st.warning("⚠️ Non-Separable Poly-Glues")
        st.warning("⚠️ Carbon-Black Pigments")
        st.warning("⚠️ PFAS Moisture Barriers")
        
    with after_col:
        st.success("### 🟢 AFTER: Vora DNA Blueprint")
        st.markdown("**New DNA Ingredients:**")
        
        # Pull ingredients from DB or use fallback
        if final_target in monster_db:
            ingredients = monster_db[final_target]['vora_fix']['recipe'].split(',')
        else:
            ingredients = ["92% Mono-Polymer Base", "5% Vora-C1 Catalyst", "3% Mineral-Anchor Nutrient"]
        
        # Display as Technical Cards
        for ingredient in ingredients:
            st.info(f"🧬 **{ingredient.strip()}**")
            
        st.caption("Industrial Instructions: Standard 'Drop-In' extrusion compatible.")

    st.divider()
    
    # Executive report
    with st.spinner("Analyzing Financial Impact..."):
        prompt = f"Executive summary for transforming {final_target} to Vora DNA focusing on EPR savings."
        response = client.chat.completions.create(model="gpt-4o", messages=[{"role": "user", "content": prompt}])
        st.markdown(response.choices[0].message.content)