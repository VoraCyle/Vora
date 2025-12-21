import streamlit as st
import google.generativeai as genai
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
import numpy as np

# --- 1. AI CONFIGURATION ---
if "GEMINI_API_KEY" in st.secrets:
    genai.configure(api_key=st.secrets["GEMINI_API_KEY"])
    model = genai.GenerativeModel('gemini-1.5-flash')
else:
    st.error("🔑 API Key Missing.")
    st.stop()

# --- 2. GLOBAL MATERIAL DATABASE ---
global_catalog = {
    "Search Global Items...": "",
    "PET (Water Bottles)": "CC1=CC=C(C=C1)C(=O)OCCO",
    "PVC (Meat Wrap)": "C=CCl",
    "Polystyrene (Food Trays)": "c1ccccc1C=C",
    "PFAS (Greaseproof Paper)": "FC(F)(C(F)(F)F)C(F)(F)F",
    "HDPE (Milk Jugs)": "CCCCCCCCCCCC",
    "BPA (Can Linings)": "CC(C1=CC=C(O)C=C1)(C2=CC=C(O)C=C2)C"
}

# --- 3. THE DUAL-PATH ANALYTICS ---
def get_dual_path_stats(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        if not mol: return None
        
        mw = Descriptors.MolWt(mol)
        has_toxins = any(a.GetSymbol() in ['Cl', 'F', 'Br', 'I'] for a in mol.GetAtoms())
        
        # Path A: Recycling Logic
        recycle_rating = 85 if (not has_toxins and mw < 150) else 30
        if "PET" in smiles or "CC1=CC" in smiles: recycle_rating = 92 # High demand for rPET
        
        # Path B: Landfill Logic
        landfill_rating = 15 if has_toxins else 45
        reaction = "Leaches Microplastics & Toxins" if has_toxins else "Persistent Mummification (Centuries)"
        
        return {
            "img": Draw.MolToImage(mol, size=(300, 300)),
            "recycle_rating": recycle_rating,
            "landfill_rating": landfill_rating,
            "reaction": reaction,
            "tox": has_toxins,
            "formula": Chem.rdMolDescriptors.CalcMolFormula(mol)
        }
    except: return None

# --- 4. THE APEX INTERFACE ---
st.set_page_config(page_title="VoraCycle Apex OS", layout="wide")
st.title("🔮 Wraith VoraCycle: Apex OS")
st.markdown("#### *The Start Line: Predicting the Endgame for Every Material*")

# SIDEBAR: PICK OR TYPE
st.sidebar.header("🔍 Material Inventory")
search_choice = st.sidebar.selectbox("Select Business Item:", list(global_catalog.keys()))
user_input = st.sidebar.text_input("🧬 Or Enter SMILES Barcode:", value=global_catalog[search_choice])

if user_input:
    data = get_dual_path_stats(user_input)
    if data:
        st.divider()
        col_img, col_metrics = st.columns([1, 2])
        
        with col_img:
            st.image(data['img'], caption=f"Molecular Identity: {data['formula']}")
        
        with col_metrics:
            c1, c2 = st.columns(2)
            c1.metric("♻️ Recycling Capability", f"{data['recycle_rating']}%", delta="High Demand" if data['recycle_rating'] > 70 else "Low Value")
            c2.metric("🚛 Landfill Safety Rating", f"{data['landfill_rating']}%", delta="Toxic Risk" if data['tox'] else "Stable", delta_color="inverse")
            
            st.warning(f"**Landfill Reaction:** {data['reaction']}")

        st.divider()

        # --- PATHWAY COMPARISON ---
        path_a, path_b = st.tabs(["🚀 PATH A: RECYCLING OUTCOME", "🌋 PATH B: LANDFILL OUTCOME"])

        with path_a:
            st.subheader("The Circularity Result")
            st.write(f"**Outcome:** {'High-Value rPCR' if data['recycle_rating'] > 70 else 'Downcycled / Rejected'}")
            st.write("""
            **Process:** Material is cleaned and shredded. If the rating is high, it maintains structural integrity 
            for up to 3-5 heat cycles before the polymer chains shorten and the plastic becomes 'brittle.'
            """)
            

        with path_b:
            st.subheader("The Reality of Disposal")
            st.write(f"**Current Reaction:** {data['reaction']}")
            st.error("🚨 CRITICAL: This item requires 'Molecular Surgery' to prevent soil contamination.")
            
            # THE "ECO-FRIENDLY BUT SAFE" SECTION
            st.info("### 🛠️ VoraCycle Re-Engineering (The Upgrade)")
            st.write("""
            To make this item **Soil-Safe** without losing **Structural Integrity** (No 'weakening'):
            
            1. **Metabolic Handles:** We introduce enzymatic 'trigger points' that remain dormant while holding food. 
            2. **Trigger-Activated:** These points only activate when exposed to soil microbes and high moisture (Landfill/Compost conditions).
            3. **Structural Shielding:** We use cross-linkers that ensure the package doesn't get 'soggy' or 'weak' during transport, maintaining the same shelf-life as traditional plastic.
            """)
            

        # FINAL DETAILED DESCRIPTION
        st.subheader("⚖️ Final Forensic Executive Summary")
        with st.spinner("AI Generating Endgame Report..."):
            prompt = (f"Explain why {user_input} ({data['formula']}) has a {data['recycle_rating']}% recycling score. "
                      f"Detail how it reacts in a landfill (currently {data['reaction']}). "
                      f"Suggest specific chemical changes to make it mineralize in soil while remaining strong enough for food safety.")
            report = model.generate_content(prompt).text
            st.write(report)
            
        
    else:
        st.warning("Please provide a valid chemical barcode.")
