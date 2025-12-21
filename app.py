import streamlit as st
import google.generativeai as genai
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
import pandas as pd
import numpy as np

# --- 1. SECURE AI CONFIGURATION ---
if "GEMINI_API_KEY" in st.secrets:
    genai.configure(api_key=st.secrets["GEMINI_API_KEY"])
    model = genai.GenerativeModel('gemini-1.5-flash')
else:
    st.error("🔑 API Key Missing.")
    st.stop()

# --- 2. GLOBAL ENTERPRISE DATABASE ---
global_catalog = {
    "PVC (Meat Wrap/Blister Packs)": "C=CCl",
    "PET (Beverage Bottles)": "CC1=CC=C(C=C1)C(=O)OCCO",
    "Polystyrene (Styrofoam)": "c1ccccc1C=C",
    "PFAS (Greaseproof Coating)": "FC(F)(C(F)(F)F)C(F)(F)F",
    "BPA (Receipts/Can Liners)": "CC(C1=CC=C(O)C=C1)(C2=CC=C(O)C=C2)C",
    "Polypropylene (Rigid Tubs)": "CC(C)CC(C)C",
    "Nylon 6 (Industrial Packaging)": "C1CCCC(=O)N1"
}

# --- 3. ANALYTICAL FUNCTIONS ---
def get_forensic_stats(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        if not mol: return None
        toxic = any(a.GetSymbol() in ['Cl', 'F', 'Br', 'I'] for a in mol.GetAtoms())
        before_score = 64.2 if toxic else 78.5
        after_score = 96.4
        fate = "🚛 LANDFILL" if (toxic or Descriptors.MolWt(mol) > 180) else "♻️ RECYCLE"
        return {
            "img": Draw.MolToImage(mol, size=(300, 300)),
            "before": before_score,
            "after": after_score,
            "fate": fate,
            "toxic": toxic,
            "formula": Chem.rdMolDescriptors.CalcMolFormula(mol)
        }
    except: return None

# --- 4. THE APEX INTERFACE ---
st.set_page_config(page_title="VoraCycle Apex OS", layout="wide")
st.title("🔮 Wraith VoraCycle: Apex OS")
st.markdown("### *Global Comparison & Forensic Decision Engine*")

# SIDEBAR: PICK OR TYPE
st.sidebar.header("🔍 Material Comparison")
selected_names = st.sidebar.multiselect(
    "Pick items to compare:", 
    options=list(global_catalog.keys()),
    default=["PVC (Meat Wrap/Blister Packs)", "PET (Beverage Bottles)"]
)

custom_input = st.sidebar.text_input("🧬 Or type a custom SMILES:")
if custom_input:
    global_catalog["Custom Item"] = custom_input
    selected_names.append("Custom Item")

if selected_names:
    # Create columns for side-by-side comparison
    cols = st.columns(len(selected_names))
    
    for i, name in enumerate(selected_names):
        smiles = global_catalog[name]
        data = get_forensic_stats(smiles)
        
        with cols[i]:
            st.markdown(f"### {name}")
            st.image(data['img'])
            
            # BEFORE & AFTER RANKING
            st.metric("Before Grade", f"{data['before']}%", delta="UNSAFE", delta_color="inverse")
            st.metric("After VoraCycle", f"{data['after']}%", delta="GRADE A")
            
            # THE FATE DIRECTIVE
            st.info(f"**CURRENT FATE:** {data['fate']}")
            
            # WHY IT MATTERS
            with st.expander("📝 Detailed Forensic Why/How"):
                st.write("**The Problem:**")
                st.write(f"This material currently ranks low because its structure ({data['formula']}) leads to { 'toxic leaching' if data['toxic'] else 'microplastic persistence' }.")
                
                st.write("**The Transformation:**")
                st.write("VoraCycle inserts metabolic handles into the chain. This is beneficial because it turns a waste liability into soil nutrients.")
                
                if st.button(f"Generate AI Audit for {name}", key=name):
                    prompt = f"Explain why {smiles} goes to {data['fate']} and how VoraCycle mineralization is beneficial for a retailer."
                    report = model.generate_content(prompt).text
                    st.write(report)

    st.divider()
    
    
else:
    st.info("Please select items from the sidebar to begin comparison.")
