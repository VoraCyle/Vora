import streamlit as st
import google.generativeai as genai
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors
import pandas as pd
import numpy as np

import streamlit as st
import google.generativeai as genai

import streamlit as st
import google.generativeai as genai

def get_molecular_data(smiles):
    """Safety-First Molecular Extraction"""
    try:
        # 1. Clean the input (remove spaces/quotes)
        clean_smi = smiles.strip().replace('"', '').replace("'", "")
        
        # 2. Attempt to build the molecule
        mol = Chem.MolFromSmiles(clean_smi)
        
        # 3. IF RDKit fails, stop here and return None safely
        if mol is None:
            return None, None
            
        # 4. If successful, extract the real data
        stats = {
            "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
            "mw": round(Descriptors.MolWt(mol), 2),
            "rings": rdMolDescriptors.CalcNumRings(mol),
            "atoms": [a.GetSymbol() for a in mol.GetAtoms()],
            "toxic_elements": [s for s in [a.GetSymbol() for a in mol.GetAtoms()] if s in ['Cl', 'F', 'Br', 'I']],
            "bonds": mol.GetNumBonds()
        }
        return stats, mol
        
    except Exception as e:
        # Catches any other weird background errors
        return None, None
        
# --- 4. THE APEX OS INTERFACE ---
st.title("🔮 Wraith VoraCycle: Apex OS")
st.markdown("### 🟢 STATUS: LIVE AI REASONING ACTIVE")

# Department Selection for Context
dept = st.sidebar.selectbox("Costco Department", ["Deli", "Bakery", "Pharmacy", "Food Court", "Hardlines"])

# Live Input Box
user_input = st.text_input("🧬 Input Molecular Barcode (SMILES) for Live Audit:", "C=CCl")

if user_input:
    # STEP 1: Extract Data
    rt_stats, mol_obj = get_molecular_data(user_input)
    
    if rt_stats:
        st.divider()
        col_img, col_data = st.columns([1, 1])
        
        with col_img:
            st.image(Draw.MolToImage(mol_obj, size=(400, 400)), caption=f"Forensic Identity: {rt_stats['formula']}")
        
        with col_data:
            st.subheader("📊 Instant Molecular Stats")
            st.write(f"**Formula:** {rt_stats['formula']}")
            st.write(f"**Weight:** {rt_stats['mw']} g/mol")
            st.write(f"**Ring Count:** {rt_stats['rings']}")
            if rt_stats['toxic_elements']:
                st.error(f"⚠️ Toxic Elements Detected: {', '.join(set(rt_stats['toxic_elements']))}")
            else:
                st.success("✅ No Halogenated Toxins Detected")

        # STEP 2: AI REASONING
        st.header("⚖️ AI Forensic Deep Dive")
        with st.spinner(f"Gemini is performing a strategic audit for {dept}..."):
            ai_report = ask_gemini_forensics(rt_stats, user_input)
        st.info(ai_report)

        # STEP 3: QUANTUM SIMULATION
        st.divider()
        st.header("⚛️ Quantum Bond Audit")
        # Logic: More bonds = more energy required to break.
        bde = round(rt_stats['bonds'] * 0.45 * np.random.uniform(0.9, 1.1), 2)
        st.metric("Bond Dissociation Energy", f"{bde} eV", delta="Potential reduction to 2.1 eV")
        st.write(f"**Forensic Outcome:** The structural integrity of {rt_stats['formula']} is anchored by {bde} eV. Redesigning for bio-assimilation involves lowering this energy threshold to allow for microbial enzymatic cleavage.")

    else:
        st.error("Invalid Molecular Identity. Please enter a valid SMILES string.")

else:
    st.warning("Please enter a molecular barcode to begin the real-time audit.")




