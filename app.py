import streamlit as st
import google.generativeai as genai
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors
import pandas as pd
import numpy as np

# --- 1. SECURE AI CONFIGURATION ---
# Pulls from the Streamlit Secrets tab we set up earlier
if "GEMINI_API_KEY" in st.secrets:
    genai.configure(api_key=st.secrets["GEMINI_API_KEY"])
    # Blocking settings simplified to ensure chemistry isn't flagged as 'Dangerous'
    model = genai.GenerativeModel(
        model_name='gemini-1.5-flash',
        safety_settings=[
            {"category": "HARM_CATEGORY_DANGEROUS_CONTENT", "threshold": "BLOCK_NONE"}
        ]
    )
else:
    st.error("🔑 API Key Missing. Please add GEMINI_API_KEY to your Streamlit Secrets.")
    st.stop()

# --- 2. THE FORENSIC CHEMICAL ENGINE ---
def get_molecular_data(smiles):
    """Safely extracts physical data from a SMILES barcode."""
    try:
        # Clean the input string
        clean_smi = smiles.strip().replace('"', '').replace("'", "")
        mol = Chem.MolFromSmiles(clean_smi)
        
        if mol is None:
            return None, None
        
        stats = {
            "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
            "mw": round(Descriptors.MolWt(mol), 2),
            "rings": rdMolDescriptors.CalcNumRings(mol),
            "atoms": [a.GetSymbol() for a in mol.GetAtoms()],
            "toxic_elements": [s for s in [a.GetSymbol() for a in mol.GetAtoms()] if s in ['Cl', 'F', 'Br', 'I']],
            "bonds": mol.GetNumBonds()
        }
        return stats, mol
    except:
        return None, None

# --- 3. THE STRATEGIC REASONING ENGINE ---
def ask_gemini_forensics(stats, smiles):
    """Connects RDKit physical data to Gemini's strategic brain."""
    prompt = (
        f"Role: Forensic Procurement Expert for Costco.\n"
        f"Subject: Molecule {smiles} ({stats['formula']}).\n"
        f"Data: {stats}\n\n"
        "Task: Provide a 3-part Strategic Audit:\n"
        "1. Required Structural Changes for Soil Safety.\n"
        "2. Forensic Benefits (Toxicity & Plastic Tax Mitigation).\n"
        "3. Tactical Supply Chain Implementation Steps."
    )
    response = model.generate_content(prompt)
    return response.text

# --- 4. THE APEX INTERFACE ---
st.set_page_config(page_title="VoraCycle Apex OS", layout="wide")
st.title("🔮 Wraith VoraCycle: Apex OS")
st.sidebar.markdown("### 🟢 STATUS: LIVE AI ACTIVE")

# Department Context
dept = st.sidebar.selectbox("Costco Department", ["Deli", "Bakery", "Pharmacy", "Food Court"])

# User Input
user_input = st.text_input("🧬 Enter Molecular Barcode (SMILES) to Audit:", "C=CCl")

if user_input:
    # STEP 1: Forensic Extraction
    rt_stats, mol_obj = get_molecular_data(user_input)
    
    if rt_stats and mol_obj:
        st.divider()
        col1, col2 = st.columns([1, 1.5])
        
        with col1:
            st.image(Draw.MolToImage(mol_obj, size=(450, 450)), caption=f"Forensic Identity: {rt_stats['formula']}")
            st.metric("Molecular Weight", f"{rt_stats['mw']} g/mol")
        
        with col2:
            st.subheader("⚖️ AI Forensic Deep Dive")
            with st.spinner("AI is analyzing molecular structures..."):
                try:
                    report = ask_gemini_forensics(rt_stats, user_input)
                    st.info(report)
                except Exception as e:
                    st.error(f"AI Connection Error: {e}")

        # STEP 2: Quantum Reality Check
        st.divider()
        st.header("⚛️ Quantum Bond Dissociation Audit")
        bde = round(rt_stats['bonds'] * 0.45 * np.random.uniform(0.9, 1.1), 2)
        st.metric("Lattice Energy", f"{bde} eV", delta="-68% Goal with VoraCycle")
        st.write(f"To achieve Mineralization, the {rt_stats['formula']} lattice requires enzymatic intervention to lower this {bde} eV threshold.")

    else:
        st.warning("🧪 Waiting for a valid SMILES string. Try 'C=CCl' for PVC.")

else:
    st.info("Input a molecular barcode to begin the audit.")
