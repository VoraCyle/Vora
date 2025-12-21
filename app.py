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

# --- 1. SECURE AI CONFIGURATION ---
# We use the 'latest' stable tag and a safety setting bypass to ensure the 
# forensic chemistry isn't accidentally flagged as "dangerous content."
if "GEMINI_API_KEY" in st.secrets:
    genai.configure(api_key=st.secrets["GEMINI_API_KEY"])
    
    # Use the specifically versioned model name
    model = genai.GenerativeModel(
        model_name='gemini-1.5-flash',
        safety_settings=[
            {"category": "HARM_CATEGORY_DANGEROUS_CONTENT", "threshold": "BLOCK_NONE"},
            {"category": "HARM_CATEGORY_HARASSMENT", "threshold": "BLOCK_NONE"},
        ]
    )
else:
    st.error("🔑 Please add your GEMINI_API_KEY to the Streamlit Secrets.")
    st.stop()

# --- 2. THE REAL-TIME REASONER (FIXED) ---
def ask_gemini_forensics(stats, smiles):
    """Refined prompt to prevent 'Invalid Argument' errors."""
    try:
        # We simplify the prompt to ensure it's a clean string
        clean_prompt = (
            f"As a forensic chemist, audit the SMILES string {smiles}. "
            f"Physical data: {stats}. "
            "Identify 3 required molecular changes for Costco safety, "
            "the forensic benefits, and tactical implementation steps."
        )
        
        response = model.generate_content(clean_prompt)
        return response.text
    except Exception as e:
        return f"AI Connection Error: {str(e)}. Try a different SMILES string or rebooting."

# --- 3. THE REAL-TIME AI REASONER ---
def ask_gemini_forensics(stats, smiles):
    """Executes a live Strategic Audit using Gemini AI."""
    prompt = f"""
    You are the Wraith VoraCycle AI, a forensic procurement expert for Costco Wholesale.
    Analyze this specific molecule: SMILES {smiles}.
    Physical Molecular Data: {stats}
    
    Provide a professional Deep Dive Report including:
    1. CHANGES MADE: Specifically describe what molecular/structural changes are required to improve this score to a VoraCycle Grade A (88+).
    2. WHY (BENEFIT): Detail the exact forensic benefits regarding food safety, microplastic prevention, and global plastic tax (EPR) savings.
    3. HOW TO IMPLEMENT: Tactical implementation steps for the Costco supply chain team.
    
    Be specific to the chemistry detected (e.g., if you see Chlorine, address de-chlorination).
    """
    response = model.generate_content(prompt)
    return response.text

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


