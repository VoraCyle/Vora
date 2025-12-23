import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from openai import OpenAI  # We only need OpenAI now!

# This connects your app to OpenAI using the secret key you'll put in Streamlit
client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])

def generate_vora_analysis(prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {"role": "system", "content": "You are the VoraCycle Arbiter, a senior forensic analyst. You provide technical, highly detailed executive reports."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.7
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Analysis Error: {str(e)}"

# --- 2. THE STRATEGIC INVENTORY --- (Unchanged)
product_inventory = {
    "Select a problematic item...": "",
    "Milk Jug (HDPE)": "CCCCCCCC",
    "Water Bottle (PET)": "CC1=CC=C(C=C1)C(=O)OCCO",
    "Aluminum Beverage Can": "[Al]",
    "Meat Wrap (PVC)": "C=CCl",
    "Deli Container (PP)": "CC(C)CC(C)C",
    "Frozen Food Bag (LDPE)": "CCCCCCCCCCCC",
    "Chip Bag (Multi-layer)": "CCCCCCCCCC.C=CC#N",
    "Styrofoam (EPS)": "C1=CC=C(C=C1)C=C",
    "Waxed Cardboard (Produce Box)": "CCCCCCCCCCCCCCCCCCCC",
}

# --- 3. THE STRATEGIC DECISION ENGINE --- (Unchanged)
def run_strategic_audit(item_name, smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        mw = Descriptors.MolWt(mol) if mol else 28.0
        toxic = any(a.GetSymbol() in ['Cl', 'F', 'Br', 'I'] for a in mol.GetAtoms()) if mol else False
        
        br = 92 if "PET" in item_name or "Aluminum" in item_name else max(10, 75 - (mw / 10))
        if toxic: br -= 30
        ar = min(99.1, br + 45) 
        
        bm = 5 if toxic else max(12, 50 - (mw / 20))
        am = 99.4 if not toxic else 94.0 
        
        if toxic or "Multi-layer" in item_name or "Styrofoam" in item_name:
            best_path = "Mineralization"
            priority = "Environmental Defense"
            reason = "Waste Trap detected. Recycling is energy-prohibitive. Mineralization ensures zero-microplastic residue."
        elif ar > 90:
            best_path = "Mechanical Recycling"
            priority = "Resource Preservation"
            reason = "High-value circularity. Resource recovery saves more energy than carbon-return."
        else:
            best_path = "Mineralization"
            priority = "Environmental Defense"
            reason = "Efficiency gap. Mineralization is the most cost-effective path for this molecular profile."

        return round(br, 1), round(ar, 1), round(bm, 1), round(am, 1), best_path, priority, reason, mw
    except:
        return None

# --- 4. THE INTERFACE ---
st.set_page_config(page_title="VoraCycle Strategic Arbiter", layout="wide")
st.title("🔮 Wraith VoraCycle: Strategic Arbiter")
st.markdown("### *Forensic Audit: Determining the Apex Move for Environment & Enterprise*")

col_search1, col_search2 = st.columns(2)
with col_search1:
    dropdown_query = st.selectbox("🧬 Preselected Problematic Items:", list(product_inventory.keys()))
with col_search2:
    manual_query = st.text_input("✍️ Manual Search (Type any item):")

query = manual_query if manual_query else dropdown_query

if query and query != "Select a problematic item...":
    smiles = product_inventory.get(query, "CCCCCCCC") 
    audit = run_strategic_audit(query, smiles)
    
    if audit:
        br, ar, bm, am, best_path, priority, reason, mw = audit
        
        st.divider()
        st.header(f"🏆 Best Strategic Path: {best_path}")
        st.info(f"**Primary Objective:** {priority} | **Logic:** {reason}")
        
        # DUAL-PATH PERFORMANCE
        st.header("📊 Comparative Forensic Performance")
        p1, p2 = st.columns(2)
        with p1:
            st.subheader("♻️ Path 1: Mechanical Recycling")
            st.metric("After VoraCycle", f"{ar}%", delta=f"Baseline: {br}%")
            st.warning(f"**Before ({br}%):** Scoring is suppressed by **Chain Scission**.")
            st.success(f"**After ({ar}%):** Optimization achieved via **Atomic Re-linkers**.")

        with p2:
            st.subheader("🌿 Path 2: Soil Mineralization")
            st.metric("After VoraCycle", f"{am}%", delta=f"Baseline: {bm}%")
            st.warning(f"**Before ({bm}%):** Scoring fails due to **Hydrophobic Locking**.")
            st.success(f"**After ({am}%):** Optimization achieved via **Metabolic Triggering**.")

        # Tabs Section (Unchanged)
        st.divider()
        st.header("⚖️ Resource Efficiency Analysis")
        t1, t2, t3 = st.tabs(["💰 Money", "⏳ Time", "🌍 Resources"])
        with t1: st.write(f"Choosing **{best_path}** eliminates Plastic Tax penalties.")
        with t2: st.write("VoraCycle reduces environmental debt from 400+ years to <180 days.")
        with t3: st.write("Ensuring structural integrity without increasing virgin plastic density.")

master_prompt = (
    f"Role: Lead Forensic Arbiter (VoraCycle).\n"
    f"Subject: Strategic Pre-Lifecycle Audit for {query}.\n\n"
    
    f"### 🎯 STRATEGIC SCORECARD\n"
    f"- CONFIDENCE IN PATH: [0-100%]\n"
    f"- CIRCULARITY POTENTIAL: [High/Med/Low]\n"
    f"- ECONOMIC RETENTION: [Value preserved vs Value lost]\n\n"

    f"### 🧪 THE 'WHY': ENVIRONMENTAL & BUSINESS SYNERGY\n"
    f"Provide a 150-word forensic justification. Explain why this specific path (Recycle or send to waste) "
    f"is the 'Highest Form of Sustainability.' Discuss how this path protects the company's "
    f"Bottom Line by avoiding the 'Waste Tax' and Energy Loss associated with recycling or mineralization. "
    f"Be technical about the material's integrity.\n\n"

    f"### ⚖️ BEST PATH JUSTIFICATION (THE SELECTION LOGIC)\n"
    f"Explain the decision-making process. Why is this better than the alternative? "
    f"If you are avoiding waste, explain the 'Carbon Avoidance' achieved. If you are picking "
    f"the company's best path, explain the 'Supply Chain Independence' it creates. "
    f"Provide 150 words of rigorous analysis.\n\n"
    
    f"### 🛡️ ENVIRONMENTAL FAILSAFE\n"
    f"If this product MUST eventually become waste, how do we ensure it returns to the "
    f"environment safely? Detail the 'Biological Cycle' or 'Non-Toxic Decomposition' plan (150 words)."
)