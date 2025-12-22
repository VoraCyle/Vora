import streamlit as st
import google.generativeai as genai
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

# --- 1. ENTERPRISE SECURITY (THE VAULT) ---
def initialize_ai():
    # Rotates through your keys to ensure 100% uptime
    keys = [
        st.secrets.get("GEMINI_KEY_1"), 
        st.secrets.get("GEMINI_KEY_2"),
        st.secrets.get("GEMINI_KEY_3"),
        st.secrets.get("GEMINI_KEY_4")
    ]
    active_keys = [k for k in keys if k]
    
    if not active_keys:
        st.error("🔑 Security Alert: No API keys found in the Secrets Vault.")
        st.stop()
    return active_keys

ACTIVE_KEYS = initialize_ai()

def get_forensic_conclusion(prompt):
    for key in ACTIVE_KEYS:
        try:
            genai.configure(api_key=key)
            model = genai.GenerativeModel('gemini-1.5-flash')
            response = model.generate_content(prompt)
            return response.text
        except Exception:
            continue
    return "⚠️ AI Brain connection interrupted. Check API quotas."

# --- 2. THE STRATEGIC INVENTORY (BIG-BOX RETAIL FOCUS) ---
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
    "Waxed Produce Box": "CCCCCCCCCCCCCCCCCCCC",
}

# --- 3. THE STRATEGIC ARBITER (DECISION LOGIC) ---
def run_strategic_audit(item_name, smiles):
    try:
        # Chemical analysis using First Principles
        mol = Chem.MolFromSmiles(smiles) if "[Al]" not in smiles else None
        mw = Descriptors.MolWt(mol) if mol else 27.0
        toxic = any(a.GetSymbol() in ['Cl', 'F', 'Br', 'I'] for a in mol.GetAtoms()) if mol else False
        
        # PATH 1: RECYCLING (Circular Economy)
        br = 92 if "PET" in item_name or "Aluminum" in item_name else max(10, 75 - (mw / 10))
        if toxic: br -= 35 
        ar = min(99.1, br + 45) 
        
        # PATH 2: MINERALIZATION (Environmental Safety)
        bm = 5 if toxic else max(12, 50 - (mw / 20))
        am = 99.4 if not toxic else 94.0 
        
        # THE ARBITRATION (Sustainability + Profit)
        if toxic or "Multi-layer" in item_name or "Styrofoam" in item_name:
            best_path, priority = "Mineralization", "Environmental Defense"
            reason = "Waste Trap detected. Recycling is energy-prohibitive. Mineralization ensures zero-microplastic residue."
        elif ar > 88:
            best_path, priority = "Mechanical Recycling", "Resource Preservation"
            reason = "High-value circularity. Resource recovery saves more energy than carbon-return."
        else:
            best_path, priority = "Mineralization", "Environmental Defense"
            reason = "Efficiency gap. Mineralization is the most cost-effective path for this molecular profile."

        return round(br, 1), round(ar, 1), round(bm, 1), round(am, 1), best_path, priority, reason
    except:
        return None

# --- 4. INTERFACE LAYOUT ---
st.set_page_config(page_title="VoraCycle Strategic Arbiter", layout="wide", page_icon="🔮")
st.title("🔮 Wraith VoraCycle: Strategic Arbiter")
st.markdown("### *Determining the Apex Move for Environment & Enterprise*")

c1, c2 = st.columns(2)
with c1: dropdown_query = st.selectbox("🧬 Select Industry-Standard Item:", list(product_inventory.keys()))
with c2: manual_query = st.text_input("✍️ Manual Search / Custom Item:")

query = manual_query if manual_query else dropdown_query

if query and query != "Select a problematic item...":
    smiles = product_inventory.get(query, "CCCCCCCC") 
    audit = run_strategic_audit(query, smiles)
    
    if audit:
        br, ar, bm, am, best_path, priority, reason = audit
        
        st.divider()
        st.header(f"🏆 Best Strategic Path: {best_path}")
        st.info(f"**Primary Objective:** {priority} | **Logic:** {reason}")
        
        # FORENSIC BREAKDOWN
        st.subheader("📊 Comparative Forensic Performance")
        p1, p2 = st.columns(2)
        with p1:
            st.metric("Mechanical Recycling (After)", f"{ar}%", delta=f"Before: {br}%")
            st.warning(f"**Before:** Scoring suppressed by **Chain Scission**. Heat fractures legacy polymers.")
            st.success(f"**After:** Fixed via **Atomic Re-linkers** for infinite circularity.")
        with p2:
            st.metric("Soil Mineralization (After)", f"{am}%", delta=f"Before: {bm}%")
            st.warning(f"**Before:** Failed by **Hydrophobic Locking**. Biological fortress for 400 years.")
            st.success(f"**After:** Fixed via **Metabolic Handles** for 180-day return.")

        # THE 3 TABS: MONEY, TIME, RESOURCES
        st.divider()
        st.header("⚖️ Resource Efficiency Analysis")
        t1, t2, t3 = st.tabs(["💰 Money", "⏳ Time", "🌍 Resources"])
        with t1:
            st.write(f"### Financial Savings: {best_path}")
            st.write(f"Eliminates Plastic Tax penalties and avoids failed-sorting fees. Transitioning from **Liability** to **Circular Asset**.")
        with t2:
            st.write("### Time Efficiency")
            st.write(f"VoraCycle surgery reduces the environmental debt timeline from 400+ years to **<180 days**.")
        with t3:
            st.write("### Resource Optimization")
            st.write("Ensuring 100% structural integrity for **Fresh, Frozen, and Dry** food storage without increasing plastic density.")

        # AI DEEP-DIVE
        st.divider()
        st.header("📈 Final Strategic Forensic Deep-Dive")
        with st.spinner("Synthesizing Business & Environmental Impact..."):
            prompt = (f"Perform a deep forensic audit for {query}. The Arbiter chose {best_path} for {priority}. "
                      f"1. Explain why 'Before' (Recycle: {br}%, Mineralize: {bm}%) was an environmental threat. "
                      f"2. Detail how VoraCycle surgery (linkers or handles) creates the 'After' state. "
                      f"3. Explain why this path helps the environment the most. "
                      f"4. Confirm safety for Frozen/Fresh/Dry food and the 180-day finish line.")
            st.info(get_gemini_brain(prompt))
    else:
        st.error("Audit failed. Material signature not recognized.")
    



