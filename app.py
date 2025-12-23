import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from openai import OpenAI  # We only need OpenAI now!

# This connects your app to OpenAI using the secret key you'll put in Streamlit
client = OpenAI(api_key=st.secrets["openai"]["OPENAI_API_KEY"])

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

# --- FINAL FORENSIC CONCLUSION (HIGH-DENSITY VERSION) ---
st.divider()
st.header("📈 Deep Forensic & Strategic Analysis")

if query: 
    with st.spinner("VoraCycle Arbiter (GPT-4o) is conducting a deep strategic audit..."):
        # The Secret Sauce: We give it a 'Persona' and 'Word Count' requirements
        master_prompt = (
            f"Role: Senior Forensic Resource Strategist.\n"
            f"Task: Conduct a high-level strategic audit for: {query}.\n\n"
            f"INSTRUCTIONS FOR DEPTH:\n"
            f"1. 💰 RESOURCE EFFICIENCY (DEEP DIVE): Provide a 200-word forensic analysis of "
            f"WHY this path is beneficial. Use terms like 'Capital Resilience', 'Opportunity Cost', "
            f"and 'EBITDA impact'. Explain the invisible waste that traditional audits miss.\n\n"
            f"2. 🛡️ BEST STRATEGIC PATH (LOGIC): Provide a 200-word technical justification for the "
            f"chosen strategy. Explain the trade-offs. Why this specific path over other options? "
            f"Analyze the risk mitigation and long-term velocity gains.\n\n"
            f"3. 🌍 SYSTEMIC OUTCOME: A final 150-word summary on asset optimization.\n\n"
            f"FORMATTING: Use ### Headers for each section. DO NOT use bullet points. "
            f"Write in full, professional, dense paragraphs. Be authoritative and decisive."
        )

        full_analysis = generate_vora_analysis(master_prompt)
        
        # Display in a professional container
        st.info(full_analysis)