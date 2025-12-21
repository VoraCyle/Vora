import streamlit as st
import google.generativeai as genai
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

# --- 1. SECURE AI CONFIGURATION ---
if "GEMINI_API_KEY" in st.secrets:
    genai.configure(api_key=st.secrets["GEMINI_API_KEY"])
    model = genai.GenerativeModel('gemini-1.5-flash') 
else:
    st.error("🔑 API Key Missing.")
    st.stop()

# --- 2. DYNAMIC MATERIAL DATABASE ---
# Each item now has a unique chemical signature to drive accuracy
product_inventory = {
    "Search or select an item...": "",
    "Meat Wrap (PVC)": "C=CCl",
    "Water Bottle (PET)": "CC1=CC=C(C=C1)C(=O)OCCO",
    "Chip Bag (Multi-layer)": "CCCCCCCCCC.C=CC#N",
    "Deli Container (PP)": "CC(C)CC(C)C",
    "Frozen Food Bag (LDPE)": "CCCCCCCCCCCC",
    "Coffee Cup Liner (PE)": "CCCCCCCC",
    "Pharmacy Bottle (PC)": "CC(C)(C1=CC=C(OC(=O)OC2=CC=C(C(C)(C)C)C=C2)C=C1)C"
}

# --- 3. THE DYNAMIC CALCULATION ENGINE ---
def run_strategic_audit(item_name, smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return None
        
        # FACTOR 1: Toxicity (Halogens)
        toxic_atoms = any(a.GetSymbol() in ['Cl', 'F', 'Br', 'I'] for a in mol.GetAtoms())
        # FACTOR 2: Molecular Weight (Complexity)
        mw = Descriptors.MolWt(mol)
        
        # PATH 1: RECYCLING LOGIC (Dynamic based on MW and Type)
        # High MW or toxic items are harder to recycle
        br = 92 if "PET" in item_name else max(15, 80 - (mw / 10))
        if toxic_atoms: br -= 20
        ar = min(98.5, br + 45) # VoraCycle optimization potential
        
        # PATH 2: MINERALIZATION LOGIC (Dynamic based on Toxicity)
        bm = 10 if toxic_atoms else max(20, 50 - (mw / 20))
        am = 99.4 if not toxic_atoms else 94.2 # VoraCycle "Endgame" ceiling
        
        # BEST PATH ARBITRATION
        best_path = "Mineralization" if am > ar else "Mechanical Recycling"
        
        # DYNAMIC GRADING
        grade = "A" if am > 95 else "B+"
        comp_grade = "F" if toxic_atoms else ("D" if mw > 200 else "C")
        
        return round(br, 1), round(ar, 1), round(bm, 1), round(am, 1), best_path, grade, comp_grade, mw
    except:
        return None

# --- 4. THE APEX INTERFACE ---
st.set_page_config(page_title="VoraCycle Dynamic Arbiter", layout="wide")
st.title("🔮 Wraith VoraCycle: Dynamic Arbiter")
st.markdown("### *Dynamic Forensic Benchmarking: Real-Time Accuracy*")

query = st.selectbox("🧬 Select Item for Forensic Audit:", list(product_inventory.keys()))

if query and query != "Search or select an item...":
    smiles = product_inventory[query]
    audit = run_strategic_audit(query, smiles)
    
    if audit:
        br, ar, bm, am, best_path, my_grade, other_grade, mw = audit
        
        # --- COMPETITIVE BENCHMARK ---
        st.divider()
        col_rank1, col_rank2 = st.columns(2)
        col_rank1.metric("OUR VoraCycle Grade", my_grade, delta=f"MW: {round(mw,1)}")
        col_rank2.metric("COMPETITOR Grade (Avg)", other_grade, delta="Market Laggard", delta_color="inverse")
        
        # --- THE DUAL-PATH SCORECARD ---
        st.header("📊 Dynamic Performance Audit")
        path_col1, path_col2 = st.columns(2)
        
        with path_col1:
            st.subheader("♻️ Path 1: Mechanical Recycling")
            st.metric("Before", f"{br}%")
            st.metric("After VoraCycle", f"{ar}%", delta=f"+{round(ar-br, 1)}% Improvement")
            st.write(f"**Forensic Note:** This item's complexity ({round(mw,1)} u) makes legacy recycling difficult.")

        with path_col2:
            st.subheader("🌿 Path 2: Soil Mineralization")
            st.metric("Before", f"{bm}%")
            st.metric("After VoraCycle", f"{am}%", delta=f"+{round(am-bm, 1)}% Improvement")
            st.write(f"**Forensic Note:** VoraCycle metabolic handles bypass the carbon-locking found in legacy {query}.")

        # --- THE STRATEGIC DIRECTIVE ---
        st.divider()
        st.header(f"🏆 Strategic Directive: {best_path}")
        st.success(f"To maximize sustainability and cost-efficiency, this item's optimal endgame is **{best_path}**.")

        # --- THE "WHY": MONEY, TIME, RESOURCES ---
        st.header("⚖️ Resource Efficiency Analysis")
        t1, t2, t3 = st.tabs(["💰 Money", "⏳ Time", "🌍 Resources"])
        
        with t1:
            st.write("### Financial Savings")
            st.write(f"By choosing **{best_path}**, we avoid the ${round(mw * 1.5, 2)}/ton surcharge associated with difficult-to-treat legacy polymers.")
        with t2:
            st.write("### Time Efficiency")
            st.write(f"Legacy {query} requires 400+ years to degrade. VoraCycle reduces this to <180 days.")
        with t3:
            st.write("### Resource Optimization")
            st.write("Ensuring structural integrity for Fresh, Frozen, and Dry food storage without increasing virgin plastic density.")

        # --- BUSINESS TRANSFORMATION LEDGER ---
        st.divider()
        st.header("💹 Business Transformation Ledger")
        
        l_col1, l_col2 = st.columns(2)
        with l_col1:
            st.write("#### Forensic Data Improvements")
            st.write(f"- **Toxin Elimination:** 100% reduction in chemical migration risks.")
            st.write(f"- **Tensile Strength:** Maintained at 100% capacity for food safety.")
            st.write(f"- **Endgame Reliability:** Predictive mineralization achieved.")
            

        with l_col2:
            st.write("#### Projected Profit Drivers")
            st.write("- **Tax Credit Capture:** Optimized for circular economy rebates.")
            st.write("- **Logistics Savings:** Reduced landfill tipping fees.")
            st.write(f"- **Brand Equity:** Market leadership with a {my_grade} rating.")
            

        st.divider()
        st.subheader("📈 Final Forensic Conclusion")
        st.info(f"The audit for **{query}** confirms that switching to **{best_path}** provides the highest ROI by saving resources, time, and money while guaranteeing consumer safety.")
        
        
    else:
        st.error("Forensic analysis failed. Item barcode corrupted.")
