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
        
        toxic_atoms = any(a.GetSymbol() in ['Cl', 'F', 'Br', 'I'] for a in mol.GetAtoms())
        mw = Descriptors.MolWt(mol)
        
        # PATH 1: RECYCLING LOGIC
        br = 92 if "PET" in item_name else max(15, 80 - (mw / 10))
        if toxic_atoms: br -= 20
        ar = min(98.5, br + 45) 
        
        # PATH 2: MINERALIZATION LOGIC
        bm = 10 if toxic_atoms else max(20, 50 - (mw / 20))
        am = 99.4 if not toxic_atoms else 94.2 
        
        # BEST PATH ARBITRATION
        best_path = "Mineralization" if am > ar else "Mechanical Recycling"
        
        grade = "A" if am > 95 else "B+"
        comp_grade = "F" if toxic_atoms else ("D" if mw > 200 else "C")
        
        return round(br, 1), round(ar, 1), round(bm, 1), round(am, 1), best_path, grade, comp_grade, mw
    except:
        return None

# --- 4. THE APEX INTERFACE ---
st.set_page_config(page_title="VoraCycle Strategic Arbiter", layout="wide")
st.title("🔮 Wraith VoraCycle: Strategic Arbiter")
st.markdown("### *Forensic Deep-Dive: Material Evolution from Start-Line to Endgame*")

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
            st.warning(f"**Legacy Baseline ({br}%):** Driven by **Polymer Fragmentation**. High-heat reprocessing shears the carbon chains. This creates 'Downcycling' where the plastic becomes too brittle for food safety, requiring virgin plastic dilution.")
            st.success(f"**VoraCycle Optimized ({ar}%):** Achieved via **Dynamic Cross-Linking**. We introduce re-extenders that 'heal' the polymer during melting, maintaining high-tensile strength for infinite food-grade cycles.")

        with path_col2:
            st.subheader("🌿 Path 2: Soil Mineralization")
            st.metric("Before", f"{bm}%")
            st.metric("After VoraCycle", f"{am}%", delta=f"+{round(am-bm, 1)}% Improvement")
            st.warning(f"**Legacy Baseline ({bm}%):** Driven by **Hydrophobic Locking**. The plastic acts as a biological fortress. Microbes cannot attach to or penetrate the surface, leading to centuries of persistence.")
            st.success(f"**VoraCycle Optimized ({am}%):** Achieved via **Metabolic Triggering**. We insert latent scission points that act as enzymatic 'beacons.' When in soil, microbes recognize these sites and digest the plastic as a nutrient source.")

        # --- THE STRATEGIC DIRECTIVE ---
        st.divider()
        st.header(f"🏆 Strategic Directive: {best_path}")
        st.success(f"Optimizing for **{best_path}** yields the highest efficiency and sustainability ROI.")

        # --- THE "WHY": MONEY, TIME, RESOURCES ---
        st.header("⚖️ Resource Efficiency Analysis")
        t1, t2, t3 = st.tabs(["💰 Money", "⏳ Time", "🌍 Resources"])
        
        with t1:
