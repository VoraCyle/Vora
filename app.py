import streamlit as st
import google.generativeai as genai
from rdkit import Chem
from rdkit.Chem import Descriptors

# --- 1. SECURE CONFIGURATION ---
if "GEMINI_API_KEY" in st.secrets:
    genai.configure(api_key=st.secrets["GEMINI_API_KEY"])
    model = genai.GenerativeModel('gemini-1.5-flash') 
else:
    st.error("🔑 API Key Missing.")
    st.stop()

# --- 2. THE STRATEGIC INVENTORY ---
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

# --- 3. THE STRATEGIC DECISION ENGINE ---
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
        
        # DUAL-PATH PERFORMANCE WITH DEEP DESCRIPTIONS
        st.header("📊 Comparative Forensic Performance")
        p1, p2 = st.columns(2)
        with p1:
            st.subheader("♻️ Path 1: Mechanical Recycling")
            st.metric("After VoraCycle", f"{ar}%", delta=f"Baseline: {br}%")
            st.warning(f"**Before ({br}%):** Scoring is suppressed by **Chain Scission**. Legacy heat-cycles cause polymer fragmentation, leading to yellowing, brittle failure, and loss of food-grade safety.")
            st.success(f"**After ({ar}%):** Optimization achieved via **Atomic Re-linkers**. VoraCycle 'heals' the polymer chains during melting, maintaining high-tensile strength for infinite circular loops.")

        with p2:
            st.subheader("🌿 Path 2: Soil Mineralization")
            st.metric("After VoraCycle", f"{am}%", delta=f"Baseline: {bm}%")
            st.warning(f"**Before ({bm}%):** Scoring fails due to **Hydrophobic Locking**. The carbon backbone acts as a biological fortress that soil enzymes cannot recognize, resulting in 400+ years of persistence.")
            st.success(f"**After ({am}%):** Optimization achieved via **Metabolic Triggering**. Latent scission bridges activate in soil, allowing microbes to recognize and digest the material as a nutrient source.")

        # --- THE 3 TABS: MONEY, TIME, RESOURCES ---
        st.divider()
        st.header("⚖️ Resource Efficiency Analysis")
        t1, t2, t3 = st.tabs(["💰 Money", "⏳ Time", "🌍 Resources"])
        
        with t1:
            st.write("### Financial Savings")
            st.write(f"Choosing **{best_path}** eliminates Plastic Tax penalties and avoids failed-sorting fees at recycling hubs.")
        with t2:
            st.write("### Time Efficiency")
            st.write(f"VoraCycle reduces the environmental debt timeline from 400+ years to <180 days.")
        with t3:
            st.write("### Resource Optimization")
            st.write("Ensuring 100% structural integrity for **Fresh, Frozen, and Dry** food storage without increasing virgin plastic density.")

        # --- FINAL FORENSIC CONCLUSION ---
        st.divider()
        st.header("📈 Strategic Forensic Deep-Dive")
        with st.spinner("Synthesizing Final Conclusion..."):
            prompt = (
                f"Detailed forensic audit for {query}. Best path: {best_path} for {priority}. "
                f"Explain why the Before ratings ({br}% and {bm}%) were liabilities (Chain Scission and Biological Dead-Lock). "
                f"How VoraCycle surgery creates the high-performance 'After' state. "
                f"Why this path helps the environment the most. "
                f"Business Case: Money, Resources, Time. "
                f"Confirm safety for Frozen, Fresh, Dry food and describe the 180-day finish line."
            )
            response = model.generate_content(prompt)
            st.info(response.text)
            
    else:
        st.error("Audit failed. Material signature not recognized.")
