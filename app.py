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
    "Select an item...": "",
    "Milk Jug (HDPE)": "CCCCCCCC",
    "Water Bottle (PET)": "CC1=CC=C(C=C1)C(=O)OCCO",
    "Aluminum Beverage Can": "[Al]",
    "Meat Wrap (PVC)": "C=CCl",
    "Deli Container (PP)": "CC(C)CC(C)C",
    "Frozen Food Bag (LDPE)": "CCCCCCCCCCCC",
    "Chip Bag (Multi-layer)": "CCCCCCCCCC.C=CC#N",
    "Styrofoam (EPS)": "C1=CC=C(C=C1)C=C",
}

# --- 3. THE STRATEGIC DECISION ENGINE ---
def run_strategic_audit(item_name, smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        mw = Descriptors.MolWt(mol) if mol else 27.0
        toxic = any(a.GetSymbol() in ['Cl', 'F', 'Br', 'I'] for a in mol.GetAtoms()) if mol else False
        
        # PATH 1: RECYCLING (Resource Preservation)
        br = 92 if "PET" in item_name or "Aluminum" in item_name else max(10, 75 - (mw / 10))
        if toxic: br -= 30
        ar = min(99.1, br + 45) 
        
        # PATH 2: MINERALIZATION (Environmental Safety)
        bm = 5 if toxic else max(12, 50 - (mw / 20))
        am = 99.4 if not toxic else 94.0 
        
        # --- THE ARBITRATION LOGIC ---
        # 1. Environment: Prevent microplastics and toxins.
        # 2. Business: Minimize tax, maximize material value.
        
        if toxic or "Multi-layer" in item_name or "Styrofoam" in item_name:
            best_path = "Mineralization"
            priority = "Environmental Defense"
            reason = "Recycling this material creates toxic runoff or high-energy waste. Mineralization ensures zero-microplastic residue."
        elif ar > 90:
            best_path = "Mechanical Recycling"
            priority = "Resource Preservation"
            reason = "High-value circularity. Keeping this material in the economy saves more energy than creating it from scratch."
        else:
            best_path = "Mineralization"
            priority = "Environmental Defense"
            reason = "Efficiency gap. Mineralization is faster and cheaper than complex sorting."

        return round(br, 1), round(ar, 1), round(bm, 1), round(am, 1), best_path, priority, reason, mw
    except:
        return None

# --- 4. THE INTERFACE ---
st.set_page_config(page_title="VoraCycle Strategic Arbiter", layout="wide")
st.title("🔮 Wraith VoraCycle: Strategic Arbiter")
st.markdown("### *Determining the Apex Move for Environment & Enterprise*")

query = st.selectbox("🧬 Forensic Material Selection:", list(product_inventory.keys()))

if query and query != "Select an item...":
    smiles = product_inventory.get(query, "CCCCCCCC") 
    audit = run_strategic_audit(query, smiles)
    
    if audit:
        br, ar, bm, am, best_path, priority, reason, mw = audit
        
        st.divider()
        st.subheader(f"🏆 Best Strategic Path: {best_path}")
        st.markdown(f"**Primary Objective:** {priority}")
        st.info(f"**Decision Logic:** {reason}")
        
        # DUAL-PATH PERFORMANCE
        p1, p2 = st.columns(2)
        with p1:
            st.metric("Mechanical Recycling (After)", f"{ar}%", delta=f"Base: {br}%")
            st.write("**Environmental Aspect:** Lowers carbon footprint by reducing virgin plastic demand.")
        with p2:
            st.metric("Soil Mineralization (After)", f"{am}%", delta=f"Base: {bm}%")
            st.write("**Environmental Aspect:** Eliminates microplastics and environmental persistence.")

        # --- FINAL FORENSIC CONCLUSION ---
        st.header("📈 Strategic Forensic Deep-Dive")
        with st.spinner("Synthesizing Business & Environmental Impact..."):
            prompt = (
                f"Perform a deep forensic audit for {query}. The Arbiter chose {best_path} for {priority}. "
                f"1. Explain why the 'Before' state (Recycle: {br}%, Mineralize: {bm}%) represents a threat to the environment and a liability for the business. "
                f"2. Describe the VoraCycle 'Molecular Surgery' performed to reach the 'After' state. "
                f"3. THOROUGHLY EXPLAIN WHY THIS PATH HELPS THE ENVIRONMENT THE MOST (e.g., zero microplastics, lower energy, or resource preservation). "
                f"4. EXPLAIN THE BUSINESS CASE: How this saves money (taxes), resources, and time while keeping food (frozen/fresh/dry) 100% safe. "
                f"5. Describe the finish line 180 days after disposal."
            )
            response = model.generate_content(prompt)
            st.info(response.text)
        
    else:
        st.error("Audit failed.")
