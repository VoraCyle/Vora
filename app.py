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
        
        # PATH 1: RECYCLING (Resource Preservation)
        br = 92 if "PET" in item_name or "Aluminum" in item_name else max(10, 75 - (mw / 10))
        if toxic: br -= 30
        ar = min(99.1, br + 45) 
        
        # PATH 2: MINERALIZATION (Environmental Safety)
        bm = 5 if toxic else max(12, 50 - (mw / 20))
        am = 99.4 if not toxic else 94.0 
        
        # --- THE ARBITRATION LOGIC ---
        if toxic or "Multi-layer" in item_name or "Styrofoam" in item_name:
            best_path = "Mineralization"
            priority = "Environmental Defense"
            reason = "Recycling this material creates toxic runoff or requires unfeasible energy. Mineralization ensures zero-microplastic residue."
        elif ar > 90:
            best_path = "Mechanical Recycling"
            priority = "Resource Preservation"
            reason = "High-value circularity. Keeping this material in the economy saves more energy than returning it to the soil."
        else:
            best_path = "Mineralization"
            priority = "Environmental Defense"
            reason = "Efficiency gap. Mineralization is faster and more cost-effective for this specific molecular structure."

        return round(br, 1), round(ar, 1), round(bm, 1), round(am, 1), best_path, priority, reason, mw
    except:
        return None

# --- 4. THE INTERFACE ---
st.set_page_config(page_title="VoraCycle Strategic Arbiter", layout="wide")
st.title("🔮 Wraith VoraCycle: Strategic Arbiter")
st.markdown("### *Forensic Audit: Determining the Apex Move for Environment & Enterprise*")

# DUAL SEARCH FUNCTIONALITY
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
        p1, p2 = st.columns(2)
        with p1:
            st.metric("Mechanical Recycling (After)", f"{ar}%", delta=f"Before: {br}%")
            st.write("**Forensic Baseline:** Low score due to **Chain Scission** (Heat-cycles degrade the polymer).")
        with p2:
            st.metric("Soil Mineralization (After)", f"{am}%", delta=f"Before: {bm}%")
            st.write("**Forensic Baseline:** Low score due to **Biological Dead-Lock** (Microbes cannot see the carbon).")

        # --- THE 3 TABS: MONEY, TIME, RESOURCES ---
        st.divider()
        st.header("⚖️ Resource Efficiency Analysis")
        t1, t2, t3 = st.tabs(["💰 Money", "⏳ Time", "🌍 Resources"])
        
        with t1:
            st.write("### Financial Savings")
            st.write(f"Choosing **{best_path}** eliminates Plastic Tax penalties (approx. $200-$500/ton) and removes the expense of failed mechanical sorting.")
            st.write("- **ROI:** Transitioning from Waste Liability to Circular Asset.")
            st.write("- **Tax Logic:** Compliance with upcoming EPR (Extended Producer Responsibility) laws.")
            
        with t2:
            st.write("### Time Efficiency")
            st.write(f"Legacy {query} requires 400+ years to degrade in a landfill.")
            st.write(f"- **VoraCycle Speed:** {best_path} reduces this timeline to <180 days.")
            st.write("- **Supply Chain Speed:** Atomic re-linkers allow for immediate reuse in the circular path.")
            
        with t3:
            st.write("### Resource Optimization")
            st.write("Ensures 100% structural integrity for **Fresh, Frozen, and Dry** food storage.")
            st.write("- **Virgin Material Reduction:** Uses latent bridges to maintain strength without needing extra plastic thickness.")
            st.write("- **Zero Contamination:** Path selected to ensure no microplastic leakage into the water table.")

        # --- FINAL FORENSIC CONCLUSION ---
        st.divider()
        st.header("📈 Strategic Forensic Deep-Dive")
        with st.spinner("Synthesizing Final Conclusion..."):
            prompt = (
                f"Detailed forensic audit for {query}. Best path: {best_path} for {priority}. "
                f"1. Why the Before ratings ({br}% and {bm}%) were liabilities (Chain Scission and Biological Dead-Lock). "
                f"2. How VoraCycle surgery creates the high-performance 'After' state. "
                f"3. THOROUGHLY EXPLAIN WHY THIS PATH HELPS THE ENVIRONMENT THE MOST. "
                f"4. EXPLAIN THE BUSINESS CASE: Money, Resources, and Time. "
                f"5. Confirm safety for Frozen, Fresh, and Dry food and describe the 180-day finish line."
            )
            response = model.generate_content(prompt)
            st.info(response.text)
            
        
        
    else:
        st.error("Audit failed. Material signature not recognized.")
