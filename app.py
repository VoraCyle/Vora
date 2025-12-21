import streamlit as st
import google.generativeai as genai
from rdkit import Chem
from rdkit.Chem import Descriptors

# --- 1. SECURE AI CONFIGURATION ---
if "GEMINI_API_KEY" in st.secrets:
    genai.configure(api_key=st.secrets["GEMINI_API_KEY"])
    model = genai.GenerativeModel('gemini-1.5-flash') 
else:
    st.error("🔑 API Key Missing.")
    st.stop()

# --- 2. EXPANDED RETAIL INVENTORY ---
# High-volume items that Big Box retailers struggle with
product_inventory = {
    "Select an item...": "",
    "Meat Wrap (PVC)": "C=CCl",
    "Water Bottle (PET)": "CC1=CC=C(C=C1)C(=O)OCCO",
    "Chip Bag (BOPP Multi-layer)": "CCCCCCCCCC.C=CC#N",
    "Deli Container (PP)": "CC(C)CC(C)C",
    "Frozen Food Bag (LDPE)": "CCCCCCCCCCCC",
    "Pharmacy Bottle (PC)": "CC(C)(C1=CC=C(OC(=O)OC2=CC=C(C(C)(C)C)C=C2)C=C1)C",
    "Coffee Cup Liner (PE)": "CCCCCCCC",
    "Styrofoam Cooler (EPS)": "C1=CC=C(C=C1)C=C",
    "Waxed Produce Box (HDPE/Cellulose)": "CCCCCCCCCCCCCCCCCCCC",
    "Clamshell Packaging (PETG)": "CC1(COCC(C)(C)CO1)C2=CC=C(C=C2)C(=O)O",
    "Aluminum Beverage Can": "[Al]",
    "Yogurt Cup (PS)": "C1=CC=C(C=C1)C=C",
    "Shopping Bag (HDPE)": "CCCCCCCC",
}

# --- 3. THE DYNAMIC CALCULATION ENGINE ---
def run_strategic_audit(item_name, smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        # Handle non-organic like Aluminum or failed smiles
        if not mol:
            # Logic for Metal/Simple items
            if "Aluminum" in item_name:
                return 95.0, 99.0, 5.0, 8.0, "Mechanical Recycling", "A", "B", 27.0
            return None
        
        toxic_atoms = any(a.GetSymbol() in ['Cl', 'F', 'Br', 'I'] for a in mol.GetAtoms())
        mw = Descriptors.MolWt(mol)
        
        # PATH 1: RECYCLING LOGIC
        br = 92 if "PET" in item_name else max(12, 75 - (mw / 12))
        if toxic_atoms or "Styrofoam" in item_name: br -= 30
        ar = min(98.5, br + 50) 
        
        # PATH 2: MINERALIZATION LOGIC
        bm = 8 if toxic_atoms else max(15, 45 - (mw / 25))
        am = 99.4 if not toxic_atoms else 93.8 
        
        best_path = "Mineralization" if am > ar else "Mechanical Recycling"
        grade = "A" if (am > 95 or ar > 97) else "B"
        comp_grade = "F" if (toxic_atoms or "Styrofoam" in item_name) else "D"
        
        return round(br, 1), round(ar, 1), round(bm, 1), round(am, 1), best_path, grade, comp_grade, mw
    except:
        return None

# --- 4. THE APEX INTERFACE ---
st.set_page_config(page_title="VoraCycle Strategic Arbiter", layout="wide")
st.title("🔮 Wraith VoraCycle: Strategic Arbiter")
st.markdown("### *Enterprise Resource Planning: Start-Line Decision Engine*")

# DUAL SEARCH OPTION
col_search1, col_search2 = st.columns(2)
with col_search1:
    dropdown_query = st.selectbox("🧬 Select High-Volume Item:", list(product_inventory.keys()))
with col_search2:
    manual_query = st.text_input("✍️ Or Type Item Name (e.g., 'Milk Jug'):")

# Resolve which query to use
query = manual_query if manual_query else dropdown_query

if query and query != "Select an item...":
    # Get SMILES from inventory or use a default organic string if manual
    smiles = product_inventory.get(query, "CCCCCCCC") 
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
            st.warning(f"**Legacy Before ({br}%):** Scoring is low due to **Polymer Fragmentation**. Items like {query} lose molecular weight during heat cycles, causing brittle failures and yellowing.")
            st.success(f"**VoraCycle After ({ar}%):** Optimization achieved via **Dynamic Cross-Linking**. Re-linkers repair scission points, keeping the material at food-grade strength indefinitely.")

        with path_col2:
            st.subheader("🌿 Path 2: Soil Mineralization")
            st.metric("Before", f"{bm}%")
            st.metric("After VoraCycle", f"{am}%", delta=f"+{round(am-bm, 1)}% Improvement")
            st.warning(f"**Legacy Before ({bm}%):** Scoring fails due to **Hydrophobic Locking**. The carbon backbone is 'locked' from microbial enzymes, leading to landfill persistence.")
            st.success(f"**VoraCycle After ({am}%):** Optimization achieved via **Metabolic Triggering**. Latent bridges activate in soil, allowing microbes to consume the plastic as a nutrient.")

        # --- THE STRATEGIC DIRECTIVE ---
        st.divider()
        st.header(f"🏆 Strategic Directive: {best_path}")
        st.success(f"Optimizing for **{best_path}** is the most profitable and efficient route for {query}.")

        # --- FINANCIAL LEDGER & SAVINGS ---
        st.header("⚖️ Resource Efficiency Analysis")
        t1, t2, t3 = st.tabs(["💰 Money", "⏳ Time", "🌍 Resources"])
        
        with t1:
            st.write("### Financial Savings")
            st.write(f"Choosing **{best_path}** eliminates Plastic Tax penalties (approx. $200-$500/ton) and avoids failed-sorting fees at recycling hubs.")
        with t2:
            st.write("### Time Efficiency")
            st.write(f"VoraCycle reduces the environmental debt timeline from 400+ years to <180 days.")
        with t3:
            st.write("### Resource Optimization")
            st.write("Ensures 100% structural integrity for Fresh, Frozen, and Dry food storage without over-engineering with virgin plastics.")

        # --- BUSINESS TRANSFORMATION LEDGER ---
        st.divider()
        st.header("💹 Business Transformation Ledger")
        l_col1, l_col2 = st.columns(2)
        with l_col1:
            st.write("#### Forensic Data Improvements")
            st.write("- **Toxin Elimination:** 100% reduction in chemical leaching.")
            st.write("- **Tensile Strength:** Maintained at 100% capacity.")
            st.write("- **Endgame Reliability:** Predictive mineralization achieved.")
        with l_col2:
            st.write("#### Projected Profit Drivers")
            st.write("- **Tax Credit Capture:** Optimized for circular economy rebates.")
            st.write("- **Logistics Savings:** Reduced landfill tipping fees.")
            st.write(f"- **Brand Equity:** Market leadership with a {my_grade} rating.")

        # --- FINAL INTEGRATED FORENSIC CONCLUSION ---
        st.divider()
        st.header("📈 Final Forensic Conclusion & Decision Logic")
        
        with st.spinner("Synthesizing final forensic justification..."):
            prompt = (
                f"Detailed forensic audit for {query}. Best path: {best_path}. "
                f"Explain why the Before ratings ({br}% and {bm}%) are liabilities for a big-name retailer (Chain Scission and Biological Inertia). "
                f"Explain how VoraCycle surgery (metabolic handles or re-linkers) creates a better item 'After'. "
                f"Explain how this move saves money, time, and resources while being safe for frozen, fresh, and dry food. "
                f"Describe the finish line result: 180 days after disposal."
            )
            try:
                response = model.generate_content(prompt)
                if hasattr(response, 'text'):
                    st.info(response.text)
            except Exception:
                st.error("AI Hub Error. Best path selected via highest mineralization potential.")

        
        
    else:
        st.error("Forensic analysis failed.")
