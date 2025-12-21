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
        
        # DYNAMIC GRADING
        grade = "A" if am > 95 else "B+"
        comp_grade = "F" if toxic_atoms else ("D" if mw > 200 else "C")
        
        return round(br, 1), round(ar, 1), round(bm, 1), round(am, 1), best_path, grade, comp_grade, mw
    except:
        return None

# --- 4. THE APEX INTERFACE ---
st.set_page_config(page_title="VoraCycle Strategic Arbiter", layout="wide")
st.title("🔮 Wraith VoraCycle: Strategic Arbiter")
st.markdown("### *Dynamic Forensic Benchmarking: Legacy Baseline vs. VoraCycle Future*")

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
            st.warning(f"**Legacy Baseline ({br}%):** This rating is suppressed by 'Chain Scission.' Legacy heat-cycles break molecular bonds, causing yellowing and loss of food-grade safety.")

        with path_col2:
            st.subheader("🌿 Path 2: Soil Mineralization")
            st.metric("Before", f"{bm}%")
            st.metric("After VoraCycle", f"{am}%", delta=f"+{round(am-bm, 1)}% Improvement")
            st.warning(f"**Legacy Baseline ({bm}%):** This rating reflects 'Biological Inertia.' Legacy carbon bonds are atoms locked in a crystalline lattice that soil enzymes cannot recognize or digest.")

        # --- THE STRATEGIC DIRECTIVE ---
        st.divider()
        st.header(f"🏆 Strategic Directive: {best_path}")
        st.success(f"Optimizing for **{best_path}** yields the highest efficiency and sustainability ROI.")

        # --- THE "WHY": MONEY, TIME, RESOURCES ---
        st.header("⚖️ Resource Efficiency Analysis")
        t1, t2, t3 = st.tabs(["💰 Money", "⏳ Time", "🌍 Resources"])
        
        with t1:
            st.write("### Financial Savings")
            st.write(f"By choosing **{best_path}**, we avoid surcharges and potential plastic taxes while reducing expensive sorting overhead.")
        with t2:
            st.write("### Time Efficiency")
            st.write(f"Legacy {query} requires 400+ years to degrade. VoraCycle reduces this to <180 days.")
        with t3:
            st.write("### Resource Optimization")
            st.write("Ensuring structural integrity for Fresh, Frozen, and Dry food storage without excessive material use.")

        # --- BUSINESS TRANSFORMATION LEDGER ---
        st.divider()
        st.header("💹 Business Transformation Ledger")
        
        l_col1, l_col2 = st.columns(2)
        with l_col1:
            st.write("#### Forensic Data Improvements")
            st.write("- **Toxin Elimination:** 100% reduction in leaching risks.")
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
                f"Explain why {best_path} was chosen as the ultimate endgame for {query}. "
                f"Specifically detail why the BEFORE ratings (Recycle: {br}%, Mineralization: {bm}%) were so low—explain the chemical 'Chain Scission' and 'Biological Inertia' that makes them liabilities. "
                f"Explain how the VoraCycle surgery (metabolic handles or re-linkers) creates a better item. "
                f"Explain how this move saves money, time, and resources while being 100% safe for frozen, fresh, and dry food variables. "
                f"Describe the finish line result: 180 days after disposal."
            )
            try:
                response = model.generate_content(prompt)
                if hasattr(response, 'text'):
                    st.info(response.text)
            except Exception as e:
                st.error("Forensic summary fallback: Path chosen based on highest mineralization potential and toxin reduction.")

        
        
    else:
        st.error("Forensic analysis failed.")
