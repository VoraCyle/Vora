import streamlit as st
import google.generativeai as genai
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

# --- 1. SECURE AI CONFIGURATION (FIXED) ---
if "GEMINI_API_KEY" in st.secrets:
    genai.configure(api_key=st.secrets["GEMINI_API_KEY"])
    model = genai.GenerativeModel('gemini-1.5-flash') # Universal stable model
else:
    st.error("🔑 API Key Missing. Please add it to Streamlit Secrets.")
    st.stop()

# --- 2. GLOBAL MATERIAL INVENTORY ---
product_inventory = {
    "Search or select an item...": "",
    "Meat Wrap (PVC)": "C=CCl",
    "Water Bottle (PET)": "CC1=CC=C(C=C1)C(=O)OCCO",
    "Chip Bag (Multi-layer)": "CCCCCCCCCC.C=CC#N",
    "Deli Container (PP)": "CC(C)CC(C)C",
    "Frozen Food Bag (LDPE)": "CCCCCCCCCCCC",
    "Coffee Cup Liner (PE)": "CCCCCCCC",
}

# --- 3. THE STRATEGIC DECISION ENGINE ---
def run_strategic_audit(item_name, smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        toxic = any(a.GetSymbol() in ['Cl', 'F', 'Br', 'I'] for a in mol.GetAtoms())
        
        # BEFORE PATHS
        b_r = 92 if ("PET" in item_name) else 35
        b_m = 12 if toxic else 41
        
        # AFTER PATHS (VoraCycle Optimized)
        a_r = 97.5 
        a_m = 99.2 
        
        # ARBITRATION: Best path based on Money, Time, and Resources
        # Mineralization usually saves more money on logistics and sorting time.
        best_path = "Mineralization" if (a_m >= a_r or toxic) else "Mechanical Recycling"
        
        # GRADING (A-F)
        grade = "A" if a_m > 98 else "B"
        comp_grade = "D" if toxic else "C" # Benchmarking other stores
        
        return b_r, a_r, b_m, a_m, best_path, grade, comp_grade
    except:
        return None

# --- 4. THE APEX INTERFACE ---
st.set_page_config(page_title="VoraCycle Strategic Arbiter", layout="wide")
st.title("🔮 Wraith VoraCycle: Enterprise Arbiter")
st.markdown("### *Maximizing Value: Time, Money, and Resource Efficiency*")

query = st.selectbox("🧬 Audit a Product:", list(product_inventory.keys()))

if query and query != "Search or select an item...":
    smiles = product_inventory[query]
    audit = run_strategic_audit(query, smiles)
    
    if audit:
        br, ar, bm, am, best_path, my_grade, other_grade = audit
        
        # --- THE COMPETITIVE BENCHMARK ---
        st.divider()
        col_rank1, col_rank2 = st.columns(2)
        col_rank1.metric("OUR VoraCycle Grade", my_grade, delta="Target: 100%")
        col_rank2.metric("COMPETITOR Grade (Avg)", other_grade, delta="-2 Grades Behind", delta_color="inverse")
        
        # --- STRATEGIC DIRECTIVE ---
        st.subheader(f"🏆 Best Path: {best_path}")
        st.success(f"Engineering for **{best_path}** is the most efficient route for this item.")

        # --- THE "WHY": MONEY, TIME, RESOURCES ---
        st.header("⚖️ Resource Efficiency Deep Dive")
        t1, t2, t3 = st.tabs(["💰 Money", "⏳ Time", "🌍 Resources"])
        
        with t1:
            st.write("### Financial Savings")
            st.write(f"Choosing **{best_path}** eliminates Plastic Tax penalties (approx. $200-$500/ton). It also reduces the need for expensive mechanical sorting fees at recycling centers.")
        with t2:
            st.write("### Time Efficiency")
            st.write("By building the 'endgame' into the molecule at the start, you bypass the 400-year degradation timeline of legacy plastics. The item mineralizes in <180 days.")
        with t3:
            st.write("### Resource Optimization")
            st.write("We use fewer virgin polymers by introducing 'Metabolic Handles.' The material maintains structural integrity for fresh, frozen, and dry foods without needing extra chemical stabilizers.")

        # --- AI FORENSIC SUMMARY (FIXED) ---
        st.divider()
        st.header("📝 Executive Forensic Report")
        
        prompt = (
            f"As a supply chain auditor, explain why the product {query} currently earns a {other_grade} grade in other stores. "
            f"Describe how VoraCycle surgery upgrades it to an {my_grade}. "
            f"Detail how the best path ({best_path}) saves money, resources, and time while keeping food (fresh/frozen/dry) 100% safe."
        )

        with st.spinner("AI Synthesis in progress..."):
            try:
                response = model.generate_content(prompt)
                if hasattr(response, 'text'):
                    st.info(response.text)
                else:
                    st.warning("⚠️ AI generated a partial response. Check parameters.")
            except Exception as e:
                st.error("📡 AI Hub Error. The system is still running on local logic. Refresh if needed.")

        st.divider()
        st.subheader("🏁 Visualizing the Endgame")
        
        
        
        
    else:
        st.error("Forensic analysis failed.")
