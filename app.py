import streamlit as st
from openai import OpenAI

# --- 1. ACCESS & CONNECTION ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except Exception as e:
    st.error("🚨 API Key Missing: Please add 'OPENAI_API_KEY' to your Streamlit Secrets.")
    st.stop()

# --- 2. THE VORA 100: STRATEGIC CATEGORIES ---
vora_100 = {
    "🚩 TOP LIABILITIES (High EPR Fines)": [
        "Multi-Layer Mylar Snack Pouches", "Black Plastic Rotisserie Trays", 
        "PVC Clamshell Electronics Packs", "PFAS-Coated Fast Food Wraps",
        "Aseptic Juice Cartons (Tetra Pak)", "BOPP Non-Recyclable Pet Food Bags"
    ],
    "⚡ QUICKEST FIXES (Immediate ROI)": [
        "LLDPE Warehouse Stretch Wrap", "BPA-Free Thermal Receipts", 
        "Mono-material Polyethylene Mailers", "Uncoated Corrugated Cardboard",
        "Natural Fiber Pallet Strapping", "Detectable Pigment Black Plastics"
    ],
    "🧬 LONG-TERM DNA SHIFTS (High Impact)": [
        "Lithium-Ion Battery Slurry", "Mixed-Textile Kirkland Apparel",
        "Integrated LED Appliance Panels", "Multi-Material Athletic Footwear",
        "EPS (Styrofoam) Cold-Chain Coolers", "Composite Construction Returns"
    ]
}

# Flatten for dropdown
dropdown_items = ["-- Select a Strategic Asset --"]
for category, items in vora_100.items():
    dropdown_items.extend(items)

# --- 3. THE ARBITER ENGINE ---
def generate_vora_analysis(prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {"role": "system", "content": "You are the VoraCycle Lead Arbiter. You specialize in Liability Removal and ROI-driven Circular Engineering."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.7 
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Analysis Error: {str(e)}"

# --- 4. USER INTERFACE ---
st.set_page_config(page_title="VoraCycle: Liability Command", layout="wide")
st.title("🛡️ VoraCycle: Strategic Liability Command")
st.markdown("### Converting Environmental Waste Liabilities into Financial Assets.")

col1, col2 = st.columns(2)
with col1:
    st.markdown("##### 📦 1. Strategic Asset Selection")
    dropdown_choice = st.selectbox("Select by Category (Liabilities vs. Quick Fixes):", dropdown_items)

with col2:
    st.markdown("##### 🔍 2. Custom SKU / DNA Audit")
    search_query = st.text_input("Search any retail product:")

final_query = search_query if search_query else (dropdown_choice if dropdown_choice != "-- Select a Strategic Asset --" else None)

if final_query:
    st.divider()
    with st.spinner(f"Calculating ROI and DNA Shift for {final_query}..."):
        
        master_prompt = f"""
        Conduct a Forensic Audit and ROI Estimate for: {final_query}.
        
        ### 📊 DUAL-PATH FORENSIC SCORECARD
        Compare 'Status Quo' vs. 'Vora Resilient Design'.
        Table: Primary Outcome, Sustainability Rating, Capital Retention %, and Environmental Toxin Risk.

        ### 💰 COST-BENEFIT ESTIMATOR (FINANCIAL WIN)
        1. **EPR Penalty Avoidance:** How much does this save the company in government waste fines?
        2. **Logistics Efficiency:** How does the DNA change reduce shipping weight or volume?
        3. **Supply Chain Savings:** How does switching to mono-materials lower raw material costs?

        ### 🧬 DNA ENGINEERING: THE "START-AT-THE-END" BLUEPRINT
        - **WHAT TO REMOVE:** Identify the specific 'System Spoilers'.
        - **WHAT TO REPLACE WITH:** Suggest the mono-material/bio-mineral fix.
        - **UNIVERSAL SAFETY:** Explain how this design stays safe even in a landfill.

        ### 📋 VORA SUPPLIER SCORECARD
        Evaluate the manufacturer. Table: Material Purity, Chemical Transparency, Circular Readiness, Vora Grade.

        ---

        ### 🏁 EXECUTIVE SUMMARY: THE LIABILITY REMOVAL VERDICT
        Focus on the 'Zero-Liability' brand advantage. Why is this rebuild a profit-center, not a cost-center?
        """

        st.markdown(generate_vora_analysis(master_prompt))

# --- 5. FOOTER ---
st.sidebar.info(f"VoraCycle v5.0.0 | Strategic Assets: {sum(len(v) for v in vora_100.values())}")