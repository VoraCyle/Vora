import streamlit as st
from openai import OpenAI

# --- 1. ACCESS & CONNECTION ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except Exception as e:
    st.error("🚨 API Key Missing: Please add 'OPENAI_API_KEY' to your Streamlit Secrets.")
    st.stop()

# --- 2. THE VORA 100 CATEGORIES ---
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

# --- 3. THE ARBITER ENGINE ---
def generate_vora_analysis(prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {"role": "system", "content": "You are the VoraCycle Lead Arbiter. You specialize in 'Path-Agnostic' DNA engineering."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.4 
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Analysis Error: {str(e)}"

# --- 4. USER INTERFACE ---
st.set_page_config(page_title="VoraCycle: Resilient Design", layout="wide")
st.title("🛡️ VoraCycle: Strategic DNA Command Center")
st.markdown("#### Engineering for 'Universal Sustainability'—Safe in Any Waste Stream.")

# Category Dropdown Logic
dropdown_items = ["-- Select a Strategic Asset --"]
for category, items in vora_100.items():
    dropdown_items.extend(items)

col1, col2 = st.columns(2)
with col1:
    st.markdown("##### 📦 1. Strategic Liability Selection")
    dropdown_choice = st.selectbox("Choose a known vulnerability:", dropdown_items)

with col2:
    st.markdown("##### 🔍 2. Custom SKU / Material Audit")
    search_query = st.text_input("Or type a custom material (e.g. 'Styrofoam'):")

# --- 5. THE LOGIC GATE (The Fix) ---
final_query = None

# If user typed something, that takes priority
if search_query.strip():
    final_query = search_query
# If user selected a REAL item from the dropdown, use that
elif dropdown_choice != "-- Select a Strategic Asset --":
    final_query = dropdown_choice

# ONLY run if final_query has a real value
if final_query:
    st.divider()
    with st.spinner(f"Simulating Forensic Path for {final_query}..."):
        
        master_prompt = f"""
        Audit the following for an Enterprise Retailer: {final_query}.

        ### 📊 DUAL-PATH FORENSIC SCORECARD
        Compare 'Status Quo Design' vs. 'Vora Resilient Design'.

        ### 🧬 DNA TRANSFORMATION: THE "PATH-AGNOSTIC" FAILSAFE
        Explain how to change the materials NOW so the item is safe in EVERY scenario:
        1. **PATH A: IF IT GOES TO WASTE (Landfill):** How do we change the DNA so it safely mineralizes?
        2. **PATH B: IF IT GOES TO RECYCLE:** How do we simplify the DNA for high-value technical recovery?

        ### 💰 COST-BENEFIT ESTIMATOR (ROI)
        - **EPR Savings:** Estimated reduction in government waste penalties.
        - **Operational ROI:** Savings from mono-material sourcing.

        ### 📋 VORA SUPPLIER SCORECARD
        Evaluate the manufacturer. Grade: Material Purity, Chemical Transparency, and Overall Vora Grade.

        ---

        ### 🏁 EXECUTIVE SUMMARY
        Provide the 'Zero-Liability' brand advantage verdict. (200 Words).
        """

        st.markdown(generate_vora_analysis(master_prompt))
else:
    # This shows when the home screen is empty
    st.info("👆 Please select an asset from the list or type a custom material above to begin the Forensic Audit.")

# --- 6. FOOTER ---
st.sidebar.info(f"VoraCycle v5.4.0 | {sum(len(v) for v in vora_100.values())} Liabilities Cataloged")