import streamlit as st
from openai import OpenAI

# --- 1. ACCESS & CONNECTION ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except Exception as e:
    st.error("🚨 API Key Missing: Please add 'OPENAI_API_KEY' to your Streamlit Secrets.")
    st.stop()

# --- 2. STRATEGIC ASSET REGISTRY (The Hard 100 Categories) ---
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

# --- 3. THE PATH-AGNOSTIC ENGINE ---
def generate_vora_analysis(prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {"role": "system", "content": """You are the VoraCycle Lead Arbiter. 
                You specialize in 'Resilient Circularity'—engineering products so they remain 
                sustainable regardless of how the consumer disposes of them. 
                Your mission is to fix the DNA so the 'End-Game' is safe by default.
                Provide forensic engineering depth: molecular recycling, bio-mineralization, and mono-material purity."""},
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

col1, col2 = st.columns(2)
with col1:
    st.markdown("##### 📦 1. Strategic Liability Selection")
    dropdown_choice = st.selectbox("Select by Category:", dropdown_items)
with col2:
    st.markdown("##### 🔍 2. Custom SKU / Material Audit")
    search_query = st.text_input("Custom SKU/Material Audit:")

final_query = search_query if search_query else (dropdown_choice if dropdown_choice != "-- Select a High-Impact Vulnerability --" else None)

if final_query:
    st.divider()
    with st.spinner("Simulating Universal Sustainability Outcomes..."):
        
        master_prompt = f"""
        Audit the following for an Enterprise Retailer: {final_query}.

        ### 📊 DUAL-PATH FORENSIC SCORECARD
        Compare 'Status Quo Design' vs. 'Vora Resilient Design'.
        - Sustainability Rating (0-10)
        - Outcome if Landfilled (Toxin vs Nutrient)
        - Outcome if Recycled (Downcycled vs Upcycled)
        - Financial Risk (High vs Zero)

        ---

        ### 🧬 DNA TRANSFORMATION: THE "PATH-AGNOSTIC" FAILSAFE
        Explain how to change the materials NOW so the item is safe in EVERY scenario:
        1. **PATH A: IF IT GOES TO WASTE (Landfill):** How do we change the chemicals/binders/DNA so it safely mineralizes as a non-toxic biological nutrient?
        2. **PATH B: IF IT GOES TO RECYCLE:** How do we simplify the DNA (e.g., Mono-materials, Detectable Pigments) to ensure it stays high-value technical nutrients?
        *Focus on removing 'Material Incompatibility' so the consumer cannot 'break' the system.* (300+ Words).

        ### 💰 COST-BENEFIT ESTIMATOR (ROI)
        - **EPR Savings:** Estimated reduction in government waste penalties.
        - **Operational ROI:** Savings from mono-material sourcing and logistics efficiency.

        ### 📋 VORA SUPPLIER SCORECARD
        Evaluate the theoretical manufacturer. Grade: Material Purity, Chemical Transparency, and Overall Vora Grade (A-F).

        ---

        ### 🏁 EXECUTIVE SUMMARY: THE SUSTAINABILITY LEADERSHIP VERDICT
        Explain why this design makes the company an industry leader. How does creating a 'consumer-proof' sustainable product eliminate future liability and maximize brand trust? (200 Words).
        """

        st.markdown(generate_vora_analysis(master_prompt))

# --- 5. FOOTER ---
st.sidebar.info(f"VoraCycle Enterprise Logic: v5.3.0 | Assets Cataloged: {sum(len(v) for v in vora_100.values())}")