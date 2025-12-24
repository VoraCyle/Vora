import streamlit as st
from openai import OpenAI

# --- 1. ACCESS & CONNECTION ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except Exception as e:
    st.error("🚨 API Key Missing: Please add 'OPENAI_API_KEY' to your Streamlit Secrets.")
    st.stop()

# --- 2. THE VORA 100 STRATEGIC REGISTRY ---
vora_100 = {
    "🥩 POULTRY & FRESH MEATS": ["MAP Poultry Trays", "Absorbent Poultry Pads", "Vacuum Wraps", "Black Meat Trays"],
    "🥬 FRESH PRODUCE & GOODS": ["Cellulose Berry Clamshells", "Bio-Produce Bags", "Waxed Boxes", "Mesh Citrus Bags"],
    "❄️ FROZEN & REFRIGERATED": ["Aqueous Frozen Bags", "Multi-Layer Meal Pouches", "Mono-PE Trays", "Poly-Ice Cream Cartons"],
    "📦 DRY GOODS & PANTRY": ["Metallized Snack Liners", "Composite Canisters", "BOPP Cereal Liners", "Multi-Wall Pet Food Bags"],
    "🚩 HIGH-RISK LIABILITIES": ["PVC Clamshells", "PFAS Wrappers", "Lithium Battery Packs", "LLDPE Stretch Wrap", "BPA Receipts"]
}

# --- 3. THE ARBITER ENGINE ---
def generate_vora_analysis(prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {"role": "system", "content": """You are the VoraCycle Chief Sustainability Officer. 
                Your specialty is 'Pre-Emptive Circularity'. You translate complex DNA engineering 
                into massive financial and ecological wins for Fortune 500 retailers. 
                Show how fixing the 'End Path' at the 'Start Line' creates a zero-waste, 
                zero-liability future."""},
                {"role": "user", "content": prompt}
            ],
            temperature=0.4 
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Analysis Error: {str(e)}"

# --- 4. USER INTERFACE ---
st.set_page_config(page_title="VoraCycle: Executive Command", layout="wide")
st.title("🛡️ VoraCycle: Strategic DNA Command Center")
st.markdown("### The Pre-Emptive Circularity Report: Future-Proofing the Enterprise.")

dropdown_items = ["-- Select a Strategic Asset --"]
for category, items in vora_100.items():
    dropdown_items.extend(items)

col1, col2 = st.columns(2)
with col1:
    dropdown_choice = st.selectbox("Select Asset for Executive Audit:", dropdown_items)
with col2:
    search_query = st.text_input("Search Custom SKU (e.g. 'Kirkland Salmon'):")

final_query = search_query if search_query else (dropdown_choice if dropdown_choice != "-- Select a Strategic Asset --" else None)

if final_query:
    st.divider()
    with st.spinner(f"Generating Executive Impact Report for {final_query}..."):
        
        master_prompt = f"""
        Generate a Pre-Emptive Impact Report for: {final_query}.

        ### 📊 EXECUTIVE IMPACT DASHBOARD (The 10-Year View)
        Compare 'Business as Usual' vs. 'Vora DNA Optimization'.
        Table: Environmental Impact (Toxins vs Nutrients), Financial Risk (High Fees vs Zero Fees), Consumer Trust (Liability vs Leader), and Supply Chain Speed.

        ### 🌍 THE ECOLOGICAL HEDGE (Healing the Planet)
        - **Pre-Emptive Deletion:** How does changing the DNA now remove the need for massive cleanup costs later?
        - **Earth-Native Stability:** Describe the transition from 'Persistent Pollution' to 'Safe Mineralization'.
        - **Resource Preservation:** How much raw material is saved by making this 100% technical-ready?

        ### 💰 THE FINANCIAL FORTRESS (Winning the Market)
        - **EPR Immunity:** Quantify the protection against 2025-2030 waste taxes and plastics bans.
        - **Operational Speed:** How does mono-material DNA simplify warehouse logistics and backhauling?
        - **The Brand Moat:** How does 'Consumer-Proof' sustainability protect the stock price and brand value?

        ### 🧬 THE DNA FAILSAFE & FOOD SAFETY
        - **Path-Agnostic Success:** Confirm the item succeeds in Waste or Recycle paths.
        - **Purity Guarantee:** Confirm the design is 100% bio-inert and safe for food contact.

        ### 🏁 THE CSO VERDICT
        A 3-sentence summary of why this specific SKU must be re-engineered immediately to lead the industry.
        """

        st.markdown(generate_vora_analysis(master_prompt))
else:
    st.info("👆 Please select an asset to generate the Strategic Impact Report.")

# --- 5. FOOTER ---
st.sidebar.info(f"VoraCycle v5.7.0 | The Future-Proof Enterprise")


