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
    "🥩 POULTRY & FRESH MEATS": ["MAP Poultry Trays", "Absorbent Poultry Pads", "Black Meat Trays"],
    "🥬 FRESH PRODUCE & GOODS": ["Cellulose Berry Clamshells", "Bio-Produce Bags", "Waxed Boxes"],
    "❄️ FROZEN & REFRIGERATED": ["Aqueous Frozen Bags", "Multi-Layer Meal Pouches", "Mono-PE Trays"],
    "📦 DRY GOODS & PANTRY": ["Metallized Snack Liners", "Composite Canisters", "BOPP Cereal Liners"],
    "🚩 HIGH-RISK LIABILITIES": ["PVC Clamshells", "PFAS Wrappers", "Lithium Battery Packs", "LLDPE Stretch Wrap"]
}

# --- 3. THE ARBITER ENGINE ---
def generate_vora_analysis(prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {"role": "system", "content": """You are the VoraCycle Lead Forensic Engineer. 
                Your primary task is to contrast current 'Toxic DNA' with 'Vora Optimized DNA'. 
                You must provide clear ratings (0-10) and explain the EXACT reason why the new DNA 
                outperforms the old one in every waste scenario."""},
                {"role": "user", "content": prompt}
            ],
            temperature=0.3 
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Analysis Error: {str(e)}"

# --- 4. USER INTERFACE ---
st.set_page_config(page_title="VoraCycle: Forensic Comparison", layout="wide")
st.title("🛡️ VoraCycle: Strategic DNA Command Center")
st.markdown("### Contrasting Current Liabilities with Pre-Emptive Circular Assets.")

dropdown_items = ["-- Select a Strategic Asset --"]
for category, items in vora_100.items():
    dropdown_items.extend(items)

col1, col2 = st.columns(2)
with col1:
    dropdown_choice = st.selectbox("Select Asset for Comparison:", dropdown_items)
with col2:
    search_query = st.text_input("Custom SKU Audit:")

final_query = search_query if search_query else (dropdown_choice if dropdown_choice != "-- Select a Strategic Asset --" else None)

if final_query:
    st.divider()
    with st.spinner(f"Running Comparative DNA Analysis for {final_query}..."):
        
        master_prompt = f"""
        Execute a side-by-side Forensic Audit for: {final_query}.

        ### 📊 THE FORENSIC COMPARISON TABLE
        Create a Markdown Table with the following columns: 
        | Feature | Status Quo (Before) | Vora Optimized (After) | Improvement Rating (1-10) |
        | :--- | :--- | :--- | :--- |
        | **Material DNA** | [Describe Current Mix] | [Describe Vora Swap] | [Score] |
        | **Waste Outcome** | [Describe Toxin/Landfill] | [Describe Mineralization] | [Score] |
        | **Recycle Value** | [Downcycle/Zero Value] | [High-Purity Asset] | [Score] |
        | **Food Safety** | [Leaching Risk] | [Bio-Inert Purity] | [Score] |
        | **Financial Risk** | [High EPR Fines] | [Zero Liability] | [Score] |

        ---

        ### 🧬 THE TRANSFORMATION LOGIC (The "Why")
        Provide a 100-word explanation of why the 'After' DNA improves upon the 'Before' DNA. 
        Focus on how removing 'Material Incompatibility' makes the product immune to consumer error.

        ### 🍎 FOOD SAFETY & PURITY ASSURANCE
        Detail why the Vora DNA is safer for direct food contact and eliminates chemical migration.

        ### 💰 PRE-EMPTIVE FINANCIAL ROI
        Estimate the savings in 'Waste Taxes' and 'Supply Chain Efficiency' gained by this switch.

        ### 🏁 FINAL VERDICT
        Why is this the ONLY path for a company like Costco to become a true sustainability leader?
        """

        st.markdown(generate_vora_analysis(master_prompt))
else:
    st.info("👆 Select an item to see the Before vs. After Forensic Comparison.")

# --- 5. FOOTER ---
st.sidebar.info("VoraCycle v5.7.0 | Comparative Forensic Engine")
