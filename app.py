import streamlit as st
from openai import OpenAI
import json

# --- 1. ACCESS & CONNECTION ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except Exception as e:
    st.error("🚨 API Key Missing.")
    st.stop()

# --- 2. THE VORA 100 STRATEGIC REGISTRY ---
vora_100 = {
    "🥩 POULTRY & FRESH MEATS": ["MAP Poultry Trays", "Absorbent Poultry Pads", "Vacuum Wraps", "Black Meat Trays"],
    "🧻 PAPER & HYGIENE WRAPS": ["Kirkland Bath Tissue Case-Wrap", "Paper Towel Overwrap"],
    "🚩 HIGH-RISK LIABILITIES": ["PVC Clamshells", "PFAS Wrappers", "Lithium Battery Packs", "LLDPE Stretch Wrap"]
}

# --- 3. THE ARBITER ENGINE ---
def generate_vora_analysis(prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {"role": "system", "content": "You are the VoraCycle CSO. You specialize in forensic DNA reconstruction of industrial materials."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.3 
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Analysis Error: {str(e)}"

# --- 4. USER INTERFACE ---
st.set_page_config(page_title="VoraCycle: Executive Command", layout="wide")
st.title("🛡️ VoraCycle: Strategic DNA Command Center")

dropdown_items = ["-- Select a Strategic Asset --"]
for category, items in vora_100.items():
    dropdown_items.extend(items)

col1, col2 = st.columns(2)
with col1:
    dropdown_choice = st.selectbox("Select Asset for Executive Audit:", dropdown_items)
with col2:
    search_query = st.text_input("Search Custom SKU:")

final_query = search_query if search_query else (dropdown_choice if dropdown_choice != "-- Select a Strategic Asset --" else None)

if final_query:
    st.divider()
    with st.spinner(f"Reconstructing DNA for {final_query}..."):
        
        # This prompt specifically asks for the "Before and After" visuals
        master_prompt = f"""
        Execute a Forensic DNA Audit for: {final_query}.

        ### 🧬 DNA TRANSFORMATION SCORECARD
        Provide a side-by-side comparison:
        - **BEFORE DNA (Status Quo):** List the toxic elements, multi-material layers, and tax liabilities.
        - **AFTER DNA (Vora Design):** List the mono-material transition, specific Vora additives, and tax exemptions.

        ### 📊 EXECUTIVE IMPACT DASHBOARD
        Compare 'Business as Usual' vs. 'Vora DNA Optimization' in a Markdown table.

        ### 💰 FINANCIAL & ECOLOGICAL VERDICT
        Explain why this specific change is a 'Financial Fortress' for the brand.
        """

        st.markdown(generate_vora_analysis(master_prompt))
else:
    st.info("👆 Please select an asset to see the DNA Transformation.")