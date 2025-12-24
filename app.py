import streamlit as st
from openai import OpenAI
import json

# --- 1. ACCESS & CONNECTION ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except Exception as e:
    st.error("🚨 API Key Missing.")
    st.stop()

# --- 2. THE VORA 100 REGISTRY ---
vora_100 = {
    "🥩 POULTRY & FRESH MEATS": ["MAP Poultry Trays", "Absorbent Poultry Pads", "Vacuum Wraps", "Black Meat Trays"],
    "🧻 PAPER & HYGIENE WRAPS": ["Kirkland Bath Tissue Case-Wrap", "Paper Towel Overwrap"],
    "🚩 HIGH-RISK LIABILITIES": ["PVC Clamshells", "PFAS Wrappers", "LLDPE Stretch Wrap"]
}

# --- 3. THE ARBITER ENGINE ---
def generate_vora_analysis(prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {"role": "system", "content": "You are the VoraCycle CSO. You create high-impact, visual material audits."},
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
    
    # --- NEW VISUAL DNA BLOCKS ---
    st.subheader(f"🧬 DNA Forensic Transformation: {final_query}")
    
    # Create the Side-by-Side columns
    before_col, after_col = st.columns(2)
    
    with before_col:
        st.error("### 🔴 BEFORE: Status Quo DNA")
        st.markdown("""
        - **Materials:** Multi-layer Laminates / Composites
        - **Path:** Linear (Landfill/Incineration)
        - **Risk:** High EPR Tax Liability
        - **DNA Markers:** PFAS, Toxic Glues, Mixed Polymers
        """)
        
    with after_col:
        st.success("### 🟢 AFTER: Vora Resilient DNA")
        st.markdown("""
        - **Materials:** Mono-Material / Vora-Infused
        - **Path:** Circular (Infinite Loop / Safe Soil)
        - **Risk:** EPR Tax Exempt
        - **DNA Markers:** Vora-C1 Catalyst, Bio-Mineralized Base
        """)

    st.divider()

    # --- FULL DETAILED REPORT BELOW ---
    with st.spinner(f"Generating Executive Impact Report..."):
        master_prompt = f"Generate a detailed Strategic Impact Report for: {final_query}. Include a Comparison Table and an Executive Verdict."
        st.markdown(generate_vora_analysis(master_prompt))
else:
    st.info("👆 Please select an asset to see the DNA Transformation.")