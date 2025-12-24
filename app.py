import streamlit as st
from openai import OpenAI
import json

# --- 1. SETUP & DATABASE ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except:
    st.error("🚨 API Key Missing in Streamlit Secrets.")

# Load the Monster Library
try:
    with open('material_library.json', 'r') as f:
        monster_db = json.load(f)
except:
    monster_db = {}

# --- 2. THE DNA ANALYZER ENGINE ---
def run_dna_analysis(query, volume=1000000):
    # Check if the query is a known 'Monster' from our Database
    if query in monster_db:
        item = monster_db[query]
        savings = item['tax_penalty'] * volume
        
        st.subheader(f"🧬 DNA Transformation: {item['name']}")
        col1, col2 = st.columns(2)
        with col1:
            st.error("🔴 **CURRENT DNA (Status Quo)**")
            st.write(f"**Material:** {item['description']}")
            st.write(f"**Annual Tax Liability:** ${savings:,.2f}")
        with col2:
            st.success("🟢 **VORA DNA (The Cure)**")
            st.write(f"**New Solution:** {item['vora_fix']['material']}")
            st.write(f"**Circular Path:** 100% Recyclable / Mono-material")
        
        with st.expander("📝 INDUSTRIAL DNA RECIPE (Send to Mixer)"):
            st.code(item['vora_fix']['recipe'], language="text")
            st.caption("Standardized instructions for third-party pellet blenders.")

    # Otherwise, use the AI for a Custom Forensic Audit
    else:
        with st.spinner("Executing Forensic DNA Reconstruction..."):
            prompt = f"Analyze the material DNA for: {query}. Explain how to transform it from a high-tax composite to a Vora-resilient mono-material."
            response = client.chat.completions.create(
                model="gpt-4o",
                messages=[{"role": "system", "content": "You are the Vora Lead Forensic Engineer."},
                          {"role": "user", "content": prompt}]
            )
            st.markdown(response.choices[0].message.content)

# --- 3. THE INTERFACE ---
st.title("🛡️ Vora: Strategic DNA Command Center")

# Search and Dropdown logic
col_a, col_b = st.columns([1, 2])
with col_a:
    selected_id = st.selectbox("Pick a known Monster:", ["-- Select --"] + list(monster_db.keys()))
with col_b:
    custom_search = st.text_input("OR Type a Custom SKU/Material:")

volume = st.sidebar.number_input("Annual Unit Volume:", value=1000000)

# Execute based on user choice
final_target = custom_search if custom_search else (selected_id if selected_id != "-- Select --" else None)

if final_target:
    st.divider()
    run_dna_analysis(final_target, volume)