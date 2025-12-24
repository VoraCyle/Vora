import streamlit as st
from openai import OpenAI
import json

# --- 1. ACCESS & DATABASE ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except:
    st.error("🚨 API Key Missing.")

# Load the Monster Library
try:
    with open('material_library.json', 'r') as f:
        monster_db = json.load(f)
except:
    monster_db = {}

# --- 2. THE COSTCO REGISTRY ---
vora_100 = {
    "🍗 THE ROTISSERIE CATEGORY": ["Multi-Layer Chicken Bags", "Absorbent Chicken Pad", "Rigid Chicken Trays"],
    "🧻 PAPER & HYGIENE WRAPS": ["Kirkland Bath Tissue Case-Wrap", "Paper Towel Overwrap"],
    "🚩 HIGH-RISK LIABILITIES": ["PVC Clamshells", "PFAS Wrappers", "Black Meat Trays"]
}

# --- 3. THE FORENSIC ENGINE ---
def run_audit(target, volume):
    # Check if it's in our HARD DATA library first
    if target in monster_db:
        item = monster_db[target]
        savings = item['tax_penalty'] * volume
        st.subheader(f"🧬 DNA Transformation: {item['name']}")
        col1, col2 = st.columns(2)
        with col1:
            st.error("🔴 **CURRENT DNA**")
            st.write(f"**Material:** {item['description']}")
            st.write(f"**Annual Liability:** ${savings:,.2f}")
        with col2:
            st.success("🟢 **VORA DNA**")
            st.write(f"**Fix:** {item['vora_fix']['material']}")
            st.write(f"**Annual Savings:** ${savings:,.2f}")
        
        with st.expander("📝 INDUSTRIAL RECIPE (Send to Mixer)"):
            st.code(item['vora_fix']['recipe'], language="text")

    # If it's not in the library, use AI for the full Forensic Analysis
    else:
        with st.spinner(f"Analyzing {target}..."):
            master_prompt = f"Execute a Forensic DNA Audit for: {target}. Compare 'Status Quo' vs 'Vora Resilient Design'. Explain Path A (Waste) and Path B (Recycle)."
            response = client.chat.completions.create(
                model="gpt-4o",
                messages=[{"role": "system", "content": "You are the Vora Lead Forensic Engineer."},
                          {"role": "user", "content": master_prompt}]
            )
            st.markdown(response.choices[0].message.content)

# --- 4. THE INTERFACE ---
st.title("🛡️ VoraCycle: Strategic DNA Command")

# Build Dropdown List
dropdown_items = ["-- Select --"]
for cat, items in vora_100.items():
    dropdown_items.extend(items)

col1, col2 = st.columns(2)
with col1:
    choice = st.selectbox("Select Warehouse SKU:", dropdown_items)
with col2:
    search = st.text_input("OR Type Custom Material:")

volume = st.sidebar.number_input("Annual Unit Volume:", value=1000000)

final_query = search if search else (choice if choice != "-- Select --" else None)

if final_query:
    st.divider()
    run_audit(final_query, volume)