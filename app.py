import streamlit as st
from openai import OpenAI

# --- 1. ACCESS & CONNECTION ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except Exception as e:
    st.error("🚨 API Key Missing: Please add 'OPENAI_API_KEY' to your Streamlit Secrets.")
    st.stop()

# --- 2. THE ASSET REGISTRY ---
problematic_items = [
    "-- Select a High-Volume Asset --",
    "Multi-Layer Pet Food Bags (Mylar/Poly)",
    "Black Plastic Food Trays (Carbon Black)",
    "LLDPE Pallet Shrink Wrap",
    "Expanded Polystyrene (EPS) Coolers",
    "Lithium-Ion Tool Batteries (Returns)",
    "Wax-Coated Produce Corrugated Boxes",
    "PVC Blister Packaging",
    "Composite Flooring Returns",
    "Textile Waste (Kirkland Clothing)",
    "Industrial Pallet Strapping (Nylon/PP)"
]

# --- 3. THE ARBITER ENGINE ---
def generate_vora_analysis(prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {"role": "system", "content": "You are the VoraCycle Lead Arbiter. You provide technical engineering blueprints and supplier grading for circular retail systems."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.7 
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Analysis Error: {str(e)}"

# --- 4. USER INTERFACE ---
st.set_page_config(page_title="VoraCycle: Enterprise DNA", layout="wide")
st.title("🛡️ VoraCycle: Strategic DNA Command Center")
st.markdown("### Solving the End-Game at the Design Stage.")

col1, col2 = st.columns(2)
with col1:
    dropdown_choice = st.selectbox("Select problematic Warehouse Stream:", problematic_items)
with col2:
    search_query = st.text_input("Custom SKU/Material Audit:")

final_query = search_query if search_query else (dropdown_choice if dropdown_choice != "-- Select a High-Volume Asset --" else None)

if final_query:
    st.divider()
    with st.spinner("Engineering Universal Sustainability Blueprint..."):
        
        master_prompt = f"""
        Audit the following for an Enterprise Retailer: {final_query}.

        ### 📊 DUAL-PATH FORENSIC SCORECARD
        Compare 'Status Quo' vs. 'Vora Resilient Design'.
        Include: Primary Outcome, Sustainability Rating, Capital Retention %, and Environmental Toxin Risk.

        ---

        ### 🧬 DNA ENGINEERING PROTOCOL: THE "WHAT & HOW"
        Provide a technical blueprint for the design team:
        1. **WHAT TO REMOVE:** List specific toxic binders, mixed-polymers, or non-separable coatings.
        2. **WHAT TO REPLACE WITH:** Suggest the specific mono-material or bio-mineral alternatives.
        3. **HOW TO ASSEMBLE:** Describe the 'Design for Disassembly' or 'Molecular Purity' process needed so that no matter where the item goes, it stays at the highest level of sustainability. (300+ Words).

        ---

        ### 📋 VORA SUPPLIER SCORECARD
        Evaluate a theoretical supplier providing this item based on the new DNA.
        Format as a Markdown Checklist/Table:
        - Material Purity (0-10)
        - Chemical Transparency (0-10)
        - Circular Recovery Readiness (0-10)
        - Overall Vora Grade (A-F)

        ---

        ### 🏁 FINAL EXECUTIVE SUMMARY: THE SUSTAINABILITY LEADERSHIP VERDICT
        How does this redesign transform a waste liability into a circular asset? Why is this the most beneficial path for the company's 10-year survival? (200 Words).
        """

        st.markdown(generate_vora_analysis(master_prompt))

# --- 5. FOOTER ---
st.sidebar.info("VoraCycle Enterprise Logic: v4.7.0")