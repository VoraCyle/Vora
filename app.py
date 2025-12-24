import streamlit as st
from openai import OpenAI

# --- 1. ACCESS & CONNECTION ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except Exception as e:
    st.error("🚨 API Key Missing: Please add 'OPENAI_API_KEY' to your Streamlit Secrets.")
    st.stop()

# --- 2. COSTCO-SCALE "HARD" ASSETS ---
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
                {"role": "system", "content": "You are the VoraCycle Lead Arbiter. You specialize in Pre-Lifecycle Forensic Design—solving the end of a product's life before it is even manufactured."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.7 
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Analysis Error: {str(e)}"

# --- 4. USER INTERFACE ---
st.set_page_config(page_title="VoraCycle: DNA Forensic Command", layout="wide")
st.title("🛡️ VoraCycle: Strategic DNA Command Center")
st.markdown("#### Engineering the End-Result at the Start: Forensic DNA Optimization")

col1, col2 = st.columns(2)
with col1:
    dropdown_choice = st.selectbox("Select problematic Warehouse Stream:", problematic_items)
with col2:
    search_query = st.text_input("Custom SKU/Material Audit:")

final_query = search_query if search_query else (dropdown_choice if dropdown_choice != "-- Select a High-Volume Asset --" else None)

if final_query:
    st.divider()
    with st.spinner("Simulating Pre-Lifecycle DNA Outcomes..."):
        
        master_prompt = f"""
        Audit the following for an Enterprise Retailer: {final_query}.

        ### 📊 DUAL-PATH FORENSIC SCORECARD
        Compare the 'Status Quo' (Default Design) vs. 'Vora Optimized' (Engineered for the End).
        Output this as a Markdown Table:
        - Primary Outcome (e.g. Landfill vs Infinite Loop)
        - Sustainability Rating (0-10)
        - Capital Retention (%)
        - Environmental Impact Score (0-10)
        - Risk/Liability Level (High/Low)

        ---

        ### 🧬 PRE-LIFECYCLE DNA: DESIGNING THE END AT THE START
        **The "Before" DNA (Current State):** Describe the toxic or non-circular elements currently in this product.
        **The "After" DNA (Vora Optimized):** Describe exactly how to re-engineer this product NOW (materials, assembly, chemicals) to ensure the best end-result. 
        
        *Focus on how changing the DNA at the start makes the final route (Refurbish or Mineralize) 100% efficient and safe.* (250+ Words).

        ### ⚖️ THE WINNING VERDICT: BEST STRATEGIC PATH
        Identify the best end-route. Provide 200 words of forensic logic explaining why this path is the most beneficial for the environment and the business P&L.

        ### 🚚 HANDLING & LOGISTICS (SOP)
        1. **Warehouse Handling:** Precise sorting/prepping instructions.
        2. **Logistics:** Transport path details and 'Backhauling' opportunities.
        3. **Facility:** Specify the target Circular Hub or Processing center.

        ---

        ### 🏁 FINAL EXECUTIVE SUMMARY: THE VORA ADVANTAGE
        Provide a 200-word concluding argument. 
        Explain how this "Start-at-the-End" strategy keeps sustainability high and eliminates the concept of waste entirely by fixing the product DNA before it starts.
        """

        st.markdown(generate_vora_analysis(master_prompt))

# --- 5. FOOTER ---
st.sidebar.info("VoraCycle Enterprise Logic: v4.5.0")
