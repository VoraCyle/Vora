import streamlit as st
from openai import OpenAI

# --- 1. ACCESS & CONNECTION ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except Exception as e:
    st.error("🚨 API Key Missing: Please add 'OPENAI_API_KEY' to your Streamlit Secrets.")
    st.stop()

# --- 2. THE ARBITER ENGINE ---
def generate_vora_analysis(prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {"role": "system", "content": "You are the VoraCycle Lead Arbiter. You specialize in Pre-Lifecycle Forensic Design. You analyze product DNA to prevent waste before it begins."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.7 
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Analysis Error: {str(e)}"

# --- 3. USER INTERFACE ---
st.set_page_config(page_title="VoraCycle DNA Arbiter", layout="wide")
st.title("🛡️ VoraCycle: Pre-Lifecycle DNA Arbiter")
st.markdown("#### Solving the waste problem at the design stage—before production even starts.")

query = st.text_input("Enter the Asset or Product Concept for DNA Audit:", 
                     placeholder="e.g., Composite Carbon Fiber Aerospace Panel")

if query:
    st.divider()
    
    with st.spinner("🧬 Mapping Product DNA and simulating environmental outcomes..."):
        # --- 4. THE DESIGN-PHASE FORENSIC PROMPT ---
        master_prompt = f"""
        Act as the VoraCycle Lead Arbiter. Conduct a 'Pre-Production Forensic Audit' for: {query}.
        
        GOAL: Optimize the product DNA before it starts to ensure the highest level of sustainability.

        ### 📊 DUAL-PATH FORENSIC SCORECARD
        Compare the 'Status Quo' (how it's usually made) vs. the 'Vora Optimized' version.
        Include these metrics in a Markdown Table:
        - Primary Outcome (Waste vs Loop)
        - Sustainability Rating (0-10)
        - Financial Value Retention (%)
        - Environmental Risk Level

        ---

        ### 🔍 THE FORENSIC COMPARISON (250+ Words)
        Provide an in-depth analysis of why the 'Vora Optimized Path' is the best strategic choice for the company. 
        Discuss how the Status Quo leads to 'Capital Leakage' and environmental liability. 
        Contrast this with the Optimized Path's ability to maintain high operational velocity and asset recovery. 
        Explain the 'Why' with rigorous, professional logic.

        ### 🧬 THE 'DNA' MODIFICATION (DESIGN FOR ENVIRONMENT)
        If the current end-result of this product leads to waste, how must we change the 
        PRODUCT DNA (materials selection, chemical binders, modularity) right now? 
        Detail exactly how to make it 'Environmentally Safe' so that even if it reaches 
        end-of-life, it acts as a biological nutrient rather than a toxin. 
        Provide 250+ words on non-toxic design and circular engineering.

        ### ⚖️ STRATEGIC REASONING & 5-YEAR OUTCOME
        Provide a 150-word final justification. Why is this path the absolute best for the 
        company's long-term survival? How does this decision create a 'Circular Moat' 
        around the business?
        """

        # --- 5. EXECUTION & DISPLAY ---
        analysis_result = generate_vora_analysis(master_prompt)
        st.markdown(analysis_result)

# --- 6. FOOTER ---
st.sidebar.info("VoraCycle DNA Logic: v3.0.0 (Pre-Production Mode)")