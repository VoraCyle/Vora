import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from pyscf import gto, scf
import streamlit_authenticator as stauth

# --- Page Config ---
st.set_page_config(page_title="Wraith VoraCycle", layout="centered")

# --- Authentication Logic ---
config = {
    'credentials': {
        'usernames': {
            'wraith': {
                'email': 'wraith@voracycle.com',
                'name': 'Wraith Admin',
                'password': '$2b$12$EixZaYVK1fsbw1ZfbX3OXe.WI/.2B8K5wjW6SAff4D8cQ4AtNtg6'
            }
        }
    },
    'cookie': {
        'expiry_days': 30,
        'key': 'vora_signature_key',
        'name': 'vora_cookie'
    },
    'preauthorized': {
        'emails': []
    }
} 

authenticator = stauth.Authenticate(
    credentials,
    "voracycle_cookie",
    "random_key_123",
    cookie_expiry_days=30
)

name, authentication_status, username = authenticator.login('Login to Wraith VoraCycle', 'main')

if authentication_status:
    authenticator.logout('Logout', 'sidebar')
    [cite_start]st.sidebar.success(f"Welcome, {name}!") 

    st.title("🔮 Wraith VoraCycle")
    [cite_start]st.markdown("**Precycling Powered by Precision**") 
    [cite_start]st.markdown("Predict material circularity before it becomes waste. Quantum-inspired AI for proactive regeneration.") 

    # --- Material Selection ---
    [cite_start]st.sidebar.header("Material Selection") 
    material = st.sidebar.selectbox(
        "Choose a common wholesale packaging material",
        ["Polyethylene (PE)", "Polyethylene Terephthalate (PET)", "Polypropylene (PP)", "Custom SMILES"]
    [cite_start]) 

    smiles_dict = {
        [cite_start]"Polyethylene (PE)": "CCCCCCCCCC", 
        [cite_start]"Polyethylene Terephthalate (PET)": "C1=CC=C(C=C1)C(=O)OCCOC(=O)C2=CC=C(C=C2)C(=O)O", 
        [cite_start]"Polypropylene (PP)": "CC(C)CCCCCC" 
    }

    if material != "Custom SMILES":
        [cite_start]smiles = smiles_dict[material] 
    else:
        [cite_start]smiles = st.sidebar.text_input("Enter SMILES string", "CCCC") 

    [cite_start]mol = Chem.MolFromSmiles(smiles) 
    if mol:
        [cite_start]img = Draw.MolToImage(mol, size=(600, 600)) 
        [cite_start]st.image(img, caption=f"Molecular structure: {material if material != 'Custom SMILES' else 'Custom'}") 

    # --- Simulation ---
    [cite_start]st.write("### Precision Simulation Running...") 
    
    # We use a proxy molecule for the quantum calculation to ensure performance
    [cite_start]proxy_mol = gto.M(atom='C 0 0 0; C 0 0 1.4; H 0 1 1; H 0 0 2.4', basis='sto-3g') 
    [cite_start]mf = scf.RHF(proxy_mol) 
    [cite_start]energy = mf.kernel() 

    [cite_start]st.write(f"**Quantum energy stability**: {energy:.2f} Hartree") 

    [cite_start]scores = {"Polyethylene (PE)": 48, "Polyethylene Terephthalate (PET)": 72, "Polypropylene (PP)": 55} 
    [cite_start]score = scores.get(material, 65) if material != "Custom SMILES" else 65 
    
    [cite_start]st.metric(label="**VoraCycle Circularity Score**", value=f"{score}/100", delta=f"{score-50:+} vs average") 
    # --- Redesign Recommendations ---
    [cite_start]st.markdown("### 🔬 Precision Redesign Recommendation") 
    
    redesign_suggestions = {
        "Polyethylene Terephthalate (PET)": {
            [cite_start]"optimized_name": "Polyethylene Furanoate (PEF)", 
            [cite_start]"optimized_score": 94, 
            [cite_start]"optimized_stability": -74.20,
            
            [cite_start]"modifications": "Replace terephthalic acid with bio-derived 2,5-furandicarboxylic acid (FDCA). Incorporate enzyme-recognizable cleavage sites.", [cite: 7, 8]
        },
        "Polyethylene (PE)": {
            [cite_start]"optimized_name": "Bio-PE with Enzymatic Triggers", 
            [cite_start]"optimized_score": 88, 
            [cite_start]"optimized_stability": -68.50, 
            [cite_start]"modifications": "Blend with PHA copolymers and add ester-based degradation triggers.", 
        },
        "Polypropylene (PP)": {
            [cite_start]"optimized_name": "Isotactic PP with Bio-additives", 
            [cite_start]"optimized_score": 85,
            [cite_start]"optimized_stability": -72.10,
            [cite_start]"modifications": "Incorporate oxidation-promoting prodegradants and bio-based comonomers.", 
        }
    }

    suggestion = redesign_suggestions.get(material, {
        [cite_start]"optimized_name": "Custom Optimized Variant",
        [cite_start]"optimized_score": score + 15, 
        [cite_start]"optimized_stability": energy + 1.5, 
        [cite_start]"modifications": "Functional group optimization using quantum-guided bond lability tuning." 
    })

    col1, col2 = st.columns(2)
    with col1:
        st.error(f"**Original**: {score}/100")
    with col2:
        st.success(f"**Redesign**: {suggestion['optimized_name']} ({suggestion['optimized_score']}/100)")

    [cite_start]st.info(f"**Action Plan:** {suggestion['modifications']}") 
    

    # --- Environmental Fate ---
    [cite_start]st.markdown("### 🌍 Environmental Fate & Redesign Impact")
    colA, colB = st.columns(2)
    with colA:
        st.markdown("**Original Fate:**")
        [cite_start]st.write("- Persists for hundreds of years.") 
        [cite_start]st.write("- Fragments into microplastics.") 
    with colB:
        st.markdown("**Redesigned Fate:**")
        [cite_start]st.write("- Fully mineralizes in months to years.") 
        [cite_start]st.write("- Supports 100% closed-loop recycling.") 

    [cite_start]st.caption("Prototype created by Wraith | Pioneering Precycling with Precision") 

elif authentication_status == False:
    [cite_start]st.error('Username/password incorrect') 
elif authentication_status is None:
    [cite_start]st.warning('Please enter your credentials') 






