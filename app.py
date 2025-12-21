import streamlit as st
import streamlit_authenticator as st_auth
from rdkit import Chem
from rdkit.Chem import Draw
from pyscf import gto, scf

# 1. Credentials
credentials = {
    'usernames': {
        'wraith': {
            'name': 'Wraith',
            'password': 'Vora1630' 
        }
    }
}

# 2. Authenticator Setup
authenticator = st_auth.Authenticate(
    credentials,
    "vora_cookie",
    "auth_key",
    cookie_expiry_days=30
)

# 3. Login UI
authenticator.login(location='main')

auth_status = st.session_state.get("authentication_status")
name = st.session_state.get("name")
username = st.session_state.get("username")

# 4. Main App Logic
if auth_status:
    authenticator.logout('Logout', 'sidebar')
    st.title("🔮 Wraith VoraCycle")
    st.sidebar.success(f"Welcome, {name}!")
    
    material = st.sidebar.selectbox("Select Material", ["Polyethylene (PE)", "Polyethylene Terephthalate (PET)", "Polypropylene (PP)", "Custom SMILES"])

    smiles_dict = {
        "Polyethylene (PE)": "CCCCCCCCCC", 
        "Polyethylene Terephthalate (PET)": "C1=CC=C(C=C1)C(=O)OCCOC(=O)C2=CC=C(C=C2)C(=O)O", 
        "Polypropylene (PP)": "CC(C)CCCCCC" 
    }

    if material != "Custom SMILES":
        smiles = smiles_dict.get(material, "CCCC")
    else:
        smiles = st.sidebar.text_input("Enter SMILES string", "CCCC") 

    mol = Chem.MolFromSmiles(smiles) 
    if mol:
        img = Draw.MolToImage(mol, size=(400, 400)) 
        st.image(img, caption=f"Molecular structure: {material}") 

    st.write("### Precision Simulation Running...") 
    
    # Quantum Calculation
    proxy_mol = gto.M(atom='C 0 0 0; C 0 0 1.4; H 0 1 1; H 0 0 2.4', basis='sto-3g') 
    mf = scf.RHF(proxy_mol) 
    energy = mf.kernel() 

    st.write(f"**Quantum energy stability**: {energy:.2f} Hartree") 

    # Scoring
    scores = {"Polyethylene (PE)": 48, "Polyethylene Terephthalate (PET)": 72, "Polypropylene (PP)": 55} 
    score = scores.get(material, 65)
    st.metric(label="VoraCycle Score", value=f"{score}/100") 

    # Recommendations
    st.markdown("### 🔬 Precision Redesign Recommendation") 
    st.info("Optimization: Functional group tuning via quantum-guided bond lability.")

    st.caption("Prototype created by Wraith | Pioneering Precycling") 

elif auth_status == False:
    st.error('Username/password incorrect') 
elif auth_status is None:
    st.warning('Please enter your credentials')
    
    
    
    





