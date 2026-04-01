import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import pandas as pd
from PIL import Image
import io

# Page Configuration
st.set_page_config(
    page_title="Erythromycin Stereochemistry Analysis",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# Custom CSS for Premium Look and Positioning
st.markdown("""
<style>
    .main {
        background-color: #0e1117;
        color: #ffffff;
    }
    .top-right-info {
        position: absolute;
        top: -60px;
        right: 10px;
        text-align: right;
        font-family: 'Inter', sans-serif;
        line-height: 1.2;
        background: rgba(255, 255, 255, 0.05);
        padding: 10px 15px;
        border-radius: 8px;
        border: 1px solid rgba(255, 255, 255, 0.1);
        z-index: 1000;
    }
    .top-right-info p {
        margin: 0;
        font-size: 0.85rem;
        color: #cccccc;
    }
    .top-right-info b {
        color: #ff4b4b;
    }
    .center-title {
        text-align: center;
        margin-top: 20px;
        margin-bottom: 30px;
    }
    .metric-card {
        background-color: #1e2130;
        padding: 20px;
        border-radius: 10px;
        border-left: 5px solid #ff4b4b;
        margin-bottom: 20px;
    }
</style>
""", unsafe_allow_html=True)

# Personal Information (Top Right)
st.markdown("""
<div class="top-right-info">
    <p>Name: <b>Kabish Sridar</b></p>
    <p>Class/Section: <b>AIML - A</b></p>
    <p>Roll No: <b>RA2511026050018</b></p>
</div>
""", unsafe_allow_html=True)

# Middle Top Header
st.markdown("<div class='center-title'><h1>Erythromycin Analysis</h1><p>Molecular Stereochemistry and Chiral Center Mapping</p></div>", unsafe_allow_html=True)

# Fixed Molecule: Erythromycin
erythromycin_smiles = "CC[C@H]1[C@@H](C)[C@@H](O)[C@H](OC(=O)[C@H](C)[C@@H](O)[C@H](C)O)[C@@H](C)[C@@H](O)[C@H](C)O[C@H]2[C@H](O)[C@@H](C)O[C@@H](O[C@@H]3[C@H](C)[C@@H](O)[C@H](C)O[C@H]3O)[C@H](O)[C@H]2O[C@@H]1C"

try:
    mol = Chem.MolFromSmiles(erythromycin_smiles)
    
    if mol:
        # 1. Visualization and Summary
        col1, col2 = st.columns([1.2, 1])
        
        with col1:
            st.subheader("2D Chemical Structure")
            img = Draw.MolToImage(mol, size=(800, 500), kekulize=True)
            st.image(img, use_container_width=True, caption="Structural representation of Erythromycin")

        # 2. Calculation
        mol_with_hs = Chem.AddHs(mol)
        Chem.AssignStereochemistry(mol_with_hs, force=True, cleanIt=True)
        chiral_centers = Chem.FindMolChiralCenters(
            mol_with_hs,
            includeUnassigned=False,
            useLegacyImplementation=False
        )
        
        with col2:
            st.subheader("Molecular Profile")
            st.markdown(f"""
            <div class="metric-card">
                <h2 style='margin:0; color:#ff4b4b;'>{len(chiral_centers)}</h2>
                <p style='margin:0; opacity:0.8;'>Stereocenters (Chiral Centers)</p>
            </div>
            """, unsafe_allow_html=True)
            
            # Descriptors
            mw = Chem.rdMolDescriptors.CalcExactMolWt(mol)
            formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
            
            st.info(f"**Chemical Formula:** {formula}")
            st.info(f"**Exact Mass:** {mw:.4f} u")
            st.success("Analysis complete. Chiral mapping generated below.")

        st.markdown("---")
        
        # 3. Detailed Results
        st.subheader("Chiral Centers Details (Atom Map)")
        if len(chiral_centers) > 0:
            data = []
            for idx, config in chiral_centers:
                atom = mol_with_hs.GetAtomWithIdx(idx)
                element = atom.GetSymbol()
                data.append({
                    "Atom Index": idx,
                    "Atom Symbol": element,
                    "Config (R/S)": config
                })
            
            df = pd.DataFrame(data)
            # Center the table slightly for aesthetic
            c1, c2, c3 = st.columns([1, 4, 1])
            with c2:
                st.dataframe(df, use_container_width=True, hide_index=True)
        else:
            st.warning("No chiral centers identified for this configuration.")

except Exception as e:
    st.error(f"Critical Error in Molecular Engine: {str(e)}")

# Minimal Footer
st.markdown("<br><br>", unsafe_allow_html=True)
st.markdown("---")
st.caption("Chemistry Project Submission | Powered by RDKit Stereochemistry Engine")
