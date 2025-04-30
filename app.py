import streamlit as st
import pubchempy as pcp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from io import BytesIO
from fpdf import FPDF
from datetime import datetime
import plotly.graph_objects as go
from scipy.constants import epsilon_0, elementary_charge, Boltzmann
import firebase_admin
from firebase_admin import credentials, firestore
import json
import os

# ----------------------------
# 1. CONFIGURAÇÕES INICIAIS
# ----------------------------
st.set_page_config(
    page_title="Simulador Profissional de Eletroforese", 
    layout="wide",
    page_icon="🔬"
)

# ----------------------------
# 2. CLASSE DE SIMULAÇÃO CIENTÍFICA
# ----------------------------
class NernstPlanckSimulator:
    def __init__(self):
        self.constantes = {
            'epsilon_0': 8.854e-12,  # F/m
            'k_B': Boltzmann,         # J/K
            'e': elementary_charge    # C
        }
    
    def calcular_mobilidade(self, carga, raio, viscosidade, permissividade_rel, temperatura, forca_ionica):
        """Modelo completo com correções de Debye-Hückel e força iônica"""
        eta = viscosidade * 1e-3  # cP para Pa·s
        epsilon = permissividade_rel * self.constantes['epsilon_0']
        
        # Termo de Stokes-Einstein
        termo_stokes = (carga * self.constantes['e']) / (6 * np.pi * eta * raio)
        
        # Correção de Debye-Hückel
        raio_debye = np.sqrt(epsilon * self.constantes['k_B'] * temperatura / 
                          (self.constantes['e']**2 * forca_ionica * 1e3 * 6.022e23))
        fator_debye = 1 / (1 + raio/raio_debye)
        
        # Efeito da temperatura
        fator_temp = np.exp(298.15/temperatura - 1)
        
        return termo_stokes * fator_debye * fator_temp * 1e8  # Unidades práticas

# ----------------------------
# 3. CONFIGURAÇÃO DO FIREBASE (OPCIONAL)
# ----------------------------
def init_firebase():
    if not firebase_admin._apps:
        # Substitua pelo seu config ou arquivo JSON
        FIREBASE_CONFIG = {
            "type": "service_account",
            # ... (seus dados de configuração)
        }
        cred = credentials.Certificate(FIREBASE_CONFIG)
        firebase_admin.initialize_app(cred)
    return firestore.client()

try:
    db = init_firebase()
    firebase_ready = True
except:
    firebase_ready = False
    st.warning("Modo offline ativado (Firebase não configurado)")

# ----------------------------
# 4. FUNÇÕES PRINCIPAIS
# ----------------------------
def buscar_molecula(nome):
    try:
        mol = pcp.get_compounds(nome, 'name')[0]
        return {
            "Nome": mol.iupac_name,
            "Fórmula": mol.molecular_formula,
            "Peso Molecular": mol.molecular_weight,
            "SMILES": mol.canonical_smiles,
            "CID": mol.cid,
            "Data": datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        }
    except Exception as e:
        st.error(f"Erro na busca: {str(e)}")
        return None

def gerar_relatorio_pdf(params, resultados):
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", size=12)
    
    # Cabeçalho
    pdf.cell(0, 10, "Relatório Científico de Eletroforese", ln=True, align='C')
    pdf.ln(10)
    
    # Parâmetros
    pdf.cell(0, 10, f"Data: {datetime.now().strftime('%d/%m/%Y %H:%M')}", ln=True)
    for key, value in params.items():
        pdf.cell(0, 10, f"{key}: {value}", ln=True)
    
    # Gráfico
    pdf.ln(10)
    pdf.cell(0, 10, "Cromatograma Simulado:", ln=True)
    pdf.image("temp_plot.png", x=10, w=180)
    
    return pdf

# ----------------------------
# 5. INTERFACE STREAMLIT
# ----------------------------
def main():
    st.title("🔬 Simulador Profissional de Eletroforese Capilar")
    
    # Barra lateral
    with st.sidebar:
        st.header("Configurações Globais")
        modo_avancado = st.checkbox("Modo Avançado", True)
        if firebase_ready:
            st.success("Firebase Conectado")
        else:
            st.warning("Modo Local Ativo")
    
    # Abas principais
    tab1, tab2, tab3 = st.tabs(["Busca Molecular", "Simulação", "Banco de Dados"])
    
    with tab1:
        st.header("🔍 Busca no PubChem")
        nome_molecula = st.text_input("Digite o nome da molécula (em inglês):", "aspirin")
        
        if nome_molecula:
            with st.spinner("Buscando no PubChem..."):
                dados = buscar_molecula(nome_molecula)
            
            if dados:
                col1, col2 = st.columns(2)
                with col1:
                    st.json(dados)
                    if st.button("Salvar no Banco"):
                        if firebase_ready:
                            db.collection("moleculas").document(str(dados["CID"])).set(dados)
                            st.success("Salvo no Firebase!")
                with col2:
                    st.image(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{dados['CID']}/PNG", 
                           caption=f"Estrutura 2D de {dados['Nome']}")
    
    with tab2:
        st.header("⚡ Simulação de Eletroforese")
        simulator = NernstPlanckSimulator()
        
        with st.expander("Parâmetros da Simulação", expanded=True):
            col1, col2 = st.columns(2)
            
            with col1:
                voltagem = st.slider("Voltagem (kV)", 5, 30, 15)
                comprimento = st.slider("Comprimento do Capilar (cm)", 10, 100, 50)
                temperatura = st.slider("Temperatura (°C)", 15, 40, 25)
                
            with col2:
                pH = st.slider("pH do Tampão", 2.0, 10.0, 7.0, 0.1)
                viscosidade = st.slider("Viscosidade (cP)", 0.8, 2.5, 1.0, 0.1)
                forca_ionica = st.slider("Força Iônica (mM)", 10, 200, 50)
        
        if st.button("Executar Simulação Completa"):
            # Cálculos científicos
            raio_estimado = 1e-9  # 1 nm como aproximação
            mobilidade = simulator.calcular_mobilidade(
                carga=2, 
                raio=raio_estimado,
                viscosidade=viscosidade,
                permissividade_rel=78.5,
                temperatura=temperatura + 273.15,
                forca_ionica=forca_ionica
            )
            
            # Simulação do cromatograma
            tempo_migracao = (comprimento * 1e-2) / (mobilidade * 1e-8 * voltagem * 1e3)
            t = np.linspace(0, tempo_migracao * 1.5, 1000)
            y = np.exp(-(t - tempo_migracao)**2 / (2 * (tempo_migracao/10)**2)) * 100
            
            # Visualização 3D
            fig_3d = go.Figure(data=[
                go.Scatter3d(
                    x=t,
                    y=[mobilidade] * len(t),
                    z=y,
                    mode='lines',
                    line=dict(width=8, color='purple')
                )
            ])
            fig_3d.update_layout(
                scene=dict(
                    xaxis_title='Tempo (s)',
                    yaxis_title='Mobilidade (10^-8 m²/Vs)',
                    zaxis_title='Intensidade'
                ),
                title='Perfil 3D da Separação'
            )
            
            # Mostrar resultados
            col1, col2 = st.columns(2)
            with col1:
                st.plotly_chart(fig_3d, use_container_width=True)
            with col2:
                st.metric("Tempo de Migração", f"{tempo_migracao:.2f} s")
                st.metric("Mobilidade Calculada", f"{mobilidade:.2e} x10⁻⁸ m²/Vs")
                
                # Gerar PDF
                pdf = gerar_relatorio_pdf(
                    params={
                        "Voltagem": f"{voltagem} kV",
                        "Comprimento": f"{comprimento} cm",
                        "Temperatura": f"{temperatura} °C",
                        "pH": pH,
                        "Força Iônica": f"{forca_ionica} mM"
                    },
                    resultados={
                        "Tempo Migração": tempo_migracao,
                        "Mobilidade": mobilidade
                    }
                )
                
                st.download_button(
                    label="📥 Baixar Relatório Completo",
                    data=pdf.output(dest='S').encode('latin1'),
                    file_name="relatorio_eletroforese.pdf",
                    mime="application/pdf"
                )
    
    with tab3:
        st.header("📊 Banco de Dados")
        if firebase_ready:
            docs = db.collection("moleculas").stream()
            dados = [doc.to_dict() for doc in docs]
            st.dataframe(pd.DataFrame(dados))
        else:
            st.warning("Conecte-se ao Firebase para ativar esta funcionalidade")

if __name__ == "__main__":
    main()
