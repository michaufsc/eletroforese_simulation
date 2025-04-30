import streamlit as st
import pubchempy as pcp
import pandas as pd
import numpy as np
from io import BytesIO
from fpdf import FPDF
from datetime import datetime
from scipy.constants import epsilon_0, elementary_charge, Boltzmann
import matplotlib.pyplot as plt
import os

# Verificação e instalação automática do Plotly (opcional)
try:
    import plotly.graph_objects as go
    PLOTLY_ENABLED = True
except ImportError:
    PLOTLY_ENABLED = False
    st.warning("Para visualizações 3D interativas, instale Plotly: `pip install plotly`")

# Configuração da página
st.set_page_config(
    page_title="Simulador Avançado de Eletroforese",
    layout="wide",
    page_icon="🔬"
)

# Classe de simulação científica
class AdvancedElectrophoresisSimulator:
    def __init__(self):
        self.constants = {
            'epsilon_0': epsilon_0,
            'e': elementary_charge,
            'k_B': Boltzmann
        }
    
    def calculate_mobility(self, charge, radius, viscosity, permittivity, temp, ionic_strength):
        """Cálculo preciso da mobilidade eletroforética"""
        eta = viscosity * 1e-3  # cP to Pa·s
        epsilon = permittivity * self.constants['epsilon_0']
        
        # Termo principal
        mobility = (charge * self.constants['e']) / (6 * np.pi * eta * radius)
        
        # Correção de Debye-Hückel
        debye_length = np.sqrt(epsilon * self.constants['k_B'] * temp / 
                             (self.constants['e']**2 * ionic_strength * 1e3 * 6.022e23))
        correction = 1 / (1 + radius/debye_length)
        
        return mobility * correction * 1e8  # Convert to practical units

# Interface principal
def main():
    st.title("🔬 Simulador Profissional de Eletroforese Capilar")
    
    # Inicializar simulador
    simulator = AdvancedElectrophoresisSimulator()
    
    # Controles na sidebar
    with st.sidebar:
        st.header("Configurações")
        advanced_mode = st.checkbox("Modo Avançado", True)
        
        if not PLOTLY_ENABLED:
            st.error("Recursos 3D desativados (Plotly não instalado)")

    # Abas principais
    tab1, tab2 = st.tabs(["Busca Molecular", "Simulação"])

    with tab1:
        st.header("🔍 Consulta ao PubChem")
        compound = st.text_input("Nome da molécula:", "caffeine")
        
        if compound:
            with st.spinner("Buscando no PubChem..."):
                try:
                    result = pcp.get_compounds(compound, 'name')[0]
                    data = {
                        "Nome": result.iupac_name,
                        "Fórmula": result.molecular_formula,
                        "Peso Molecular": result.molecular_weight,
                        "CID": result.cid
                    }
                    
                    col1, col2 = st.columns(2)
                    with col1:
                        st.json(data)
                    with col2:
                        st.image(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{result.cid}/PNG",
                                caption=f"Estrutura de {result.iupac_name}")
                
                except Exception as e:
                    st.error(f"Erro na busca: {str(e)}")

    with tab2:
        st.header("⚡ Simulação de Eletroforese")
        
        with st.expander("Parâmetros de Controle", expanded=True):
            col1, col2 = st.columns(2)
            
            with col1:
                voltage = st.slider("Voltagem (kV)", 5, 30, 15)
                length = st.slider("Comprimento do capilar (cm)", 10, 100, 50)
                temp = st.slider("Temperatura (°C)", 15, 40, 25)
            
            with col2:
                ph = st.slider("pH do tampão", 2.0, 10.0, 7.0, 0.1)
                viscosity = st.slider("Viscosidade (cP)", 0.8, 2.5, 1.0, 0.1)
                ionic_str = st.slider("Força iônica (mM)", 10, 200, 50)
        
        # Cálculos científicos
        if st.button("Executar Simulação", type="primary"):
            with st.spinner("Calculando..."):
                try:
                    # Parâmetros de exemplo (substitua pelos seus valores reais)
                    charge = -1 if ph > 7 else 1
                    radius = 1e-9  # 1 nm
                    
                    mobility = simulator.calculate_mobility(
                        charge=charge,
                        radius=radius,
                        viscosity=viscosity,
                        permittivity=78.5,
                        temp=temp + 273.15,
                        ionic_strength=ionic_str
                    )
                    
                    migration_time = (length * 1e-2) / (mobility * 1e-8 * voltage * 1e3)
                    
                    # Gerar dados do cromatograma
                    t = np.linspace(0, migration_time * 2, 1000)
                    signal = np.exp(-(t - migration_time)**2 / (0.1 * migration_time**2)) * 100
                    
                    # Visualização
                    st.success(f"Tempo de migração: {migration_time:.2f} s | Mobilidade: {mobility:.2e} m²/Vs")
                    
                    if PLOTLY_ENABLED:
                        # Gráfico 3D interativo
                        fig = go.Figure(data=[
                            go.Scatter3d(
                                x=t,
                                y=[mobility] * len(t),
                                z=signal,
                                mode='lines',
                                line=dict(width=8, color='#FF6B6B')
                        ])
                        fig.update_layout(
                            scene=dict(
                                xaxis_title='Tempo (s)',
                                yaxis_title='Mobilidade',
                                zaxis_title='Intensidade'
                            ),
                            title='Perfil 3D da Separação',
                            margin=dict(l=0, r=0, b=0, t=30)
                        )
                        st.plotly_chart(fig, use_container_width=True)
                    else:
                        # Fallback para matplotlib
                        fig, ax = plt.subplots()
                        ax.plot(t, signal, color='purple')
                        ax.set_xlabel('Tempo (s)')
                        ax.set_ylabel('Intensidade')
                        st.pyplot(fig)
                
                except Exception as e:
                    st.error(f"Erro na simulação: {str(e)}")

if __name__ == "__main__":
    main()
