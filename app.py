import streamlit as st
import pubchempy as pcp
import pandas as pd
import numpy as np
from io import BytesIO
from fpdf import FPDF
from datetime import datetime
import os

# Configuração da página
st.set_page_config(
    page_title="Simulador de Eletroforese Capilar",
    layout="wide",
    page_icon="🔬"
)

# Verificação de dependências gráficas
try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_ENABLED = True
except ImportError:
    MATPLOTLIB_ENABLED = False
    st.warning("Visualizações avançadas desativadas (Matplotlib não instalado)")

# Classe de simulação científica
class EletroforeseSimulator:
    def __init__(self):
        self.constantes = {
            'fator_conversao': 1e5  # Fator para unidades práticas
        }
    
    def calcular_mobilidade(self, carga, massa, pH, temperatura=25):
        """Calcula a mobilidade eletroforética simplificada"""
        # Correção para temperatura (fictício para exemplo)
        fator_temp = 1 + 0.02 * (temperatura - 25)
        
        # Cálculo básico da mobilidade
        return (carga / massa) * (1 + (pH - 7) * 0.1) * self.constantes['fator_conversao'] * fator_temp

# Interface principal
def main():
    st.title("🔬 Simulador de Eletroforese Capilar")
    
    # Inicializar simulador
    simulator = EletroforeseSimulator()
    
    # Controles na sidebar
    with st.sidebar:
        st.header("Configurações")
        modo_avancado = st.checkbox("Usar parâmetros avançados", False)
    
    # Abas principais
    tab1, tab2 = st.tabs(["Busca Molecular", "Simulação"])

    with tab1:
        st.header("🔍 Consulta ao PubChem")
        nome_molecula = st.text_input("Nome da molécula (em inglês):", "aspirin")
        
        if nome_molecula:
            with st.spinner("Buscando no PubChem..."):
                try:
                    resultado = pcp.get_compounds(nome_molecula, 'name')[0]
                    dados = {
                        "Nome": resultado.iupac_name,
                        "Fórmula": resultado.molecular_formula,
                        "Peso Molecular": resultado.molecular_weight,
                        "CID": resultado.cid
                    }
                    
                    col1, col2 = st.columns(2)
                    with col1:
                        st.json(dados)
                    with col2:
                        st.image(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{resultado.cid}/PNG",
                                caption=f"Estrutura de {resultado.iupac_name}")
                
                except Exception as e:
                    st.error(f"Erro na busca: {str(e)}")

    with tab2:
        st.header("⚡ Simulação de Eletroforese")
        
        # Parâmetros básicos
        col1, col2 = st.columns(2)
        with col1:
            voltagem = st.slider("Voltagem (kV)", 5, 30, 15)
            comprimento = st.slider("Comprimento do capilar (cm)", 10, 100, 50)
        with col2:
            pH = st.slider("pH do tampão", 2.0, 10.0, 7.0, 0.1)
            temperatura = st.slider("Temperatura (°C)", 15, 40, 25)
        
        # Parâmetros avançados
        if modo_avancado:
            with st.expander("Parâmetros Avançados"):
                viscosidade = st.slider("Viscosidade (cP)", 0.8, 2.5, 1.0, 0.1)
                forca_ionica = st.slider("Força iônica (mM)", 10, 200, 50)
        
        # Controle de simulação
        if st.button("Executar Simulação", type="primary"):
            with st.spinner("Calculando..."):
                try:
                    # Exemplo com valores padrão
                    massa = 180.16  # Massa molecular da aspirina (exemplo)
                    carga = -1 if pH > 7 else 1  # Carga simplificada
                    
                    mobilidade = simulator.calcular_mobilidade(
                        carga=carga,
                        massa=massa,
                        pH=pH,
                        temperatura=temperatura
                    )
                    
                    tempo_migracao = (comprimento * 1e-2) / (mobilidade * voltagem * 1e3)
                    
                    # Gerar dados do cromatograma
                    t = np.linspace(0, tempo_migracao * 2, 100)
                    sinal = np.exp(-(t - tempo_migracao)**2 / (0.2 * tempo_migracao**2)) * 100
                    
                    # Visualização
                    st.success(f"Tempo de migração estimado: {tempo_migracao:.2f} segundos")
                    
                    if MATPLOTLIB_ENABLED:
                        fig, ax = plt.subplots()
                        ax.plot(t, sinal, color='purple', linewidth=2)
                        ax.set_xlabel('Tempo (s)')
                        ax.set_ylabel('Intensidade')
                        ax.grid(True, linestyle='--', alpha=0.6)
                        st.pyplot(fig)
                    else:
                        st.line_chart(
                            pd.DataFrame({
                                'Tempo': t,
                                'Intensidade': sinal
                            }).set_index('Tempo')
                        )
                    
                    # Relatório em PDF
                    pdf = FPDF()
                    pdf.add_page()
                    pdf.set_font("Arial", size=12)
                    pdf.cell(0, 10, "Relatório de Simulação", ln=True, align='C')
                    pdf.cell(0, 10, f"Molécula: {nome_molecula}", ln=True)
                    pdf.cell(0, 10, f"Tempo de migração: {tempo_migracao:.2f} s", ln=True)
                    
                    img_path = "temp_plot.png"
                    if MATPLOTLIB_ENABLED:
                        fig.savefig(img_path, bbox_inches='tight')
                        pdf.image(img_path, x=10, y=40, w=180)
                        os.remove(img_path)
                    
                    st.download_button(
                        label="📥 Baixar Relatório (PDF)",
                        data=pdf.output(dest='S').encode('latin1'),
                        file_name="relatorio_eletroforese.pdf",
                        mime="application/pdf"
                    )
                
                except Exception as e:
                    st.error(f"Erro na simulação: {str(e)}")

if __name__ == "__main__":
    main()
