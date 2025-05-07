import streamlit as st
import pubchempy as pcp
import pandas as pd
import numpy as np
from fpdf import FPDF
from datetime import datetime
import matplotlib.pyplot as plt

# Configuração da página
st.set_page_config(
    page_title="CE Simulator PRO",
    layout="wide",
    page_icon="🧪"
)

# Banco de dados de moléculas
BANCO_MOLECULAS = {
    "Ácido Gálico": {
        "massa": 170.12,
        "carga": -1,
        "raio_hidro": 3.5,
        "lambda_max": 270
    },
    "Quercetina": {
        "massa": 302.23,
        "carga": -1,
        "raio_hidro": 4.2,
        "lambda_max": 370
    },
    "Cafeína": {
        "massa": 194.19,
        "carga": 0,
        "raio_hidro": 3.8,
        "lambda_max": 273
    },
    "Ácido Ascórbico": {
        "massa": 176.12,
        "carga": -1,
        "raio_hidro": 3.6,
        "lambda_max": 265
    }
}

# Classe de simulação científica aprimorada
class EletroforeseSimulator:
    def __init__(self):
        self.constantes = {
            'epsilon': 78.5,          # Constante dielétrica da água
            'viscosidade_agua': 0.89, # cP a 25°C
            'faraday': 96485,         # C/mol
            'R': 8.314                # J/(mol·K)
        }

    def calcular_mobilidade(self, composto, pH, temperatura):
        """Calcula a mobilidade considerando propriedades do composto"""
        T = temperatura + 273.15
        carga = self.calcular_carga_efetiva(composto['carga'], pH)
        
        # Cálculo da mobilidade eletroforética
        mu_ef = (carga * self.constantes['faraday']) / (
            6 * np.pi * self.constantes['viscosidade_agua'] * 
            composto['raio_hidro'] * 1e-10
        )
        
        # Cálculo da mobilidade eletroosmótica (EOF)
        zeta = -0.04 * (pH - 3)  # Modelo simplificado
        mu_eo = (self.constantes['epsilon'] * zeta * 1e-3) / (
            4 * np.pi * self.constantes['viscosidade_agua'] * 1e-3
        )
        
        return (mu_ef + mu_eo) * 1e8  # Em unidades práticas

    def calcular_carga_efetiva(self, carga_base, pH):
        """Modelo simplificado de ionização"""
        return carga_base * (1 / (1 + 10**(pH - 7)))

# Interface principal
def main():
    st.title("🧪 CE Simulator PRO")
    simulator = EletroforeseSimulator()

    # Sidebar com configurações
    with st.sidebar:
        st.header("⚙ Configurações")
        voltagem = st.slider("Voltagem (kV)", 5, 30, 15)
        comprimento = st.slider("Comprimento do capilar (cm)", 10, 100, 50)
        pH = st.slider("pH do tampão", 2.0, 10.0, 7.0, 0.1)
        temperatura = st.slider("Temperatura (°C)", 15, 40, 25)
        lambda_detecao = st.slider("Comprimento de onda (nm)", 200, 400, 280)

    # Abas principais
    tab1, tab2, tab3 = st.tabs(["Banco de Moléculas", "Simulação", "Relatório"])

    with tab1:
        st.header("📚 Banco de Moléculas")
        
        # Seleção de moléculas do banco
        selecionadas = st.multiselect(
            "Selecione moléculas para simulação:",
            list(BANCO_MOLECULAS.keys()),
            default=["Ácido Gálico", "Quercetina"]
        )

        # Busca no PubChem
        st.subheader("🔍 Adicionar Nova Molécula")
        nova_molecula = st.text_input("Nome da molécula (em inglês):", "aspirin")
        
        if st.button("Buscar no PubChem"):
            try:
                compound = pcp.get_compounds(nova_molecula, 'name')[0]
                info = {
                    "Nome": compound.iupac_name,
                    "Fórmula": compound.molecular_formula,
                    "Massa": compound.molecular_weight
                }
                st.success("Molécula encontrada!")
                st.json(info)
                
                if st.button("Adicionar ao Banco Temporário"):
                    BANCO_MOLECULAS[compound.iupac_name] = {
                        "massa": compound.molecular_weight,
                        "carga": -1,  # Valor padrão
                        "raio_hidro": 4.0,
                        "lambda_max": 270
                    }
            except Exception as e:
                st.error(f"Erro na busca: {str(e)}")

    with tab2:
        st.header("⚡ Simulação")
        
        if st.button("Executar Simulação", type="primary"):
            if not selecionadas:
                st.warning("Selecione pelo menos uma molécula!")
                return

            resultados = []
            tempos = []
            sinais = []

            for molecula in selecionadas:
                props = BANCO_MOLECULAS[molecula]
                
                # Calcular mobilidade
                mu = simulator.calcular_mobilidade(props, pH, temperatura)
                
                # Calcular tempo de migração
                tempo = (comprimento * 1e-2) / (mu * voltagem * 1e3)
                tempos.append(tempo)
                
                # Gerar pico
                t = np.linspace(0, max(tempos)*1.5, 1000)
                sinal = np.exp(-(t - tempo)**2 / (0.1 * tempo**2)) * 100
                sinais.append(sinal)
                
                # Coletar resultados
                resultados.append({
                    "Molécula": molecula,
                    "Mobilidade": mu,
                    "Tempo": tempo,
                    "Intensidade": np.max(sinal)
                })

            # Plotar resultados
            fig, ax = plt.subplots(figsize=(10, 6))
            for sinal, tempo, molecula in zip(sinais, tempos, selecionadas):
                ax.plot(t, sinal, label=f"{molecula} ({tempo:.2f}s)")
            
            ax.set_xlabel("Tempo (s)")
            ax.set_ylabel("Intensidade (UA)")
            ax.set_title("Eletroferograma Simulado")
            ax.legend()
            ax.grid(True)
            st.pyplot(fig)

            # Exibir tabela de resultados
            df_resultados = pd.DataFrame(resultados)
            st.dataframe(df_resultados.style.format({"Mobilidade": "{:.2e}", "Tempo": "{:.2f}"}))

    with tab3:
        st.header("📊 Relatório")
        if 'resultados' in locals():
            pdf = FPDF()
            pdf.add_page()
            pdf.set_font("Arial", size=12)
            
            # Cabeçalho
            pdf.cell(0, 10, "Relatório de Simulação CE", ln=True, align='C')
            pdf.cell(0, 10, f"Data: {datetime.now().strftime('%d/%m/%Y %H:%M')}", ln=True)
            
            # Parâmetros
            pdf.cell(0, 10, "Parâmetros da Simulação:", ln=True)
            pdf.cell(0, 10, f"- Voltagem: {voltagem} kV", ln=True)
            pdf.cell(0, 10, f"- pH: {pH}", ln=True)
            pdf.cell(0, 10, f"- Temperatura: {temperatura} °C", ln=True)
            
            # Resultados
            pdf.cell(0, 10, "Resultados:", ln=True)
            for idx, res in enumerate(resultados, 1):
                pdf.cell(0, 10, f"{idx}. {res['Molécula']}: {res['Tempo']:.2f} s", ln=True)
            
            # Salvar PDF
            st.download_button(
                label="⬇️ Baixar Relatório Completo",
                data=pdf.output(dest='S').encode('latin1'),
                file_name="relatorio_ce.pdf",
                mime="application/pdf"
            )
        else:
            st.info("Execute uma simulação primeiro para gerar o relatório")

if __name__ == "__main__":
    main()
