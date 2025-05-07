import streamlit as st
import pubchempy as pcp
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from io import BytesIO
from fpdf import FPDF

st.set_page_config(page_title="Consulta e Simulação - PubChem + Eletroforese", layout="centered")
st.title("🔬 Consulta de Moléculas + Simulação de Eletroforese")

# Criar banco local se não existir
ARQUIVO_CSV = "banco_moleculas.csv"
if not os.path.exists(ARQUIVO_CSV):
    df_vazio = pd.DataFrame(columns=["Nome", "Fórmula", "Peso Molecular", "SMILES", "CID"])
    df_vazio.to_csv(ARQUIVO_CSV, index=False)

# Função para buscar no PubChem
def buscar_molecula(nome):
    try:
        mol = pcp.get_compounds(nome, 'name')[0]
        dados = {
            "Nome": mol.iupac_name,
            "Fórmula": mol.molecular_formula,
            "Peso Molecular": mol.molecular_weight,
            "SMILES": mol.canonical_smiles,
            "CID": mol.cid
        }
        return dados
    except IndexError:
        return None

# Entrada do usuário
nome_molecula = st.text_input("Digite o nome da molécula:")

if nome_molecula:
    dados = buscar_molecula(nome_molecula)
    if dados:
        st.success("Molécula encontrada!")
        st.dataframe(pd.DataFrame([dados]))

        # Mostrar imagem da estrutura
        st.image(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{dados['CID']}/PNG", caption="Estrutura da Molécula")

        # Botão para salvar no banco local
        if st.button("Salvar no banco local"):
            banco_df = pd.read_csv(ARQUIVO_CSV)
            if dados['CID'] not in banco_df['CID'].values:
                banco_df = pd.concat([banco_df, pd.DataFrame([dados])], ignore_index=True)
                banco_df.to_csv(ARQUIVO_CSV, index=False)
                st.success("Molécula salva no banco!")
            else:
                st.warning("Essa molécula já está no banco.")
    else:
        st.error("Molécula não encontrada no PubChem.")

# Mostrar banco de moléculas local
st.subheader("📁 Banco de Moléculas Local")
banco_df = pd.read_csv(ARQUIVO_CSV)
st.dataframe(banco_df, use_container_width=True)

# ---------------------- Simulação de Eletroforese ----------------------
st.subheader("⚡ Simulação de Eletroforese Capilar")
st.markdown("Simule tempos de migração com base em massa molar e carga relativa.")

if not banco_df.empty:
    selecionadas = st.multiselect("Escolha moléculas para simular:", banco_df["Nome"].tolist())
    voltagem = st.slider("Voltagem aplicada (kV):", 5, 30, 15)
    comprimento_capilar = st.slider("Comprimento do capilar (cm):", 10, 100, 50)
    pH = st.slider("pH da solução tampão:", 2.0, 10.0, 7.0, step=0.1)
    ruido = st.checkbox("Adicionar ruído ao cromatograma", value=True)

    if selecionadas:
        tempo_base = comprimento_capilar / (voltagem * 1e3)
        tempos = []
        intensidades = []
        massas = []

        for nome in selecionadas:
            linha = banco_df[banco_df["Nome"] == nome].iloc[0]
            massa = linha["Peso Molecular"]
            carga_simulada = -1 if massa > 120 else 1
            mobilidade = (carga_simulada / massa) * (1 + (pH - 7) * 0.1) * 1e5
            tempo_migracao = comprimento_capilar / (mobilidade * voltagem)
            intensidade = np.exp(-massa / 300) * 100
            tempos.append((nome, tempo_migracao))
            intensidades.append(intensidade)
            massas.append(massa)

        tempos_ordenados = sorted(zip(tempos, intensidades, massas), key=lambda x: x[0][1])
        t = np.linspace(0, max([x[0][1] for x in tempos_ordenados]) + 5, 1000)
        y = np.zeros_like(t)

        for (nome, tempo), intensidade, massa in tempos_ordenados:
            largura = 0.5 + massa / 500
            pico = intensidade * np.exp(-((t - tempo)**2) / (2 * largura**2))
            y += pico

        if ruido:
            y += np.random.normal(0, 0.5, size=len(y))

        fig, ax = plt.subplots()
        ax.plot(t, y, color='purple')
        ax.set_xlabel("Tempo (s)")
        ax.set_ylabel("Intensidade (u.a.)")
        ax.set_title("Cromatograma Simulado de Eletroforese")
        st.pyplot(fig)

        # Exportar PDF
        buffer = BytesIO()
        fig.savefig(buffer, format="png")
        buffer.seek(0)

        class PDF(FPDF):
            def header(self):
                self.set_font("Arial", "B", 12)
                self.cell(0, 10, "Relatório de Simulação de Eletroforese", ln=True, align="C")

        if st.button("📄 Exportar PDF da Simulação"):
            pdf = PDF()
            pdf.add_page()
            pdf.set_font("Arial", size=10)

            pdf.cell(0, 10, f"Moléculas Simuladas: {', '.join([x[0][0] for x in tempos_ordenados])}", ln=True)
            pdf.cell(0, 10, f"Voltagem: {voltagem} kV | Comprimento do capilar: {comprimento_capilar} cm | pH: {pH}", ln=True)

            img_path = "cromatograma_temp.png"
            with open(img_path, "wb") as f:
                f.write(buffer.read())
            pdf.image(img_path, x=10, y=40, w=180)
            pdf.output("simulacao_eletroforese.pdf")

            with open("simulacao_eletroforese.pdf", "rb") as f:
                st.download_button("📥 Baixar PDF", f, file_name="simulacao_eletroforese.pdf")
    else:
        st.info("Selecione moléculas para simular.")
else:
    st.warning("Adicione moléculas ao banco para usar a simulação.")
