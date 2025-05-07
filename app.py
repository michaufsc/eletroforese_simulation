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

# Funções para gerenciamento do banco de dados
def carregar_banco():
    """Carrega o banco de moléculas do arquivo CSV"""
    if not os.path.exists("banco_moleculas.csv"):
        return pd.DataFrame(columns=["Nome", "Fórmula", "Peso Molecular", "CID", "Data"])
    return pd.read_csv("banco_moleculas.csv")

def salvar_banco(df):
    """Salva o banco de moléculas no arquivo CSV"""
    df.to_csv("banco_moleculas.csv", index=False)

# Interface de inserção de moléculas
def interface_insercao():
    st.header("🧪 Inserir Moléculas no Banco")
    
    col1, col2 = st.columns([2, 3])
    
    with col1:
        with st.form("form_insercao"):
            nome_molecula = st.text_input("Nome da molécula (em inglês):", 
                                        "aspirin",
                                        help="Ex: caffeine, glucose, dopamine")
            
            opcoes_busca = st.radio("Tipo de busca:", 
                                   ["Automática (PubChem)", "Manual"],
                                   horizontal=True)
            
            if opcoes_busca == "Manual":
                formula = st.text_input("Fórmula molecular:")
                peso = st.number_input("Peso molecular (g/mol):", min_value=0.0)
                cid = st.text_input("CID (opcional):")
            else:
                formula = peso = cid = None
            
            if st.form_submit_button("Adicionar ao Banco"):
                with st.spinner("Processando..."):
                    try:
                        banco_df = carregar_banco()
                        
                        if opcoes_busca == "Automática (PubChem)":
                            resultado = pcp.get_compounds(nome_molecula, 'name')
                            if not resultado:
                                st.error("Molécula não encontrada no PubChem!")
                                return
                            
                            mol = resultado[0]
                            novo_registro = {
                                "Nome": mol.iupac_name,
                                "Fórmula": mol.molecular_formula,
                                "Peso Molecular": mol.molecular_weight,
                                "CID": mol.cid,
                                "Data": datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                            }
                        else:
                            if not all([formula, peso]):
                                st.error("Preencha todos os campos obrigatórios!")
                                return
                            
                            novo_registro = {
                                "Nome": nome_molecula,
                                "Fórmula": formula,
                                "Peso Molecular": float(peso),
                                "CID": cid if cid else "N/A",
                                "Data": datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                            }
                        
                        # Verifica se já existe
                        if not banco_df.empty and str(novo_registro["CID"]) in banco_df["CID"].values:
                            st.warning("Esta molécula já está no banco!")
                        else:
                            banco_df = pd.concat([banco_df, pd.DataFrame([novo_registro])], ignore_index=True)
                            salvar_banco(banco_df)
                            st.success("Molécula adicionada com sucesso!")
                            st.rerun()
                    
                    except Exception as e:
                        st.error(f"Erro ao adicionar molécula: {str(e)}")
    
    with col2:
        st.subheader("Banco de Moléculas Atual")
        banco_df = carregar_banco()
        
        if not banco_df.empty:
            # Formatação do DataFrame para exibição
            df_exibicao = banco_df.copy()
            df_exibicao["Peso Molecular"] = df_exibicao["Peso Molecular"].round(2)
            
            st.dataframe(
                df_exibicao.sort_values("Data", ascending=False),
                column_config={
                    "Peso Molecular": st.column_config.NumberColumn(
                        format="%.2f g/mol"
                    ),
                    "Data": st.column_config.DatetimeColumn(
                        format="DD/MM/YYYY HH:mm"
                    )
                },
                use_container_width=True,
                hide_index=True
            )
            
            # Opções de gerenciamento
            with st.expander("Gerenciar Banco"):
                cid_remover = st.selectbox(
                    "Selecione para remover (por CID):",
                    banco_df["CID"].unique()
                )
                
                if st.button("Remover Molécula", type="secondary"):
                    banco_df = banco_df[banco_df["CID"] != cid_remover]
                    salvar_banco(banco_df)
                    st.success(f"Molécula CID {cid_remover} removida!")
                    st.rerun()
                
                if st.download_button(
                    "Exportar Banco (CSV)",
                    data=banco_df.to_csv(index=False).encode('utf-8'),
                    file_name="banco_moleculas.csv",
                    mime="text/csv"
                ):
                    st.toast("Banco exportado com sucesso!", icon="✅")
        else:
            st.info("Nenhuma molécula cadastrada ainda.")

# ... (restante do código da simulação mantido igual)

def main():
    st.title("🔬 Simulador de Eletroforese Capilar")
    
    tab1, tab2 = st.tabs(["Gerenciamento de Moléculas", "Simulação"])
    
    with tab1:
        interface_insercao()
    
    with tab2:
        # ... (código da simulação anterior)
        pass

if __name__ == "__main__":
    main()
