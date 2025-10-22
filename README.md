# 🧬 Species-Specific Sequence Finder

A **Streamlit web app** for discovering *species-specific genomic fragments*  
and computing marker coverage across target genomes.

**Live Demo (Render):** _Deploying soon on Render.com_

---

## ⚙️ Features
- Upload **target** and **non-target** genomes (.fna / .fasta)
- Automatic **sliding-window slicing** and **BLAST-based filtering**
- Extract **species-specific markers**
- Compute **marker coverage** across target genomes
- Download results (FASTA / CSV / TSV)

---

## 🚀 How to run locally
```bash
git clone https://github.com/jinying0705/specific-sequence-finder.git
cd specific-sequence-finder
pip install -r requirements.txt
streamlit run app.py
