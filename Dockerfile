FROM python:3.10
ENV USER=username
ENV HOME=/home/$USER
RUN useradd -m -u 1000 $USER
WORKDIR $HOME/kg
COPY .streamlit $HOME/kg/.streamlit
COPY requirements.txt $HOME/kg/requirements.txt
COPY pages $HOME/kg/pages/
COPY AMR-KG_Database.py $HOME/kg/AMR-KG_Database.py
COPY data $HOME/kg/data/
COPY amrkg_chemspace.html $HOME/kg/amrkg_chemspace.html
COPY start-script.sh $HOME/kg/start-script.sh
RUN mkdir $HOME/kg/models
RUN apt-get update && apt-get install curl -y
RUN curl https://zenodo.org/api/records/13868088/files/ecfp4_rf.pkl/content -o $HOME/kg/models/ecfp4_rf.pkl
RUN curl https://zenodo.org/api/records/13868088/files/erg_rf.pkl/content -o $HOME/kg/models/erg_rf.pkl
RUN curl https://zenodo.org/api/records/13868088/files/maccs_rf.pkl/content -o $HOME/kg/models/maccs_rf.pkl
RUN curl https://zenodo.org/api/records/13868088/files/mhfp6_rf.pkl/content -o $HOME/kg/models/mhfp6_rf.pkl
RUN curl https://zenodo.org/api/records/13868088/files/rdkit_rf.pkl/content -o $HOME/kg/models/rdkit_rf.pkl
RUN pip install --no-cache-dir -r requirements.txt \
    && chmod +x start-script.sh \
    && chown -R $USER:$USER $HOME \
    && rm -rf /var/lib/apt/lists/*
USER $USER
EXPOSE 8501
HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health
ENTRYPOINT ["./start-script.sh"]
