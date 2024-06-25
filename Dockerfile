# Select base image (can be ubuntu, python, shiny etc)
FROM python:3.10-slim

# Create user name and home directory variables. 
# The variables are later used as $USER and $HOME. 
ENV USER=username
ENV HOME=/home/$USER

# Add user to system
RUN useradd -m -u 1000 $USER

# Set working directory (this is where the code should go)
WORKDIR $HOME/kg

# Update system and install dependencies.
RUN apt-get update && apt-get install --no-install-recommends -y \
    build-essential \
    software-properties-common

# Copy code and start script (this will place the files in home/username/)
COPY requirements.txt $HOME/kg/requirements.txt
COPY pages $HOME/kg/pages/
COPY AMR-KG_Database.py $HOME/kg/AMR-KG_Database.py
COPY data $HOME/kg/data/
COPY amrkg_chemspace.html $HOME/kg/amrkg_chemspace.html
COPY start-script.sh $HOME/kg/start-script.sh
RUN mkdir $HOME/kg/models

RUN apt-get install curl -y
RUN curl  https://zenodo.org/api/records/3407840/files/chem_phys_rf.pkl -o $HOME/kg/models/chem_phys_rf.pkl -L
RUN curl  https://zenodo.org/api/records/3407840/files/ecfp_rf.pkl -o $HOME/kg/models/ecfp_rf.pkl -L
RUN curl  https://zenodo.org/api/records/3407840/files/erg_rf.pkl -o $HOME/kg/models/erg_rf.pkl -L
RUN curl  https://zenodo.org/api/records/3407840/files/maccs_rf.pkl -o $HOME/kg/models/maccs_rf.pkl -L
RUN curl  https://zenodo.org/api/records/3407840/files/mhfp6_rf.pkl -o $HOME/kg/models/mhfp6_rf.pkl -L
RUN curl  https://zenodo.org/api/records/3407840/files/rdkit_rf.pkl -o $HOME/kg/models/rdkit_rf.pkl -L

RUN pip install --no-cache-dir -r requirements.txt \
    && chmod +x start-script.sh \
    && chown -R $USER:$USER $HOME \
    && rm -rf /var/lib/apt/lists/*

USER $USER
EXPOSE 8501

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

ENTRYPOINT ["./start-script.sh"]