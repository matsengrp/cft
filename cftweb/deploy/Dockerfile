# A Docker container for running the cftweb web server app.
#
# Build with:
# 	docker build --tag cftweb .
# Run with:
# 	docker run -it --rm -p 5000:5000 -v $PWD/output:/home/cftweb/data cftweb

FROM continuumio/miniconda

# Update aptitude with new repo
RUN apt-get update

# Install software and github private repo permissions
RUN apt-get install -y git
RUN apt-get install -y openssh-server

RUN conda config --add channels conda-forge 
RUN conda config --add channels bioconda
RUN conda install -y pandas biopython flask nestly pyqt
RUN conda install -y pip
RUN pip install seqmagick
RUN pip install ete3 scons flask-breadcrumbs


RUN useradd --create-home --shell /bin/bash cftweb
# add ssh key for acessing private github repo.
# from http://stackoverflow.com/a/25977750/1135316
# repo-key was created with,
#     ssh-keygen -q -t rsa -N '' -f repo-key
# then add repo-key.pub to your repository deployment keys.
# On GitHub, go to [your repository] -> Settings -> Deploy keys


RUN mkdir /home/cftweb/.ssh
ADD repo-key /home/cftweb/.ssh/id_rsa
RUN touch /home/cftweb/.ssh/known_hosts
RUN ssh-keyscan -T 60 github.com >> /home/cftweb/.ssh/known_hosts

RUN chown -R cftweb.cftweb /home/cftweb
RUN chmod 700 /home/cftweb/.ssh
RUN chmod 600 /home/cftweb/.ssh/id_rsa

USER cftweb
WORKDIR /home/cftweb


RUN git clone git@github.com:matsengrp/cft.git

WORKDIR /home/cftweb/cft/cftweb
CMD python -m cftweb --file ../../data/metadata.json 
