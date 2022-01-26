FROM continuumio/miniconda3:latest

LABEL org.opencontainers.image.title="baghera"
LABEL org.opencontainers.image.description="BAyesian Gene-level HERitability Analysis."
LABEL org.opencontainers.image.url="https://github.com/stracquadaniolab/baghera"
LABEL org.opencontainers.image.documentation="https://github.com/stracquadaniolab/baghera"
LABEL org.opencontainers.image.source="https://github.com/stracquadaniolab/baghera"
LABEL org.opencontainers.image.revision="2.1.12"
LABEL org.opencontainers.image.vendor="stracquadaniolab"
LABEL org.opencontainers.image.authors="Viola Fanfani, Giovanni Stracquadanio"

######################################################################
# update sources and add tini
######################################################################
RUN apt-get update --fix-missing \
    && apt-get clean \
    && apt-get install tini \
    && apt-get autoremove \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# create directory layout
RUN mkdir -p /opt/baghera && mkdir -p /data
COPY . /opt/baghera

######################################################################
# adding requirements using conda
######################################################################
RUN conda env update -f /opt/baghera/docker/environment.yml \
    && conda clean --all

######################################################################
# installing python packages
######################################################################
RUN pip --no-cache-dir install --compile /opt/baghera

ENTRYPOINT ["/usr/bin/tini", "--", "baghera-tool"]

