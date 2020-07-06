FROM continuumio/miniconda3:latest

LABEL org.stracquadaniolab.name="baghera"
LABEL org.stracquadaniolab.description="BAyesian Gene-level HERitability Analysis"
LABEL org.stracquadaniolab.maintainer="Viola Fanfani"
LABEL org.stracquadaniolab.maintainer="Giovanni Stracquadanio"
LABEL org.stracquadaniolab.version="$VERSION"

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

