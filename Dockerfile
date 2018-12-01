FROM nfcore/base
LABEL authors="Qi Zhao" \
      description="Docker image containing all requirements for nf-core/multiexseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-multiexseq-1.0dev/bin:$PATH
