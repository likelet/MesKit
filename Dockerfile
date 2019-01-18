FROM nfcore/base
LABEL authors="Qi Zhao" \
      description="Docker image containing all requirements for nf-core/multiexseq pipeline"

COPY environment1.yml environment2.yml ./

ENV PATH /opt/conda/envs/multiexseq_facets/bin:$PATH
RUN conda env create -f /environment1.yml -n multiexseq_facets && conda clean -a

ENV PATH /opt/conda/envs/multiexseq_freec/bin:$PATH
RUN conda env create -f /environment2.yml -n multiexseq_freec && conda clean -a

ENV PATH /opt/conda/envs/nf-core-multiexseq-1.0dev/bin:$PATH

CMD ["/bin/bash"]