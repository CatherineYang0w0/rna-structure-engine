ARG BASE_IMAGE=mambaorg/micromamba:1.5.10
FROM ${BASE_IMAGE}

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean -a -y

ENV MPLBACKEND=Agg
ENV PATH=/opt/conda/bin:$PATH

WORKDIR /work
CMD ["python", "-m", "engine.run_rnastructure", "--case", "xist", "--test", "test-1", "--scheme", "rnastr-deigan"]
