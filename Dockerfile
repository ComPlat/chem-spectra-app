# Private variables, not passed in from outside, but helpful for this file
ARG BASE=ubuntu:24.04

# Stage 1: the base image
FROM ${BASE} AS base

# set timezone
ARG TZ=Europe/Berlin
RUN ln -s /usr/share/zoneinfo/${TZ} /etc/localtime

# locales
ENV LANG=en_US.UTF-8 \
    LANGUAGE=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8
RUN echo "LANG=${LANG}" >/etc/locale.conf && \
    echo "LC_ALL=${LANG}" >>/etc/locale.conf && \
    echo "${LANG} UTF-8" >/etc/locale.gen

# install system packages
RUN apt-get -y update && apt-get -y upgrade && \
    apt-get install -y --no-install-recommends --autoremove --fix-missing locales && \
    apt-get clean && \
    locale-gen ${LC_ALL}

# install yq
ADD https://github.com/mikefarah/yq/releases/latest/download/yq_linux_amd64 /bin/yq
RUN chmod +x /bin/yq

# Stage 2: the spectra image, libfreetype6-dev is for matplotlib
FROM base AS spectra
RUN apt-get install -y --no-install-recommends --autoremove --fix-missing git ca-certificates git curl bzip2 \
    gcc g++ libxrender1 libxext-dev pkg-config libfreetype6-dev libgtk2.0-dev

RUN git clone --single-branch --branch f-pbf-deploy --depth=1 https://github.com/ComPlat/chem-spectra-app /app && \
    rm -rf /app/.git

ADD https://raw.githubusercontent.com/ptrxyz/chemotion/e6af03a3fa25c2a830d2e98fd08552b624a77e30/spectra/additives/spectra_config.py /app/instance/config.py
ADD https://raw.githubusercontent.com/ptrxyz/chemotion/e6af03a3fa25c2a830d2e98fd08552b624a77e30/spectra/additives/fake-docker.py /bin/docker
RUN chmod +x /bin/docker

RUN curl -sL https://micro.mamba.pm/api/micromamba/linux-64/latest --output micromamba.tar.bz2 && \
    tar -xf micromamba.tar.bz2 bin && rm micromamba.tar.bz2 && \
    echo "PATH=/bin/micromamba:$PATH" >> /root/.profile

ENV MAMBA_ROOT_PREFIX="/root/micromamba"

RUN /bin/micromamba shell init --shell bash --root-prefix /root/micromamba 

RUN micromamba create -f /app/environment.yml -c conda-forge

RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    micromamba clean --all --yes

RUN mkdir -p /shared /app/instance && \
    ln -s /shared /app/chem_spectra/tmp

WORKDIR /app

# Stage 3: finalize the image
FROM spectra AS app

ENV MAMBA_ROOT_PREFIX="/root/micromamba" \
    FLASK_ENV=production \
    FLASK_DEBUG=0 \
    MSC_HOST=msconvert \
    MSC_PORT=4000 \
    MSC_VALIDATE=true \
    SPECTRA_PORT=4000

EXPOSE 4000

ADD https://github.com/krallin/tini/releases/latest/download/tini /tini
RUN chmod +x /tini

WORKDIR "/app"

ENTRYPOINT ["/tini", "--", "micromamba", "run", "-n", "chemSpec"]
CMD ["gunicorn", "--timeout", "600", "-w", "4", "-b", "0.0.0.0:4000", "server:app"]

# HEALTHCHECK --interval=5s --timeout=3s --start-period=30s --retries=3 \
#     CMD curl --fail http://localhost:4000/ping || exit 1
