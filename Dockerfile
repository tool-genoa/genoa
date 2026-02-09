FROM ubuntu:24.04

RUN apt-get update && apt-get install -y --no-install-recommends \
    vim \
    git \
    build-essential \
    python3 \
    python3-pip \
    python3-dev \
    swig \
    libopenbabel-dev \
    openbabel \
    imagemagick \
    graphviz \
    gfortran \
    libnetcdff-dev \
    gawk \
    cmake \
    wget \
    unzip \
    curl \
    zlib1g-dev \
    libcurl4-openssl-dev \
    m4 \
    scons \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# Fix openbabel include path
RUN if [ ! -d "/usr/local/include/openbabel3" ]; then \
    ln -s /usr/include/openbabel3 /usr/local/include/openbabel3; \
    fi

# GENOA
WORKDIR /app/genoa
COPY ./genoa/requirements.txt .
RUN pip3 install --no-cache-dir --break-system-packages -r requirements.txt
COPY ./genoa/. .

# Enable genoa
RUN pip install -e . --break-system-packages

# Box model: BOXMODEL4GECKO
WORKDIR /app/genoa/vendor/boxmodel4gecko-v-1-0-genoa
RUN ./build.sh --arch docker.gnu

# Box model: SSH-aerosol
## RUN git clone https://github.com/sshaerosol/ssh-aerosol.git /app/genoa/vendor/ssh-aerosol
# WORKDIR /app/genoa/vendor/ssh-aerosol/tools/install_docker
# RUN chmod +x setup_ssh-aerosol.sh
# RUN ./setup_ssh-aerosol.sh > setup.log 2>&1

# WORKDIR /app/genoa/vendor/ssh-aerosol
# RUN ./compile > compile.log 2>&1

# Work directory for runtime
WORKDIR /app/genoa
