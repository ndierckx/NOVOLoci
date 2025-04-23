# 1. Base image
FROM ubuntu:24.10

# 2. Metadata
LABEL maintainer="ndierckx@github" \
      version="0.1" \
      description="NOVOLoci: BLAST+ • MAFFT • Perl MCE & Parallel::ForkManager"

# 3. Environment
ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    PATH=/usr/local/bin:/usr/bin:$PATH

# 4. Install system packages & Perl deps
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      ncbi-blast+ \
      mafft \
      cpanminus \
      libparallel-forkmanager-perl \
      libmce-perl \
      build-essential \
      curl && \
    cpanm --notest MCE Parallel::ForkManager && \
    rm -rf /var/lib/apt/lists/*

# 5. Copy your code in
COPY . /opt/novoloci

# 6. Make scripts executable
RUN chmod -R +x /opt/novoloci

# 7. Set working dir
WORKDIR /opt/novoloci

ENV PATH=/opt/novoloci:$PATH

# 8. Default entrypoint
ENTRYPOINT ["perl", "NOVOLoci0.1.pl"]
