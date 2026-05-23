FROM ubuntu:22.04

# Some metadata
LABEL base_image="ubuntu:22.04"
LABEL version="1"
LABEL software="msgfplus"
LABEL software.version="2026.04.15"
LABEL about.summary="MS-GF+ performs peptide identification by scoring MS/MS spectra against peptides derived from a protein sequence database."
LABEL about.home="https://github.com/bigbio/msgfplus"
LABEL about.documentation="https://github.com/bigbio/msgfplus"
LABEL about.license_file="https://github.com/bigbio/msgfplus/blob/dev/LICENSE.txt"
LABEL about.tags="Proteomics"
LABEL maintainer="Yasset Perez-Riverol <ypriverol@gmail.com>"

ENV DEBIAN_FRONTEND=noninteractive

# Install OpenJDK 17 JDK and Maven for building
RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-17-jdk-headless \
    maven \
    && rm -rf /var/lib/apt/lists/*

# Build MS-GF+ from source
WORKDIR /build
COPY pom.xml .
COPY src/ src/
RUN mvn -B package -DskipTests -q && \
    cp target/MSGFPlus.jar /usr/local/bin/MSGFPlus.jar && \
    rm -rf /build /root/.m2

# Remove Maven (no longer needed at runtime)
RUN apt-get remove -y maven && apt-get autoremove -y && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set appropriate permissions
RUN chmod +x /usr/local/bin/MSGFPlus.jar

WORKDIR /data/
