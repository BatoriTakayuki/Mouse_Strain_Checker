FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive PYTHONUNBUFFERED=1 PYTHONPATH=/app STRAIN_DATA=/data
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-venv hisat2 samtools fastp ca-certificates \
    && rm -rf /var/lib/apt/lists/*
RUN python3 -c "import hashlib,urllib.request; p='/tmp/sratoolkit.tar.gz'; urllib.request.urlretrieve('https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.4.1/sratoolkit.3.4.1-ubuntu64.tar.gz',p); assert hashlib.file_digest(open(p,'rb'),'md5').hexdigest() == 'ec6e9056a2bfebcf23c6cd6e02951ef2', 'SRA Toolkit checksum mismatch'" \
    && tar -xzf /tmp/sratoolkit.tar.gz -C /opt \
    && mv /opt/sratoolkit.3.4.1-ubuntu64 /opt/sratoolkit \
    && rm /tmp/sratoolkit.tar.gz
RUN python3 -m venv /opt/venv
ENV PATH=/opt/venv/bin:/opt/sratoolkit/bin:$PATH
RUN prefetch --version && fasterq-dump --version && vdb-validate --version && fastp --version
COPY requirements-web.txt /app/requirements-web.txt
RUN pip install --no-cache-dir -r /app/requirements-web.txt
COPY mouse_strain /app/mouse_strain
COPY tests /app/tests
COPY start.bat start.ps1 /app/
WORKDIR /data
EXPOSE 8501
ENTRYPOINT ["python", "-m", "mouse_strain"]
CMD ["web", "--host", "0.0.0.0"]
