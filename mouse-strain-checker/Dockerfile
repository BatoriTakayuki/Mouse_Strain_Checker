FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive PYTHONUNBUFFERED=1 PYTHONPATH=/app STRAIN_DATA=/data
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-venv hisat2 samtools sra-toolkit ca-certificates \
    && rm -rf /var/lib/apt/lists/*
RUN python3 -m venv /opt/venv
ENV PATH=/opt/venv/bin:$PATH
COPY requirements-web.txt /app/requirements-web.txt
RUN pip install --no-cache-dir -r /app/requirements-web.txt
COPY mouse_strain /app/mouse_strain
COPY tests /app/tests
WORKDIR /data
EXPOSE 8501
ENTRYPOINT ["python", "-m", "mouse_strain"]
CMD ["web", "--host", "0.0.0.0"]
