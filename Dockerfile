FROM python:3.10-slim

LABEL maintainer="Priyal Tripathi <priyaltripathi2910@gmail.com>"
LABEL description="BioSeq — A bioinformatics sequence analysis toolkit"

ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

RUN useradd --create-home appuser

WORKDIR /bio_seq_v1

COPY pyproject.toml README.md ./
COPY bio_seq_v1/ bio_seq_v1/

RUN pip install --no-cache-dir .

USER appuser

ENTRYPOINT ["bioseq"]