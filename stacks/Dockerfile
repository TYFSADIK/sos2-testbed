# BASE can be any Python 3.12 image; locally a prebuilt one saves the downloads
ARG BASE=python:3.12-slim
FROM ${BASE}
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt pytest==8.3.3
COPY sos sos
COPY sensor sensor
COPY services services
COPY tools tools
COPY tests tests
ENV PYTHONUNBUFFERED=1 PYTHONPATH=/app
