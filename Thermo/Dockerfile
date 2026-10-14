FROM python:3.14-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_ROOT_USER_ACTION=ignore

WORKDIR /app

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY gunicorn.conf.py .
COPY StateFunctionCalculator StateFunctionCalculator
COPY EfficiencyCalculator EfficiencyCalculator
COPY web_app web_app

# Some sources are owner-only (600/700) on disk; make the tree readable by the
# unprivileged user the server runs as.
RUN chmod -R a+rX /app && useradd --create-home app
USER app

CMD ["gunicorn", "--config", "gunicorn.conf.py"]
