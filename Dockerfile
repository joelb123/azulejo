FROM python:3.8-slim AS builder

WORKDIR /project

RUN pip install poetry==1.1.4

RUN apt update && apt install -y --no-install-recommends \
  g++ \
  gcc \
  libc6-dev \
  make \
  patch \
  && rm -rf /var/lib/apt/lists/*
COPY pyproject.toml poetry.lock ./
RUN POETRY_VIRTUALENVS_IN_PROJECT=true poetry install --no-root --no-dev 
# https://github.com/mbello/memory-tempfile/pull/8
RUN sed -i'' 's/mp\[8\]/mp[mp.index("-",6)+1]/' $(find . -name memory_tempfile.py)
COPY . .
RUN POETRY_VIRTUALENVS_IN_PROJECT=true poetry install --no-dev
RUN POETRY_VIRTUALENVS_IN_PROJECT=true poetry run azulejo install dagchainer_tool
RUN POETRY_VIRTUALENVS_IN_PROJECT=true echo y | poetry run azulejo install usearch

FROM python:3.8-slim AS final
RUN apt update && apt install -y --no-install-recommends \
  muscle \
  ncbi-blast+ \
  && rm -rf /var/lib/apt/lists/*

COPY --from=builder /project /project
ENV VIRTUAL_ENV=/project/.venv
ENV PATH=${VIRTUAL_ENV}/bin:${PATH}

ENTRYPOINT ["azulejo"]
