VERSION 0.6
FROM python:3.10-slim
WORKDIR /code

deps:
    RUN pip install wheel
    COPY requirements.txt ./
    RUN pip wheel -r requirements.txt --wheel-dir=wheels

build:
    FROM +deps
    COPY dggrid4py dggrid4py
    SAVE ARTIFACT dggrid4py /dggrid4py
    SAVE ARTIFACT wheels /wheels

docker:
    COPY +build/dggrid4py dggrid4py
    COPY +build/wheels wheels
    COPY requirements.txt ./
    RUN pip install --no-index --find-links=wheels -r requirements.txt
    # ENTRYPOINT ["python3", "./dggrid4py/hello.py"]
    SAVE IMAGE dggrid4py:0.2.6

