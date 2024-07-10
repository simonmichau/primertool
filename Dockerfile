FROM python:3.8-slim
# Uncomment the following lines and insert the proxy URL if you are behind a proxy
#ENV HTTP_PROXY="http://[USER]:[PASSWORD]@[PROXY_IP]:8080/"
#ENV HTTPS_PROXY="http://[USER]:[PASSWORD]@[PROXY_IP]:8080/"
#ENV NO_PROXY=localhost,[LOCAL_IP],0.0.0.0

# Setup the working directory for the container
WORKDIR /primertool

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    software-properties-common \
    git \
    && rm -rf /var/lib/apt/lists/*

#RUN git clone https://github.com/streamlit/streamlit-example.git .

# Copy the requirements file to the container
COPY ./requirements.txt ./

RUN pip3 install -r requirements.txt

# Copy the rest of the application code to the container
COPY ./ ./

EXPOSE 8501

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

ENTRYPOINT ["streamlit", "run", "streamlit_main.py", "--server.port=8501", "--server.address=0.0.0.0"]

