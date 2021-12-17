FROM debian:bullseye-slim 

RUN apt-get update && apt-get install -y \
    curl \
    zip \
    python3 \
    python3-pip \
    git

# install AWS parallelcluster, add path to root
RUN pip3 install "aws-parallelcluster==2.10.1" --upgrade --user
ENV PATH "$PATH:/root/.local/bin"

# install aws commandline v2
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
  unzip awscliv2.zip && \
  ./aws/install

ENTRYPOINT ["tail", "-f", "/dev/null"]