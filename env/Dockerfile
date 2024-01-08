FROM ubuntu:20.04

ARG TARGETPLATFORM

RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata
RUN apt-get install --no-install-recommends -y \
    curl \
    zip \
    unzip \
    python3 \
    python3-pip \
    git \
    vim \
    ssh \
    jq

# install AWS parallelcluster, add path to root
RUN pip3 install "aws-parallelcluster==3.8.0" --upgrade --user
ENV PATH "$PATH:/root/.local/bin"

ENV NODE_VERSION 18.13.0

# install node 
RUN curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.38.0/install.sh | bash && \
  . /root/.nvm/nvm.sh && \
  nvm install $NODE_VERSION && \
  nvm alias default $NODE_VERSION && \
  nvm use default

ENV NVM_DIR /root/.nvm 
ENV NODE_PATH $NVM_DIR/v$NODE_VERSION/lib/node_modules
ENV PATH $NVM_DIR/v$NODE_VERSION/bin:$PATH

# install aws commandline v2

# install the correct version of AWS CLI based on the architecture
RUN if [ "$TARGETPLATFORM" = "linux/amd64" ]; then \
      curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"; \
    elif [ "$TARGETPLATFORM" = "linux/arm64" ]; then \
      curl "https://awscli.amazonaws.com/awscli-exe-linux-aarch64.zip" -o "awscliv2.zip"; \
    else \
      echo "Unsupported platform: $TARGETPLATFORM"; \
      exit 1; \
    fi

RUN unzip awscliv2.zip && ./aws/install

# cd /tmp
RUN apt install -y wget && \
  wget https://releases.hashicorp.com/packer/1.8.5/packer_1.8.5_linux_arm64.zip && \
  unzip packer* && \
  mkdir /opt/packer && \
  mv packer /opt/packer/ && \
  ln -s /opt/packer/packer /usr/local/bin/packer 


ENTRYPOINT ["tail", "-f", "/dev/null"]
