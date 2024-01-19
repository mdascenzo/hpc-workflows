FROM ubuntu:20.04

# set default values for build args
ARG ID=rnaseq

# copy build and test scripts to container
COPY aws/ami/builds/${ID}/install.sh /install.sh

# create directory for test scripts and copy to container
RUN mkdir /tests
COPY tests/${ID} /tests/${ID}

# modify permissions and run install script
RUN chmod +x /install.sh
RUN /install.sh

# modify permissions and run test script
RUN chmod +x /tests/${ID}/test.sh
WORKDIR /tests/${ID}
ENTRYPOINT ["./test.sh"]