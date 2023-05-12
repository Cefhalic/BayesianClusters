FROM ubuntu:focal


ENV TZ=Europe/London
ENV DEBIAN_FRONTEND=noninteractive
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone




# curl is required to download Anaconda.
RUN apt-get update && apt-get install -y \
    curl \
	libgl1-mesa-glx \
	libegl1-mesa  \
	libxrandr2 \
	libxss1 \
	libxcursor1 \
	libxcomposite1 \
	libasound2 \
	libxi6 \
	libxtst6 \
	python3-dev \
	pip \
	build-essential \
    gedit \
	vim \
	git \
	wget \
	apt-transport-https \
	pkg-config \
	software-properties-common \
	&& apt-get clean


RUN cd /tmp && curl -O https://repo.anaconda.com/archive/Anaconda3-2023.03-1-Linux-x86_64.sh

ENV CONDA_PATH=/opt/anaconda3/Anaconda3-2023.03-1-Linux-x86_64.sh
ENV ENVIRONMENT_NAME=main
SHELL ["/bin/bash", "-c"]
RUN chmod +x /tmp/Anaconda3-2023.03-1-Linux-x86_64.sh
RUN mkdir /root/.conda
RUN bash -c "/tmp/Anaconda3-2023.03-1-Linux-x86_64.sh -b -p ${CONDA_PATH}"

# Initializes Conda for bash shell interaction.
RUN ${CONDA_PATH}/bin/conda init bash

# Upgrade Conda to the latest version
RUN ${CONDA_PATH}/bin/conda update -n base -c defaults conda -y

# Create the work environment and setup its activation on start.
RUN ${CONDA_PATH}/bin/conda create --name ${ENVIRONMENT_NAME} -y
RUN echo conda activate ${ENVIRONMENT_NAME} >> /root/.bashrc

# set conda-forge channel and packages - matplot - cxx-compiler - boost - gsl- make & pip packages
COPY ./environment.yml /tmp/
RUN . ${CONDA_PATH}/bin/activate ${ENVIRONMENT_NAME} \
  && conda env update --file /tmp/environment.yml --prune
  
RUN git clone https://github.com/Cefhalic/BayesianClusters.git Bayesian
RUN cd Bayesian && \
    git checkout PythonBindings

COPY ./Makefile ./Bayesian/
RUN . ${CONDA_PATH}/bin/activate ${ENVIRONMENT_NAME} && \
     cd Bayesian && \
     make




# sudo docker build -t bayes_conda .
# sudo docker run -ti -v "$(pwd):/home" -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix bayes_conda
# source ~/.bashrc

# sudo docker run -it rootproject/root  .

# sudo docker run -ti -v "$(pwd):/home" -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix bayesian/andrew
# sudo docker system prune --all --force

# https://docs.anaconda.com/free/anaconda/install/linux/
