FROM ubuntu:20.04
LABEL maintainer="Ignacio Vizzo <ivizzo@uni-bonn.de>"

# setup environment
ENV TERM xterm
ENV DEBIAN_FRONTEND=noninteractive

# install essentials
RUN apt-get update && apt-get install --no-install-recommends -y \
    build-essential \
    cmake \
    git \
    pybind11-dev \
    python3 \
    python3-dev \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install Open3D dependencies: https://github.com/intel-isl/Open3D/issues/3388
RUN apt-get update && apt-get install --no-install-recommends -y \
    libgl1 \
    libgomp1 \
    libusb-1.0-0 \
    && rm -rf /var/lib/apt/lists/*

# Install python dependencies
RUN pip3 install click numpy open3d
