# 1. Start with a standard Ubuntu Linux environment
FROM ubuntu:22.04

# 2. Prevent interactive prompts from stalling the build
ENV DEBIAN_FRONTEND=noninteractive

# 3. Install the standard C++ build tools, CMake, Git, and Python 3
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    python3 \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

# 4. Set the working directory inside the container to /app
WORKDIR /app

# 5. Copy your entire project into the container
COPY . /app

# 6. Install Python dependencies in a virtual environment
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN pip3 install --no-cache-dir -r tools/requirements.txt

# 7. Build the C++ startracker engine
RUN mkdir -p build && cd build && cmake .. && make

# 8. Make the test script executable
RUN chmod +x test_phase1.sh

# 9. By default, open a bash terminal so you can poke around
CMD ["/bin/bash"]
