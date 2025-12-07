Flye Installation
=================

Availability
------------

Flye is available for Linux and MacOS platforms.

Bioconda Releases
-----------------

You can get the latest stable release through Bioconda:

    conda install flye

Alternatively, you can get a release version from the github releases page


Docker Installation
-------------------

Flye is also available as a Docker image, which provides a complete, isolated environment
with all dependencies pre-installed. This is the easiest way to get started.

### Using Pre-built Docker Image (CURRENTLY NOT AVAILABLE)

Once available on Docker Hub: 

    docker pull fenderglass/flye:latest

### Building Docker Image from Source

To build the Docker image from the Flye repository:

    git clone https://github.com/fenderglass/Flye
    cd Flye
    docker build -t flye:2.9.6 .

After building or pulling the image, you can run Flye with:

    docker run --rm -v $(pwd):/data flye:2.9.6 --help

See the [Usage documentation](USAGE.md#docker-usage) for more details on running Flye with Docker.


Building Requirements
---------------------

* Python 2.7 or 3.5+ (with setuptools package installed)
* C++ compiler with C++11 support (GCC 4.8+ / Clang 3.3+ / Apple Clang 5.0+)
* GNU make
* Git
* Core OS development headers (zlib, ...)


Local building (without installation)
-------------------------------------

You may use the package locally without system installation.
To get and compile the latest git version, run:

    git clone https://github.com/fenderglass/Flye
    cd Flye
    make

Then, Flye will be available as:

    python bin/flye

Building on ARM architecture
----------------------------

In case building Flye on ARM architecture fails, you might need to modify minimap2 compilation parameters:
https://github.com/fenderglass/Flye/issues/386

Installing from source
----------------------

To install the Flye package into your system, run:

    git clone https://github.com/fenderglass/Flye
	cd Flye
    python setup.py install

Depending on your OS, you might need to add
```--user``` or ```--prefix``` options to the 
install command for the local installation.

After installation, Flye could be invoked via:

    flye

Optionally, run some tests to ensure that installation was successful:

    python flye/tests/test_toy.py
