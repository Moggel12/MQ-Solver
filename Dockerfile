FROM ubuntu:22.04

ENV TZ=Europe/Copenhagen
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt update
RUN apt upgrade -y
RUN apt install build-essential valgrind git cmake sagemath python3 python3-pip texlive-full latexmk -y
RUN pip3 install compiledb
