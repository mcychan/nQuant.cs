FROM ubuntu:24.04
WORKDIR /tmp
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies, add the official backports PPA, and install the SDK
RUN apt-get update && apt-get install -y \
    software-properties-common \
    && add-apt-repository -y ppa:dotnet/backports \
    && add-apt-repository universe \
    && apt-get update && apt-get install -y \
    dotnet-sdk-10.0 libgdiplus \
    && rm -rf /var/lib/apt/lists/*
ADD . /tmp/nQuant.cs
WORKDIR /tmp/nQuant.cs
RUN dotnet build -c Release nQuant.Master/nQuant.Master.csproj -f net6.0 -r linux-x64
RUN dotnet build -c Release nQuant.Console/nQuant.Console.csproj -f net6.0 -r linux-x64
RUN dotnet publish -c Release nQuant.Console/nQuant.Console.csproj -f net6.0 -o ../build -r linux-x64 --self-contained
RUN cp -R samples /tmp/build/
WORKDIR /tmp/build
# docker system prune -a
# docker build -t nquantcs .
# docker run -it nquantcs bash
# docker cp <containerId>:/file/path/within/container /host/path/target
# docker cp foo.txt <containerId>:/foo.txt
