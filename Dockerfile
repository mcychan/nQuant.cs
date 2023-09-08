FROM ubuntu:22.04
WORKDIR /tmp
RUN apt update -y
RUN apt install -y ca-certificates dotnet6 libgdiplus
ADD . /tmp/nQuant.cs
WORKDIR /tmp/nQuant.cs
RUN dotnet build -c Release nQuant.Master/nQuant.Master.csproj -r linux-x64
RUN dotnet build -c Release nQuant.Console/nQuant.Console.csproj -r linux-x64
RUN dotnet publish -c Release nQuant.Console/nQuant.Console.csproj -o ../build -r linux-x64 --self-contained
RUN cp -R samples /tmp/build/
WORKDIR /tmp/build
# docker system prune -a
# docker build -t nquantcs .
# docker run -it nquantcs bash
# docker cp <containerId>:/file/path/within/container /host/path/target
# docker cp foo.txt <containerId>:/foo.txt
