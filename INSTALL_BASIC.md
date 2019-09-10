This page is only for debug purpose.

For the latest installation, please refer to Anaconda &  Docker website.

## Anaconda for Linux

Please refer to `https://www.anaconda.com/`.

```
$ cd /tmp
$ curl -O https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
$ bash Anaconda3-2019.10-Linux-x86_64.sh
```

add `export PATH=~/anaconda3/bin:$PATH` to `~/.bashrc`


## Docker for Ubuntu

Please refer to `https://docs.docker.com/install/`.

```
$ sudo apt-get update
$ sudo apt-get install \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg-agent \
    software-properties-common
$ curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
$ sudo add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
   $(lsb_release -cs) \
   stable"
```

```
$ sudo apt-get update
$ sudo apt-get install docker-ce docker-ce-cli containerd.io

$ sudo docker run hello-world
```

assume the username is: `ubuntu`.

```
$ sudo usermod -a -G docker ubuntu
```
