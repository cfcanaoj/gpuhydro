# How to use github

To download the source code, you need an account of [github](https://github.com/).

## Make ssh key
After you login one of cfca servers, e.g., `g00.cfca.nao.ac.jp`, first you need to make ssh keys. You have to determine the password.
    
    ssh-keygen -t rsa -f .ssh/id_rsa_git
    
Then two files are created. `.ssh/id_rsa_git` is the private key. `.ssh/id_rsa_git.pub` is the public key. Submit the public key in [github site](https://github.com/).

## Setup

Please edit `~/.ssh/config` and add the following setting.
    
    Host github.com
         IdentityFile ~/.ssh/id_rsa_git

## Download
You can download the file with the following command.

    git clone git@github.com:cfcanaoj/gpuhydro gpuhydro
    
