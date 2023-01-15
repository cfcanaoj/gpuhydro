# How to use github

To download the source code, you need an account of [github](https://github.com/).

## Make ssh key
After you login one of cfca servers, e.g., `g00.cfca.nao.ac.jp`, first you need to make ssh keys.
    
    Host github.com
         IdentityFile ~/.ssh/id_rsa_git
	 

Please edit `~/.ssh/config`.
    
    Host github.com
         IdentityFile ~/.ssh/id_rsa_git


After you login the server, `g00.cfca.nao.ac.jp`, follow the instruction. 
    cd /gwork0/<username>
    git clone git@github.com:cfcanaoj/gpuhydro gpuhydro
    cd gpuhydro/KHf90openaccc
    
