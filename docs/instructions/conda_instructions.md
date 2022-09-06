# Anaconda Installation

Linux environment is needed to install and use the [Anaconda package manager](https://www.anaconda.com/products/distribution).  

`wget` links below are only for example.  One should check for updated distributions.

If using a Mac download the Mac distribution, Mac 64-Bit Command Line installer

```
wget https://repo.anaconda.com/archive/Anaconda3-2022.05-MacOSX-x86_64.sh
```

If using WSL, download the Linux 64-Bit Installer

```
wget https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh
```

Install Anaconda using the downloaded file.

```
bash ./Anaconda3-2022.05-*-x86_64.sh
```

Press `Enter` to review agreement.  Exit agreement, `q`.  Accept terms, `yes`.  Press enter to install in default home directory.  After installation agree to `conda init`.

Close and reopen terminal.

# Anaconda Environment

Do not install packages in base.  Instead make an environment.  

Summary of commands to [manage environments](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)

Create new environment
```
conda create --name myenv
```
```
conda activate myenv
```
For vSNP3 see [README](../../README.md)

[Additional Tools](../../docs/instructions/additional_tools.md)