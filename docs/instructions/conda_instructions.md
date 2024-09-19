# Miniconda Installation

A Linux environment is needed to install and use [Miniconda](https://docs.conda.io/en/latest/miniconda.html), a minimal installer for conda.

`wget` links below are for example. Always check for the latest distributions on the official Miniconda website.

If using a Mac, download the Mac distribution:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
```

If using WSL or Linux, download the Linux 64-Bit Installer:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Install Miniconda using the downloaded file:

```
bash Miniconda3-latest-*-x86_64.sh
```

Follow the prompts:
1. Press `Enter` to review the license agreement.
2. Press `q` to exit the agreement view.
3. Type `yes` to accept the terms.
4. Press `Enter` to confirm the default installation location or enter a custom path.
5. When asked if you wish to initialize Miniconda3, type `yes`.

Close and reopen your terminal for the changes to take effect.

# Conda Environment

It's best practice not to install packages in the base environment. Instead, create a new environment for your project.

Summary of commands to [manage environments](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html):

Create a new environment:
```
conda create --name myenv
```

Activate the environment:
```
conda activate myenv
```

For vSNP3, please refer to the [README](../../README.md).

For information on additional tools, see [Additional Tools](../../docs/instructions/additional_tools.md).

Remember, Miniconda provides a minimal set of packages. If you need additional packages, you can install them using `conda install` within your activated environment.
