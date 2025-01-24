from setuptools import setup, find_packages
import os
import subprocess
import sys

# Function to check if R is installed
def check_r_installed():
    try:
        subprocess.check_call(['R', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except subprocess.CalledProcessError:
        return False

# Function to check if RevBayes is installed
def check_revbayes_installed():
    try:
        subprocess.check_call(['rb', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except subprocess.CalledProcessError:
        return False

# Function to check R packages
def check_r_packages():
    """Check for required R packages."""
    r_script = """
        missing_pkgs <- setdiff(c('ggplot2', 'RevGadgets', 'coda'), installed.packages()[,'Package'])
        if (length(missing_pkgs) > 0) {
            message('The following R packages are missing: ', paste(missing_pkgs, collapse=', '))
            message('Please install them in R with: install.packages(c("ggplot2", "RevGadgets", "coda"))')
        } else {
            message('All required R packages are installed.')
        }
    """
    try:
        subprocess.run(["R", "-e", r_script], check=True)
    except subprocess.CalledProcessError:
        return False
        



if not check_r_installed():
    print("R is not installed. Please install R from: https://cran.r-project.org/")
if not check_revbayes_installed():
    print("RevBayes is not installed. Please install RevBayes from: https://revbayes.github.io/")
if not check_r_packages():
    print("Could not verify R packages. Please ensure 'ggplot2', 'RevGadgets', and 'coda' are installed in R.")


setup(
    name="traderpros",
    version="0.1.0",
    author="Guilherme Azevedo",
    description="A CLI tool for generating RevBayes scripts and run a Trait Dependend Rates Protracted Speciation",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/ghfazevedo/traderpros",
    packages=["traderpros"],
    package_dir={"traderpros": "src/"},  
    include_package_data=True,
    install_requires=[
        'numpy',
        'dendropy',
        'seaborn',
        'matplotlib'
    ],
    entry_points={
        'console_scripts': [
            'traderpros = traderpros.traderpros:main',
            'conspecific_binary=traderpros.conspecific_binary:main',
            'conspecific_probs=traderpros.conspecific_probs:main',
            'sum_and_plot=traderpros.sum_and_plot:main',
            'traderpros_sim=traderpros.traderpros_sim:main'
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
    ],
    python_requires='>=3.6',
)

