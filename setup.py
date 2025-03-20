import os
from setuptools import setup, find_packages

# Safely read README.md to avoid errors if it's missing
def read_long_description():
    try:
        with open("README.md", "r", encoding="utf-8") as f:
            return f.read()
    except FileNotFoundError:
        return "MetaTrimX: A CLI tool for automated metabarcoding data processing."

setup(
    name="MetaTrimX",
    version="1.0.0",
    author="Subramaniam Vijayakumar",
    author_email="subramanyamvkumar@gmail.com",
    description="A CLI tool for automated metabarcoding data processing",
    long_description=read_long_description(),
    long_description_content_type="text/markdown",
    url="https://github.com/SUBRAMANIAM96/MetaTrimX",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    python_requires=">=3.6",
    install_requires=[
        "cutadapt>=4.0",
        "vsearch>=2.21"  # Add vsearch dependency if needed
    ],
    extras_require={
        "dev": ["pytest", "black", "flake8"],
        "docs": ["sphinx", "myst-parser"]
    },
    entry_points={
        "console_scripts": [
            "metatrimx=metatrimx.cli:main"
        ]
    },
    include_package_data=True,
    package_data={
        "metatrimx": ["data/*", "trim_pipeline.sh"]  # Ensure trim_pipeline.sh is included
    },
    scripts=["metatrimx/core/metatrimx.sh"],  # Ensure metatrimx.sh is installed as an executable
    license="MIT",
    keywords="metabarcoding, bioinformatics, DNA sequencing, pipeline"
)

# Ensure metatrimx.sh and trim_pipeline.sh are executable
script_paths = ["metatrimx/core/metatrimx.sh", "metatrimx/trim_pipeline.sh"]
for script in script_paths:
    if os.path.exists(script):
        os.chmod(script, 0o755)
