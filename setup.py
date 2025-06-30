from setuptools import setup, find_packages

setup(
    name="modelomics",
    version="0.0.0-beta",
    description="Tools and models for protein design all in one place",
    author="A. Clay Richard",
    author_email="alanclayrichard@gmail.com",
    url="https://github.com/alanclayrichard/modelomics",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.10",
    install_requires=[
        "esm>=3.2.0",
        "biopython>=1.85",
    ],
)