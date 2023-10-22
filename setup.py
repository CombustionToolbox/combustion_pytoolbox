import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name = "Combustion PyToolbox",
    version = "0.0.2",
    author = "Alberto Cuadra-Lara",
    author_email = "acuadra@ing.uc3m.es",
    description = "A Python Thermochemical Tool",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/AlbertoCuadra/ThermochemicalCode_Python",
    packages=setuptools.find_packages(),
    install_requires  = ['numpy', 'scipy', 'pandas', 'palettable', 'seaborn'], # List all your dependencies inside the list
    license = 'GPL-3.0 License'
)