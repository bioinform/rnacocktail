from setuptools import setup, find_packages

version = "Unknown"
for line in open("src/_version.py"):
    if line.startswith("__version__"):
        version = line.strip().split("=")[1].strip().replace('"', '')

print version
setup(
    name='RNACocktail Pipeline',
    version=version,
    description='RNACocktail: A comprehensive framework for accurate and efficient RNA-Seq analysis',
    author='Roche Sequencing Solutions, Inc',
    author_email='bina.rd@roche.com',
    url='https://github.com/bioinform/rnacocktail',
    packages=find_packages(),
    install_requires=["pysam", "pybedtools"],
    scripts=['scripts/run_rnacocktail.py','scripts/hisat2_jun2bed.py',
             'scripts/gpd2gtf.py']
)
