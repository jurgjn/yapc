import setuptools
import yapc

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name = "yapc",
    version = yapc.__version__,
    author = "jurgjn",
    author_email = "jurgjn@users.noreply.github.com",
    description = "Yapc is a (yet another) peak caller for genomic high-throughput sequencing data",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/jurgjn/yapc",
    packages = setuptools.find_packages(),
    scripts = ['bin/yapc'],
    classifiers = [
        "Development Status :: 2 - Pre-Alpha"
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    install_requires = [
        'pandas',
        'pybigwig>=0.3',
        #'idr>=2.0', #TODO cannot find (the actual) IDR on pip...
    ]
)
