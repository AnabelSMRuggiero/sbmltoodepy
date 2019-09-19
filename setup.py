import setuptools

with open("README.md", "r", encoding = "utf8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sbmltoodepy",
    version="1.0.3.post1",
    author="Steve M Ruggiero, Ashlee N Ford Versypt",
    author_email="SteveMRuggiero@gmail.com",
    description="A package that creates a python implementation of an SBML model.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/SMRuggiero/sbmltoodepy",
    packages=setuptools.find_packages(),
    package_data = {
        'sbmltoodepy' : ['sbml_files/*.xml']
    },
    install_requires=[
            'python-libsbml',
            'numpy',
            'scipy'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
)