import setuptools

with open("README.md", "r", encoding = "utf8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sbmltoodepy_SMRuggiero",
    version="0.0.2",
    author="SMR",
    author_email="SteveMRuggiero@gmail.com",
    description="A small example package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    package_data = {
        'sbmltoodepy' : ['sbml_files/*.xml']
    },
    install_requires=[
            'python-libsbml',
            'numpy'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
)