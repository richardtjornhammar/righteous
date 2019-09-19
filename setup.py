import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name = "righteous-fa",
    version = "1.2.0",
    author = "Richard Tjörnhammar",
    author_email = "richard.tjornhammar@gmail.com",
    description = "Righteous Pathway Factor Analysis",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/richardtjornhammar/impetuous",
    packages = setuptools.find_packages('src'),
    package_dir = {'righteous':'src/righteous','quantification':'src/quantification','convert':'src/convert','examples':'src/examples'},
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
)
