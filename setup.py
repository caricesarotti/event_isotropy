from setuptools import setup, find_packages

extras_require = {"test": ["matplotlib", "prettytable"]}

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="eventIsotropy",
    version="0.0.1",
    author="Cari Cesarotti",
    description="a robust measure of event isotropy at colliders ",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    url="https://github.com/caricesarotti/event_isotropy",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    python_requires=">=3.6",
    install_requires=["POT", "astropy-healpix"],
    extras_require=extras_require,
)
