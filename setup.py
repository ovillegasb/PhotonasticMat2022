
from setuptools import setup, find_packages

with open('requirements.txt') as f:
    install_requires = f.read().strip().split('\n')

setup(
    name="StampTools",
    version="0.1.0",
    packages=find_packages(),
    install_requires=install_requires,
    entry_points={
        "console_scripts": [
            "stamptools = stamptools.__main__:main"
        ]
    },
    package_data={
        'stamptools': [
            '*'
        ],
    },
    author="Orlando Villegas",
    author_email="ovillegas.bello0317@gmail.com",
    description='Suite of tools for working with Stamp4 output and input files',
)
