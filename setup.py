import setuptools

with open("README.md", "r", encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name="GS-PRACTICE",
    version="0.1.0",
    author="shirotak",
    author_email="tshiro@kuhp.kyoto-u.ac.jp",
    description="Tumor genomic subtyping using mutational signatures",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/shirotak/GS-PRACTICE",
    project_urls={
        "Bug Tracker": "https://github.com/shirotak/GS-PRACTICE/issues",
    },
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    include_package_data=True,
    entry_points={
        'console_scripts': ['gs-practice = gspractice.run_gspractice:runAsStandalone',
            ]
    },
    python_requires=">=3.7",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    )
