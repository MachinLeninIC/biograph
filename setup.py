import os
import subprocess
from setuptools import setup
from distutils.command.build import build as DistutilsBuild
from setuptools.command.install import install as SetuptoolsInstall



class BioGraphBuild(DistutilsBuild):
    """Build BioGraph dependencies"""
    description = "Build biograph dependencies"
    user_options = []
    def initialize_options(self):
        DistutilsBuild.initialize_options(self)

    def finalize_options(self):
        DistutilsBuild.finalize_options(self)


    def git_clone(self, url, name):
        if os.path.exists(url.split("/")[-1]):
            print("Folder for {} already exists, skipping...".format(name))
            return
        errcode = subprocess.call([
            "git", "clone", url
        ])
        if errcode != 0:
            print("Error downloading {}!".format(name))
            raise Exception("Could not clone {} from {}".format(name, url))

    def check_python_version(self):
        result = subprocess.run(['python', '--version'], stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT)

        if "Python 2" in str(result.stdout):
            return "python3"
        return "python"

    def run(self):
        # TODO: future -- maybe use directly requirements.txt?
        # https://pip.readthedocs.io/en/1.1/requirements.html#requirements-file-format
        # if the package has a `setup.py develop` then it should work
        self.git_clone("https://github.com/weizhongli/cdhit", "CDHit")
        #self.git_clone("https://github.com/Eigenstate/vmd-python", "VMD")
        self.git_clone("https://github.com/getcontacts/getcontacts", "GetContacts")

        #subprocess.call(["make", "zlib=no"], cwd=os.path.join(os.getcwd(), "cdhit"))
        #subprocess.call(["python3", "setup.py", "build"], cwd=os.path.join(os.getcwd(), "vmd-python"))
        #subprocess.call(["python3", "setup.py", "install"], cwd=os.path.join(os.getcwd(), "vmd-python"))
        # vmd =>netcdf, tk

        # perseus?

        #requirements.txt
        DistutilsBuild.run(self)

class BioGraphInstall(SetuptoolsInstall):
    def initialize_options(self):
        SetuptoolsInstall.initialize_options(self)

    def finalize_options(self):
        SetuptoolsInstall.finalize_options(self)

    def run(self):
        #self.run_command("build")
        SetuptoolsInstall.run(self)

setup(
    name="biograph",
    version="0.1",
    description="Represent proteins as graphs for machine learning",
    author="Leonardo Cordoba, Leandro Lombardi, Sebastian Prillo, Joaquin Torre Zaffaroni",
    author_email="joaquintorrezaffaroni@gmail.com",
    packages=["biograph"],
    package_data={
        "biograph": [
            "*"
        ]
    },
    cmdclass={
        "build": BioGraphBuild,
        "install": BioGraphInstall
    },
    install_requires = [
        "networkx>=2.3",
        "pandas>=0.24.2",
        "numpy>=1.17.2",
        "scipy>=1.3.0",
        "requests>=2.18.4",
        "matplotlib>=3.1.0",
        "biopython>=0.1.0",
        "flatten_dict>=0.2.0",
        "tk",
        "scikit-learn",
        "pytest",
        "seaborn",
        "cython"
    ]
)