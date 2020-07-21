import os
import subprocess
from setuptools import setup, find_packages
from distutils.command.build import build as DistutilsBuild
from setuptools.command.install import install as SetuptoolsInstall



class BioGraphBuild(DistutilsBuild):
    """Build BioGraph dependencies"""
    description = "Build biograph dependencies"
    user_options = DistutilsBuild.user_options + [("ignore-vmd=", None, "whether to ignore VMD during build (yes/no, default no)")]
    def initialize_options(self):
        self.ignore_vmd = None
        DistutilsBuild.initialize_options(self)

    def finalize_options(self):
        self.ignore_vmd = self.ignore_vmd == "yes"
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
        self.git_clone("https://github.com/getcontacts/getcontacts", "GetContacts")


        subprocess.call(["make", "zlib=no"], cwd=os.path.join(os.getcwd(), "cdhit"))
        subprocess.call(["cp", "-r", "cdhit", "biograph"])

        if not self.ignore_vmd:
            self.git_clone("https://github.com/Eigenstate/vmd-python", "VMD")
            subprocess.call(["python3", "setup.py", "build"], cwd=os.path.join(os.getcwd(), "vmd-python"))
            subprocess.call(["python3", "setup.py", "install"], cwd=os.path.join(os.getcwd(), "vmd-python"))

        # TODO: compile perseus
        DistutilsBuild.run(self)

class BioGraphInstall(SetuptoolsInstall):
    def initialize_options(self):
        SetuptoolsInstall.initialize_options(self)


    def finalize_options(self):
        SetuptoolsInstall.finalize_options(self)

    def run(self):
        SetuptoolsInstall.run(self)

setup(
    name="biograph",
    version="0.3",
    description="Represent proteins as graphs for machine learning",
    author="Joaquin Torre Zaffaroni, Leandro Lombardi, Leonardo Cordoba",
    author_email="joaquintorrezaffaroni@gmail.com",
    scripts=["scripts/biograph_pymol_path"],
    packages=["biograph", "getcontacts", "getcontacts.contact_calc", "getcontacts.Applications"],
    package_data={
        "biograph": [
            "perseus/*",
            "cdhit/*"
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
        "cython",
        "gudhi"
    ]
)
