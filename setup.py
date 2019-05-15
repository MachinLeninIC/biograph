from setuptools import setup

setup(name = "pyprot",
      packages=["pyprot"],
      package_data={
          'pyprot': [
              '*',
              'perseus/*'
          ]}
      )
