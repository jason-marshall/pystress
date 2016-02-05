from setuptools import setup

setup(name='PyStress',
      version='0.1',
      description='Manipulation and plotting of stress states',
      url='https://github.com/jason-marshall/pystress',
      author='Jason Marshall',
      author_email='jason.p.marshall@gmail.com',
      license='MIT',
      packages=['pystress'],
      install_requires=[
          'numpy',
          'matplotlib'
      ],
      zip_safe=False)
