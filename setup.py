from setuptools import setup


setup(name = 'duneggd',
      version = '0.0',
      description = 'DUNE ND Geometries generated by GGD',
      author = 'Jose Palomino and Guang Yang',
      author_email = gyang@nngroup.physics.sunysb.edu',
      license = 'GPLv2',
      url = 'https://github.com/gyang9/GuJo',
      package_dir = {'':'python'},
      packages = ['duneggd', 'duneggd.test'],
      install_requires = [l for l in open("requirements.txt").readlines() if l.strip()],
      # implicitly depends on ROOT
              
  )

