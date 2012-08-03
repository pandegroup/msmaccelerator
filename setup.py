from distutils.core import setup
import glob

setup(name='msmaccelerator',
      version = '0.1',
      packages = ['msmaccelerator', 'msmaccelerator.scripts'],
      package_dir = {'msmaccelerator':'lib',
                     'msmaccelerator.scripts': 'scripts'},
      scripts = glob.glob('scripts/*'),
)