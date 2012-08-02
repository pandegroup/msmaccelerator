from distutils.core import setup
import glob

setup(name='msmaccelerator',
      version = '0.1',
      packages = ['msmaccelerator'],
      package_dir = {'msmaccelerator':'lib'},
      scripts = glob.glob('scripts/'),
)
