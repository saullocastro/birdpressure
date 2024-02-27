import os
import inspect
import subprocess
from setuptools import setup, find_packages


is_released = True
version = '0.0.1'


def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout=subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        git_revision = out.strip().decode('ascii')
    except OSError:
        git_revision = "Unknown"

    return git_revision


def get_version_info(version, is_released):
    fullversion = version
    if not is_released:
        git_revision = git_version()
        fullversion += '.dev0+' + git_revision[:7]
    return fullversion


def write_version_py(version, is_released, filename='birdpressure/version.py'):
    fullversion = get_version_info(version, is_released)
    with open("./birdpressure/version.py", "wb") as f:
        f.write(b'__version__ = "%s"\n' % fullversion.encode())
    return fullversion


# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    setupdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    return open(os.path.join(setupdir, fname)).read()


#_____________________________________________________________________________

install_requires = [
        "numpy",
        "sympy",
        "matplotlib",
        ]

#Trove classifiers
CLASSIFIERS = """\

Development Status :: 1 - Planning
Intended Audience :: Education
Intended Audience :: Science/Research
Intended Audience :: Developers
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Mathematics
Topic :: Education
Topic :: Software Development
Topic :: Software Development :: Libraries :: Python Modules
Operating System :: Microsoft :: Windows
Operating System :: Unix
Operating System :: POSIX :: BSD
Programming Language :: Python :: 3.8
Programming Language :: Python :: 3.9
Programming Language :: Python :: 3.10
Programming Language :: Python :: 3.11
License :: OSI Approved :: BSD License

"""

fullversion = write_version_py(version, is_released)

data_files = [('', [
        'README.md',
        'LICENSE',
        ])]

keywords = [
    'crashworthiness',
    'bird strike',
    'pressure',
    'preliminary design',
]

s = setup(
    name = "birdpressure",
    version = fullversion,
    author = "Jonas Bertholdt, Saullo G. P. Castro",
    author_email = "jonas.bertholdt@gmail.com, S.G.P.Castro@tudelft.nl",
    description = ("Calculate peak pressure due to bird impact"),
    license = "3-Clause BSD",
    keywords = keywords,
    url = "https://github.com/saullocastro/birdpressure",
    packages=find_packages(),
    data_files=data_files,
    long_description=read('README.md'),
    long_description_content_type = 'text/markdown',
    classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
    install_requires=install_requires,
    include_package_data=True,
)

