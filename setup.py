from setuptools import setup, find_packages

setup(
    name='tcwindprofile',
    version='2.0.3',
    packages=find_packages(),
    install_requires=[],
    author='Dan Chavas',
    author_email='drchavas@gmail.com',
    description='Create a fast and robust radial profile of the tropical cyclone rotating wind and pressure field from inputs Vmax, R34kt, latitude, translation speed, and environmental pressure',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/drchavas/tcwindprofile',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    license='MIT',
#    license_files=['LICENSE'],
    python_requires='>=3.6',
)
