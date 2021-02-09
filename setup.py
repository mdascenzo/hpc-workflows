from setuptools import setup, find_packages
setup(
    name="workflows",
    version="0.0.1",
    packages=find_packages(),
    scripts=['src/bin/create_config.py'],
    install_requires=['natsort', 'ruamel.yaml'],
    # todo: possibly include workflows with install
    # data_files=[
    #    ('workflows/rnaseq', ['src/rnaseq/rnaseq.smk']),
    # ],
    # include_package_data=True
)
