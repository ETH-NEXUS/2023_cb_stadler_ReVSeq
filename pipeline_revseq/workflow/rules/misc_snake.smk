import os.path
import sys
import inspect
import copy
import pandas as pd
from snakemake.utils import validate

# import sample map and retrieve sample names
sample_map = pd.read_table(config["inputOutput"]["sample_map"], header=0)
samples = sample_map.set_index("sample", drop=False)
sample_ids = sample_map["sample"].tolist()
validate(samples, "../schema/sample_map.schema.yaml")
lanes = sample_map.set_index("lane", drop=False)
lane_ids = sample_map["lane"].tolist()


#########################################
###  Check config file for missing values
#########################################

fail_instantly = False

# define error object and message
class Error(object):
    def __init__(self, key, name):
        self.__key = key
        self.__name = name

    def __add__(self, other):
        return self

    def __call__(self, wildcards=None):
        sys.exit(
            """
            ===============================================
            You have not specified '{}' for '{}'
            ===============================================
            """.format(
                self.__key, self.__name
            )
        )

    def __getitem__(self, value):
        return Error(key=self.__key, name=self.__name)


# define config object
class Config(object):
    def __init__(self, kwargs, name="Config"):
        self.__name = name
        self.__members = {}
        for (key, value) in kwargs.items():
            if isinstance(value, dict):
                self.__members[key] = Config(kwargs=value, name=key)
            else:
                self.__members[key] = value

    def __getitem__(self, key):
        if key in self.__members:
            return self.__members[key]
        else:
            if fail_instantly:
                sys.exit(
                    """
                    ===============================================
                    You have not specified '{}' for '{}'
                    ===============================================
                    """.format(
                        key, self.__name
                    )
                )
            else:
                return Error(key=key, name=self.__name)


# check with the above class definitions if the config file contains all necessary values
config = Config(config)
