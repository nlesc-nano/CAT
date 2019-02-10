""" A module designed for sanitizing and interpreting the input file. """

__all__ = ['get_job_settings', 'lower_dict_keys']

import time
import itertools
from os.path import join

from ..misc import get_time


def lower_dict_keys(dic):
    """ Turn all keys in a dictionary, or class derived from dictionary, to lowercase.
    Dictionaries are searched recursivly. """
    if isinstance(dic, dict):
        for key in dic:
            # Check if a key is lowercase; turn to lowercase if not
            try:
                key_lower = key.lower()
                if key != key_lower:
                    dic[key_lower] = dic[key]
                    del dic[key]
            except AttributeError:
                pass

            # Check if a value is a dictionary or a class derived from dictionary
            if isinstance(dic[key_lower], dict):
                dic[key_lower] = lower_dict_keys(dic[key_lower])

    return dic


def get_job_settings(arg_dict, jobs=1):
    """
    """
    if isinstance(arg_dict, bool):
        ret = [None for i in range(jobs*2)]
        print(get_time() + 'No user-specified jobs & settings found for qd_dissociate, \
              switching to defaults')
    else:
        try:
            ret = [arg_dict[item] for item in arg_dict]
            len_ret = len(ret)
        except TypeError:
            raise TypeError('Only booleans, dictiories or dictionary derived objects are \
                            valid when defining jobs')

        # Pad with <None> ret if is smaller than 2 * *jobs*
        if len_ret < jobs*2:
            for i in range(jobs*2 - len_ret):
                ret.append(None)
            print(get_time() + 'No jobs & settings have been specified found for the \
                  last ' + str(jobs - len_ret/2) + ' jobs, switching to defaults')
        # Pop entries from ret if it is larger than 2 * *jobs*
        elif len_ret > jobs*2:
            ret = ret[0:jobs*2]
            print(get_time() + str(len_ret / 2) + ' jobs have been specified while the \
                  argument only support ' + str(jobs) + ', the last ' + str(len_ret/2 - jobs)
                  + ' jobs and their settings will be ignored')

    return ret
