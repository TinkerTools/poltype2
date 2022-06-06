# This file is part of BINANA, released under the Apache 2.0 License. See
# LICENSE.md or go to https://opensource.org/licenses/Apache-2.0 for full
# details. Copyright 2021 Jacob D. Durrant.

from binana._utils.shim import _set_default


def hashtable_entry_add_one(hashtable, key, toadd=None):
    # note that dictionaries (hashtables) are passed by reference in
    # python

    # This is to keep track of the different kinds of interactions (counts), I
    # think.

    toadd = _set_default(toadd, 1)
    hashtable[key] = hashtable[key] + toadd if key in hashtable else toadd


def list_alphebetize_and_combine(list_obj):
    list_obj.sort()
    return "_".join(list_obj)
