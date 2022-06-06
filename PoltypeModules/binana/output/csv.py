# This file is part of BINANA, released under the Apache 2.0 License. See
# LICENSE.md or go to https://opensource.org/licenses/Apache-2.0 for full
# details. Copyright 2021 Jacob D. Durrant.

def _sanitize(s):
    return s.replace(",", "-").replace('"', "-")


def _recurse(pre_commas, data, csv):
    keys = list(data.keys())
    keys.sort()
    for key in keys:
        key = _sanitize(key)
        val = data[key]

        csv += pre_commas + key

        # __pragma__ ('skip')
        try:
            single_types = [int, float, str, unicode]
        except:
            single_types = [int, float, str]
        # __pragma__ ('noskip')

        """?
        single_types = [int, float, str]
        ?"""

        if type(val) in single_types:
            csv += "," + _sanitize(str(val)) + "\n"
        elif type(val) is list:
            if len(val) == 0:
                csv += ",none\n"
            elif type(val[0]) is dict:
                for i, item in enumerate(val):
                    csv += "\n"
                    if len(val) > 1:
                        csv += pre_commas + "," + key + "." + str(i + 1) + "\n"
                        csv = _recurse(pre_commas + ",,", item, csv)
                    else:
                        csv = _recurse(pre_commas + ",", item, csv)
            else:
                csv += "\n"
        elif type(val) is dict:
            csv += "\n"
            csv = _recurse(pre_commas + ",", val, csv)

    return csv


def collect(data):
    """Collects all the characterized interactions between the protein and
    ligand into a CSV-formatted string.

    Args:
        data (dict): A dictionary containing information about all the
            interactions. The output of 
            :py:func:`~binana.output.dictionary.collect`

    Returns:
        str: A CSV-formatted string containing the same information present in
        the input dictionary.
    """

    csv = _recurse("", data, "")

    while "\n\n" in csv:
        csv = csv.replace("\n\n", "\n")

    return csv
